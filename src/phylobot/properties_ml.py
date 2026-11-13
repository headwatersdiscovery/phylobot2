"""Neural property predictors: disorder and ESM-2 likelihood.

These need torch and model weights, which only the `phylobot-ml` image carries
(`docker/ml.Dockerfile`, the `ml` dependency group), so nothing here is imported until a
property is actually computed: `phylobot.properties` registers them with thin wrappers
and the web tier, the tools image and the tests never import torch. Every function takes
a *batch* of sequences, because one model call over many sequences is where the time
goes; `properties.draw_values` hands over each ancestor's distinct draws at once.

Both read the protein, not its alignment: gaps are dropped (`properties.residues`) and
so is anything outside the twenty standard amino acids (metapredict refuses U or X;
reconstructed ancestors never carry them, but a template or an AltAll might).

Measured 2026-10-21, one CPU thread, per 600-residue sequence: metapredict about 7 ms;
ESM-2 35M about 0.4 s native and 1.9 s under the amd64 emulation Docker Desktop uses on
Apple Silicon.
"""

import argparse
import math
import os
from collections.abc import Sequence
from functools import cache
from typing import Any

STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
GAP = "-"

#: metapredict V3's own defaults for calling IDRs (its `predict_disorder` signature), kept
#: here so the number we publish is reproducible from these values and the version alone.
DISORDER_VERSION = "V3"

#: ESM-2's context is 1024 tokens, two of them the start and end markers. Longer proteins
#: (unmasked ancestors exceed this: Ime2's MAFFT ancestors are ~1050 residues) are scored
#: in overlapping windows, each residue taken from the window where it sits furthest
#: from an edge, so every residue is scored once and with as much context as fits.
ESM_WINDOW = 1022
#: One window per forward pass. Batching buys nothing here — measured 2026-10-21, a batch
#: of 8 was no faster per sequence than a batch of 1 on CPU — and it costs a great deal of
#: memory, because attention is quadratic in length: a batch of 8 full windows peaked at
#: 2.15 GB, over the 2 GB a pipeline task declares. So: one at a time.
ESM_BATCH = 1


def protein(sequence: str) -> str:
    """The residues a predictor sees: no gaps, standard amino acids only."""
    return "".join(c for c in sequence.upper() if c in STANDARD)


def disorder_fraction(sequences: Sequence[str]) -> list[float]:
    """The share of each sequence's residues inside an IDR called by metapredict V3.

    metapredict's IDR calls (minimum size 12, gaps of up to 10 closed) rather than the
    raw per-residue scores against a cutoff: a lone residue above 0.5 in a folded domain
    is not a disordered region, and the domain calls are what metapredict's authors
    recommend reading. 0 for a sequence with no residues.
    """
    import metapredict

    cleaned = [protein(s) for s in sequences]
    present = [s for s in cleaned if s]
    by_sequence: dict[str, float] = {}
    if present:
        distinct = list(dict.fromkeys(present))
        calls = metapredict.predict_disorder(
            distinct, version=DISORDER_VERSION, return_domains=True, silence_warnings=True
        )
        for sequence, call in zip(distinct, calls, strict=True):
            inside = sum(end - start for start, end in call.disordered_domain_boundaries)
            by_sequence[sequence] = inside / len(sequence)
    return [by_sequence.get(s, 0.0) for s in cleaned]


@cache
def _esm() -> tuple[Any, Any]:
    """The model and its alphabet, loaded once per process from the baked-in weights."""
    import esm
    import torch

    torch.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    # fair-esm's checkpoints carry their training arguments as an argparse.Namespace,
    # which torch >= 2.6 refuses to unpickle unless it is declared safe.
    torch.serialization.add_safe_globals([argparse.Namespace])
    path = os.environ.get("PHYLOBOT_ESM_MODEL", "/opt/esm/esm2_t12_35M_UR50D.pt")
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(path)
    model.eval()
    return model, alphabet


def windows(length: int, size: int = ESM_WINDOW) -> list[tuple[int, int, int, int]]:
    """How a sequence of `length` is split for scoring: (start, end, keep_from, keep_to)
    per window, where [keep_from, keep_to) are the residues taken from that window. The
    kept ranges tile [0, length) exactly once."""
    if length <= size:
        return [(0, length, 0, length)]
    stride = size // 2
    starts = list(range(0, length - size, stride)) + [length - size]
    out = []
    for i, start in enumerate(starts):
        end = start + size
        # Hand over to the next window halfway through the overlap.
        keep_from = 0 if i == 0 else (starts[i - 1] + size + start) // 2
        keep_to = length if i == len(starts) - 1 else (end + starts[i + 1]) // 2
        out.append((start, end, keep_from, keep_to))
    return out


def esm_log_likelihood(sequences: Sequence[str]) -> list[float]:
    """ESM-2's mean per-residue log-probability of each sequence, from one unmasked pass.

    One forward pass with nothing masked, reading off log p(x_i) at every position: a
    sequence costs one pass instead of the L passes of a true pseudo-log-likelihood, and
    on Ime2 it separated real sequences from 10%-mutated and shuffled ones more cleanly
    than a 60-site masked estimate did (the probe report). It is an *intensive* measure
    (per residue), so a longer draw is not penalised for its length. Natural logarithm;
    higher is more plausible. NaN for a sequence with no residues.
    """
    import torch

    model, alphabet = _esm()
    converter = alphabet.get_batch_converter()
    cleaned = [protein(s) for s in sequences]
    pieces: list[tuple[str, int, int]] = []  # (window's residues, keep_from, keep_to)
    owners: list[int] = []
    distinct = list(dict.fromkeys(s for s in cleaned if s))
    for index, sequence in enumerate(distinct):
        for start, end, keep_from, keep_to in windows(len(sequence)):
            pieces.append((sequence[start:end], keep_from - start, keep_to - start))
            owners.append(index)
    totals = [0.0] * len(distinct)
    for first in range(0, len(pieces), ESM_BATCH):
        batch = pieces[first : first + ESM_BATCH]
        _, _, tokens = converter([(str(i), p[0]) for i, p in enumerate(batch)])
        with torch.inference_mode():
            logits = model(tokens)["logits"]
        scores = torch.log_softmax(logits, dim=-1).gather(-1, tokens.unsqueeze(-1)).squeeze(-1)
        for offset, (_, keep_from, keep_to) in enumerate(batch):
            # Token 0 is the start marker, so residue i is token i + 1.
            kept = scores[offset, keep_from + 1 : keep_to + 1]
            totals[owners[first + offset]] += float(kept.sum())
    means = {s: totals[i] / len(s) for i, s in enumerate(distinct)}
    return [means.get(s, math.nan) for s in cleaned]
