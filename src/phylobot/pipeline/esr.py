"""The ESR job's pipeline steps (`phylobot-esr`; engine in `phylobot.esr`).

Three steps, so every held-out tip is its own Nextflow task (cached on resume, run side
by side):

- `fit <staged> <out>`: checks the staged files (an alignment, a tree, a model and
  settings.json; uploaded ones are hostile), renames the sequences t1, t2, ... so that
  no name of theirs reaches IQ-TREE, fits the model on the tree (its branch lengths
  first, for an uploaded tree), and saves the fit with the IDs of the tips to hold out
  (tips.txt).
- `tip <fit> <id> <out>`: one held-out tip at each level asked for.
- `collect <fit> <out> <tip folders...>`: every tip's rows with the real names put back
  (esr_tips.tsv, esr_sites.tsv) and the job page's figures (esr_summary.json: per level,
  accuracy against mean posterior, and the reliability bins).
"""

import argparse
import csv
import json
import math
import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from phylobot import esr
from phylobot.pipeline.export import quote
from phylobot.seqio import read_fasta, write_fasta
from phylobot.trees import load_tree

SETTINGS = "settings.json"
TIPS = "tips.txt"
NAMES = "names.tsv"
SUMMARY = "esr_summary.json"
BINS = 10  # reliability bins of the ML residue's posterior, 0-0.1 ... 0.9-1
MIN_SEQUENCES = 4  # a held-out tip leaves at least three, so the tree still has an edge
SAFE_ID = re.compile(r"^t[0-9]{1,6}$")
TIP_ID = re.compile(r"(?<=[(,])t[0-9]{1,6}(?=[:,)])")
MODEL = re.compile(r"^[A-Za-z0-9.]+(\+[A-Za-z0-9]+)*$")
GAPS = set("-?.")


def _settings(staged: Path) -> dict[str, Any]:
    data = json.loads((staged / SETTINGS).read_text())
    if not isinstance(data, dict):
        raise ValueError("settings.json must hold an object")
    levels = data.get("levels")
    if not isinstance(levels, list) or not levels or set(levels) - set(esr.LEVELS):
        raise ValueError(f"levels must be some of {list(esr.LEVELS)}")
    tips = data.get("tips")
    if not isinstance(tips, list) or not tips or not all(isinstance(t, str) for t in tips):
        raise ValueError("tips must be a list of sequence names")
    seed = data.get("seed", 1)
    if not isinstance(seed, int) or isinstance(seed, bool) or not 1 <= seed < 2**31:
        raise ValueError("seed must be a positive integer")
    return {"levels": [lvl for lvl in esr.LEVELS if lvl in levels], "tips": tips, "seed": seed,
            "fit_lengths": bool(data.get("fit_lengths"))}  # fmt: skip


def check(alignment: dict[str, str], newick: str, model: str) -> None:
    """What ESR needs of its inputs, said as the web form says it (ValueError)."""
    if len(alignment) < MIN_SEQUENCES:
        raise ValueError(f"ESR needs at least {MIN_SEQUENCES} sequences")
    lengths = {len(s) for s in alignment.values()}
    if len(lengths) != 1:
        raise ValueError("the sequences are not aligned (their lengths differ)")
    residues = set("".join(alignment.values()).upper()) - GAPS
    if residues <= set("ACGTUN"):
        raise ValueError("ESR handles protein alignments only")
    if bad := sorted(residues - set(esr.AMINO_ACIDS) - set("XBZ")):
        raise ValueError(f"the alignment holds characters that are not amino acids: {''.join(bad)}")
    if not MODEL.match(model) or not esr.supported(model):
        raise ValueError(f"ESR cannot fix the parameters of model {model!r}")
    try:
        tree = load_tree(newick)
    except Exception as error:  # dendropy's parse errors
        if "Duplicate taxon" in str(error):
            raise ValueError("the tree names a sequence twice") from error
        raise ValueError("the tree is not a Newick tree that can be read") from error
    leaves = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    if set(leaves) != set(alignment):
        missing = sorted(set(alignment) - set(leaves))[:3]
        extra = sorted(set(leaves) - set(alignment))[:3]
        raise ValueError(
            "the tree and the alignment name different sequences"
            + (f"; not in the tree: {', '.join(missing)}" if missing else "")
            + (f"; not in the alignment: {', '.join(extra)}" if extra else "")
        )
    if set(esr.DUMMIES) & set(alignment):
        raise ValueError(
            f"sequence names reserved for ESR: {sorted(set(esr.DUMMIES) & set(alignment))}"
        )


def binary(newick: str) -> bool:
    """Whether every node has two children (a trifurcating top node allowed): E1 must
    prune a tip and merge the edges it leaves."""
    tree = load_tree(newick)
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        allowed = (2, 3) if node is tree.seed_node else (2,)
        if len(node.child_nodes()) not in allowed:
            return False
    return True


def rename(alignment: dict[str, str], newick: str) -> tuple[dict[str, str], str, dict[str, str]]:
    """The alignment and tree with sequences named t1, t2, ... (alignment order), and
    {id: name}."""
    ids = {name: f"t{i}" for i, name in enumerate(alignment, start=1)}
    tree = load_tree(newick)
    for leaf in tree.leaf_node_iter():
        leaf.taxon.label = ids[leaf.taxon.label]
    renamed = esr.plain(tree.as_string(schema="newick", suppress_rooting=True))
    return {ids[k]: v for k, v in alignment.items()}, renamed, {v: k for k, v in ids.items()}


def fit(staged: Path, out: Path) -> None:
    settings = _settings(staged)
    alignment = read_fasta(staged / "alignment.fasta")
    newick = (staged / "tree.nwk").read_text()
    model = (staged / "model.txt").read_text().strip()
    check(alignment, newick, model)
    if unknown := [t for t in settings["tips"] if t not in alignment]:
        raise ValueError(f"not in the alignment: {', '.join(unknown[:3])}")
    if "e1" in settings["levels"] and not binary(newick):
        raise ValueError("E1 needs a binary tree (it prunes each tip): this one has polytomies")
    renamed, tree, names = rename(alignment, newick)
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="esr-fit-") as tmp:
        given = Path(tmp) / "given"
        given.mkdir()
        write_fasta(given / "alignment.fasta", renamed)
        (given / "tree.nwk").write_text(tree)
        (given / "model.txt").write_text(model + "\n")
        work = Path(tmp) / "work"
        work.mkdir()
        run = esr.local_runner(settings["seed"])
        esr.Combination(given, work, run, fit_lengths=settings["fit_lengths"]).save(out)
    ids = {name: i for i, name in names.items()}
    (out / TIPS).write_text("".join(f"{ids[t]}\n" for t in dict.fromkeys(settings["tips"])))
    (out / NAMES).write_text("id\tname\n" + "".join(f"{i}\t{n}\n" for i, n in names.items()))
    kept = {"levels": settings["levels"], "seed": settings["seed"],
            "fit_lengths": settings["fit_lengths"]}  # fmt: skip
    (out / SETTINGS).write_text(json.dumps(kept, indent=2) + "\n")


def tip(fitted: Path, tip_id: str, out: Path) -> None:
    if not SAFE_ID.match(tip_id):
        raise ValueError(f"not a tip ID: {tip_id!r}")
    settings = json.loads((fitted / SETTINGS).read_text())
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="esr-tip-") as tmp:
        combination = esr.Combination.load(fitted, Path(tmp), esr.local_runner(settings["seed"]))
        if tip_id not in combination.alignment:
            raise ValueError(f"{tip_id} is not in the alignment")
        index = int(tip_id[1:])
        results = [getattr(combination, level)(index, tip_id) for level in settings["levels"]]
    esr.write_tsv(out / "esr_tips.tsv", esr.TIP_FIELDS, (r.summary() for r in results))
    esr.write_tsv(out / "esr_sites.tsv", esr.SITE_FIELDS, (s for r in results for s in r.sites))


def _read(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def summarise(sites: list[dict[str, str]], tips: list[dict[str, str]]) -> dict[str, Any]:
    """Per level: sites and tips scored, accuracy (the ML residue is the true one),
    mean posterior of the ML residue, mean log score, and the reliability bins."""
    levels: dict[str, Any] = {}
    for level in esr.LEVELS:
        rows = [s for s in sites if s["level"] == level]
        if not rows:
            continue
        bins: list[dict[str, Any]] = [
            {"low": i / BINS, "high": (i + 1) / BINS, "sites": 0, "correct": 0, "pp": 0.0}
            for i in range(BINS)
        ]
        correct = 0
        pp = logs = 0.0
        for row in rows:
            p = float(row["pp_ml"])
            hit = row["true"] == row["ml"]
            b = bins[min(int(p * BINS), BINS - 1)]
            b["sites"] += 1
            b["correct"] += hit
            b["pp"] += p
            correct += hit
            pp += p
            logs += math.log(max(float(row["pp_true"]), esr.LOG_FLOOR))
        for b in bins:
            n = b["sites"]
            b["accuracy"] = b["correct"] / n if n else None
            b["pp"] = b["pp"] / n if n else None
        n = len(rows)
        levels[level] = {
            "tips": sum(t["level"] == level for t in tips),
            "sites": n,
            "accuracy": correct / n,
            "mean_pp": pp / n,
            "log_score": logs / n,
            "bins": bins,
        }
    return {"levels": levels}


def collect(fitted: Path, out: Path, folders: Sequence[Path]) -> None:
    names = {row["id"]: row["name"] for row in _read(fitted / NAMES)}
    order = {line.strip(): i for i, line in enumerate((fitted / TIPS).read_text().split())}
    tips: list[dict[str, str]] = []
    sites: list[dict[str, str]] = []
    for folder in folders:
        tips += _read(folder / "esr_tips.tsv")
        sites += _read(folder / "esr_sites.tsv")
    rank = {lvl: i for i, lvl in enumerate(esr.LEVELS)}
    tips.sort(key=lambda r: (order.get(r["tip"], len(order)), rank[r["level"]]))
    sites.sort(key=lambda r: (order.get(r["tip"], len(order)), rank[r["level"]], int(r["site"])))
    summary = summarise(sites, tips)
    for row in (*tips, *sites):
        row["tip"] = names[row["tip"]]
    out.mkdir(parents=True, exist_ok=True)
    esr.write_tsv(out / "esr_tips.tsv", esr.TIP_FIELDS, tips)
    esr.write_tsv(out / "esr_sites.tsv", esr.SITE_FIELDS, sites)
    fit_record = json.loads((fitted / "fit.json").read_text())
    settings = json.loads((fitted / SETTINGS).read_text())
    summary |= {"fit": fit_record, "settings": settings}
    (out / SUMMARY).write_text(json.dumps(summary, indent=2) + "\n")
    # Written with the IDs, then each swapped for its quoted name, as EXPORT writes trees:
    # load_tree reads an underscore as itself, so a name must never be left unquoted.
    tree = esr.plain((fitted / "tree.nwk").read_text())
    (out / "tree.nwk").write_text(TIP_ID.sub(lambda m: quote(names[m.group(0)]), tree))


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="phylobot-esr", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("fit")
    p.add_argument("staged", type=Path)
    p.add_argument("out", type=Path)
    p = sub.add_parser("tip")
    p.add_argument("fit", type=Path)
    p.add_argument("tip")
    p.add_argument("out", type=Path)
    p = sub.add_parser("collect")
    p.add_argument("fit", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("folders", type=Path, nargs="*")
    args = parser.parse_args(argv)
    if args.step == "fit":
        fit(args.staged, args.out)
    elif args.step == "tip":
        tip(args.fit, args.tip, args.out)
    else:
        collect(args.fit, args.out, args.folders)


if __name__ == "__main__":
    main()
