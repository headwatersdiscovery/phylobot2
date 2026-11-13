"""The alignment viewer: per-column metrics ("tracks"), row order, and exports.

The page draws the grid and the tracks in a canvas from compact JSON (pilot jobs reach
500 sequences x 5000 columns, far too many cells for HTML); this module computes what
it draws. Every metric is computed on the amino-acid alignment's extant sequences
(ancestors shown as rows never enter them), one value per column, None where it is
undefined (e.g. an all-gap column).

Conservation follows Capra & Singh (2007, Bioinformatics 23:1875): Jensen-Shannon
divergence from the BLOSUM62 background distribution, Henikoff sequence weights, the
score scaled by the weighted fraction of non-gaps, columns over 30% gaps left unscored,
and each score averaged half-and-half with its window's (3 columns either side).
"""

import math
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from phylobot.trees import load_tree

GAP = "-"
AMINO = "ARNDCQEGHILKMFPSTWYV"  # BLOSUM62's order

BLOSUM62_ROWS = """
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62 = {
    (row[0], AMINO[j]): int(value)
    for row in (line.split() for line in BLOSUM62_ROWS.strip().splitlines())
    for j, value in enumerate(row[1:])
}
# BLOSUM62's background amino-acid distribution, as Capra & Singh use it (AMINO order).
BACKGROUND = dict(zip(AMINO, (
    0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062,
    0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072,
), strict=True))  # fmt: skip
JSD_PSEUDOCOUNT = 1e-6
JSD_GAP_CUTOFF = 0.3  # columns with more gaps than this are not scored
JSD_WINDOW = 3
JSD_WINDOW_WEIGHT = 0.5

# Kyte & Doolittle (1982) hydropathy.
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip
HYDROPATHY_WINDOW = 9
# Chou & Fasman (1978) conformational parameters: helix, strand, turn (1 = indifferent).
CHOU_FASMAN = {
    "A": (1.42, 0.83, 0.66), "R": (0.98, 0.93, 0.95), "N": (0.67, 0.89, 1.56),
    "D": (1.01, 0.54, 1.46), "C": (0.70, 1.19, 1.19), "Q": (1.11, 1.10, 0.98),
    "E": (1.51, 0.37, 0.74), "G": (0.57, 0.75, 1.56), "H": (1.00, 0.87, 0.95),
    "I": (1.08, 1.60, 0.47), "L": (1.21, 1.30, 0.59), "K": (1.14, 0.74, 1.01),
    "M": (1.45, 1.05, 0.60), "F": (1.13, 1.38, 0.60), "P": (0.57, 0.55, 1.52),
    "S": (0.77, 0.75, 1.43), "T": (0.83, 1.19, 0.96), "W": (1.08, 1.37, 0.96),
    "Y": (0.69, 1.47, 1.14), "V": (1.06, 1.70, 0.50),
}  # fmt: skip
CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0}  # at neutral pH; His counted as 0


def parse_fasta(text: str) -> dict[str, str]:
    """{name: aligned row}, upper case; names end at the first whitespace."""
    records: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0] if line[1:].strip() else ""
            records[name] = ""
        elif name is not None:
            records[name] += line.strip().upper().replace(".", GAP)
    return records


def columns(rows: Sequence[str]) -> list[tuple[str, ...]]:
    width = max((len(r) for r in rows), default=0)
    return list(zip(*(r.ljust(width, GAP) for r in rows), strict=True))


def henikoff_weights(cols: Sequence[Sequence[str]], n: int) -> list[float]:
    """Position-based sequence weights (Henikoff & Henikoff 1994), summing to 1; gaps
    count as a symbol, as in Capra & Singh's code."""
    weights = [0.0] * n
    for col in cols:
        counts = Counter(col)
        kinds = len(counts)
        for i, symbol in enumerate(col):
            weights[i] += 1.0 / (kinds * counts[symbol])
    total = sum(weights)
    return [w / total for w in weights] if total else [1.0 / n] * n


def jensen_shannon(col: Sequence[str], weights: Sequence[float]) -> float | None:
    """One column's conservation (0-1) before windowing; None if over the gap cutoff."""
    gapped = sum(w for symbol, w in zip(col, weights, strict=True) if symbol not in BACKGROUND)
    if sum(1 for s in col if s not in BACKGROUND) / len(col) > JSD_GAP_CUTOFF:
        return None
    freq = dict.fromkeys(AMINO, JSD_PSEUDOCOUNT)
    for symbol, w in zip(col, weights, strict=True):
        if symbol in freq:
            freq[symbol] += w
    total = sum(freq.values())
    divergence = 0.0
    for aa in AMINO:
        p, q = freq[aa] / total, BACKGROUND[aa]
        r = (p + q) / 2
        divergence += p * math.log2(p / r) + q * math.log2(q / r)
    return divergence / 2 * (1 - gapped)


def windowed(scores: Sequence[float | None], window: int = JSD_WINDOW) -> list[float | None]:
    """Capra & Singh's window: each score averaged with the mean of its neighbours
    (scored columns within `window` either side), weight JSD_WINDOW_WEIGHT."""
    if window <= 0:
        return list(scores)
    out: list[float | None] = []
    for i, score in enumerate(scores):
        if score is None:
            out.append(None)
            continue
        near = [
            s for j in range(max(0, i - window), min(len(scores), i + window + 1))
            if j != i and (s := scores[j]) is not None
        ]  # fmt: skip
        if near:
            score = (1 - JSD_WINDOW_WEIGHT) * score + JSD_WINDOW_WEIGHT * sum(near) / len(near)
        out.append(score)
    return out


def entropy(counts: Counter[str], states: int) -> float | None:
    """Shannon entropy of the counts, normalised by log(states): 0 = one state, 1 = all
    equally frequent."""
    total = sum(counts.values())
    if not total:
        return None
    h = -sum(c / total * math.log(c / total) for c in counts.values() if c)
    return h / math.log(states)


def mean(values: Iterable[float]) -> float | None:
    values = list(values)
    return sum(values) / len(values) if values else None


def sum_of_pairs(residues: Sequence[str]) -> float | None:
    """Mean BLOSUM62 score over every pair of residues in the column (gaps left out)."""
    counts = Counter(r for r in residues if r in BACKGROUND)
    n = sum(counts.values())
    if n < 2:
        return None
    total = 0
    kinds = list(counts)
    for i, a in enumerate(kinds):
        total += counts[a] * (counts[a] - 1) // 2 * BLOSUM62[a, a]
        for b in kinds[i + 1 :]:
            total += counts[a] * counts[b] * BLOSUM62[a, b]
    return total / (n * (n - 1) / 2)


def profile(row: str, window: int = HYDROPATHY_WINDOW) -> list[float | None]:
    """A sequence's sliding-window Kyte-Doolittle profile, placed on its alignment
    columns (None at its gaps): the mean over `window` residues centred on each."""
    positions = [i for i, r in enumerate(row) if r in KYTE_DOOLITTLE]
    values = [KYTE_DOOLITTLE[row[i]] for i in positions]
    half = window // 2
    out: list[float | None] = [None] * len(row)
    for k, column in enumerate(positions):
        near = values[max(0, k - half) : k + half + 1]
        out[column] = sum(near) / len(near)
    return out


@dataclass
class Series:
    label: str
    values: list[float | None]


@dataclass
class Track:
    """One row of plots above the alignment: a metric per column (several series are
    drawn as lines in one plot, e.g. each gap method's calls)."""

    key: str
    label: str
    description: str
    series: list[Series]
    low: float | None = None  # fixed axis range; None: from the data
    high: float | None = None
    text: list[str] = field(default_factory=list)  # a letter per column (consensus)

    def json(self) -> dict[str, Any]:
        def rounded(values: list[float | None]) -> list[float | None]:
            return [None if v is None else round(v, 4) for v in values]

        return {
            "key": self.key, "label": self.label, "description": self.description,
            "series": [{"label": s.label, "values": rounded(s.values)} for s in self.series],
            "low": self.low, "high": self.high, "text": "".join(self.text),
        }  # fmt: skip


@dataclass(frozen=True)
class TrackKind:
    label: str
    group: str
    description: str


# Every track the viewer offers, in menu order. Evolution and selection tracks need a
# combination's (or a follow-up job's) files and are offered only where they exist.
TRACKS = {
    "jsd": TrackKind("Conservation (Jensen-Shannon)", "Composition and conservation",
                     "Jensen-Shannon divergence from BLOSUM62's background (Capra & Singh "
                     "2007): 1 = conserved; weighted by the non-gap fraction, columns over "
                     "30% gaps unscored, window of 3."),
    "entropy": TrackKind("Shannon entropy", "Composition and conservation",
                         "Of the residues (gaps left out), normalised: 0 = one residue, 1 = "
                         "all 20 equally frequent."),
    "entropy_gaps": TrackKind("Shannon entropy, gaps as a state", "Composition and conservation",
                              "As Shannon entropy, with the gap as a 21st state."),
    "consensus": TrackKind("Consensus", "Composition and conservation",
                           "The most frequent residue and its share of the sequences."),
    "occupancy": TrackKind("Occupancy", "Composition and conservation",
                           "Share of sequences with a residue (not a gap)."),
    "distinct": TrackKind("Distinct residues", "Composition and conservation",
                          "Number of different amino acids in the column."),
    "blosum": TrackKind("BLOSUM62 sum of pairs", "Composition and conservation",
                        "Mean BLOSUM62 score over all pairs of residues in the column."),
    "hydropathy": TrackKind("Hydropathy (mean)", "Physico-chemical",
                            "Mean Kyte-Doolittle hydropathy of the column's residues "
                            "(positive = hydrophobic)."),
    "hydropathy_profile": TrackKind("Hydropathy profile of the reference", "Physico-chemical",
                                    "Kyte-Doolittle, mean over 9 residues of the reference "
                                    "sequence (the first sequence, or the chosen one)."),
    "helix": TrackKind("Helix propensity", "Physico-chemical",
                       "Mean Chou-Fasman helix parameter (above 1 favours helix)."),
    "strand": TrackKind("Strand propensity", "Physico-chemical",
                        "Mean Chou-Fasman strand parameter (above 1 favours strand)."),
    "turn": TrackKind("Turn propensity", "Physico-chemical",
                      "Mean Chou-Fasman turn parameter (above 1 favours turns)."),
    "charge": TrackKind("Charge", "Physico-chemical",
                        "Mean charge of the column's residues at neutral pH (K, R +1; D, E "
                        "-1)."),
    "rate": TrackKind("Site rate", "Evolution on the tree",
                      "IQ-TREE's empirical Bayes rate per column (mean 1), on the chosen "
                      "tree and model; needs rate heterogeneity (+G, +R)."),
    "site_lh": TrackKind("Site log-likelihood", "Evolution on the tree",
                         "Each column's log-likelihood on the chosen tree and model."),
    "changes": TrackKind("Substitutions on the tree", "Evolution on the tree",
                         "Branches whose two ends have different residues in the "
                         "reconstruction (most probable states; gaps left out)."),
    "ancestor_pp": TrackKind("Ancestors' certainty", "Evolution on the tree",
                             "Mean posterior probability of the ancestors' most probable "
                             "residues (IQ-TREE's reconstruction)."),
    "gap_calls": TrackKind("Gaps in ancestors", "Evolution on the tree",
                           "Share of ancestral nodes with a gap, per gap method."),
    "branch_site": TrackKind("Branch-site tests: P(positive)", "Selection",
                             "Each branch-site test's Bayes empirical Bayes probability of "
                             "positive selection on its foreground (codeml Model A)."),
}  # fmt: skip
DEFAULT_TRACKS = ("jsd",)  # owner decision 2026-11-11 (was jsd and entropy)
COMPOSITION = [k for k, kind in TRACKS.items() if kind.group in
               ("Composition and conservation", "Physico-chemical")]  # fmt: skip


def composition_tracks(rows: Sequence[str], reference: str = "") -> dict[str, Track]:
    """Every alignment-only track (composition, conservation, physico-chemical) from
    the extant rows; `reference` is the row for the hydropathy profile."""
    cols = columns(rows)
    n = len(rows)
    weights = henikoff_weights(cols, n) if n else []
    values: dict[str, list[float | None]] = {key: [] for key in COMPOSITION}
    consensus: list[str] = []
    jsd_raw: list[float | None] = []
    for col in cols:
        residues = [s for s in col if s in BACKGROUND]
        counts = Counter(residues)
        top = counts.most_common(1)
        consensus.append(top[0][0] if top else GAP)
        values["consensus"].append(top[0][1] / n if top else None)
        values["occupancy"].append(len(residues) / n)
        values["distinct"].append(float(len(counts)))
        values["entropy"].append(entropy(counts, len(AMINO)) if residues else None)
        with_gaps = counts + Counter({GAP: n - len(residues)})
        values["entropy_gaps"].append(entropy(with_gaps, len(AMINO) + 1))
        values["blosum"].append(sum_of_pairs(residues))
        values["hydropathy"].append(mean(KYTE_DOOLITTLE[r] for r in residues))
        for k, key in enumerate(("helix", "strand", "turn")):
            values[key].append(mean(CHOU_FASMAN[r][k] for r in residues))
        values["charge"].append(mean(CHARGE.get(r, 0.0) for r in residues))
        jsd_raw.append(jensen_shannon(col, weights))
    values["jsd"] = windowed(jsd_raw)
    values["hydropathy_profile"] = profile(reference or (rows[0] if rows else ""))
    ranges = {"jsd": (0.0, 1.0), "entropy": (0.0, 1.0), "entropy_gaps": (0.0, 1.0),
              "consensus": (0.0, 1.0), "occupancy": (0.0, 1.0), "distinct": (0.0, None),
              "helix": (0.0, 2.0), "strand": (0.0, 2.0), "turn": (0.0, 2.0),
              "hydropathy": (-4.5, 4.5), "hydropathy_profile": (-4.5, 4.5),
              "charge": (-1.0, 1.0)}  # fmt: skip
    tracks = {}
    for key in COMPOSITION:
        low, high = ranges.get(key, (None, None))
        kind = TRACKS[key]
        tracks[key] = Track(key, kind.label, kind.description,
                            [Series(kind.label, values[key])], low, high)  # fmt: skip
    tracks["consensus"].text = consensus
    return tracks


def tsv_column(tsv: str, column: str, width: int) -> list[float | None]:
    """A per-site TSV's column (site 1-based) as one value per alignment column."""
    lines = tsv.splitlines()
    if not lines:
        return [None] * width
    header = lines[0].split("\t")
    site, wanted = header.index("site"), header.index(column)
    out: list[float | None] = [None] * width
    for line in lines[1:]:
        fields = line.split("\t")
        i = int(fields[site]) - 1
        if 0 <= i < width:
            try:
                out[i] = float(fields[wanted])
            except ValueError:
                pass
    return out


def mean_pp(ancestors_pp: str, width: int) -> list[float | None]:
    """Mean posterior probability of the ML state per site, over every node
    (ancestors_pp.tsv: node, site, state, pp)."""
    total, count = [0.0] * width, [0] * width
    for line in ancestors_pp.splitlines()[1:]:
        _, site, _, pp = line.split("\t")
        i = int(site) - 1
        if 0 <= i < width:
            total[i] += float(pp)
            count[i] += 1
    return [t / c if c else None for t, c in zip(total, count, strict=True)]


def substitutions(newick: str, sequences: dict[str, str], width: int) -> list[float | None]:
    """Per column, the branches whose two ends differ (both residues; a gap at either
    end is not a substitution). `sequences` has every leaf's and node's row by name;
    branches whose ends lack a row are skipped."""
    tree = load_tree(newick)
    counts = [0.0] * width
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is None:
            continue
        child_name = node.taxon.label if node.taxon is not None else node.label
        a, b = sequences.get(parent.label or ""), sequences.get(child_name or "")
        if a is None or b is None:
            continue
        for i, (x, y) in enumerate(zip(a[:width], b[:width], strict=False)):
            if x != y and x != GAP and y != GAP:
                counts[i] += 1
    return list(counts)


def gap_shares(versions: dict[str, dict[str, str]], width: int) -> list[Series]:
    """Per gap version ({label: {node: row}}), the share of nodes with a gap per column."""
    out = []
    for label, rows in versions.items():
        values: list[float | None] = [None] * width
        if rows:
            n = len(rows)
            values = [sum(1 for r in rows.values() if i < len(r) and r[i] == GAP) / n
                      for i in range(width)]  # fmt: skip
        out.append(Series(label, values))
    return out


def tree_order(newick: str, outgroup: Iterable[str] = (), internal: bool = False) -> list[str]:
    """Names top to bottom as the tree viewer draws them: ladderized (smaller clades
    first), the outgroup's clade last. With `internal`, each internal node's label is
    listed just above its first descendant (preorder), for ancestors as rows."""
    tree = load_tree(newick)
    tree.ladderize(ascending=True)
    outgroup = set(outgroup)
    root = tree.seed_node
    kids = root.child_nodes()
    if outgroup:

        def holds_outgroup(child: Any) -> bool:
            return any(leaf.taxon.label in outgroup for leaf in child.leaf_iter())

        for child in [k for k in kids if holds_outgroup(k)]:
            root.remove_child(child)
            root.add_child(child)
    names = []
    for node in tree.preorder_node_iter():
        if node.taxon is not None:
            names.append(node.taxon.label)
        elif internal and node.label:
            names.append(node.label)
    return names


# Clustal's conservation line: '*' identical, ':' one strong group, '.' one weak group.
STRONG_GROUPS = ("STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW")
WEAK_GROUPS = ("CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK",
               "NEQHRK", "FVLIM", "HFY")  # fmt: skip


def conservation_mark(col: Sequence[str]) -> str:
    if GAP in col or not col:
        return " "
    kinds = set(col)
    if len(kinds) == 1:
        return "*"
    if any(kinds <= set(g) for g in STRONG_GROUPS):
        return ":"
    if any(kinds <= set(g) for g in WEAK_GROUPS):
        return "."
    return " "


def clustal(rows: dict[str, str], block: int = 60, protein: bool = True) -> str:
    """CLUSTAL W (.aln) text: blocks of `block` columns, conservation line for proteins."""
    names = list(rows)
    pad = max((len(n) for n in names), default=0) + 6
    width = max((len(r) for r in rows.values()), default=0)
    lines = ["CLUSTAL W multiple sequence alignment (PhyloBot)", "", ""]
    for start in range(0, width, block):
        chunk = {n: rows[n][start : start + block] for n in names}
        lines += [f"{n.ljust(pad)}{chunk[n]}" for n in names]
        if protein:
            marks = "".join(conservation_mark(c) for c in columns(list(chunk.values())))
            lines.append(" " * pad + marks)
        lines.append("")
    return "\n".join(lines) + "\n"


def tracks_tsv(tracks: Sequence[Track], width: int) -> str:
    """Every track's values, one row per alignment column."""
    header = ["column"]
    cells: list[list[str]] = [[str(i + 1)] for i in range(width)]
    for track in tracks:
        if track.text:
            header.append(f"{track.label}: residue")
            for i in range(width):
                cells[i].append(track.text[i] if i < len(track.text) else "")
        for series in track.series:
            name = track.label if len(track.series) == 1 else f"{track.label}: {series.label}"
            header.append(name)
            for i in range(width):
                v = series.values[i] if i < len(series.values) else None
                cells[i].append("" if v is None else f"{v:.6g}")
    return "\n".join("\t".join(r) for r in [header, *cells]) + "\n"
