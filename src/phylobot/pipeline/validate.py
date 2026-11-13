"""Validate user inputs and rename everything to safe IDs.

Treats all input as hostile. Only safe IDs (sequences s0001..., provided sources
provided1...) reach downstream tools; user names and labels are kept in names.tsv
and sources.tsv and restored at export.

    phylobot-validate manifest.json

The manifest (written by the workflow, never by the user) describes the staged files:
input.fasta (if sequences were given), alignment1.fasta..., tree1.nwk..., each with
its original file name and optional user label, plus the limits to enforce, the
sequence type ("protein", the default, or "codon" for coding DNA) and the genetic code
(NCBI number, default 1).

For coding DNA, sanitized.fasta holds the translations (every amino-acid step runs on
them unchanged) and sanitized.cds.fasta the coding sequences, without final stop
codons; a provided alignment must be a codon alignment (gaps in whole codons), written
translated to provided/alignments/ and as given to provided/codon_alignments/.

Writes sanitized.fasta, names.tsv, sources.tsv, validation.json, outgroup.txt (the
outgroup's safe IDs, one per line; empty if none was given),
modelfinder_args.txt (ModelFinder's candidate-set flags; empty = IQ-TREE's default set),
codon_modelfinder_args.txt (the same for codon ModelFinder),
and provided/alignments/providedN.fasta, provided/trees/providedN.nwk.
An alignment entry may carry a `key`, an aligner's key (`ALIGNMENT_KEYS`): an alignment
another job made with that aligner, inherited by a job that extends it. It is checked
like any other and keeps that key (provided/alignments/mafft.fasta), so the extension's
results name it as a job that ran the aligner itself would; unkeyed ones are numbered
provided1... among themselves.
With "species" in the manifest (reconciliation), species.tsv and species_tree.nwk are
staged too and checked (`phylobot.species.checks`), writing species.safe.tsv,
species_tree.safe.nwk and species_names.tsv.
Exits 1 with a message on invalid input.
"""

import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import dendropy

from phylobot import branch_site, catalog, codeml, codons
from phylobot.species import checks as species_checks
from phylobot.trees import canonical_split, newick_problem

AMINO = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS = set("BZJX")
GAPS = set("-.")
NUCLEOTIDE = set("ACGTUN")
SEQTYPES = ("protein", "codon")
MIN_SEQUENCES = 4
MAX_LABEL_LENGTH = 100
MAX_TREE_BYTES = 5_000_000
MAX_LISTED = 10  # names shown per error message
MAX_AMBIGUOUS_LISTED = 500  # ambiguous codons itemised in validation.json (all counted)
AMBIGUITY_WARNING = "ambiguous bases"  # how its warning starts (the web page shows its own)
# Aligners' keys, as main.nf's --aligners takes them: the keys an inherited alignment keeps.
ALIGNMENT_KEYS = ("mafft", "muscle", "clustalo", "prank", "prank_codon")


class InputError(Exception):
    pass


def listing(names: set[str] | list[str]) -> str:
    shown = sorted(names)[:MAX_LISTED]
    more = len(names) - len(shown)
    return ", ".join(repr(n) for n in shown) + (f" and {more} more" if more > 0 else "")


def parse_fasta(text: str, what: str, unique: bool = True) -> list[tuple[str, str, str]]:
    """Return (name, description, sequence text) records; with `unique` (the default),
    a name used twice is refused, saying whether the two records are the same sequence."""
    records: list[tuple[str, str, list[str]]] = []
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            if not header:
                raise InputError(f"{what}, line {number}: empty sequence name")
            name, _, description = header.partition(" ")
            if any(ord(c) < 32 for c in name):
                raise InputError(f"{what}: sequence name {name!r} contains control characters")
            records.append((name, description.strip(), []))
        elif not records:
            raise InputError(f"{what}, line {number}: sequence data before the first '>' header")
        else:
            records[-1][2].append(line)
    out = [(n, d, "".join("".join(parts).split())) for n, d, parts in records]
    seen: dict[str, str] = {}
    for name, _, seq in out if unique else []:
        residues = seq.upper().replace("-", "").replace(".", "")
        if name in seen:
            same = seen[name] == residues
            raise InputError(
                f"{what}: duplicate sequence name {name!r} "
                + ("(the same sequence twice: remove one)" if same else
                   "(two different sequences: rename one)")
            )  # fmt: skip
        seen[name] = residues
    return out


def read_text(path: Path, what: str) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise InputError(f"{what} is not UTF-8 text") from exc


def check_residues(name: str, seq: str, what: str, allowed: set[str]) -> None:
    bad = sorted(set(seq) - allowed)
    if bad:
        raise InputError(f"{what}: sequence {name!r} contains invalid characters: {''.join(bad)}")


def read_sequences(path: Path, limits: dict[str, int]) -> list[tuple[str, str, str]]:
    what = "sequences"
    records = parse_fasta(read_text(path, what), what)
    cleaned = []
    for name, description, raw in records:
        seq = raw.upper().replace("-", "").replace(".", "").rstrip("*")
        if not seq:
            raise InputError(f"{what}: sequence {name!r} is empty")
        check_residues(name, seq, what, AMINO | AMBIGUOUS)
        cleaned.append((name, description, seq))
    return cleaned


def coding_record(
    name: str, raw: str, code: codons.GeneticCode, what: str, notes: dict[str, list[str]]
) -> tuple[str, str]:
    """(coding sequence, protein) for one ungapped coding sequence; notes collects
    remarks (final stop removed) by kind, for one aggregated warning. IUPAC ambiguity
    codes are accepted (`codons.translate_codon`; `ambiguity_report` lists them)."""
    seq = raw.upper().replace("-", "").replace(".", "")
    if not seq:
        raise InputError(f"{what}: sequence {name!r} is empty")
    bad = sorted(set(seq) - codons.DNA)
    if bad:
        raise InputError(
            f"{what}: sequence {name!r} contains characters that are not nucleotide (IUPAC) "
            f"codes: {''.join(bad)}"
        )
    try:
        cds, protein, remarks = codons.check_cds(seq, code)
    except codons.CodonError as exc:
        raise InputError(f"{what}: sequence {name!r}: {exc}") from None
    for remark in remarks:
        notes.setdefault(remark, []).append(name)
    return cds, protein


def ambiguity_report(cds: dict[str, str], code: codons.GeneticCode) -> dict[str, Any] | None:
    """validation.json's "ambiguous_codons" for coding sequences ({user name: coding
    sequence}): counts, the first MAX_AMBIGUOUS_LISTED codons, and a fingerprint the New
    Job page's confirmation must match (the same sequences, codons and genetic code).
    None if no codon has an ambiguity code."""
    found = [(name, a) for name, seq in cds.items() for a in codons.ambiguous_codons(seq, code)]
    if not found:
        return None
    key = json.dumps([code.ncbi, [(n, a.number, a.codon) for n, a in found]])
    return {
        "codons": len(found),
        "sequences": len({n for n, _ in found}),
        "unknown": sum(a.unknown for _, a in found),
        "translated": sum(a.translation != "X" for _, a in found),
        "fingerprint": hashlib.sha256(key.encode()).hexdigest()[:16],
        "listed": [
            {
                "sequence": n,
                "number": a.number,
                "codon": a.codon,
                "aminos": a.aminos,
                "translation": a.translation,
            }
            for n, a in found[:MAX_AMBIGUOUS_LISTED]
        ],
    }


def check_looks_like_dna(texts: list[str], what: str) -> None:
    letters = "".join(texts).upper().replace("-", "").replace(".", "")
    if letters and sum(c in NUCLEOTIDE for c in letters) / len(letters) < 0.9:
        raise InputError(
            f"{what} look like protein, not coding DNA; choose protein input "
            "(pipeline: --seqtype protein)"
        )


def read_coding_sequences(
    path: Path, code: codons.GeneticCode, notes: dict[str, list[str]]
) -> tuple[list[tuple[str, str, str]], dict[str, str]]:
    """(records with each sequence's translation, {name: coding sequence})."""
    what = "sequences"
    records = parse_fasta(read_text(path, what), what)
    check_looks_like_dna([raw for _, _, raw in records], what)
    translated, cds = [], {}
    for name, description, raw in records:
        cds[name], protein = coding_record(name, raw, code, what, notes)
        translated.append((name, description, protein))
    return translated, cds


def codon_alignment(
    rows: dict[str, str], code: codons.GeneticCode, what: str, notes: dict[str, list[str]]
) -> tuple[dict[str, str], dict[str, str]]:
    """({name: codon row}, {name: protein row}) for a provided codon alignment. A final
    stop codon is removed from each row, then columns left all gaps are dropped."""
    check_looks_like_dna(list(rows.values()), what)
    trimmed = {}
    for name, row in rows.items():
        last = len(row.rstrip("-")) // 3 * 3 - 3
        if last >= 0:
            if codons.translate_codon(row[last : last + 3], code) == "*":
                row = row[:last] + "---" + row[last + 3 :]
                notes.setdefault("final stop codon removed", []).append(name)
        trimmed[name] = row
    length = len(next(iter(trimmed.values())))
    if length % 3:
        raise InputError(f"{what}: aligned length {length} is not a multiple of 3")
    keep = [
        i for i in range(0, length, 3) if any(row[i : i + 3] != "---" for row in trimmed.values())
    ]
    codon_rows = {name: "".join(row[i : i + 3] for i in keep) for name, row in trimmed.items()}
    protein_rows = {}
    for name, row in codon_rows.items():
        try:
            protein = codons.translate_row(row, code)
        except codons.CodonError as exc:
            raise InputError(f"{what}: sequence {name!r}: {exc}") from None
        if (stop := protein.find("*")) >= 0:
            raise InputError(
                f"{what}: sequence {name!r}: internal stop codon {row[3 * stop : 3 * stop + 3]} "
                f"in codon column {stop + 1} (genetic code {code.ncbi}, {code.name})"
            )
        protein_rows[name] = protein
    return codon_rows, protein_rows


def detect_seqtype(fasta: str) -> str:
    """'codon' if the sequences look like DNA (90% A/C/G/T/U/N), else 'protein'; the
    New Job page's guess, which the user can override."""
    letters = [
        c for line in fasta.splitlines() if not line.startswith(">") for c in line.upper()
        if c.isalpha()
    ]  # fmt: skip
    if letters and sum(c in NUCLEOTIDE for c in letters) / len(letters) >= 0.9:
        return "codon"
    return "protein"


def check_counts(
    records: list[tuple[str, str, str]], limits: dict[str, int], coding: bool = False
) -> None:
    if len(records) < MIN_SEQUENCES:
        raise InputError(f"need at least {MIN_SEQUENCES} sequences, found {len(records)}")
    if len(records) > limits["max_sequences"]:
        raise InputError(
            f"at most {limits['max_sequences']} sequences allowed, found {len(records)}"
        )
    unit = "codons" if coding else "residues"
    for name, _, seq in records:
        if len(seq) > limits["max_length"]:
            raise InputError(
                f"sequence {name!r} is {len(seq)} {unit}; limit is {limits['max_length']}"
            )
    if coding:
        return  # the records are translations
    residues = "".join(seq for _, _, seq in records)
    if sum(c in NUCLEOTIDE for c in residues) / len(residues) > 0.9:
        raise InputError(
            "input looks like nucleotide sequence; for coding DNA choose coding-DNA input "
            "(pipeline: --seqtype codon)"
        )


def read_alignment(
    path: Path, what: str, warnings: list[str], coding: bool = False
) -> tuple[dict[str, str], dict[str, str]]:
    """Return ({name: aligned row}, {name: description}); rows uppercased, '.' -> '-'."""
    records = parse_fasta(read_text(path, what), what)
    rows: dict[str, str] = {}
    lowercase = 0
    for name, _, raw in records:
        lowercase += sum(c.islower() for c in raw) if not coding else 0
        row = raw.upper().replace(".", "-")
        check_residues(name, row, what, (codons.DNA if coding else AMINO | AMBIGUOUS) | {"-"})
        rows[name] = row
    lengths = {len(row) for row in rows.values()}
    if len(lengths) > 1:
        raise InputError(f"{what}: rows have different lengths ({sorted(lengths)}); not aligned")
    if lowercase:
        warnings.append(
            f"{what}: {lowercase} lowercase residues treated as uppercase "
            "(some aligners use lowercase for unaligned regions)"
        )
    return rows, {name: description for name, description, _ in records}


def check_same_names(names: set[str], expected: set[str], what: str) -> None:
    missing, extra = expected - names, names - expected
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing {listing(missing)}")
        if extra:
            parts.append(f"not in the sequences: {listing(extra)}")
        raise InputError(f"{what}: names do not match the sequences; " + "; ".join(parts))


def check_newick_syntax(text: str, what: str) -> None:
    """Reject malformed Newick that parsers tend to repair silently (`trees.newick_problem`)."""
    if problem := newick_problem(text):
        raise InputError(f"{what}: {problem}")


def check_branch_lengths(tree: dendropy.Tree, what: str) -> None:
    """A tree whose branch lengths are kept as given (--tree_branch_lengths keep): every
    branch but the root's needs one, and none may be negative."""
    missing = negative = 0
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        length = node.edge.length
        if length is None:
            missing += 1
        elif length < 0:
            negative += 1
    if missing:
        raise InputError(
            f"{what}: {missing} branch{'es' * (missing > 1)} without a length, "
            "and its branch lengths are to be kept as given"
        )
    if negative:
        raise InputError(f"{what}: {negative} negative branch length{'s' * (negative > 1)}")


def read_tree(
    path: Path,
    what: str,
    names: set[str],
    marks: list[tuple[dendropy.Node, str]] | None = None,
) -> dendropy.Tree:
    """A provided tree, checked against the sequence names. With `marks` (a list to
    fill), codeml's foreground marks (#1, $1) are allowed and collected as (node, scope)."""
    if path.stat().st_size > MAX_TREE_BYTES:
        raise InputError(f"{what} is larger than {MAX_TREE_BYTES // 1_000_000} MB")
    text = read_text(path, what)
    if marks is not None:
        try:
            text = branch_site.encode_marks(text)
        except branch_site.MarkError as exc:
            raise InputError(f"{what}: {exc}") from None
    depth = deepest = 0
    for c in text:  # a tree over n leaves nests at most n deep; reject pathological input
        depth += (c == "(") - (c == ")")
        deepest = max(deepest, depth)
    if deepest > len(names) + 1:
        raise InputError(f"{what}: nesting is deeper than possible for {len(names)} sequences")
    check_newick_syntax(text, what)
    try:
        # Leaf names stay plain node labels (no taxa), so every name check below is ours.
        tree = dendropy.Tree.get(
            data=text,
            schema="newick",
            preserve_underscores=True,
            suppress_leaf_node_taxa=True,
        )
    except Exception as exc:
        raise InputError(f"{what} is not valid Newick: {exc}") from exc
    if marks is not None:
        for node in tree.preorder_node_iter():
            node.label, scope = branch_site.strip_mark(node.label)
            if scope is not None:
                marks.append((node, scope))
    labels = [leaf.label or "" for leaf in tree.leaf_node_iter()]
    if "" in labels:
        raise InputError(f"{what}: a leaf has no name")
    duplicates = {label for label in labels if labels.count(label) > 1}
    if duplicates:
        raise InputError(f"{what}: duplicate leaf names {listing(duplicates)}")
    check_same_names(set(labels), names, what)
    return tree


def foreground_tree(
    tree: dendropy.Tree, newick: str, marks: list[tuple[dendropy.Node, str]], warnings: list[str]
) -> str:
    """The marked tree for codeml (leaves already renamed to safe IDs): each mark's
    branches, as unrooted splits, marked #1 on the unrooted tree."""
    foreground: set[branch_site.Split] = set()
    try:
        for node, scope in marks:
            foreground |= branch_site.splits_of(tree, node, scope)
        marked = branch_site.mark_tree(newick, foreground)
    except branch_site.MarkError as exc:
        raise InputError(f"tree 1: {exc}") from None
    root = tree.seed_node.child_nodes()
    everyone = branch_site.leaves_below(tree.seed_node)
    if len(root) == 2:
        sides = {canonical_split(branch_site.leaves_below(c), everyone) for c in root}
        if sides & foreground:
            warnings.append(
                "tree 1 is rooted: codeml analyses it unrooted, so the two branches at "
                "the root are one foreground branch"
            )
    return marked


def clean_label(label: Any, fallback: str, what: str) -> str:
    text = str(label).strip() if label not in (None, "") else fallback
    if any(ord(c) < 32 or ord(c) == 127 for c in text):
        raise InputError(f"{what}: label contains control characters")
    if len(text) > MAX_LABEL_LENGTH:
        raise InputError(f"{what}: label is longer than {MAX_LABEL_LENGTH} characters")
    return text


def alignment_keys(entries: list[dict[str, Any]], coding: bool) -> list[str]:
    """Each provided alignment's key: its own `key` (an inherited alignment), else
    provided1... counted among the unkeyed ones."""
    keys: list[str] = []
    for index, entry in enumerate(entries, 1):
        key = entry.get("key")
        if key in (None, ""):
            keys.append(f"provided{sum(1 for k in keys if k.startswith('provided')) + 1}")
            continue
        if key not in ALIGNMENT_KEYS:
            raise InputError(f"alignment {index}: unknown aligner key {key!r}")
        if key == "prank_codon" and not coding:
            raise InputError(f"alignment {index}: a prank_codon alignment needs coding DNA")
        if key in keys:
            raise InputError(f"alignment {index}: aligner key {key!r} is used twice")
        keys.append(key)
    return keys


def validate(manifest_path: Path, workdir: Path = Path(".")) -> None:
    """Validate the staged files in `workdir` (default: the current directory, as in a
    pipeline task) and write the outputs there."""
    manifest = json.loads(manifest_path.read_text())
    limits = {key: int(manifest[key]) for key in ("max_sequences", "max_length", "max_provided")}
    alignments: list[dict[str, Any]] = manifest.get("alignments", [])
    trees: list[dict[str, Any]] = manifest.get("trees", [])
    own = [a for a in alignments if a.get("key") in (None, "")]  # not inherited
    for kind, entries in (("alignments", own), ("trees", trees)):
        if len(entries) > limits["max_provided"]:
            raise InputError(f"at most {limits['max_provided']} provided {kind} allowed")
    if not manifest.get("sequences") and not alignments:
        raise InputError("no sequences: give an input FASTA or a provided alignment")
    models: list[str] = manifest.get("models") or []
    # ModelFinder's candidates; absent (manifests before they had their own list): models.
    candidates: list[str] = (
        manifest["candidates"] if "candidates" in manifest else catalog.candidate_models(models)
    ) or []
    for model in models + candidates:
        if not catalog.is_model(model):
            raise InputError(f"unknown substitution model {model!r}")
    if mixtures := [m for m in candidates if catalog.is_mixture(m)]:
        raise InputError(
            f"site-profile mixtures are not ModelFinder candidates: {listing(mixtures)}"
        )
    codon_models: list[str] = manifest.get("codon_models") or []
    codon_candidates: list[str] = manifest.get("codon_candidates") or []
    for model in codon_models + codon_candidates:
        if not catalog.is_codon_model(model):
            raise InputError(f"unknown codon model {model!r}")
    codon_modelfinder = bool(manifest.get("codon_modelfinder"))
    codeml_models: list[str] = manifest.get("codeml_models") or []
    # The branch-site test (a follow-up job): one codon alignment and one tree with
    # codeml's foreground marks; nothing else runs.
    branch_site_test = bool(manifest.get("branch_site"))
    if branch_site_test and (len(alignments) != 1 or len(trees) != 1):
        raise InputError("the branch-site test needs one codon alignment and one marked tree")
    unknown_codeml = sorted(set(codeml_models) - set(codeml.MODELS))
    if unknown_codeml:
        raise InputError(f"unknown codeml models {listing(unknown_codeml)}")
    raxml = bool(manifest.get("raxml"))
    if raxml and (mixtures := [m for m in models if catalog.is_mixture(m)]):
        raise InputError(
            f"RAxML-NG has no site-profile mixtures (ticked: {', '.join(mixtures)}); untick "
            "those models or RAxML-NG"
        )
    if raxml and (missing := catalog.raxml_missing(list(dict.fromkeys(models + candidates)))):
        raise InputError(
            f"RAxML-NG does not have the {', '.join(sorted({m.split('+')[0] for m in missing}))} "
            f"matrix (ticked: {', '.join(missing)}); untick those models or RAxML-NG"
        )
    modelfinder_args = catalog.modelfinder_args(candidates, raxml=raxml)
    seqtype = str(manifest.get("seqtype") or "protein")
    if seqtype not in SEQTYPES:
        raise InputError(f"unknown sequence type {seqtype!r} (choose {' or '.join(SEQTYPES)})")
    coding = seqtype == "codon"
    keys = alignment_keys(alignments, coding)
    if (codon_models or codon_modelfinder or codeml_models or branch_site_test) and not coding:
        raise InputError("codon models need coding DNA")

    try:
        code = codons.genetic_code(manifest.get("genetic_code") or 1)
    except codons.CodonError as exc:
        raise InputError(str(exc)) from None
    if (codeml_models or branch_site_test) and code.paml_icode is None:
        raise InputError(
            f"PAML does not have genetic code {code.ncbi} ({code.name}): untick codeml or "
            "choose a code PAML offers"
        )

    warnings: list[str] = []
    notes: dict[str, list[str]] = {}  # coding-sequence remarks, by kind
    provided_rows: list[dict[str, str]] = []  # protein rows (translated, for coding DNA)
    provided_codon_rows: list[dict[str, str]] = []
    for index, _ in enumerate(alignments, 1):
        path = workdir / f"alignment{index}.fasta"
        what = f"alignment {index}"
        rows, _ = read_alignment(path, what, warnings, coding)
        if coding:
            codon_rows, rows = codon_alignment(rows, code, what, notes)
            provided_codon_rows.append(codon_rows)
        provided_rows.append(rows)

    cds: dict[str, str] = {}  # coding sequences, for coding DNA
    if manifest.get("sequences"):
        if coding:
            records, cds = read_coding_sequences(workdir / "input.fasta", code, notes)
        else:
            records = read_sequences(workdir / "input.fasta", limits)
    else:
        # Sequences come from the first provided alignment, with gaps removed.
        first = read_alignment(workdir / "alignment1.fasta", "alignment 1", [], coding)
        if coding:
            first_rows = provided_codon_rows[0]
            records = []
            for n in first_rows:
                cds[n], protein = coding_record(
                    n, first_rows[n], code, "alignment 1", {}
                )  # its notes were collected with the alignment
                records.append((n, first[1][n], protein))
        else:
            records = [(n, first[1][n], row.replace("-", "")) for n, row in first[0].items()]
        for name, _, seq in records:
            if not seq:
                raise InputError(f"alignment 1: sequence {name!r} is all gaps")
    check_counts(records, limits, coding)
    sequences = {name: seq for name, _, seq in records}
    names = set(sequences)

    compared = provided_codon_rows if coding else provided_rows
    expected = cds if coding else sequences
    for index, rows in enumerate(compared, 1):
        what = f"alignment {index}"
        check_same_names(set(rows), names, what)
        differ = {n for n, row in rows.items() if row.replace("-", "") != expected[n]}
        if differ:
            raise InputError(
                f"{what}: with gaps removed, these rows differ from the sequences: "
                f"{listing(differ)}"
            )
    for kind, which in notes.items():
        unique = sorted(set(which))
        warnings.append(
            f"{kind}: {len(unique)} sequence{'s' * (len(unique) > 1)} ({listing(unique)})"
        )
    ambiguity = ambiguity_report(cds, code) if coding else None
    if ambiguity:
        many = ambiguity["codons"] > 1
        warnings.append(
            f"{AMBIGUITY_WARNING}: {ambiguity['codons']} codon{'s' * many} in "
            f"{ambiguity['sequences']} sequence{'s' * (ambiguity['sequences'] > 1)} "
            "(translated as the amino acid every reading gives, else X; codeml weighs each "
            "possible codon, IQ-TREE's codon models treat the codon as unknown)"
        )

    safe = {name: f"s{index:04d}" for index, (name, _, _) in enumerate(records, 1)}
    outgroup: list[str] = [str(name) for name in manifest.get("outgroup") or []]
    unknown = sorted(set(outgroup) - names)
    if unknown:
        raise InputError(f"outgroup sequences not among the sequences: {listing(unknown)}")
    if outgroup and len(names) - len(set(outgroup)) < 3:
        raise InputError("the outgroup must leave at least 3 ingroup sequences")
    (workdir / "outgroup.txt").write_text("".join(f"{safe[n]}\n" for n in sorted(set(outgroup))))
    species_report = None
    if manifest.get("species"):  # reconciliation: species.tsv and species_tree.nwk staged
        try:
            prepared = species_checks.prepare(
                read_text(workdir / "species.tsv", "species table"),
                read_text(workdir / "species_tree.nwk", "species tree"),
                safe,
            )
        except species_checks.SpeciesError as exc:
            raise InputError(str(exc)) from None
        (workdir / "species.safe.tsv").write_text(prepared.mapping)
        (workdir / "species_tree.safe.nwk").write_text(prepared.tree)
        (workdir / "species_names.tsv").write_text(prepared.names)
        species_report = prepared.report
    marks: list[tuple[dendropy.Node, str]] | None = [] if branch_site_test else None
    parsed_trees = [
        read_tree(workdir / f"tree{index}.nwk", f"tree {index}", names, marks)
        for index, _ in enumerate(trees, 1)
    ]
    if marks is not None and not marks:
        raise InputError("tree 1: no foreground branch is marked (add #1 after a node)")
    for index, (entry, tree) in enumerate(zip(trees, parsed_trees, strict=True), 1):
        if entry.get("keep_lengths"):
            check_branch_lengths(tree, f"tree {index}")

    with open(workdir / "sanitized.fasta", "w") as fasta, open(workdir / "names.tsv", "w") as out:
        out.write("id\tname\tdescription\n")
        for name, description, seq in records:
            fasta.write(f">{safe[name]}\n{seq}\n")
            out.write(f"{safe[name]}\t{name}\t{description}\n")
    if coding:
        with open(workdir / "sanitized.cds.fasta", "w") as fasta:
            for name, _, _ in records:
                fasta.write(f">{safe[name]}\n{cds[name]}\n")

    sources = ["key\tkind\tlabel\tfilename"]
    (workdir / "provided/alignments").mkdir(parents=True, exist_ok=True)
    (workdir / "provided/trees").mkdir(parents=True, exist_ok=True)
    seen_labels: dict[str, set[str]] = {"alignment": set(), "tree": set()}
    for kind, entries in (("alignment", alignments), ("tree", trees)):
        for index, entry in enumerate(entries, 1):
            what = f"{kind} {index}"
            filename = clean_label(entry.get("filename"), f"{kind}{index}", what)
            label = clean_label(entry.get("label"), filename, what)
            if label.casefold() in seen_labels[kind]:
                raise InputError(f"{what}: label {label!r} is used twice")
            seen_labels[kind].add(label.casefold())
            key = keys[index - 1] if kind == "alignment" else f"provided{index}"
            sources.append(f"{key}\t{kind}\t{label}\t{filename}")

    for key, rows in zip(keys, provided_rows, strict=True):
        with open(workdir / f"provided/alignments/{key}.fasta", "w") as handle:
            for name, _, _ in records:  # same order as the sequences
                handle.write(f">{safe[name]}\n{rows[name]}\n")
    if coding:
        (workdir / "provided/codon_alignments").mkdir(parents=True, exist_ok=True)
    for key, rows in zip(keys if coding else [], provided_codon_rows, strict=True):
        with open(workdir / f"provided/codon_alignments/{key}.fasta", "w") as handle:
            for name, _, _ in records:
                handle.write(f">{safe[name]}\n{rows[name]}\n")
    for index, tree in enumerate(parsed_trees, 1):
        for node in tree.preorder_node_iter():
            node.label = safe[node.label] if node.is_leaf() else None
        newick = tree.as_string(
            schema="newick",
            suppress_rooting=True,
            suppress_leaf_node_labels=False,
            suppress_internal_node_labels=True,
        ).strip()
        if marks:  # codeml's tree: unrooted, with #1 on each foreground branch
            newick = foreground_tree(tree, newick, marks, warnings)
        (workdir / f"provided/trees/provided{index}.nwk").write_text(newick + "\n")

    (workdir / "sources.tsv").write_text("\n".join(sources) + "\n")
    (workdir / "modelfinder_args.txt").write_text(" ".join(modelfinder_args) + "\n")
    if codeml_models:
        assert code.paml_icode is not None
        codeml.write_controls(workdir / "codeml", codeml_models, code.paml_icode)
    if branch_site_test:
        assert code.paml_icode is not None
        branch_site.write_controls(workdir / "codeml", code.paml_icode)
    codon_args = catalog.codon_modelfinder_args(codon_candidates)
    (workdir / "codon_modelfinder_args.txt").write_text(" ".join(codon_args) + "\n")
    by_sequence: dict[str, list[str]] = {}
    for name, _, seq in records:
        by_sequence.setdefault(cds[name] if coding else seq, []).append(name)
    report = {
        "seqtype": seqtype,
        **({"genetic_code": {"ncbi": code.ncbi, "name": code.name}} if coding else {}),
        "sequences": len(records),
        "min_length": min(len(s) for s in sequences.values()),
        "max_length": max(len(s) for s in sequences.values()),
        "identical_sequence_groups": [g for g in by_sequence.values() if len(g) > 1],
        "provided_alignments": len(alignments),
        "provided_trees": len(trees),
        "outgroup": sorted(set(outgroup)),
        "modelfinder_candidates": (
            catalog.modelfinder_candidates(candidates).models()
            if not catalog.is_default_set(candidates, raxml=raxml)
            else "IQ-TREE default"
            + (" without " + ", ".join(sorted(catalog.RAXML_MISSING)) if raxml else "")
        ),
        **({"species": species_report} if species_report is not None else {}),
        **({"ambiguous_codons": ambiguity} if ambiguity else {}),
        "warnings": warnings,
    }
    (workdir / "validation.json").write_text(json.dumps(report, indent=2) + "\n")


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    try:
        validate(Path(sys.argv[1]))
    except InputError as exc:
        print(f"Input error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
