"""What a library hands to the bench: the protein FASTA, the manifest, the
names table and the methods paragraph.

Pure text in, pure text out. The view gathers the members and what each came from; nothing
here reads storage or the database, so every document can be checked against a handful of
members in a unit test.

An export is a *frozen* thing: its files and the spec that produced them are kept under
`libraries/<id>/exports/<n>/` so a set that was ordered can be read back months later, and
restored into the editable library. The design note's reasoning: a library is an input to
somebody else's pipeline, so export fidelity matters more than anything on the page.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from phylobot import library

FASTA = "members.fasta"
MANIFEST = "manifest.tsv"
NAMES = "names.tsv"
METHODS = "methods.txt"
SPEC = "spec.json"
FILES = (FASTA, MANIFEST, NAMES, METHODS, SPEC)
LABELS = {
    FASTA: "Protein FASTA",
    MANIFEST: "Manifest (TSV)",
    NAMES: "Names (TSV)",
    METHODS: "Methods paragraph",
    SPEC: "The design, as saved",
}
WRAP = 60  # FASTA line length, as every other tool in the pipeline writes it


@dataclass(frozen=True)
class Source:
    """Where one member came from, named as a person would say it."""

    combination: str  # "mafft/iqtree/modelfinder"
    label: str  # "MAFFT x IQ-TREE x ModelFinder"
    node: str
    version: str  # the ancestor viewer's key
    version_label: str  # "Parsimony (shared presence)"
    short_label: str = ""  # "Parsimony (shared)", for prose rather than a table cell
    weight: float = 0.0


@dataclass(frozen=True)
class Row:
    """One member as the export writes it: the sequence, its name, and its provenance."""

    name: str
    member: library.Member
    source: Source
    template: str  # its ancestor's ML row, for residue numbering and the difference count
    selected: bool
    mean_pp: float | None = None

    @property
    def swaps(self) -> list[tuple[str, int, str]]:
        return library.parse_swaps(self.member.variant)

    @property
    def sites(self) -> str:
        """The swapped sites as "K300R (residue 287)": the alignment column names the
        member, the residue number is where it sits in the protein being ordered."""
        parts = []
        for ml, column, alt in self.swaps:
            residue = library.residue_number(self.template, column, self.member.width)
            where = f"residue {residue}" if residue else "a gap in this ancestor"
            parts.append(f"{ml}{column}{alt} ({where})")
        return "; ".join(parts)

    @property
    def description(self) -> str:
        variant = library.variant_label(self.member.variant, self.member.kind)
        where = f"{self.source.node} · {self.source.label} · {self.source.version_label}"
        return f"{variant} of {where}"


def wrapped(sequence: str, width: int = WRAP) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width)) or ""


def fasta(rows: Sequence[Row]) -> str:
    """The proteins to order, ungapped: what goes to a synthesis vendor."""
    out = []
    for row in rows:
        out.append(f">{row.name} {row.description}\n{wrapped(row.member.unaligned)}")
    return "\n".join(out) + "\n" if out else ""


MANIFEST_COLUMNS = (
    "name", "variant", "kind", "selected", "length", "support", "weight", "mean_pp",
    "changes_from_ml", "sites", "combination", "node", "reconstruction",
)  # fmt: skip


def manifest(rows: Sequence[Row]) -> str:
    """One line per member, with every number the page showed and the provenance behind
    it. Tab-separated, as every other table the pipeline writes."""
    lines = ["\t".join(MANIFEST_COLUMNS)]
    for row in rows:
        member = row.member
        comparable = len(row.template) == len(member.sequence)
        changes = len(member.differences(row.template)) if comparable else ""
        lines.append(
            "\t".join(
                str(value)
                for value in (
                    row.name,
                    member.variant,
                    member.kind,
                    "yes" if row.selected else "no",
                    member.length,
                    member.support,
                    f"{member.weight:.4f}",
                    "" if row.mean_pp is None else f"{row.mean_pp:.4f}",
                    changes,
                    row.sites,
                    row.source.combination,
                    row.source.node,
                    row.source.version_label,
                )
            )
        )
    return "\n".join(lines) + "\n"


def names(rows: Sequence[Row]) -> str:
    """Name to description: the table a vendor's order form and a lab notebook share, so
    a short name on a tube can be read back to the ancestor it came from."""
    lines = ["name\tdescription"]
    lines += [f"{row.name}\t{row.description}" for row in rows]
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class Methods:
    """Everything the paragraph has to state. The view fills it in; the wording lives in
    `paragraph` so it can be read and checked in one place."""

    library_name: str
    code: str
    job_name: str
    design: library.Design
    sources: Sequence[Source]
    members: int
    selected: int
    ambiguous: int
    coverage: str
    seed: int
    version: str


CITATIONS = (
    "Hanson-Smith V, Johnson A (2016) PhyloBot: A Web Portal for Automated Phylogenetics, "
    "Ancestral Sequence Reconstruction, and Exploration of Mutational Trajectories. "
    "PLoS Comput Biol 12(7): e1004976. https://doi.org/10.1371/journal.pcbi.1004976",
    "Eick GN, Bridgham JT, Anderson DP, Harms MJ, Thornton JW (2017) Robustness of "
    "Reconstructed Ancestral Protein Functions to Statistical Uncertainty. "
    "Mol Biol Evol 34(2): 247-261. https://doi.org/10.1093/molbev/msw223",
)


def paragraph(m: Methods) -> str:
    """The methods paragraph, to paste into a manuscript.

    It states the caveats as well as the settings (the design note is firm about this):
    draws come from IQ-TREE's *marginal* per-site posteriors, with no joint distribution
    available, and single-site alternates say nothing about epistasis. Those two sentences
    are what make the library defensible rather than merely impressive.
    """
    axes = m.design.axes
    where = (
        "; ".join(
            f"{s.node} from {s.label} ({s.short_label or s.version_label})" for s in m.sources
        )
        or "no ancestors"
    )
    variants = ["the maximum-likelihood sequence of each ancestor"]
    if axes.alt_all:
        variants.append(
            f"an AltAll control carrying the second-most-likely state at every site whose "
            f"alternative reached a posterior probability of {axes.threshold:.2f}"
        )
    if axes.samples:
        variants.append(
            f"{axes.samples} sequence(s) drawn per ancestor from the per-site posteriors"
        )
    if axes.max_substitutions:
        variants.append(
            f"every combination of up to {axes.max_substitutions} single-site alternate(s) "
            f"at those sites"
        )
    strategy = {
        "supported": "those that the greatest number of reconstructions agreed on",
        "spanning": "a maximally dissimilar subset, to span the disagreement",
        "per_node": "equal shares per ancestor",
    }[m.design.strategy]
    return (
        "\n\n".join(
            [
                f"Ancestral library {m.library_name} ({m.code}), from PhyloBot job "
                f"“{m.job_name}”. Ancestors were reconstructed at: {where}.",
                f"Candidate sequences comprised {', '.join(variants)}. Sites were counted as "
                f"uncertain at a posterior probability of {axes.threshold:.2f}; the collected "
                f"ancestors had {m.ambiguous} such site(s). Deduplicating on the unaligned "
                f"protein gave {m.members} distinct sequence(s), of which {m.selected} were "
                f"selected for a budget of {m.design.budget} by taking {strategy}. "
                f"{m.coverage}",
                "Caveats. Sequences drawn from the posteriors use IQ-TREE's marginal per-site "
                "posterior probabilities; no joint distribution over sites is available, so "
                "such draws are not samples from the joint posterior and can combine states "
                "that are individually plausible but jointly unlikely. Single-site alternates "
                "vary one site at a time and therefore do not address epistasis between "
                "uncertain sites.",
                f"Reproducibility. Random seed {m.seed}; PhyloBot {m.version}. The manifest "
                f"accompanying this library names every sequence, the ancestor and "
                f"reconstruction it came from, and the sites at which it differs.",
                "References. " + " ".join(CITATIONS),
            ]
        )
        + "\n"
    )


def documents(rows: Sequence[Row], methods: Methods, spec: str) -> dict[str, str]:
    """Every file one export writes, by name.

    The FASTA and the names table hold the *selected* set — what goes to a vendor — while
    the manifest lists every member the design produced, with a column saying which were
    taken: the deliverable should not quietly include sequences the budget cut, and the
    record should not quietly lose the ones it did.
    """
    taken = [row for row in rows if row.selected]
    return {
        FASTA: fasta(taken),
        MANIFEST: manifest(rows),
        NAMES: names(taken),
        METHODS: paragraph(methods),
        SPEC: spec,
    }


def safe_name(name: str) -> str:
    """A download's filename from the library's name: letters, digits, dash, underscore."""
    kept = "".join(c if c.isalnum() or c in "-_" else "-" for c in name).strip("-")
    return kept[:60] or "library"


def listing(exports: Sequence[Any]) -> list[dict[str, Any]]:
    """The exports a page lists, newest first, each with its files."""
    return [
        {
            "export": export,
            "files": [(name, LABELS[name]) for name in FILES],
        }
        for export in exports
    ]
