"""Example sequence sets offered on the New Job page ("Load examples").

Each entry ships one or two FASTA files under ``examples/``: the full homolog set
and, for the larger families, a ``.subset`` for a quick run. The list is the nine
paper-1 families, plus the bundled Ime2 smoke-test set and the lysin coding-DNA set.

Citations are the dataset's source paper where one exists; curated sets we assembled
for PhyloBot 2 (no prior resurrection study) carry no DOI. ``build`` in the sibling
``examples/README.md`` records where each FASTA came from — a change here must keep the
shipped files in step.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Example:
    key: str
    label: str  # short name shown in the drawer
    blurb: str  # one sentence for the "About these datasets" popover
    count: int  # sequences in the full set
    seqtype: str = "protein"  # "protein" or "codon" (coding DNA)
    citation: str = ""  # human-readable reference, "" when none exists
    doi: str = ""  # DOI without the resolver prefix
    url: str = ""  # a link to use when there is no DOI
    # An optional second reference (e.g. the dataset's paper plus the family's origin).
    citation2: str = ""
    doi2: str = ""
    url2: str = ""
    subset_count: int = 0  # sequences in the subset (0 = no subset shipped)
    outgroup: tuple[str, ...] = field(default_factory=tuple)

    @property
    def full_file(self) -> str:
        return f"examples/{self.key}.fasta"

    @property
    def subset_file(self) -> str:
        return f"examples/{self.key}.subset.fasta"

    @property
    def has_subset(self) -> bool:
        return self.subset_count > 0

    @property
    def coding(self) -> bool:
        return self.seqtype == "codon"

    @property
    def link(self) -> str:
        if self.doi:
            return f"https://doi.org/{self.doi}"
        return self.url

    @property
    def references(self) -> list[tuple[str, str]]:
        """Each reference as (citation, link); link is "" when there is none."""
        refs: list[tuple[str, str]] = []
        if self.citation:
            refs.append((self.citation, self.link))
        if self.citation2:
            link2 = f"https://doi.org/{self.doi2}" if self.doi2 else self.url2
            refs.append((self.citation2, link2))
        return refs


# Listed smallest to largest (sequences in the full set): the order New Job shows them in.
EXAMPLES: tuple[Example, ...] = tuple(
    sorted(
        (
            Example(
                key="ime2",
                label="Ime2 small test",
                blurb="Seven fungal Ime2 meiosis kinases — a small example that runs quickly.",
                count=7,
                outgroup=("Yarrowia.lipolytica.IME2",),
            ),
            Example(
                key="gfp",
                label="GFP",
                blurb="Coral green/red fluorescent proteins (coding DNA).",
                count=109,
                seqtype="codon",
                citation="Alieva et al. 2008, PLoS ONE 3:e2680",
                doi="10.1371/journal.pone.0002680",
                subset_count=20,
            ),
            Example(
                key="gba1",
                label="GBA1",
                blurb="Lysosomal glucocerebrosidase (GBA1) orthologues, curated for PhyloBot 2.",
                count=131,
                subset_count=20,
            ),
            Example(
                key="maguk",
                label="MAGUK",
                blurb="MAGUK scaffolds (DLG/PSD-95), the family with the GK enzyme-to-scaffold "
                "functional switch.",
                count=135,
                citation="Anderson et al. 2016, eLife 5:e10147",
                doi="10.7554/eLife.10147",
                subset_count=20,
            ),
            Example(
                key="pink1",
                label="PINK1",
                blurb="PINK1 mitophagy-kinase orthologues sampled across animals, curated for "
                "PhyloBot 2.",
                count=105,
                subset_count=20,
            ),
            Example(
                key="rubisco",
                label="Rubisco",
                blurb="Rubisco large subunits and Rubisco-like proteins across Forms I–IV.",
                count=81,
                citation="Kacar et al. 2017, Geobiology 15:628–640",
                doi="10.1111/gbi.12243",
                subset_count=20,
            ),
            Example(
                key="snca",
                label="Synuclein",
                blurb="Vertebrate α/β/γ-synucleins (46 species); α-synuclein is the methods "
                "control. Curated for PhyloBot 2.",
                count=115,
                subset_count=20,
            ),
            Example(
                key="src_abl",
                label="Src/Abl",
                blurb="Src and Abl tyrosine-kinase domains — a published extant-sequence "
                "reconstruction benchmark.",
                count=76,
                citation="Sennett & Theobald 2024, J. Mol. Evol. 92:181–206",
                doi="10.1007/s00239-024-10162-3",
                subset_count=20,
            ),
            Example(
                key="sr238",
                label="Steroid Receptors",
                blurb=(
                    "Steroid hormone receptors (AR, ER, GR, MR, PR). Contains DNA-binding domain, "
                    "hinge, and ligand-binding domain."
                ),
                count=238,
                citation="Eick et al. 2012, PLoS Genet. 8:e1003072",
                doi="10.1371/journal.pgen.1003072",
                citation2="Bridgham, Carroll & Thornton 2006, Science 312:97–101",
                doi2="10.1126/science.1123348",
                subset_count=20,
            ),
            Example(
                key="vatpase",
                label="V-ATPase",
                blurb="Fungal V-ATPase ring subunits 3, 11 and 16.",
                count=140,
                citation="Finnigan et al. 2012, Nature 481:360–364",
                doi="10.1038/nature10724",
                subset_count=20,
            ),
            Example(
                key="lysin",
                label="Lysin",
                blurb="Abalone sperm lysin from 25 species (coding DNA) — a classic "
                "positive-selection codon set.",
                count=25,
                seqtype="codon",
                citation="Lee, Ota & Vacquier 1995, Mol. Biol. Evol. 12:231–238",
                doi="10.1093/oxfordjournals.molbev.a040200",
                subset_count=12,
            ),
        ),
        key=lambda e: e.count,
    )
)

BY_KEY: dict[str, Example] = {e.key: e for e in EXAMPLES}


def resource(key: str, variant: str = "full") -> str:
    """The packaged FASTA path for one example and variant ("full" or "subset").

    Falls back to the full set when a subset is asked for but none was shipped.
    """
    example = BY_KEY[key]
    if variant == "subset" and example.has_subset:
        return example.subset_file
    return example.full_file
