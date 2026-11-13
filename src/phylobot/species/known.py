"""Species PhyloBot knows offline, from UniProt's controlled vocabulary of species.

UniProt's speclist.txt (CC BY 4.0, UniProt Consortium; bundled as
bundled/uniprot_species.tsv.gz by scripts/fetch-uniprot-species.py, release on its first
line) names about 28,000 organisms with their NCBI taxids, common names and synonyms.
From it, without asking anyone:

- `code`: UniProt species codes (IME2_YEAST, Q6C4M5_YARLI); virtual codes (9SACH,
  Saccharomycetaceae) name a group above species;
- `species`: a binomial ("Homo sapiens", or a synonym such as "Candida glabrata") to its
  official name and taxid;
- `common`: a common name ("human", "baker's yeast"), when it is unambiguous;
- `by_taxid`: a taxid's name;
- `abbreviation`: "S.cerevisiae", "H_sapiens" (genus initial + epithet) and four-letter
  lab codes ("Scer", "Hsap"), only when one species fits: the curated model organisms
  first, since "Scer" fits ten species in UniProt's list but means only one in a lab.

Everything is built once per process from one 700 KB file (`table`).
"""

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import cache
from importlib.resources import files


@dataclass(frozen=True)
class Known:
    name: str  # official scientific name, strain notes dropped
    taxid: int
    group: bool = False  # a UniProt virtual code: a taxon above species


# Organisms people abbreviate (the usual model and reference species): a four-letter
# code or an initial means one of these before any other species in UniProt's list.
MODEL_ORGANISMS = (
    "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio", "Xenopus laevis",
    "Xenopus tropicalis", "Gallus gallus", "Bos taurus", "Sus scrofa", "Canis lupus",
    "Macaca mulatta", "Pan troglodytes", "Oryzias latipes", "Takifugu rubripes",
    "Ciona intestinalis", "Strongylocentrotus purpuratus", "Nematostella vectensis",
    "Drosophila melanogaster", "Anopheles gambiae", "Apis mellifera",
    "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Schizosaccharomyces pombe",
    "Candida albicans", "Kluyveromyces lactis", "Yarrowia lipolytica",
    "Neurospora crassa", "Aspergillus nidulans", "Aspergillus fumigatus",
    "Cryptococcus neoformans", "Dictyostelium discoideum", "Arabidopsis thaliana",
    "Oryza sativa", "Zea mays", "Chlamydomonas reinhardtii", "Physcomitrium patens",
    "Plasmodium falciparum", "Trypanosoma brucei", "Toxoplasma gondii",
    "Tetrahymena thermophila", "Escherichia coli", "Bacillus subtilis",
    "Mycobacterium tuberculosis", "Pseudomonas aeruginosa", "Staphylococcus aureus",
    "Methanocaldococcus jannaschii", "Haloferax volcanii",
)  # fmt: skip

STRAIN = re.compile(r"\s*\(.*?\)")


def binomial(name: str) -> str | None:
    """'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)' -> 'Saccharomyces
    cerevisiae'; None for names that are not Genus species."""
    words = STRAIN.sub(" ", name).split()
    if (
        len(words) >= 2
        and re.fullmatch(r"[A-Z][a-z]+", words[0])
        and re.fullmatch(r"[a-z][a-z-]+", words[1])
    ):
        return f"{words[0]} {words[1]}"
    return None


@dataclass
class Tables:
    release: str
    codes: dict[str, Known]
    species: dict[str, Known]  # lower-case binomial or synonym -> species
    common: dict[str, list[Known]]  # lower-case common name -> every species with it
    taxids: dict[int, str]
    genera: set[str]
    initials: dict[tuple[str, str], list[Known]]  # (genus initial, epithet) -> species
    fours: dict[str, list[Known]]  # "scer" -> species


@cache
def table() -> Tables:
    data = files("phylobot.species").joinpath("bundled", "uniprot_species.tsv.gz").read_bytes()
    lines = gzip.decompress(data).decode().splitlines()
    codes: dict[str, Known] = {}
    species: dict[str, Known] = {}
    common: dict[str, dict[str, Known]] = defaultdict(dict)  # name -> {species: Known}
    taxids: dict[int, str] = {}
    for line in lines[1:]:
        code, kind, taxid_text, name, commons, synonyms = (line.split("\t") + ["", ""])[:6]
        taxid, group = int(taxid_text), kind == "group"
        short = binomial(name)
        known = Known(STRAIN.sub("", name).strip() if group or not short else short, taxid, group)
        codes[code] = known
        taxids.setdefault(taxid, STRAIN.sub("", name).strip())
        if group or short is None:
            continue
        # one entry per species: the species-level one (its name is exactly the
        # binomial) wins over its strains
        if short.lower() not in species or name == short:
            species[short.lower()] = Known(short, taxid)
        for synonym in filter(None, synonyms.split("|")):
            if (other := binomial(synonym)) and other.lower() not in species:
                species[other.lower()] = Known(short, taxid)
        for name_ in filter(None, commons.split("|")):
            holders = common[name_.strip().lower()]
            if name == short or short not in holders:  # the species' own taxid over strains'
                holders[short] = Known(short, taxid)
    official = {k.name for k in species.values()}
    genera = {name.split()[0] for name in official}
    initials: dict[tuple[str, str], list[Known]] = defaultdict(list)
    fours: dict[str, list[Known]] = defaultdict(list)
    for name in sorted(official):
        known = species[name.lower()]
        genus, epithet = name.split()
        initials[(genus[0], epithet)].append(known)
        fours[(genus[0] + epithet[:3]).lower()].append(known)
    by_common = {name: list(holders.values()) for name, holders in common.items()}
    return Tables(
        lines[0].lstrip("# "), codes, species, by_common, taxids, genera, dict(initials),
        dict(fours),
    )  # fmt: skip


def release() -> str:
    return table().release


def code(value: str) -> Known | None:
    return table().codes.get(value)


def species(name: str) -> Known | None:
    return table().species.get(" ".join(name.split()).lower())


def common(name: str) -> Known | None:
    """A common name's species, when one species (or one model organism) has it."""
    return one(table().common.get(" ".join(name.split()).lower(), []))


def by_taxid(taxid: int) -> str | None:
    return table().taxids.get(taxid)


def is_genus(word: str) -> bool:
    return word in table().genera


def one(candidates: list[Known]) -> Known | None:
    """The single model organism among candidates, else the only candidate."""
    models = [k for k in candidates if k.name in MODEL_ORGANISMS]
    if len(models) == 1:
        return models[0]
    return candidates[0] if len(candidates) == 1 else None


def abbreviation(text: str) -> Known | None:
    """'S.cerevisiae', 'S_cerevisiae', 'S. cerevisiae' or a four-letter lab code ('Scer',
    'HSAP'), when exactly one species fits (model organisms first)."""
    if match := re.fullmatch(r"([A-Z])[._ ]{1,2}([a-z][a-z-]{2,})", text.strip()):
        return one(table().initials.get((match[1], match[2]), []))
    if re.fullmatch(r"[A-Za-z]{4}", text):  # lab codes: model organisms only
        models = [k for k in table().fours.get(text.lower(), []) if k.name in MODEL_ORGANISMS]
        return models[0] if len(models) == 1 else None
    return None
