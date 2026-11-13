"""The species a FASTA header names, if it names one, read offline.

Detection is a suggestion, never the truth: the user reviews
and corrects every row before a job uses it. Recognised, most reliable first:

1. Annotations: UniProtKB `OS=Saccharomyces cerevisiae OX=4932`; UniRef `Tax=...
   TaxID=...` (a cluster's common taxon, which may be above species); GenBank
   `[organism=...]`; NCBI's `... [Saccharomyces cerevisiae]` at the end.
2. Identifiers that imply their species (`known` and `IDS`): UniProt entry names
   (`sp|P32581|IME2_YEAST`, `Q6C4M5_YARLI`, by UniProt's species code); taxid-prefixed
   ids from STRING, eggNOG and OrthoDB (`9606.ENSP00000269305`, `9606_0:001c7b`);
   Ensembl ids (`ENSMUSP...`, by their species prefix); model-organism ids (`AT1G01010`,
   `YAL001C`, `FBpp0070000`, `WBGene00000001`, `PF3D7_0100100`, ...).
3. A binomial UniProt knows, anywhere in the name or description, in any case
   (`IME2_Saccharomyces_cerevisiae`, `homo_sapiens`); else the naming conventions at
   the start of a name: `Genus.species.GENE` (the Ime2 example),
   `Genus_species_Accession` (PhyloBot's homolog export), `Genus_species_subspecies`;
   else a bracketed binomial mid-description (`[Homo sapiens] isoform 2`).
4. Guesses, labelled as such: an abbreviation leading the name (`S.cerevisiae_IME2`,
   `H_sapiens`, or a model organism's four-letter code: `Scer_IME2`, `Hsap`), taken
   only when exactly one species fits.

A lowercase third word in a name may be a subspecies or a lowercase gene name; the
resolver falls back to the first two words when three find nothing.
"""

import re
from dataclasses import dataclass

from phylobot.species import known

UNIPROT = re.compile(r"(?:^|\s)OS=(?P<taxon>.+?)\s+OX=(?P<taxid>\d{1,10})(?:\s|$)")
UNIREF = re.compile(r"(?:^|\s)Tax=(?P<taxon>.+?)\s+TaxID=(?P<taxid>\d{1,10})(?:\s|$)")
NCBI = re.compile(r"\[(?P<taxon>[A-Z][^\[\]]{2,120})\]\s*$")
ORGANISM = re.compile(r"\[organism=(?P<taxon>[^\[\]]{3,120})\]")
BRACKET = re.compile(r"\[(?P<taxon>[A-Z][a-z]+ [a-z][a-z.-]+(?: [^\[\]]{1,80})?)\]")
ENTRY_NAME = re.compile(
    r"^(?:(?:sp|tr)\|[A-Z0-9]{6,10}\|)?[A-Z0-9]{1,10}_(?P<code>[A-Z0-9]{1,5})(?=$|[| ])"
)
TAXID_PREFIX = re.compile(r"^(?P<taxid>\d{1,7})(?:\.(?=[A-Za-z0-9_]{3})|_\d+:)")
NAME = re.compile(
    r"^(?P<genus>[A-Z][a-z]{1,40})[._ ](?P<species>[a-z][a-z-]{1,40})"
    r"(?:[._ ](?P<third>[a-z][a-z-]{1,40}))?(?=[._ |]|$)"
)
TOKENS = re.compile(r"[^A-Za-z]+")
LEADING = re.compile(r"^([A-Z][._]{1,2}[a-z][a-z-]{2,}|[A-Za-z]{4})(?=[._| -]|$)")
MAX_TAXON = 120

# Ensembl stable ids: ENS + a species prefix + P/G/T/E + 11 digits (no prefix: human).
ENSEMBL = re.compile(r"(?<![A-Z])ENS(?P<prefix>[A-Z]{0,5}?)[PGTE]\d{11}(?:\.\d+)?")
ENSEMBL_SPECIES = {
    "": "Homo sapiens", "MUS": "Mus musculus", "RNO": "Rattus norvegicus",
    "DAR": "Danio rerio", "GAL": "Gallus gallus", "XET": "Xenopus tropicalis",
    "BTA": "Bos taurus", "SSC": "Sus scrofa", "CAF": "Canis lupus familiaris",
    "FCA": "Felis catus", "ECA": "Equus caballus", "OAR": "Ovis aries",
    "OCU": "Oryctolagus cuniculus", "MMU": "Macaca mulatta", "PTR": "Pan troglodytes",
    "GGO": "Gorilla gorilla", "PPY": "Pongo abelii", "CJA": "Callithrix jacchus",
    "OAN": "Ornithorhynchus anatinus", "MOD": "Monodelphis domestica",
    "TGU": "Taeniopygia guttata", "MGA": "Meleagris gallopavo", "ACA": "Anolis carolinensis",
    "TRU": "Takifugu rubripes", "ORL": "Oryzias latipes", "GAC": "Gasterosteus aculeatus",
    "CIN": "Ciona intestinalis", "LAF": "Loxodonta africana", "CPO": "Cavia porcellus",
}  # fmt: skip
# Identifiers only one organism uses: (pattern, species, taxid).
IDS = [
    (re.compile(r"\bAT[1-5CM]G\d{5}(?:\.\d+)?\b"), "Arabidopsis thaliana", 3702),
    (re.compile(r"^Y[A-P][LR]\d{3}[WC](?:-[A-Z])?\b"), "Saccharomyces cerevisiae", 4932),
    (re.compile(r"^SP[ABC]C\d{1,4}[A-Z]?\d*\.\d{2}c?\b"), "Schizosaccharomyces pombe", 4896),
    (re.compile(r"\bFB(?:pp|gn|tr)\d{7}\b"), "Drosophila melanogaster", 7227),
    (re.compile(r"\bWBGene\d{8}\b"), "Caenorhabditis elegans", 6239),
    (re.compile(r"\bDDB_G\d{7}\b"), "Dictyostelium discoideum", 44689),
    (re.compile(r"\bZDB-GENE-\d{6}-\d+\b"), "Danio rerio", 7955),
    (re.compile(r"\bPF3D7_\d{7}\b"), "Plasmodium falciparum", 5833),
    (re.compile(r"\bLOC_Os\d{2}g\d{5}\b|\bOs\d{2}g\d{7}\b"), "Oryza sativa", 4530),
    (re.compile(r"\bGRMZM\dG\d{6}\b|\bZm\d{5}d\d{6}\b"), "Zea mays", 4577),
    (re.compile(r"\bHGNC:\d+\b"), "Homo sapiens", 9606),
    (re.compile(r"\bMGI:\d+\b"), "Mus musculus", 10090),
    (re.compile(r"\bRGD:\d+\b"), "Rattus norvegicus", 10116),
]


@dataclass(frozen=True)
class Detected:
    """The name to look up, the NCBI taxid if the header gave one, and how it was
    found ("guess" sources are shown to the user as guesses)."""

    name: str
    taxid: int | None
    source: str  # see species_web.SOURCES
    note: str = ""


def clean_taxon(text: str) -> str:
    """One line of printable text, single spaces, strain notes in brackets dropped."""
    text = re.sub(r"\(.*?\)", " ", text)
    text = re.sub(r"[^\x20-\x7e]", " ", text)
    return " ".join(text.split())[:MAX_TAXON]


def annotated(description: str) -> Detected | None:
    for pattern, source in ((UNIPROT, "uniprot"), (UNIREF, "uniref")):
        if (match := pattern.search(description)) and (taxon := clean_taxon(match["taxon"])):
            return Detected(taxon, int(match["taxid"]), source)
    for pattern, source in ((ORGANISM, "genbank"), (NCBI, "ncbi")):
        if (match := pattern.search(description)) and (taxon := clean_taxon(match["taxon"])):
            return Detected(taxon, None, source)
    return None


def identifier(name: str) -> Detected | None:
    if (match := ENTRY_NAME.match(name)) and (code := known.code(match["code"])):
        # a strain's code (YEAST: S288c) names its species; take the species' own taxid
        own = None if code.group else known.species(code.name)
        return Detected(code.name, own.taxid if own else code.taxid, "uniprot code")
    if match := TAXID_PREFIX.match(name):
        taxid = int(match["taxid"])
        if taxid > 0:
            return Detected(known.by_taxid(taxid) or "", taxid, "taxid prefix")
    if (match := ENSEMBL.search(name)) and match["prefix"] in ENSEMBL_SPECIES:
        species = ENSEMBL_SPECIES[match["prefix"]]
        found = known.species(species)
        return Detected(species, found.taxid if found else None, "ensembl id")
    for pattern, species, taxid in IDS:
        if pattern.search(name):
            return Detected(species, taxid, "model id")
    return None


def binomial(*texts: str) -> Detected | None:
    """The first pair of adjacent words, anywhere, that UniProt knows as a species,
    written as Genus species but otherwise as found (Open Tree resolves synonyms)."""
    for text in texts:
        words = [w for w in TOKENS.split(text) if w]
        for first, second in zip(words, words[1:], strict=False):
            if known.species(f"{first} {second}") is not None:
                return Detected(f"{first.capitalize()} {second.lower()}", None, "name")
    return None


def guess(name: str) -> Detected | None:
    if (match := LEADING.match(name)) and (found := known.abbreviation(match[1])):
        return Detected(found.name, found.taxid, "guess", f"guessed from {match[1]!r}: check it")
    return None


def convention(name: str) -> Detected | None:
    """Genus.species.GENE, Genus_species_Accession, Genus_species_subspecies."""
    if match := NAME.match(name):
        words = [match["genus"], match["species"]]
        if match["third"] and match["species"] not in ("sp", "cf", "aff"):
            words.append(match["third"])
        return Detected(" ".join(words), None, "name")
    return None


def detect(name: str, description: str = "") -> Detected | None:
    if found := annotated(description) or identifier(name):
        return found
    start = convention(name)
    if start is not None and known.is_genus(start.name.split()[0]):
        return start
    if found := binomial(name, description):
        return found
    if start is not None:  # a genus UniProt's list lacks: maybe real, the resolver decides
        return start
    if (match := BRACKET.search(description)) and (taxon := clean_taxon(match["taxon"])):
        return Detected(taxon, None, "ncbi")
    return guess(name)
