"""Genetic codes, translation and codon alignments, for coding-DNA input.

A coding sequence is translated with one of the genetic codes below; the protein feeds
the amino-acid aligners, and each protein alignment is back-translated to a codon
alignment by putting every sequence's own codons back in place of its residues. Every
step checks that codons and residues still correspond, so a tool that silently changed
a sequence (PRANK replaces internal stop codons with NNN) fails the job.

IUPAC ambiguity codes are accepted (`IUPAC`). A codon with one translates to the amino
acid every reading gives (CTN and CTR are Leu), else to X (`ambiguous_codons` lists them
for the user). Downstream tools differ (measured 2026-10-18 on IQ-TREE 3.1.3, PAML
4.10.10, PRANK 251117): codeml sums over a codon's readings, IQ-TREE's codon models treat
any codon with an ambiguous base as unknown, and PRANK -codon replaces it with NNN
(`restore_codons` puts the user's codon back).

The codes are those IQ-TREE (`-st CODONn`, NCBI numbering) or PAML (`icode`) offers;
neither accepts a custom code. Tables are NCBI's (gc.prt, 2026-10-01), codons in TCAG
order. PAML's tables were compared with NCBI's: its icode 0-10 are NCBI 1 (and 11, the
same amino acids), 2, 3, 4, 5, 6, 9, 10, 12, 13 and 15; its icode 11 is not an NCBI code
and is not offered.
"""

import itertools
from dataclasses import dataclass

BASES = "TCAG"
CODONS = ["".join(c) for c in itertools.product(BASES, repeat=3)]
# IUPAC nucleotide codes and the bases each stands for.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}  # fmt: skip
DNA = set(IUPAC)


class CodonError(ValueError):
    pass


@dataclass(frozen=True)
class GeneticCode:
    ncbi: int
    name: str
    amino: str  # one letter per codon in CODONS order; '*' = stop
    iqtree: bool  # IQ-TREE has -st CODON<ncbi>
    paml_icode: int | None  # codeml's icode, if PAML has this code

    @property
    def table(self) -> dict[str, str]:
        return dict(zip(CODONS, self.amino, strict=True))

    @property
    def tools(self) -> list[str]:
        return ["IQ-TREE"] * self.iqtree + ["PAML"] * (self.paml_icode is not None)


GENETIC_CODES: tuple[GeneticCode, ...] = (
    GeneticCode(1, "Standard",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=0),
    GeneticCode(2, "Vertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=1),
    GeneticCode(3, "Yeast Mitochondrial",
                "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=2),
    GeneticCode(4, "Mold, Protozoan and Coelenterate Mitochondrial, Mycoplasma, Spiroplasma",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=3),
    GeneticCode(5, "Invertebrate Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=4),
    GeneticCode(6, "Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear",
                "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=5),
    GeneticCode(9, "Echinoderm Mitochondrial, Flatworm Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=6),
    GeneticCode(10, "Euplotid Nuclear",
                "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=7),
    GeneticCode(11, "Bacterial, Archaeal and Plant Plastid",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=0),
    GeneticCode(12, "Alternative Yeast Nuclear",
                "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=8),
    GeneticCode(13, "Ascidian Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=9),
    GeneticCode(14, "Alternative Flatworm Mitochondrial",
                "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(15, "Blepharisma Macronuclear",
                "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=False, paml_icode=10),
    GeneticCode(16, "Chlorophycean Mitochondrial",
                "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(21, "Trematode Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(22, "Scenedesmus obliquus Mitochondrial",
                "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(23, "Thraustochytrium Mitochondrial",
                "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(24, "Rhabdopleuridae Mitochondrial",
                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
    GeneticCode(25, "Candidate Division SR1 and Gracilibacteria",
                "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                iqtree=True, paml_icode=None),
)  # fmt: skip
CODES = {code.ncbi: code for code in GENETIC_CODES}
STANDARD = CODES[1]


def genetic_code(ncbi: int | str) -> GeneticCode:
    try:
        return CODES[int(ncbi)]
    except KeyError, ValueError:
        raise CodonError(f"unknown genetic code {ncbi!r} (offered: {sorted(CODES)})") from None


def is_ambiguous(codon: str) -> bool:
    """Whether a codon has an ambiguity code (any base but A, C, G or T)."""
    return any(base not in "ACGT" for base in codon)


def readings(codon: str) -> list[str]:
    """Every unambiguous codon an (uppercase) codon may be, e.g. CTR -> CTA, CTG."""
    return ["".join(bases) for bases in itertools.product(*(IUPAC[b] for b in codon))]


def possible_aminos(codon: str, code: GeneticCode) -> str:
    """The amino acids (and '*' for a stop) a codon's readings encode, sorted."""
    return "".join(sorted({code.table[c] for c in readings(codon)}))


def compatible(a: str, b: str) -> bool:
    """Whether two codons may be the same codon (they share a reading)."""
    return all(set(IUPAC[x]) & set(IUPAC[y]) for x, y in zip(a, b, strict=True))


def translate_codon(codon: str, code: GeneticCode) -> str:
    """One codon's amino acid, '*' for a stop. With an ambiguity code: the amino acid
    (or stop) every reading gives, else X."""
    if not is_ambiguous(codon):
        return code.table[codon]
    aminos = possible_aminos(codon, code)
    return aminos if len(aminos) == 1 else "X"


def translate(cds: str, code: GeneticCode) -> str:
    """An ungapped coding sequence's protein (its length must be whole codons)."""
    if len(cds) % 3:
        raise CodonError(f"length {len(cds)} is not a multiple of 3")
    return "".join(translate_codon(cds[i : i + 3], code) for i in range(0, len(cds), 3))


def check_cds(cds: str, code: GeneticCode) -> tuple[str, str, list[str]]:
    """(coding sequence without a final stop codon, its protein, notes) for an uppercase
    ungapped sequence of IUPAC codes; raises CodonError for a length that is not
    whole codons (a frameshift or truncation) or an internal stop codon."""
    if len(cds) % 3:
        raise CodonError(
            f"length {len(cds)} is not a multiple of 3 (a frameshift or a truncated sequence?)"
        )
    notes = []
    protein = translate(cds, code)
    if protein.endswith("*"):
        cds, protein = cds[:-3], protein[:-1]
        notes.append("final stop codon removed")
    if (stop := protein.find("*")) >= 0:
        raise CodonError(
            f"internal stop codon {cds[3 * stop : 3 * stop + 3]} at codon {stop + 1} "
            f"(genetic code {code.ncbi}, {code.name})"
        )
    if not protein:
        raise CodonError("no codons besides a stop codon")
    return cds, protein, notes


@dataclass(frozen=True)
class AmbiguousCodon:
    """A codon of a coding sequence with an ambiguity code."""

    number: int  # 1-based, in the ungapped coding sequence
    codon: str
    aminos: str  # what its readings encode ('*' = stop), sorted
    translation: str  # what it is translated as: the one amino acid, or X

    @property
    def unknown(self) -> bool:
        """Every base is N: nothing is known about the codon."""
        return self.codon == "NNN"


def ambiguous_codons(cds: str, code: GeneticCode) -> list[AmbiguousCodon]:
    """The codons of an ungapped coding sequence that have an ambiguity code."""
    found = []
    for i in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[i : i + 3]
        if is_ambiguous(codon):
            aminos = possible_aminos(codon, code)
            found.append(AmbiguousCodon(i // 3 + 1, codon, aminos, translate_codon(codon, code)))
    return found


def restore_codons(row: str, cds: str) -> str:
    """An aligned codon row with the sequence's own codons put back where the aligner
    wrote NNN in place of a codon with an ambiguity code (PRANK -codon does); refuses
    any other change to the sequence."""
    if len(row) % 3:
        raise CodonError(f"aligned length {len(row)} is not a multiple of 3")
    out, used = [], 0
    for i in range(0, len(row), 3):
        codon = row[i : i + 3]
        if codon == "---":
            out.append(codon)
            continue
        own = cds[3 * used : 3 * used + 3]
        if codon != own and not (codon == "NNN" and is_ambiguous(own)):
            raise CodonError("the aligned sequence differs from the input")
        out.append(own)
        used += 1
    if 3 * used != len(cds):
        raise CodonError("the aligned sequence differs from the input")
    return "".join(out)


def translate_row(row: str, code: GeneticCode) -> str:
    """An aligned codon row's protein row ('---' -> '-'); gaps must be whole codons."""
    if len(row) % 3:
        raise CodonError(f"aligned length {len(row)} is not a multiple of 3")
    out = []
    for i in range(0, len(row), 3):
        codon = row[i : i + 3]
        if codon == "---":
            out.append("-")
        elif "-" in codon:
            raise CodonError(f"gap inside codon {i // 3 + 1} ({codon}); gaps must be whole codons")
        else:
            out.append(translate_codon(codon, code))
    return "".join(out)


def backtranslate(protein_row: str, cds: str, code: GeneticCode) -> str:
    """The codon row for an aligned protein row: each residue replaced by the sequence's
    own next codon, each gap by '---'. Refuses any residue its codon does not encode,
    and leftover or missing codons."""
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    out, used = [], 0
    for column, residue in enumerate(protein_row, 1):
        if residue == "-":
            out.append("---")
            continue
        if used == len(codons):
            raise CodonError(f"more residues than codons (column {column})")
        codon = codons[used]
        if translate_codon(codon, code) != residue.upper():
            raise CodonError(
                f"column {column}: residue {residue} but codon {used + 1} is {codon} "
                f"({translate_codon(codon, code)})"
            )
        out.append(codon)
        used += 1
    if used != len(codons):
        raise CodonError(f"{len(codons) - used} codons left over after the last residue")
    return "".join(out)


# Codon posteriors kept per ancestral site (ancestors_codon_posteriors.tsv): codons at or
# above this, as the viewer shows tiles (ancestors.MIN_SHOWN).
CODON_KEPT = 0.01


def codon_lines(node: str, site: int, probs: dict[str, float]) -> list[str]:
    """ancestors_codon_posteriors.tsv rows for one site: node, site, codon, p, most
    probable first, codons below CODON_KEPT left out (the best codon always kept)."""
    ranked = sorted(probs.items(), key=lambda cp: (-cp[1], cp[0]))
    kept = [cp for i, cp in enumerate(ranked) if i == 0 or cp[1] >= CODON_KEPT]
    return [f"{node}\t{site}\t{codon}\t{p:.5f}\n" for codon, p in kept]


# IQ-TREE's order of protein states (its .state columns p_A, p_R, ...).
IQTREE_AMINO_ORDER = "ARNDCQEGHILKMFPSTWYV"


def amino_posteriors(
    codon_probs: dict[str, float], code: GeneticCode
) -> tuple[dict[str, float], str]:
    """(posterior of each amino acid, IQ-TREE's order; the best codon for the best amino
    acid) at one ancestral site, from its codon posteriors. The amino acid comes first
    (synonymous codons summed): the most probable codon can encode a less probable
    amino acid (20 of 3,105 node-sites on lysin), and the coding ancestor must translate
    to the protein one."""
    amino = dict.fromkeys(IQTREE_AMINO_ORDER, 0.0)
    for codon, p in codon_probs.items():
        amino[code.table[codon]] += p
    best = max(amino, key=lambda a: amino[a])
    codon = max((c for c in codon_probs if code.table[c] == best), key=lambda c: codon_probs[c])
    return amino, codon
