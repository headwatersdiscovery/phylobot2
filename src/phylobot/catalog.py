"""The substitution models PhyloBot offers: protein matrices x rate/frequency variants.

Users pick models only from this catalog (a grid: matrices down, variants across), so a
model string is always a known-good IQ-TREE model and safe to use as a path segment and
command argument. Every matrix and variant here was checked against IQ-TREE 3.1.3.

Two ways to use ticked models (both may be enabled in one job):
- fixed: each ticked model is its own combination;
- ModelFinder: one combination that chooses among candidates. IQ-TREE only accepts a
  candidate set as matrices x rate types x frequency types (-mset/-mrate/-mfreq), so the
  candidates are the cross product of what the ticked models use, which can be more than
  the ticked cells. With +R, ModelFinder also tries several category counts.

Site-profile mixtures (MIXTURE_VARIANTS) are fixed models only: never ModelFinder candidates,
because IQ-TREE's -mset cross product and its shortcuts do not cover them, and each is a slow
fit. Every step that fits one passes -mwopt (`pipeline/modules/mixtures.nf`): without it IQ-TREE
keeps C20's published weights whatever the data (on real PINK1, 501 BIC worse). +F beside a
mixture adds a class with the alignment's own frequencies, which some real families need
(without it, C20 can fall behind a plain model in BIC).
"""

from dataclasses import dataclass

MODELFINDER = "modelfinder"  # combination key for ModelFinder model selection


@dataclass(frozen=True)
class Matrix:
    name: str
    group: str


@dataclass(frozen=True)
class Variant:
    suffix: str  # appended to the matrix name: "LG" + "+G4"
    rate: str  # ModelFinder -mrate code
    freq: str  # ModelFinder -mfreq code
    label: str


MATRICES: tuple[Matrix, ...] = (
    *(Matrix(name, "General") for name in (
        "LG", "WAG", "JTT", "JTTDCMut", "Dayhoff", "DCMut", "VT", "PMB", "Blosum62",
    )),
    *(Matrix(name, "Q (estimated from large databases)") for name in (
        "Q.pfam", "Q.insect", "Q.yeast", "Q.plant", "Q.mammal", "Q.bird",
    )),
    *(Matrix(name, "Organelle and virus") for name in (
        "mtREV", "mtART", "mtZOA", "mtMAM", "mtInv", "mtMet", "mtVer",
        "cpREV", "rtREV", "HIVb", "HIVw", "FLU", "FLAVI",
    )),
)  # fmt: skip

VARIANTS: tuple[Variant, ...] = (
    Variant("", "E", "FU", "equal rates"),
    Variant("+I", "I", "FU", "+ invariant sites"),
    Variant("+G4", "G", "FU", "+ gamma (4 categories)"),
    Variant("+I+G4", "I+G", "FU", "+ invariant + gamma"),
    Variant("+R4", "R", "FU", "+ FreeRate (4 categories)"),
    Variant("+I+R4", "I+R", "FU", "+ invariant + FreeRate"),
    Variant("+F", "E", "F", "+ empirical frequencies"),
    Variant("+F+G4", "G", "F", "+ empirical frequencies + gamma"),
    Variant("+F+R4", "R", "F", "+ empirical frequencies + FreeRate"),
)

# Fixed models only (see the module docstring). Their rate and freq codes are for
# `workload.heterogeneous`: no ModelFinder set ever holds them.
MIXTURE_VARIANTS: tuple[Variant, ...] = (
    Variant("+C20+G4", "G", "FU", "+ C20 site profiles + gamma"),
    Variant("+C20+F+G4", "G", "F", "+ C20 site profiles + a class of the alignment's own "
                                   "frequencies + gamma"),
)  # fmt: skip

# Matrices RAxML-NG 2.0.3 lacks. Every other catalog model string (and ModelFinder's other
# +R category counts) means the same in RAxML-NG: checked by fitting a fixed tree in both
# programs (log-likelihoods within 0.1). A RAxML-NG job can't use these, as a fixed model
# or ModelFinder candidate: VALIDATE refuses them if ticked, and leaves them out of
# ModelFinder's default candidates.
RAXML_MISSING = frozenset({"JTTDCMut"})

_MATRIX_NAMES = {m.name for m in MATRICES}
# The grid's rate and frequency types, in order: also IQ-TREE's defaults for ModelFinder
# (checked against a default run's modelfinder.tsv).
ALL_RATES = list(dict.fromkeys(v.rate for v in VARIANTS))
ALL_FREQS = list(dict.fromkeys(v.freq for v in VARIANTS))
_BY_SUFFIX = {v.suffix: v for v in (*VARIANTS, *MIXTURE_VARIANTS)}


def split(model: str) -> tuple[str, Variant]:
    """("LG", Variant("+G4", ...)) for "LG+G4"; raises ValueError if not in the catalog."""
    matrix, plus, rest = model.partition("+")
    variant = _BY_SUFFIX.get(plus + rest)
    if matrix not in _MATRIX_NAMES or variant is None:
        raise ValueError(f"unknown model {model!r}")
    return matrix, variant


def is_model(model: str) -> bool:
    try:
        split(model)
    except ValueError:
        return False
    return True


def all_models() -> list[str]:
    """ModelFinder's grid: every model that can be a candidate."""
    return [m.name + v.suffix for m in MATRICES for v in VARIANTS]


def all_fixed_models() -> list[str]:
    """The full-search grid: every model, mixtures too."""
    return all_models() + [m.name + v.suffix for m in MATRICES for v in MIXTURE_VARIANTS]


def is_mixture(model: str) -> bool:
    """Whether a catalog model is a site-profile mixture (see the module docstring)."""
    return is_model(model) and split(model)[1] in MIXTURE_VARIANTS


def candidate_models(models: list[str]) -> list[str]:
    """These models without the ones ModelFinder can't test: what a job with no candidate
    list of its own falls back to."""
    return [model for model in models if not is_mixture(model)]


@dataclass(frozen=True)
class Candidates:
    matrices: list[str]
    rates: list[str]
    freqs: list[str]

    @property
    def count(self) -> int:
        """Models ModelFinder will compare (+R types count once; it varies categories)."""
        return len(self.matrices) * len(self.rates) * len(self.freqs)

    @property
    def shortcuts(self) -> bool:
        """Whether IQ-TREE may use its own shortcuts: every rate and frequency type is
        ticked, so -mrate and -mfreq are left out (they are IQ-TREE's defaults), and it
        tries the other rate types only on its best matrices. Giving the full lists
        explicitly turns the shortcuts off: 27 matrices on Ime2 took 101 s instead of
        21 s, for the same best model (LG+F+G4)."""
        return self.rates == ALL_RATES and self.freqs == ALL_FREQS

    @property
    def fitted(self) -> int:
        """Models ModelFinder fits without shortcuts: +R and +I+R each try 2-10 rate
        categories (IQ-TREE's "test up to N models")."""
        rates = sum(9 if "R" in rate else 1 for rate in self.rates)
        return len(self.matrices) * rates * len(self.freqs)

    def args(self) -> list[str]:
        if self.shortcuts:
            return ["-mset", ",".join(self.matrices)]
        return [
            "-mset", ",".join(self.matrices),
            "-mrate", ",".join(self.rates),
            "-mfreq", ",".join(self.freqs),
        ]  # fmt: skip

    def models(self) -> list[str]:
        """The candidate models, as catalog strings."""
        return [
            m + v.suffix
            for m in self.matrices
            for v in VARIANTS
            if v.rate in self.rates and v.freq in self.freqs
        ]


def raxml_missing(models: list[str]) -> list[str]:
    """The models RAxML-NG can't run: their matrix is in RAXML_MISSING, or they are mixtures
    ("Wrong model specification: +C20+G4")."""
    return [m for m in models if split(m)[0] in RAXML_MISSING or is_mixture(m)]


def is_default_set(models: list[str], *, raxml: bool = False) -> bool:
    """Whether ModelFinder should use IQ-TREE's default candidates: no model ticked, or
    the whole grid (with RAxML-NG, minus the matrices it lacks). The grid's axes are
    IQ-TREE's default ones (every matrix here; rates E, I, G, I+G, R, I+R; frequencies FU,
    F: checked against a default run's modelfinder.tsv), and without flags IQ-TREE
    searches that set with its own shortcuts (Candidates.shortcuts)."""
    if not models:
        return True
    candidates = modelfinder_candidates(models)
    matrices = [m.name for m in MATRICES if not (raxml and m.name in RAXML_MISSING)]
    return candidates.matrices == matrices and candidates.shortcuts


def modelfinder_args(models: list[str], *, raxml: bool = False) -> list[str]:
    """ModelFinder's candidate flags: the cross product of the ticked models, or IQ-TREE's
    default set ([]) when none or all are ticked (is_default_set), minus the matrices
    RAxML-NG lacks if the job uses RAxML-NG (IQ-TREE's default rate and frequency types
    are kept)."""
    if not is_default_set(models, raxml=raxml):
        return modelfinder_candidates(models).args()
    if raxml:
        return ["-mset", ",".join(m.name for m in MATRICES if m.name not in RAXML_MISSING)]
    return []


def modelfinder_candidates(models: list[str]) -> Candidates:
    """The cross product ModelFinder will test for these ticked models, in catalog order."""
    parts = [split(model) for model in models]
    if not parts:
        raise ValueError("no models selected")
    if any(variant in MIXTURE_VARIANTS for _, variant in parts):
        raise ValueError("site-profile mixtures are not ModelFinder candidates")
    matrices = {matrix for matrix, _ in parts}
    rates = {variant.rate for _, variant in parts}
    freqs = {variant.freq for _, variant in parts}
    return Candidates(
        matrices=[m.name for m in MATRICES if m.name in matrices],
        rates=[r for r in dict.fromkeys(v.rate for v in VARIANTS) if r in rates],
        freqs=[f for f in dict.fromkeys(v.freq for v in VARIANTS) if f in freqs],
    )


# Codon models (coding DNA only), full IQ-TREE strings like the protein ones.
# Frequencies are always explicit: a bare "MG" silently means MG+F3X4, and MG/MGK refuse
# empirical codon frequencies (+F). Every model here was fitted with IQ-TREE 3.1.3 on
# lysin (2026-10-01). Two ways to use them, as for protein models: ModelFinder chooses
# among candidates (CODON_MODELFINDER combinations), and each ticked model gets its own
# combination. ModelFinder's codon candidates are a cross product too (-mset/-mfreq/
# -mrate; checked: MG,GY,ECMK07 x F1X4,F3X4 x E,G tested exactly those 12 models).

CODON_MODELFINDER = "modelfinder_codon"  # combination key for codon model selection


@dataclass(frozen=True)
class CodonMatrix:
    name: str
    group: str
    label: str
    freqs: tuple[str, ...]  # the frequency types IQ-TREE accepts with it


@dataclass(frozen=True)
class CodonVariant:
    suffix: str  # appended to the matrix name
    freq: str  # ModelFinder -mfreq code
    rate: str  # ModelFinder -mrate code
    label: str


CODON_MATRICES: tuple[CodonMatrix, ...] = (
    CodonMatrix("MG", "Mechanistic", "Muse-Gaut 1994: dN/dS", ("F1X4", "F3X4")),
    CodonMatrix("MGK", "Mechanistic", "MG + transition/transversion ratio", ("F1X4", "F3X4")),
    CodonMatrix("GY0K", "Mechanistic", "Goldman-Yang 1994: dN/dS", ("F1X4", "F3X4", "F")),
    CodonMatrix(
        "GY", "Mechanistic", "GY + transition/transversion ratio (codeml M0)",
        ("F1X4", "F3X4", "F"),
    ),
    CodonMatrix("ECMK07", "Empirical", "Kosiol et al. 2007 (KOSI07)", ("F1X4", "F3X4", "F")),
    CodonMatrix(
        "ECMrest", "Empirical", "Kosiol et al. 2007, one nucleotide change at a time",
        ("F1X4", "F3X4", "F"),
    ),
    CodonMatrix("ECMS05", "Empirical", "Schneider et al. 2005 (SCHN05)", ("F1X4", "F3X4", "F")),
)  # fmt: skip
_RATES = (("", "E", ""), ("+G4", "G", " + gamma"), ("+I+G4", "I+G", " + invariant + gamma"),
          ("+R4", "R", " + FreeRate"))  # fmt: skip
_FREQS = (("F3X4", "F3X4"), ("F1X4", "F1X4"), ("F", "empirical codon frequencies"))
CODON_VARIANTS: tuple[CodonVariant, ...] = tuple(
    CodonVariant(f"+{freq}{suffix}", freq, rate, f"{label}{rate_label}")
    for freq, label in _FREQS
    for suffix, rate, rate_label in _RATES
)
# IQ-TREE's own default codon candidates (its ModelFinder without flags; read from a
# run's report): these matrices (KOSI07 = ECMK07, SCHN05 = ECMS05) with F1X4, F3X4, F
# (and the empirical models' own frequencies) and rates E, I, G, I+G, R, I+R, searched
# with its shortcuts (38 of up to 330 models fitted on lysin, 76 s).
IQTREE_DEFAULT_CODON_MATRICES = ("MG", "MGK", "GY", "ECMK07", "ECMS05")
# Suggested for a full search under one codon model (GY+F3X4 is codeml's M0).
DEFAULT_CODON_MODELS = ["GY+F3X4+G4"]


def all_codon_models() -> list[str]:
    return [m.name + v.suffix for m in CODON_MATRICES for v in CODON_VARIANTS if v.freq in m.freqs]


def is_codon_model(model: str) -> bool:
    return model in _CODON_MODELS


def is_modelfinder(key: str) -> bool:
    """Whether a combination's model key is a ModelFinder choice (protein or codon)."""
    return key in (MODELFINDER, CODON_MODELFINDER)


_CODON_MODELS = frozenset(all_codon_models())
_CODON_SPLIT = {
    m.name + v.suffix: (m.name, v)
    for m in CODON_MATRICES
    for v in CODON_VARIANTS
    if v.freq in m.freqs
}


def codon_candidates(models: list[str]) -> Candidates:
    """The cross product codon ModelFinder will test for these ticked codon models."""
    parts = [_CODON_SPLIT[m] for m in models if m in _CODON_SPLIT]
    if not parts:
        raise ValueError("no codon models selected")
    matrices = {matrix for matrix, _ in parts}
    rates = {variant.rate for _, variant in parts}
    freqs = {variant.freq for _, variant in parts}
    return Candidates(
        matrices=[m.name for m in CODON_MATRICES if m.name in matrices],
        rates=[r for r in dict.fromkeys(v.rate for v in CODON_VARIANTS) if r in rates],
        freqs=[f for f in dict.fromkeys(v.freq for v in CODON_VARIANTS) if f in freqs],
    )


def is_default_codon_set(models: list[str]) -> bool:
    """No codon candidate ticked, or the whole grid: IQ-TREE's own default codon set."""
    return not models or set(models) >= _CODON_MODELS


def codon_modelfinder_args(models: list[str]) -> list[str]:
    """Codon ModelFinder's flags: [] (IQ-TREE's default set) or the ticked cross product,
    always with explicit rates and frequencies (codon data has no shortcut equivalent of
    the protein grid's)."""
    if is_default_codon_set(models):
        return []
    candidates = codon_candidates(models)
    return [
        "-mset", ",".join(candidates.matrices),
        "-mrate", ",".join(candidates.rates),
        "-mfreq", ",".join(candidates.freqs),
    ]  # fmt: skip


def codon_candidate_count(models: list[str]) -> int:
    """Models codon ModelFinder fits: ~38 for IQ-TREE's default set (with shortcuts), or
    the cross product (each +R type tries several category counts)."""
    if is_default_codon_set(models):
        return 38
    return codon_candidates(models).fitted


def label(key: str) -> str:
    """Display label for a model combination key."""
    return {MODELFINDER: "ModelFinder", CODON_MODELFINDER: "ModelFinder (codon)"}.get(key, key)
