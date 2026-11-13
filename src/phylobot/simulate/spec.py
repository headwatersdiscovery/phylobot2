"""Simulation specs: what to simulate, checked before anything runs.

A spec is the golden `spec.toml` schema, read from TOML or JSON into frozen dataclasses:

    description = "..."
    [tree]        file = "...", optional prune = [...], then outgroup = "..." or
                  midpoint = true (roots the tree there; AliSim evolves from the root)
                  -- or a random tree: random = "yh" | "u" | "cat" | "bal" | "bd",
                  taxa = n, lengths = [min, mean, max], and birth/death rates for "bd"
    [root]        file = "...", name = "..." (a sequence in that file)
                  -- or a random root: length = sites (amino acids, or codons)
    [model]       substitution = "LG+G4{0.747}" (an IQ-TREE model string), optional
                  seqtype = "CODON1"; or, for codeml's M3, seqtype, kappa, f3x4 and
                  [model.site_classes] proportions = [...], omegas = [...]
    [indels]      rates = "ins,del", sizes = "GEO{3},GEO{3}" (AliSim's syntax)
    [family]      duplication and loss rates (per copy per unit branch length) on the
                  [tree], read as a rooted species tree; optional origin (a stem above
                  its root), min_genes, and chosen events: duplications = [{clade =
                  [species...], at = 0.5}], losses = [{clade = [...], copy = "1b",
                  at = 0.5}] (families.py)
    [simulation]  seeds = [...] (the first is the golden quick tier), optional image
                  (the pinned IQ-TREE container), names = "keep" | "neutral"

Users will reach this through the web form, so every value that
ends up on AliSim's command line is checked against a strict pattern or bound here; the
commands are argument lists, never shell strings. File names in a spec are resolved by
the caller (`input_files`), never inside the pipeline: there the tree and root come in
as staged files.
"""

import json
import math
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

SHAPES = {
    "yh": "Yule-Harding",
    "u": "uniform",
    "cat": "caterpillar",
    "bal": "balanced",
    "bd": "birth-death",
}
MAX_TAXA = 5000
MAX_ROOT_LENGTH = 100_000
MAX_SEEDS = 1000
MAX_FORCED = 50

# IQ-TREE model strings: names, numbers and their separators only.
MODEL_PATTERN = re.compile(r"[A-Za-z0-9+{}.,/_-]{1,2000}")
SEQTYPE_PATTERN = re.compile(r"CODON[0-9]{1,2}")
RATES_PATTERN = re.compile(r"[0-9]*\.?[0-9]+(,[0-9]*\.?[0-9]+)?")
SIZE = r"(GEO|NB|POW|LAV)\{[0-9.]+(/[0-9.]+){0,2}\}"
SIZES_PATTERN = re.compile(rf"{SIZE}(,{SIZE})?")
NAME_PATTERN = re.compile(r"[A-Za-z0-9_.|-]{1,200}")
COPY_PATTERN = re.compile(r"[1-9][0-9]?[ab]")
IMAGE_PATTERN = re.compile(r"[a-z0-9./_-]+:[A-Za-z0-9._-]+")


class SpecError(ValueError):
    """A spec that cannot be simulated, with a message for the person who wrote it."""


@dataclass(frozen=True)
class FileTree:
    file: str = ""  # resolved by the caller; empty when the tree is staged
    prune: tuple[str, ...] = ()
    outgroup: str | None = None
    midpoint: bool = False


@dataclass(frozen=True)
class RandomTree:
    shape: str
    taxa: int
    lengths: tuple[float, float, float]  # min, mean, max branch length
    birth: float | None = None
    death: float | None = None

    @property
    def argument(self) -> str:
        """AliSim's -t value, e.g. RANDOM{yh/20} or RANDOM{bd{0.1/0.05}/20}."""
        shape = f"bd{{{self.birth}/{self.death}}}" if self.shape == "bd" else self.shape
        return f"RANDOM{{{shape}/{self.taxa}}}"


@dataclass(frozen=True)
class FileRoot:
    name: str
    file: str = ""


@dataclass(frozen=True)
class RandomRoot:
    length: int  # sites: amino acids, or codons


@dataclass(frozen=True)
class SiteClasses:
    proportions: tuple[float, ...]
    omegas: tuple[float, ...]


@dataclass(frozen=True)
class Model:
    substitution: str | None = None
    seqtype: str | None = None  # None: amino acids; CODONn: codons, NCBI code n
    kappa: float | None = None
    f3x4: tuple[float, ...] | None = None
    site_classes: SiteClasses | None = None

    @property
    def codon(self) -> bool:
        return self.seqtype is not None


@dataclass(frozen=True)
class Indels:
    rates: str
    sizes: str


@dataclass(frozen=True)
class ForcedDuplication:
    clade: tuple[str, ...]  # the species below the branch
    at: float = 0.5  # fraction of the branch's length, from its top


@dataclass(frozen=True)
class ForcedLoss:
    clade: tuple[str, ...]
    copy: str  # "1b": copy b of forced duplication 1, and everything descended from it
    at: float = 0.5


@dataclass(frozen=True)
class Family:
    duplication: float  # per gene copy per unit of species-tree branch length
    loss: float
    origin: float = 0.0  # length of a stem above the species tree's root
    min_genes: int = 4  # redraw until at least this many copies survive
    duplications: tuple[ForcedDuplication, ...] = ()
    losses: tuple[ForcedLoss, ...] = ()


@dataclass(frozen=True)
class Spec:
    tree: FileTree | RandomTree
    root: FileRoot | RandomRoot
    model: Model
    seeds: tuple[int, ...]
    indels: Indels | None = None
    family: Family | None = None
    image: str | None = None
    names: str = "keep"  # or "neutral": sim001, sim002, ... (names.tsv maps them back)
    description: str = ""
    extra: dict[str, Any] = field(default_factory=dict, compare=False)


def load(path: Path) -> Spec:
    """A spec from a .toml or .json file."""
    text = path.read_text()
    data = tomllib.loads(text) if path.suffix == ".toml" else json.loads(text)
    return parse(data)


def input_files(spec: Spec, base: Path) -> tuple[Path | None, Path | None]:
    """The tree and root files a spec names, relative to its directory."""
    tree = base / spec.tree.file if isinstance(spec.tree, FileTree) and spec.tree.file else None
    root = base / spec.root.file if isinstance(spec.root, FileRoot) and spec.root.file else None
    return tree, root


def parse(data: dict[str, Any]) -> Spec:
    if not isinstance(data, dict):
        raise SpecError("a spec is a table of sections")
    model = _model(_section(data, "model"))
    tree = _tree(_section(data, "tree"))
    root = _root(_section(data, "root"))
    indels = _indels(data["indels"]) if "indels" in data else None
    family = _family(data["family"]) if "family" in data else None
    simulation = _section(data, "simulation")
    seeds = simulation.get("seeds")
    if not isinstance(seeds, list) or not seeds or len(seeds) > MAX_SEEDS:
        raise SpecError(f"simulation.seeds: a list of 1 to {MAX_SEEDS} seeds")
    if not all(_is_int(s) and 0 < s < 2**31 for s in seeds) or len(set(seeds)) != len(seeds):
        raise SpecError("simulation.seeds: distinct whole numbers from 1 to 2^31 - 1")
    image = simulation.get("image")
    pinned = isinstance(image, str) and IMAGE_PATTERN.fullmatch(image)
    if image is not None and (not pinned or image.endswith(":latest")):
        raise SpecError("simulation.image: a container pinned to an exact tag")
    names = simulation.get("names", "keep")
    if names not in ("keep", "neutral"):
        raise SpecError('simulation.names: "keep" or "neutral"')

    if model.site_classes is not None:
        if indels is not None:
            raise SpecError("site classes are simulated class by class: indels are not supported")
        if not isinstance(root, FileRoot):
            raise SpecError("site classes need a root sequence file")
    if family is not None and not isinstance(tree, FileTree):
        raise SpecError("[family] grows on a species tree: the tree must be a file")
    if family is not None and family.min_genes > MAX_TAXA:
        raise SpecError(f"family.min_genes: at most {MAX_TAXA}")
    if isinstance(tree, RandomTree) and tree.taxa > MAX_TAXA:
        raise SpecError(f"tree.taxa: at most {MAX_TAXA}")
    description = data.get("description", "")
    if not isinstance(description, str):
        raise SpecError("description: text")
    known = {"description", "tree", "root", "model", "indels", "family", "simulation"}
    return Spec(
        tree=tree,
        root=root,
        model=model,
        seeds=tuple(seeds),
        indels=indels,
        family=family,
        image=image,
        names=names,
        description=description,
        extra={k: v for k, v in data.items() if k not in known},
    )


def _section(data: dict[str, Any], name: str) -> dict[str, Any]:
    section = data.get(name)
    if not isinstance(section, dict):
        raise SpecError(f"[{name}] is missing")
    return section


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _number(value: Any, label: str, low: float = 0.0, high: float = math.inf) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise SpecError(f"{label}: a number")
    number = float(value)
    if not (low < number <= high) or not math.isfinite(number):
        raise SpecError(f"{label}: must be above {low:g} and at most {high:g}")
    return value  # as written: it is formatted into AliSim's arguments


def _name(value: Any, label: str) -> str:
    if not isinstance(value, str) or not NAME_PATTERN.fullmatch(value):
        raise SpecError(f"{label}: letters, digits and _ . | - only")
    return value


def _file(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or "\0" in value:
        raise SpecError(f"{label}: a file name")
    return value


def _tree(section: dict[str, Any]) -> FileTree | RandomTree:
    if "random" in section:
        shape = section["random"]
        if shape not in SHAPES:
            raise SpecError(f"tree.random: one of {', '.join(SHAPES)}")
        taxa = section.get("taxa")
        if not isinstance(taxa, int) or isinstance(taxa, bool) or taxa < 3:
            raise SpecError("tree.taxa: a whole number, at least 3")
        lengths = section.get("lengths")
        if not isinstance(lengths, list) or len(lengths) != 3:
            raise SpecError("tree.lengths: [min, mean, max] branch length")
        low, mean, high = (_number(x, "tree.lengths", high=10.0) for x in lengths)
        if not low <= mean <= high or low == high:
            raise SpecError("tree.lengths: min <= mean <= max, with min below max")
        birth: float | None = None
        death: float | None = None
        if shape == "bd":
            birth = _number(section.get("birth"), "tree.birth", high=100.0)
            death = section.get("death", 0)
            if death != 0:
                death = _number(death, "tree.death", high=100.0)
            if not isinstance(death, int | float) or death >= birth:
                raise SpecError("tree.death: below the birth rate")
        return RandomTree(shape, taxa, (low, mean, high), birth, death)
    prune = section.get("prune", [])
    if not isinstance(prune, list) or not all(isinstance(p, str) for p in prune):
        raise SpecError("tree.prune: a list of leaf names")
    outgroup = section.get("outgroup")
    if outgroup is not None and not isinstance(outgroup, str):
        raise SpecError("tree.outgroup: a leaf name")
    midpoint = section.get("midpoint", False)
    if not isinstance(midpoint, bool) or (midpoint and outgroup):
        raise SpecError("tree: root on an outgroup or at the midpoint, not both")
    file = _file(section["file"], "tree.file") if "file" in section else ""
    return FileTree(file, tuple(prune), outgroup, midpoint)


def _root(section: dict[str, Any]) -> FileRoot | RandomRoot:
    if "length" in section:
        length = section["length"]
        if not _is_int(length) or not 0 < length <= MAX_ROOT_LENGTH:
            raise SpecError(f"root.length: a whole number from 1 to {MAX_ROOT_LENGTH}")
        return RandomRoot(length)
    file = _file(section["file"], "root.file") if "file" in section else ""
    return FileRoot(_name(section.get("name"), "root.name"), file)


def _model(section: dict[str, Any]) -> Model:
    seqtype = section.get("seqtype")
    if seqtype is not None and (
        not isinstance(seqtype, str) or not SEQTYPE_PATTERN.fullmatch(seqtype)
    ):
        raise SpecError("model.seqtype: CODONn (n = NCBI genetic code), or absent for proteins")
    if "site_classes" in section:
        if seqtype is None:
            raise SpecError("model.site_classes: codon specs only")
        classes = section["site_classes"]
        proportions, omegas = classes.get("proportions"), classes.get("omegas")
        if not isinstance(proportions, list) or not isinstance(omegas, list):
            raise SpecError("model.site_classes: proportions and omegas lists")
        if len(proportions) != len(omegas) or not 2 <= len(omegas) <= 10:
            raise SpecError("model.site_classes: 2 to 10 classes, as many proportions as omegas")
        for p in proportions:
            _number(p, "model.site_classes.proportions", high=1.0)
        for w in omegas:
            _number(w, "model.site_classes.omegas", high=1000.0)
        if abs(sum(proportions) - 1) > 1e-3:
            raise SpecError("model.site_classes.proportions: must sum to 1")
        kappa = _number(section.get("kappa"), "model.kappa", high=1000.0)
        f3x4 = section.get("f3x4")
        if not isinstance(f3x4, list) or len(f3x4) != 12:
            raise SpecError("model.f3x4: 12 frequencies (A C G T at each codon position)")
        for f in f3x4:
            _number(f, "model.f3x4", high=1.0)
        if "substitution" in section:
            raise SpecError("model: site classes replace the substitution string")
        if seqtype != "CODON1":
            raise SpecError("model.site_classes: the standard genetic code (CODON1) only")
        return Model(
            seqtype=seqtype,
            kappa=kappa,
            f3x4=tuple(f3x4),
            site_classes=SiteClasses(tuple(proportions), tuple(omegas)),
        )
    substitution = section.get("substitution")
    if not isinstance(substitution, str) or not MODEL_PATTERN.fullmatch(substitution):
        raise SpecError("model.substitution: an IQ-TREE model string, e.g. LG+G4{0.75}")
    return Model(substitution=substitution, seqtype=seqtype)


def _indels(section: Any) -> Indels:
    if not isinstance(section, dict):
        raise SpecError("[indels]: rates and sizes")
    rates, sizes = section.get("rates"), section.get("sizes")
    if not isinstance(rates, str) or not RATES_PATTERN.fullmatch(rates):
        raise SpecError('indels.rates: "insertion,deletion" rates, e.g. "0.03,0.03"')
    if any(float(r) > 10 for r in rates.split(",")):
        raise SpecError("indels.rates: at most 10 (relative to substitutions)")
    if not isinstance(sizes, str) or not SIZES_PATTERN.fullmatch(sizes):
        raise SpecError('indels.sizes: AliSim size distributions, e.g. "GEO{3},GEO{3}"')
    return Indels(rates, sizes)


def _rate(value: Any, label: str) -> float:
    """A rate of zero or more (zero turns the event off)."""
    if value == 0 and not isinstance(value, bool):
        return 0.0
    return float(_number(value, label, high=100.0))


def _fraction(value: Any, label: str) -> float:
    number = float(_number(value, label, high=1.0))
    if number >= 1:
        raise SpecError(f"{label}: a fraction of the branch, above 0 and below 1")
    return number


def _clade(value: Any, label: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not value or len(value) > MAX_TAXA:
        raise SpecError(f"{label}: a list of the species below the branch")
    clade = tuple(_name(v, label) for v in value)
    if len(set(clade)) != len(clade):
        raise SpecError(f"{label}: each species once")
    return clade


def _family(section: Any) -> Family:
    if not isinstance(section, dict):
        raise SpecError("[family]: duplication and loss rates")
    unknown = set(section) - {
        "duplication", "loss", "origin", "min_genes", "duplications", "losses"
    }  # fmt: skip
    if unknown:
        raise SpecError(f"[family]: unknown {', '.join(sorted(unknown))}")
    duplication = _rate(section.get("duplication"), "family.duplication")
    loss = _rate(section.get("loss"), "family.loss")
    origin = section.get("origin", 0)
    origin = (
        0.0
        if origin == 0 and not isinstance(origin, bool)
        else float(_number(origin, "family.origin", high=10.0))
    )
    min_genes = section.get("min_genes", 4)
    if not _is_int(min_genes) or min_genes < 3:
        raise SpecError("family.min_genes: a whole number, at least 3")

    def listed(key: str) -> list[dict[str, Any]]:
        items = section.get(key, [])
        if not isinstance(items, list) or len(items) > MAX_FORCED:
            raise SpecError(f"family.{key}: a list of at most {MAX_FORCED}")
        if not all(isinstance(item, dict) for item in items):
            raise SpecError(f"family.{key}: each a table")
        return items

    duplications = tuple(
        ForcedDuplication(
            _clade(d.get("clade"), f"family.duplications[{i}].clade"),
            _fraction(d.get("at", 0.5), f"family.duplications[{i}].at"),
        )
        for i, d in enumerate(listed("duplications"), 1)
    )
    losses = []
    for i, d in enumerate(listed("losses"), 1):
        copy = d.get("copy")
        if not isinstance(copy, str) or not COPY_PATTERN.fullmatch(copy):
            raise SpecError(f'family.losses[{i}].copy: e.g. "1b", copy b of duplication 1')
        if int(copy[:-1]) > len(duplications):
            raise SpecError(f"family.losses[{i}].copy: there is no duplication {copy[:-1]}")
        losses.append(
            ForcedLoss(
                _clade(d.get("clade"), f"family.losses[{i}].clade"),
                copy,
                _fraction(d.get("at", 0.5), f"family.losses[{i}].at"),
            )
        )
    return Family(duplication, loss, origin, min_genes, duplications, tuple(losses))
