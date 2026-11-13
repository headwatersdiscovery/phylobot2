"""How much computing a job asks for, before it runs: the New Job page's work estimate.

A work unit is one reference CPU-minute: a minute of one CPU on the machine where the
costs below were measured (the owner's 8-core laptop, amd64 tools under emulation,
2026-09-30, each step timed alone). Units mean the same on
every installation, so they compare jobs and could later be priced; `seconds_per_unit`
converts them to this server's time, from its own jobs' `cpu_seconds` and the units
estimated when they were submitted (`Job.work_units`).

Each step's cost is c x (N/20)^a x (L/400)^b for N sequences of mean length L (gaps
removed): the exponents come from simulated sweeps (N 10-80, L 200-1600), the levels
`c` from real data (Ime2 and UniRef90 hits), which are harder than simulated sequences.
`tasks` lists one entry per pipeline task, the same tasks `stages.expected_kinds` counts
(tested). Within a step, costs depend on the options that measurably change them: rate
heterogeneity, branch support, and ModelFinder's candidates.
"""

from dataclasses import dataclass
from statistics import median
from typing import Any

from phylobot import catalog, properties
from phylobot.executors import (
    alignment_count,
    as_list,
    branch_site,
    branch_site_runs,
    candidate_models,
    codeml_combinations,
    codeml_runs,
    coding,
    codon_amino_trees,
    codon_search,
    codon_tree_count,
    combination_count,
    esr,
    esr_tips,
    modelfinder_runs,
    planned_combinations,
    properties_only,
    property_targets,
    provided,
    reconciles,
    simulation,
    simulation_seeds,
    uncertainty_only,
)
from phylobot.models import Job


@dataclass(frozen=True)
class Size:
    sequences: int
    length: float  # mean residues per sequence, gaps removed


def size_record(size: Size) -> dict[str, float]:
    """Job.input_size: the size an estimate used, so it can be redone per stage later."""
    return {"sequences": size.sequences, "length": round(size.length, 1)}


def recorded_size(record: dict[str, Any]) -> Size | None:
    try:
        return Size(int(record["sequences"]), float(record["length"]))
    except KeyError, TypeError, ValueError:
        return None


REFERENCE_SIZE = Size(20, 400)
IME2_SIZE = Size(7, 600)  # the bundled example; assumed for jobs without an estimate


@dataclass(frozen=True)
class Cost:
    seconds: float  # at REFERENCE_SIZE
    a: float = 0.0  # exponent of sequence count
    b: float = 0.0  # exponent of length

    def at(self, size: Size) -> float:
        n = size.sequences / REFERENCE_SIZE.sequences
        length = size.length / REFERENCE_SIZE.length
        return float(self.seconds * n**self.a * length**self.b)


# Fitted to those timings (exponents from the simulated sweeps, levels scaled by
# real data: real proteins took 1.1-4x the simulated time at the same size, PRANK 4x).
# Steps of a second or two (Nextflow rounds to whole seconds) have hand-set exponents
# from the same points; the smallest are fixed. Container start-up is not counted.
STEPS = {
    "VALIDATE": Cost(0.7),
    "SUMMARY": Cost(0.5),
    # Parses every combination's tree (~0.01 s each at the reference size; unmeasured).
    "ANCESTOR_MAP": Cost(0.5),
    "ALIGN_MAFFT": Cost(3.0, 1.3, 0.9),
    "ALIGN_MUSCLE": Cost(9.1, 2.0, 1.9),
    "ALIGN_CLUSTALO": Cost(0.8, 1.2, 1.0),
    "ALIGN_PRANK": Cost(88.0, 1.33, 1.85),
    # PRANK -codon on coding DNA: 4x protein PRANK on the same sequences (one measurement:
    # lysin, 25 x 130 codons, 24 s vs 6 s, 2026-10-01); protein PRANK's exponents.
    "ALIGN_PRANK_CODON": Cost(4 * 88.0, 1.33, 1.85),
    "CODON_ALIGNMENT": Cost(0.3),
    "CODEML_EXPORT": Cost(2.0),
    "BRANCH_SITE_EXPORT": Cost(0.5),
    "EXPORT_ALIGNMENT": Cost(0.5),
    "EXPORT_TREE": Cost(0.5),
    "EXPORT": Cost(0.65, 0.2, 0.2),
    "ROOT_TREE": Cost(0.5),
    # Tier-0 properties: pure Python over params.property_draws draws per node, so it
    # scales with sequence count (nodes) and length (per draw). Priced from the core
    # properties at 1000 draws on Ime2; to be remeasured on larger families.
    "PROPERTIES": Cost(3.0, 0.5, 0.5),
    # Uncertainty: numpy over cells x columns per ancestor, so priced per
    # combination, scaling with nodes (sequences) and columns (length). One measurement,
    # 9 Ime2 combinations in 0.14 s on the laptop's host Python (2026-10-22), plus start-up
    # and reading: provisional, **owner review** with the rest of the work estimate.
    "UNCERTAINTY": Cost(0.2, 1.0, 1.0),
    "GAPS_ML": Cost(0.5),
    "GAPS_PRANK": Cost(0.5),
    "GAPS_GRASP": Cost(4.0, 0.45, 0.5),
    # Simulations (measured 2026-10-14): AliSim takes under 0.3 s up
    # to 300 x 1000 with indels, 2.7 s at 1000 x 1000 (sequences x sites); priced at
    # its size, the two Python steps fixed.
    # Reconciliation: GeneRax EVAL (fixed topology) took about a second on Ime2 (7
    # genes, 2026-10-16); assumed linear in sequences until measured on larger families.
    "SPECIES_TREE": Cost(0.5),
    "RECONCILE": Cost(1.0, 1.0, 0.0),
    "RECONCILE_EXPORT": Cost(0.5),
    # GeneRax's SPR search took 0.5-1.6x IQ-TREE's search on the same family
    # (2026-10-16): priced as an IQ-TREE search, its
    # fit and the AU test (a fit per tree source, RELL resampling) as IQ-TREE fits (below).
    "SIM_PREPARE": Cost(0.5),
    "SIM_ALISIM": Cost(0.1, 1.0, 1.0),
    "SIM_COLLECT": Cost(0.5, 0.5, 0.5),
    # ESR's scoring: reads every tip's rows (pure Python).
    "ESR_COLLECT": Cost(0.5, 0.5, 0.5),
}
# Model-dependent steps: without / with rate heterogeneity (+G4, +R4; ModelFinder's pick
# is assumed to have it, as it nearly always does).
TREE_IQTREE = {False: Cost(21.7, 1.57, 0.75), True: Cost(74.9, 1.69, 0.81)}
TREE_RAXML = {False: Cost(9.5, 1.93, 0.99), True: Cost(27.7, 1.93, 1.01)}
FIT_IQTREE = {False: Cost(0.5, 0.7, 0.5), True: Cost(1.2, 1.1, 0.5)}
ANCESTRAL_IQTREE = {False: Cost(0.4, 0.7, 0.3), True: Cost(1.0, 1.0, 0.3)}


@dataclass(frozen=True)
class CostFactor:
    search: float
    fit: float
    ancestors: float


# Codon models against protein models on the same data (lysin, 25 x 135 codons,
# MAFFT, 2026-10-01): IQ-TREE search with UFBoot GY+F3X4+G4 8 min 54 s vs ModelFinder's
# pick (Q.PFAM+G4) 1 min 6 s; fit to a provided tree 18 s vs 1 s; reconstruction 4-18 s
# vs 1 s. One dataset: refit when there are more.
CODON_FACTOR = CostFactor(search=8.0, fit=15.0, ancestors=10.0)
# A site-profile mixture against the same model with gamma alone, -mwopt on. Measured on SNCA
# (115 x 138; fit and ancestors 11-12 s against 1 s) and Ime2 (7 x 1052; a search 3 min
# against 7 s; SNCA 23 min against 98 s). Provisional, from two
# families: **owner review** with the rest of the work estimate.
MIXTURE_FACTOR = CostFactor(search=20, fit=12, ancestors=12)
# Codon ModelFinder with IQ-TREE's default set on lysin (25 x 135 codons, its own
# starting tree, 38 models fitted): 76 s alone, 202 s inside a pipeline run beside other
# steps (the figure used); a ticked subset costs in proportion to the models it fits (12
# models: 24 s alone). Scaled to the reference size with the protein ModelFinder's
# exponents (an assumption until measured on more data).
CODON_MODELFINDER_DEFAULT = Cost(202 / ((25 / 20) ** 1.46 * (135 / 400) ** 0.72), 1.46, 0.72)
# codeml (PAML 4.10.10), one run with ancestors (verbose 2) per model and starting omega,
# on lysin (25 x 135 codons, 2026-10-01; six runs side by side): M0 145 s, M1a 121 s,
# M2a 237 s, M3 198 s, M7 479 s, M8 852 s; free ratios (FR, beside M0 only, 2026-10-02)
# 1056 s. Assumed linear in sequences and length until measured on more data (FR, with
# a parameter per branch, likely grows faster with sequences).
_LYSIN = (25 / 20) * (135 / 400)
CODEML = {
    key: Cost(seconds / _LYSIN, 1.0, 1.0)
    for key, seconds in {"M0": 145, "M1a": 121, "M2a": 237, "M3": 198, "M7": 479,
                         "M8": 852, "FR": 1056}.items()
}  # fmt: skip
# The branch-site test (verbose 2, no ancestors) on lysin with one foreground branch
# (2026-10-02): Model A 389 s, its null 186 s.
BRANCH_SITE = {
    key: Cost(seconds / _LYSIN, 1.0, 1.0) for key, seconds in {"A": 389, "A_null": 186}.items()
}
# Branch support: no measurable cost at up to 40 sequences (IQ-TREE +G4 search, simulated
# N=20: 49 s without support, 50 s with SH-aLRT, 48 s with SH-aLRT and UFBoot).
UFBOOT_PER_1000 = 0.0  # extra fraction of a tree search per 1000 replicates
ALRT_PER_1000 = 0.0
# ModelFinder: with IQ-TREE's shortcuts (every rate and frequency type ticked,
# catalog.Candidates.shortcuts) per matrix: the default set (28 matrices) grew as
# N^1.46 L^0.72, real data 2.3x simulated (63 s at N=20, L=400). Without shortcuts, per
# fitted model (Candidates.fitted): too few points for its own exponents, so the same
# scaling with every model given explicitly costing 5x the default set (measured 2.3-5.6x;
# real Ime2 4.7x, 97 s).
MODELFINDER_PER_MATRIX = Cost(63.4 / 28, 1.46, 0.72)
MODELFINDER_PER_MODEL = Cost(5 * 63.4 / 1232, 1.46, 0.72)


def modelfinder_seconds(ticked: list[str], raxml: bool, size: Size) -> float:
    """ModelFinder's cost for these ticked candidates (none: IQ-TREE's default)."""
    if ticked:
        candidates = catalog.modelfinder_candidates(ticked)
    else:  # every matrix (catalog.modelfinder_args), all types
        candidates = catalog.modelfinder_candidates(catalog.all_models())
    matrices = [m for m in candidates.matrices if not (raxml and m in catalog.RAXML_MISSING)]
    if candidates.shortcuts:
        return MODELFINDER_PER_MATRIX.at(size) * len(matrices)
    return MODELFINDER_PER_MODEL.at(size) * candidates.fitted


SECONDS_PER_UNIT = 60.0

# Display groups, in pipeline order; ModelFinder is its own group because it is often
# the largest single cost and the user controls it directly.
GROUPS = {
    "VALIDATE": "setup",
    "SUMMARY": "setup",
    "ANCESTOR_MAP": "setup",
    "EXPORT_ALIGNMENT": "align",
    "MODELFINDER": "modelfinder",
    "MODELFINDER_CODON": "modelfinder",
    "CODEML": "codeml",
    "CODEML_EXPORT": "codeml",
    "BRANCH_SITE_EXPORT": "codeml",
    "TREE_IQTREE": "trees",
    "TREE_RAXML": "trees",
    "FIT_IQTREE": "trees",
    "FIT_CODON_AA": "trees",
    "EXPORT_TREE": "trees",
    "SPECIES_TREE": "trees",
    "RECONCILE": "trees",
    "RECONCILE_EXPORT": "trees",
    "TREE_GENERAX": "trees",
    "FIT_GENERAX": "trees",
    "TOPOLOGY_TEST": "trees",
    # Measuring properties is its own group, not ancestor work: it is the one step whose
    # cost the user sets directly (draws x properties), and with a neural predictor it can
    # dwarf the reconstruction it follows -- about 1.5 h per ancestor per gap version for
    # ESM-2 at 1000 draws, against seconds for the pure-Python properties. Counted under
    # "ancestors" it would be invisible inside the step it is cheapest to mistake for.
    "PROPERTIES": "properties",
    "SIM_PREPARE": "simulation",
    "SIM_ALISIM": "simulation",
    "SIM_COLLECT": "simulation",
    # ESR reconstructs extant sequences as ancestors: ancestor work.
    "ESR_FIT": "ancestors",
    "ESR_TIP": "ancestors",
    "ESR_COLLECT": "ancestors",
}
GROUP_LABELS = {
    "setup": "Checks",
    "align": "Alignments",
    "modelfinder": "ModelFinder",
    "trees": "Trees",
    "ancestors": "Ancestors",
    "properties": "Ancestral properties",
    "codeml": "codeml",
    "simulation": "Simulation",
}


def group(kind: str) -> str:
    if kind in GROUPS:
        return GROUPS[kind]
    return "align" if kind.startswith("ALIGN_") else "ancestors"


def heterogeneous(model: str) -> bool:
    """Whether a model (catalog string, or "modelfinder") has rate heterogeneity."""
    if catalog.is_modelfinder(model):
        return True
    if catalog.is_codon_model(model):
        return "+G" in model or "+R" in model
    return catalog.split(model)[1].rate not in ("E", "I")


def count(value: Any) -> int:
    try:
        return max(0, int(value or 0))
    except TypeError, ValueError:
        return 0


#: The gap versions PROPERTIES draws from: the unmasked reconstruction plus each masking
#: method run (PRANK's and GRASP's carry no posteriors; `pipeline.properties`).
MASKED_GAP_METHODS = ("parsimony", "fitch", "ml")


def neural_property_seconds(
    params: dict[str, Any], length: float, ancestors: int, versions: int
) -> float:
    """What the neural predictors add to one PROPERTIES task, in reference CPU-seconds.

    `PROPERTIES`' own cost covers the pure-Python properties, which are cheap at any
    size. A neural predictor costs per draw (ESM-2: about a second per 400 residues on
    the reference laptop), so it is priced from the registry's per-property figure x
    draws x length x ancestors x gap versions — the same figure the ancestor viewer's
    checklist shows (`properties_run.cost_of`). An estimate from single measurements
    (there is no cost model yet): **owner review**, since work units are
    what billing will price.
    """
    keys = [k for k in as_list(params.get("properties")) if k in properties.REGISTRY]
    neural = [properties.REGISTRY[k] for k in keys if properties.REGISTRY[k].image == properties.ML]
    if not neural:
        return 0.0
    draws = count(params.get("property_draws")) or properties.DEFAULT_DRAWS
    scale = (draws / 1000) * (max(length, 1.0) / 400) * max(ancestors, 1) * max(versions, 1)
    return sum(p.cost for p in neural) * scale


def esr_tasks(params: dict[str, Any], size: Size) -> list[tuple[str, float]]:
    """An ESR job: ESR_FIT fits the model on the fixed tree twice (the fit and
    its check) and, for an uploaded tree, its branch lengths first; each ESR_TIP is a
    reconstruction per level, E1 with a refit before it. Priced as the pipeline's own
    IQ-TREE fits and reconstructions, which the ESR spike's timings match at 100 tips
    (E0 0.9 s, E1 3.2 s a tip under emulation, LG+G4); a richer model costs more (E1
    about 60 s a tip under LG+FO+G12 on 76 tips). Provisional: **owner review** with the
    rest of the work estimate."""
    given = params.get("esr_of")
    plan: dict[str, Any] = given if isinstance(given, dict) else {}
    model = str(plan.get("model") or "modelfinder")
    rates = heterogeneous(model) if catalog.is_model(model) else True
    fit, ancestral = FIT_IQTREE[rates].at(size), ANCESTRAL_IQTREE[rates].at(size)
    levels = [lvl for lvl in plan.get("levels") or ["e0"] if lvl in ("e0", "e1")]
    tip = sum(ancestral + (fit if lvl == "e1" else 0.0) for lvl in levels)
    first = fit * (3 if plan.get("fit_lengths") else 2)
    return [("ESR_FIT", first), *[("ESR_TIP", tip)] * esr_tips(params),
            ("ESR_COLLECT", STEPS["ESR_COLLECT"].at(size))]  # fmt: skip


def tasks(params: dict[str, Any], size: Size = IME2_SIZE) -> list[tuple[str, float]]:
    """(kind, reference CPU-seconds) for each task a job with these params runs on
    sequences of this size."""

    def step(kind: str) -> float:
        return STEPS[kind].at(size)

    if simulation(params):
        steps = [(kind, step(kind)) for kind in ("SIM_PREPARE", "SIM_ALISIM", "SIM_COLLECT")]
        return steps * len(simulation_seeds(params))
    if esr(params):
        return esr_tasks(params, size)
    if branch_site(params):
        runs = [("CODEML", BRANCH_SITE[key].at(size)) for key in branch_site_runs()]
        export = ("BRANCH_SITE_EXPORT", step("BRANCH_SITE_EXPORT"))
        return [("VALIDATE", step("VALIDATE")), *runs, export]
    if uncertainty_only(params):
        plan = params.get("uncertainty_of")
        combos = plan.get("combinations") if isinstance(plan, dict) else None
        n = combos if isinstance(combos, int) and combos > 0 else 1
        return [("UNCERTAINTY", 1.0 + n * step("UNCERTAINTY"))]
    if properties_only(params):
        # One PROPERTIES per combination, measuring the named ancestors rather than
        # every one — so this overstates a narrow scope, which the panel's own estimate
        # (phylobot.properties_run) does not.
        plan = params.get("properties_of")
        combos = plan.get("combinations") if isinstance(plan, dict) else None
        versions = 2 if as_list(params.get("gap_methods")) else 1
        if isinstance(combos, dict) and combos:
            return [
                (
                    "PROPERTIES",
                    step("PROPERTIES")
                    + neural_property_seconds(params, size.length, len(nodes), versions),
                )
                for nodes in combos.values()
            ]
        extra = neural_property_seconds(params, size.length, 1, versions)
        return [("PROPERTIES", step("PROPERTIES") + extra)] * property_targets(params)
    aligners = [a for a in as_list(params.get("aligners", "mafft")) if a != "none"]
    out = [("VALIDATE", step("VALIDATE")), ("SUMMARY", step("SUMMARY"))]
    out += [(f"ALIGN_{a.upper()}", step(f"ALIGN_{a.upper()}")) for a in aligners]
    if coding(params):
        out += [("CODON_ALIGNMENT", step("CODON_ALIGNMENT"))] * len(aligners)
    alignments = alignment_count(params)
    out += [("EXPORT_ALIGNMENT", step("EXPORT_ALIGNMENT"))] * alignments
    if not combination_count(params):
        return out
    out.append(("ANCESTOR_MAP", step("ANCESTOR_MAP")))

    fixed = [m for m in as_list(params.get("models")) if catalog.is_model(m)]
    models = [catalog.MODELFINDER] if params.get("modelfinder", True) else []
    models += fixed if params.get("fixed_models") else []
    builders = [t for t in as_list(params.get("tree_builders", "iqtree")) if t != "none"]
    planned = planned_combinations(params)  # an extension: only what it runs
    if params.get("modelfinder", True):
        ticked = [m for m in candidate_models(params) if catalog.is_model(m)]
        cost = modelfinder_seconds(ticked, "raxml" in builders, size)
        picks = alignments if planned is None else modelfinder_runs(params, planned)
        out += [("MODELFINDER", cost)] * picks
    if coding(params) and params.get("codon_modelfinder"):
        ticked = [m for m in as_list(params.get("codon_candidates")) if catalog.is_codon_model(m)]
        share = min(1.0, catalog.codon_candidate_count(ticked) / 38)
        out += [("MODELFINDER_CODON", CODON_MODELFINDER_DEFAULT.at(size) * share)] * alignments

    support = (
        1
        + UFBOOT_PER_1000 * count(params.get("bootstrap")) / 1000
        + ALRT_PER_1000 * count(params.get("alrt")) / 1000
    )
    ancestors = params.get("stop_after", "ancestors") == "ancestors"
    per_combination = [("EXPORT", step("EXPORT"))]
    if as_list(params.get("outgroup")):
        per_combination.append(("ROOT_TREE", step("ROOT_TREE")))
    for method in as_list(params.get("gap_methods")):
        if (kind := f"GAPS_{method.upper()}") in STEPS:
            per_combination.append((kind, step(kind)))
    if as_list(params.get("properties")):
        versions = 1 + sum(m in MASKED_GAP_METHODS for m in as_list(params.get("gap_methods")))
        # An unrooted tree of N sequences has N - 2 ancestral nodes.
        extra = neural_property_seconds(params, size.length, size.sequences - 2, versions)
        per_combination.append(("PROPERTIES", step("PROPERTIES") + extra))
    if planned is not None:  # amino-acid models only (extend.refusal)
        for _aligner, tree, model in planned:
            hetero = heterogeneous(model)
            scale = MIXTURE_FACTOR if catalog.is_mixture(model) else CostFactor(1, 1, 1)
            if tree == "iqtree":
                out.append(("TREE_IQTREE", TREE_IQTREE[hetero].at(size) * support * scale.search))
            else:  # RAxML-NG's search, or a provided tree: then IQ-TREE's fit
                if tree == "raxml":
                    out.append(("TREE_RAXML", TREE_RAXML[hetero].at(size)))
                out.append(("FIT_IQTREE", FIT_IQTREE[hetero].at(size) * scale.fit))
            if ancestors:
                anc = ANCESTRAL_IQTREE[hetero].at(size) * scale.ancestors
                out += [("ANCESTRAL_IQTREE", anc), *per_combination]
            else:
                out.append(("EXPORT_TREE", step("EXPORT_TREE")))
        return out
    codon_models = as_list(params.get("codon_models")) if coding(params) else []
    if coding(params) and params.get("codon_modelfinder"):
        codon_models = [catalog.CODON_MODELFINDER, *codon_models]
    for model in models + codon_models:
        hetero = heterogeneous(model)
        codon = model in codon_models
        scale = CODON_FACTOR if codon else CostFactor(1, 1, 1)
        if catalog.is_mixture(model):
            scale = MIXTURE_FACTOR
        raxml = "raxml" in builders and not codon  # RAxML-NG has no codon models
        fits = provided(params, "tree") + raxml
        sources = codon_tree_count(params) if codon else len(builders) + provided(params, "tree")
        one = []
        if codon_search(params) if codon else "iqtree" in builders:
            one.append(("TREE_IQTREE", TREE_IQTREE[hetero].at(size) * support * scale.search))
        if raxml:
            one.append(("TREE_RAXML", TREE_RAXML[hetero].at(size)))
        one += [("FIT_IQTREE", FIT_IQTREE[hetero].at(size) * scale.fit)] * fits
        if "generax" in builders and not codon:  # GeneRax: amino-acid models only
            fit = FIT_IQTREE[hetero].at(size) * scale.fit
            one += [
                ("TREE_GENERAX", TREE_IQTREE[hetero].at(size)),
                ("FIT_GENERAX", fit),
                ("TOPOLOGY_TEST", fit * (len(builders) + provided(params, "tree"))),
            ]
        if codon:  # fitted on amino-acid trees: a fit like a provided tree's
            fit = FIT_IQTREE[hetero].at(size) * scale.fit
            one += [("FIT_CODON_AA", fit)] * codon_amino_trees(params)
        if ancestors:
            anc = ("ANCESTRAL_IQTREE", ANCESTRAL_IQTREE[hetero].at(size) * scale.ancestors)
            one += [anc, *per_combination] * sources
        else:
            one += [("EXPORT_TREE", step("EXPORT_TREE"))] * sources
        out += one * alignments
    if chosen := codeml_combinations(params):
        runs = [("CODEML", CODEML[key].at(size)) for key in codeml_runs(params)]
        out += (runs + [("CODEML_EXPORT", step("CODEML_EXPORT"))]) * chosen
    if reconciles(params):
        out.append(("SPECIES_TREE", step("SPECIES_TREE")))
        per = [("RECONCILE", step("RECONCILE")), ("RECONCILE_EXPORT", step("RECONCILE_EXPORT"))]
        out += per * combination_count(params)
    return out


def reference_seconds(params: dict[str, Any], size: Size = IME2_SIZE) -> float:
    return sum(cost for _, cost in tasks(params, size))


@dataclass
class Estimate:
    units: float
    groups: list[tuple[str, str, float]]  # (key, label, units), pipeline order, non-empty
    minutes: float | None  # on this server, steps one after another; None without history


def estimate(
    params: dict[str, Any], size: Size = IME2_SIZE, seconds_per_unit: float | None = None
) -> Estimate:
    totals = dict.fromkeys(GROUP_LABELS, 0.0)
    for kind, cost in tasks(params, size):
        totals[group(kind)] += cost / SECONDS_PER_UNIT
    units = sum(totals.values())
    return Estimate(
        units=units,
        groups=[(key, GROUP_LABELS[key], value) for key, value in totals.items() if value],
        minutes=None if seconds_per_unit is None else units * seconds_per_unit / 60,
    )


HISTORY = 50  # most recent finished jobs to calibrate from
MIN_HISTORY = 3


def seconds_per_unit(history: int = HISTORY) -> float | None:
    """This server's CPU-seconds per work unit (see calibration); None until MIN_HISTORY
    successful jobs exist."""
    value, _ = calibration(history)
    return value


def calibration(history: int = HISTORY) -> tuple[float | None, int]:
    """(CPU-seconds per work unit on this server, jobs it is learned from): the median
    ratio of recorded work to the estimate over recent successful jobs, from every user
    (only the number is shown). Only jobs that stored their estimate count once there
    are MIN_HISTORY of them; until then older jobs are estimated now, at the Ime2
    example's size (their real size is not recorded)."""
    recent = (
        Job.objects.filter(status=Job.Status.SUCCEEDED, cpu_seconds__gt=0)
        .order_by("-finished_at")
        .values_list("params", "cpu_seconds", "work_units")[:history]
    )
    stored: list[float] = []
    guessed: list[float] = []
    for params, used, units in recent:
        if units is None:
            guess = reference_seconds(params or {}) / SECONDS_PER_UNIT
            if used and guess > 0:
                guessed.append(used / guess)
        elif used and units > 0:
            stored.append(used / units)
    ratios = stored if len(stored) >= MIN_HISTORY else stored + guessed
    return (median(ratios) if len(ratios) >= MIN_HISTORY else None), len(ratios)


@dataclass
class CostRow:
    label: str
    variant: str  # model class or "", e.g. "+G4/+R4"
    cost: Cost
    per: str  # what one cost is for: "task", "matrix", "model"


def reference_costs() -> list[CostRow]:
    """Every cost the estimate uses, for the Work estimate page, in pipeline order."""
    from phylobot.stages import LABELS

    rows = [CostRow(LABELS.get(kind, kind), "", cost, "task") for kind, cost in STEPS.items()
            if kind.startswith("ALIGN_")]  # fmt: skip
    rows += [
        CostRow("ModelFinder", "every rate and frequency type (shortcuts)",
                MODELFINDER_PER_MATRIX, "matrix"),
        CostRow("ModelFinder", "some types unticked (no shortcuts)",
                MODELFINDER_PER_MODEL, "fitted model"),
    ]  # fmt: skip
    by_model = {"TREE_IQTREE": TREE_IQTREE, "TREE_RAXML": TREE_RAXML,
                "FIT_IQTREE": FIT_IQTREE, "ANCESTRAL_IQTREE": ANCESTRAL_IQTREE}  # fmt: skip
    for kind, costs in by_model.items():
        rows += [CostRow(LABELS[kind], "plain or +I" if not hetero else "+G4 or +R4",
                         costs[hetero], "task") for hetero in (False, True)]  # fmt: skip
    rows += [CostRow(LABELS.get(kind, kind), "", cost, "task") for kind, cost in STEPS.items()
             if not kind.startswith("ALIGN_")]  # fmt: skip
    return rows


# Example jobs for the Work estimate page: the New Job page's defaults (ModelFinder and
# one IQ-TREE search per alignment, SH-aLRT, Fitch gap placement) at several sizes.
EXAMPLE_PARAMS: dict[str, Any] = {
    "aligners": ["mafft"], "tree_builders": ["iqtree"], "modelfinder": True, "alrt": 1000,
    "gap_methods": ["fitch"], "outgroup": ["x"],
}  # fmt: skip
EXAMPLE_SIZES = [Size(7, 600), Size(20, 400), Size(60, 600), Size(200, 500), Size(500, 600)]


def size_of(fasta: bytes, codons: bool = False) -> Size | None:
    """Sequence count and mean ungapped length of a FASTA (for coding DNA in codons, the
    length of the translations every cost is measured on); None if it can't be read."""
    from phylobot.pipeline.validate import InputError, parse_fasta

    try:
        records = parse_fasta(fasta.decode("utf-8"), "sequences")
    except InputError, UnicodeDecodeError:
        return None
    lengths = [sum(c.isalpha() for c in sequence) for _, _, sequence in records]
    if not lengths:
        return None
    mean = sum(lengths) / len(lengths)
    return Size(len(lengths), mean / 3 if codons else mean)
