"""The job page's "All files": every result file in compact tables, one per group (the
job, each alignment, each combination, logs), each file with a one-line description.
A job with more files than INLINE_FILES shows folders instead, each with its own page
(`folder_listing`), down to a level small enough for the tables.

DESCRIPTIONS is the one place files are described, by name within their group; a test
checks that every file the pipeline writes has one.
"""

import re
from dataclasses import dataclass
from typing import Any

from phylobot import stages

# Job-level files (results/<name>).
JOB = {
    "summary.json": "What ran: every combination, its labels and fit (this page reads it).",
    "model_comparison.tsv": "Fit statistics of every combination, ranked by BIC within each "
    "alignment and data type (amino acid, codon).",
    "validation.json": "Input checks: number, type and lengths of the sequences, warnings.",
    "versions.tsv": "Every program the job ran, and its version.",
    "ancestor_map.tsv": "Every combination's ancestors by clade: rows with the same clade are "
    "the same ancestor (the same sequences below it on the outgroup-rooted tree); with the "
    "same parts, the same subtrees meet there too.",
    "sources.tsv": "Your uploaded alignments and trees, with their labels.",
}
# Each aligner's files (results/<aligner>/<name>).
BRANCH_SITE = {
    "test.tsv": "The likelihood-ratio test: 2ΔlnL, degrees of freedom, p (χ², 1 df) and p "
    "under the 50:50 mixture of 0 and χ²₁.",
    "models.tsv": "Model A and its null (ω₂ = 1): log-likelihood, parameters, kappa, site "
    "class proportions and dN/dS, every start's log-likelihood.",
    "beb.tsv": "Model A's Bayes empirical Bayes sites: posterior probability of positive "
    "selection on the foreground (codeml lists those above 0.5).",
    "sites.tsv": "Every site: Model A's BEB posteriors of the four site classes (0, 1, 2a, "
    "2b) and of positive selection on the foreground (2a + 2b).",
    "foreground.nwk": "The tree codeml analysed, unrooted, with #1 on each foreground branch "
    "(original sequence names).",
}
BRANCH_SITE_CODEML = {
    "codeml/README.txt": "How to rerun codeml on these files.",
    "codeml/aln.fasta": "The codon alignment codeml read (PhyloBot's sequence IDs).",
    "codeml/tree.nwk": "The tree codeml read, #1 on the foreground (PhyloBot's IDs).",
    "codeml/names.tsv": "PhyloBot's sequence IDs and the original names.",
}


def describe_branch_site(name: str) -> str:
    if name in BRANCH_SITE or name in BRANCH_SITE_CODEML:
        return BRANCH_SITE.get(name) or BRANCH_SITE_CODEML[name]
    run, _, kind = name.removeprefix("codeml/").rpartition(".")
    model = "the null (ω₂ = 1)" if run.startswith("A_null") else "Model A"
    start = run.rpartition(".start")[2]
    what = {"ctl": "control file", "mlc": "main output",
            "rst": "per-site class probabilities (rst)"}.get(kind, "file")  # fmt: skip
    return f"codeml's {what}: {model}, start {start}." if name.startswith("codeml/") else ""


ALIGNMENT = {
    "alignment.fasta": "The alignment (for coding DNA, of the translated proteins).",
    "alignment.cds.fasta": "The codon alignment (coding DNA): each residue's own codon, gaps "
    "in whole codons.",
}
# A combination's files (results/<aligner>/<tree>/<model>/<name>).
COMBINATION = {
    "alignment.fasta": "The alignment this combination used.",
    "alignment.cds.fasta": "Its codon alignment (codon models run on it).",
    "tree.nwk": "The tree: branch lengths in substitutions per site, ancestral nodes named "
    "Node1, Node2, ... as in the ancestor files.",
    "tree.rooted.nwk": "The tree rooted on the outgroup.",
    "rooting.json": "How the tree was rooted: the outgroup, and whether it forms one group.",
    "supports.tsv": "Branch support per node (SH-aLRT, UFBoot).",
    "site_rates.tsv": "Each alignment column's rate (IQ-TREE's posterior mean, with rate "
    "heterogeneity only) and most probable rate category.",
    "site_lh.tsv": "Each alignment column's log-likelihood on the final tree.",
    "model.txt": "The substitution model used (ModelFinder's pick, or the model you chose).",
    "fit.tsv": "Fit of the final tree: log-likelihood, parameters, AIC, AICc, BIC.",
    "modelfinder.tsv": "ModelFinder's ranking of every candidate model (on its starting tree).",
    "ancestors.fasta": "Ancestral sequences: the most probable residue at every alignment "
    "column (no gaps).",
    "ancestors_pp.tsv": "Per node and site: the most probable state and its posterior probability.",
    "ancestors_posteriors.tsv": "Per node and site: the posterior probability of every amino acid.",
    "ancestors.cds.fasta": "Coding ancestors: at each site, the best codon for the best "
    "amino acid.",
    "ancestors_codon_pp.tsv": "Per node and site: that codon and its posterior probability.",
    "ancestors_codon_posteriors.tsv": "Per node and site: every codon with posterior "
    "probability 0.01 or more.",
    "gaps-ml.tsv": "Per node and site: the probability of a gap (ML gap placement).",
    "grasp_model.txt": "The model GRASP used (the combination's matrix if GRASP has it, else LG).",
    "codeml/codeml.tsv": "PAML codeml's site models: log-likelihood, parameters, dN/dS "
    "estimates, and each starting omega's log-likelihood.",
    "codeml/selection.tsv": "Tests for positive selection (M1a vs M2a, M7 vs M8): 2ΔlnL, df, "
    "p-value.",
    "codeml/beb.tsv": "Sites under positive selection: posterior P(dN/dS > 1) by Bayes "
    "empirical Bayes (M2a, M8).",
    "properties/draws.parquet": "Every property's value per posterior draw, for every "
    "ancestor and gap version (Parquet: one row per version, node and draw).",
    "properties/summaries.json": "Each ancestor's property summaries: mean, sd, interval, "
    "f(ML) and its shift, per gap version.",
    "properties/index.json": "Which properties were measured, on which gap versions, with "
    "how many draws.",
}
POSTERIOR_INDEX = re.compile(r"(.+)\.index\.tsv$")

GAPS = {
    "parsimony": "columns judged absent masked as gaps by parsimony (shared presence)",
    "fitch": "columns judged absent masked as gaps by Fitch parsimony",
    "ml": "columns judged absent masked as gaps by maximum likelihood (P(gap) > 0.5)",
    "prank": "PRANK's own reconstruction (its residues and gaps), on our nodes",
    "grasp": "GRASP's own reconstruction (indel-aware), on our nodes",
}
LOGS = {
    "nextflow.log": "Nextflow's own log (for debugging the pipeline).",
    "pipeline.log": "The pipeline's console output.",
    "params.json": "The parameters the pipeline ran with.",
}
REPORT = {
    "report.html": "Nextflow's execution report: time and memory per step.",
    "trace.tsv": "One row per task: status, start and end, CPU and memory.",
}
GAP_FILE = re.compile(r"^ancestors\.gaps-(\w+)(\.cds)?\.fasta$")
CODEML_FILE = re.compile(r"^codeml/(M\w+|FR)/(.+)$")
TASK_LOG = re.compile(r"^tasks/(.+)\.(log|err)$")


def describe_combination_file(name: str) -> str:
    """A combination's file (path inside the combination), described."""
    if name in COMBINATION:
        return COMBINATION[name]
    if match := POSTERIOR_INDEX.match(name):
        whole = match.group(1).rsplit("/", 1)[-1] + ".tsv"
        return (
            f"Where each node's rows are in {whole} (byte offset and length), so the "
            "viewers read one ancestor instead of the whole table."
        )
    if match := GAP_FILE.match(name):
        method, coding = match.group(1), match.group(2)
        text = GAPS.get(method, f"gap placement by {method}")
        what = "Coding ancestors" if coding else "Ancestors"
        if method in ("prank", "grasp"):
            return f"{text[0].upper()}{text[1:]}."
        return f"{what} with {text}."
    if (match := CODEML_FILE.match(name)) and match.group(2) in COMBINATION:
        text = COMBINATION[match.group(2)]
        return f"codeml {match.group(1)}: {text[0].lower()}{text[1:]}"
    return ""


def describe_log(name: str) -> str:
    """A file in logs/ (path inside logs/), described."""
    if name in LOGS:
        return LOGS[name]
    if match := TASK_LOG.match(name):
        task = match.group(1)
        # Step names contain underscores (ALIGN_MAFFT, CODEML_EXPORT): the longest known one.
        steps = sorted(stages.LABELS, key=len, reverse=True)
        step = next((s for s in steps if task == s or task.startswith(s + "_")), None)
        if step is None:
            return ""
        tag = task[len(step) + 1 :]
        stream = "output" if match.group(2) == "log" else "errors (stderr)"
        where = f" ({tag})" if tag and tag != "1" else ""
        return f"Step “{stages.LABELS[step]}”{where}: its {stream}."
    return ""


# A simulation's truth bundle per replicate (phylobot.simulate).
SIMULATION = {
    "input.fasta": "The simulated sequences, unaligned: what a job on this replicate is given.",
    "true_alignment.fasta": "The true alignment of the simulated sequences.",
    "true_ancestors.fasta": "Every ancestor's true sequence, aligned to the true alignment.",
    "true_tree.nwk": "The tree the sequences evolved on, ancestors named as in "
    "true_ancestors.fasta.",
    "names.tsv": "Each simulated sequence's name and the random tree's own label for it.",
    "species_tree.nwk": "The species tree the gene family grew on, internal nodes S1, S2, ...",
    "species.tsv": "Each simulated gene and the species it belongs to.",
    "true_events.tsv": "Every duplication, loss and speciation in the family's history, with "
    "the ancestor it became.",
    "family.json": "The family's counts: genes, duplications (and how many are visible), losses.",
    "provenance.json": "How it was made: AliSim's version and exact commands, the seed, and "
    "checksums of every file.",
}


# An ESR job's results (results/esr/<name>; phylobot.pipeline.esr).
ESR = {
    "esr_tips.tsv": "Each held-out sequence at each level: its pendant branch length, the "
    "sites scored, accuracy (the most probable residue is the real one), mean posterior of "
    "the most probable residue, mean log score.",
    "esr_sites.tsv": "Every scored site of every held-out sequence: the real residue, the "
    "most probable one, their posteriors and the whole posterior over the 20 amino acids.",
    "esr_summary.json": "The figures this page shows: per level, accuracy and mean posterior "
    "over every site, the reliability bins, and the fitted model.",
    "tree.nwk": "The tree ESR ran on, with its branch lengths (fitted on the alignment, for "
    "an uploaded tree).",
}


@dataclass
class FileTable:
    title: str
    rows: list[tuple[str, str, str]]  # (result path, name shown, description)


def file_tables(keys: list[str], summary: dict[str, Any] | None = None) -> list[FileTable]:
    """The result files in tables: the job's, each alignment's, each combination's (with
    codeml's under it), the logs, and Nextflow's report."""
    summary = summary or {}
    combos = {c["key"]: c.get("label", c["key"]) for c in summary.get("combinations", [])}
    aligners = {a["key"]: a.get("label", a["key"]) for a in summary.get("alignments", [])}
    groups: dict[str, FileTable] = {}

    def add(group: str, title: str, key: str, name: str, description: str) -> None:
        groups.setdefault(group, FileTable(title, [])).rows.append((key, name, description))

    for key in sorted(keys):
        parts = key.split("/")
        if parts[0] == "logs":
            name = "/".join(parts[1:])
            add("logs", "Logs", key, name, describe_log(name))
        elif parts[0] == "pipeline_info":
            add("report", "Pipeline report", key, parts[-1], REPORT.get(parts[-1], ""))
        elif len(parts) == 1:
            add("job", "Job", key, key, JOB.get(key, ""))
        elif parts[0] == "replicates" and len(parts) == 3:
            add(f"r:{int(parts[1]):012d}", f"Replicate (seed {parts[1]})", key, parts[2],
                SIMULATION.get(parts[2], ""))  # fmt: skip
        elif parts[0] == "esr" and len(parts) == 2:
            add("esr", "ESR", key, parts[1], ESR.get(parts[1], ""))
        elif parts[0] == "branch_site":
            name = "/".join(parts[1:])
            add("branch_site", "Branch-site test", key, name, describe_branch_site(name))
        elif len(parts) == 2:
            title = f"Alignment: {aligners.get(parts[0], parts[0])}"
            add(f"a:{parts[0]}", title, key, parts[1], ALIGNMENT.get(parts[1], ""))
        else:
            combo, name = "/".join(parts[:3]), "/".join(parts[3:])
            add(f"c:{combo}", combos.get(combo, combo), key, name, describe_combination_file(name))
    order = ["job", "esr", "branch_site", *sorted(g for g in groups if g.startswith("r:"))]
    order += sorted(g for g in groups if g.startswith("a:"))
    order += sorted(g for g in groups if g.startswith("c:")) + ["logs", "report"]
    return [groups[g] for g in order if g in groups]


# More files than this and a page lists folders, not files (the 285-combination GFP job
# has 524,069; a combination with properties on 107 nodes, about 1,800).
INLINE_FILES = 400


@dataclass
class Folder:
    path: str  # relative to results/, ending in "/"
    name: str  # its last part
    files: int
    size: int
    combinations: int  # combinations at or below it


@dataclass
class FolderListing:
    folder: str  # "" for results/ itself
    total_files: int
    total_size: int
    folders: list[Folder]  # empty when the files are shown as tables
    tables: list[FileTable]  # the files directly in `folder`, or every file under it


def folder_listing(
    sizes: dict[str, int], folder: str = "", summary: dict[str, Any] | None = None
) -> FolderListing:
    """The files under `folder` (relative to results/, "" or ending in "/"): as tables if
    there are at most INLINE_FILES, else its subfolders with counts and sizes and the
    files directly in it as tables."""
    summary = summary or {}
    under = {k: n for k, n in sizes.items() if k.startswith(folder)}
    listing = FolderListing(folder, len(under), sum(under.values()), [], [])
    if len(under) <= INLINE_FILES:
        listing.tables = file_tables(list(under), summary)
        return listing
    combos = {c["key"]: c.get("label", c["key"]) for c in summary.get("combinations", [])}
    here: list[str] = []
    children: dict[str, Folder] = {}
    for key, size in under.items():
        rest = key[len(folder) :]
        if "/" not in rest:
            here.append(key)
            continue
        path = folder + rest.split("/", 1)[0] + "/"
        child = children.get(path)
        if child is None:
            count = sum(1 for c in combos if (c + "/").startswith(path))
            child = children[path] = Folder(path, path[len(folder) : -1], 0, 0, count)
        child.files += 1
        child.size += size
    listing.folders = [children[p] for p in sorted(children)]
    listing.tables = file_tables(here, summary)
    return listing


def breadcrumbs(folder: str) -> list[tuple[str, str]]:
    """(name, path) for each folder above and including `folder`, results/ first."""
    crumbs = [("results", "")]
    parts = [p for p in folder.split("/") if p]
    for n, part in enumerate(parts, 1):
        crumbs.append((part, "/".join(parts[:n]) + "/"))
    return crumbs
