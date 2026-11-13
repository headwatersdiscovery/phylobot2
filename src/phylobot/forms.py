import json
import re
import tempfile
import uuid
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from django import forms

from phylobot import (
    alignments,
    branch_site,
    catalog,
    codons,
    extend,
    inputs,
    machine,
    site_settings,
    species_web,
    storage,
    workload,
)
from phylobot import codeml as paml  # (SubmitForm has a field named codeml)
from phylobot import properties as properties_mod  # (SubmitForm has a field named properties)
from phylobot.homologs import filters, remote, sources
from phylobot.models import Alignment, Job, Limits
from phylobot.pipeline.validate import (
    AMBIGUITY_WARNING,
    MIN_SEQUENCES,
    InputError,
    detect_seqtype,
    validate,
)
from phylobot.species import table as species_table
from phylobot.species.polytomies import ARBITRARY, GUIDED

# Tools with a pipeline module; keep in step with parseTools in pipeline/main.nf (tested).
ALIGNERS = [
    ("mafft", "MAFFT"),
    ("muscle", "MUSCLE 5"),
    ("clustalo", "Clustal Omega"),
    ("prank", "PRANK (phylogeny-aware; slow)"),
    ("prank_codon", "PRANK with its codon model (aligns the DNA itself; slower)"),
]
CODON_ALIGNERS = {"prank_codon"}  # coding DNA only
SEQTYPES = [("auto", "Detect"), ("protein", "Protein"), ("codon", "Coding DNA")]
# How coding DNA is analysed (owner requirement: codon models, amino-acid models on the
# translation, or both, side by side; both by default, as first decided).
ANALYSES = [
    ("both", "Both, side by side"),
    ("codon", "Codon models"),
    ("amino", "Amino-acid models (on the translation)"),
]
AVAILABLE_ANALYSES = {"amino", "codon", "both"}


def analysis_flags(seqtype: str, analyses: str) -> tuple[bool, bool]:
    """(amino-acid models run, codon models run) for a job: protein jobs use amino-acid
    models only; coding DNA uses what the user chose."""
    if seqtype != "codon":
        return True, False
    choice = analyses if analyses in AVAILABLE_ANALYSES else "both"
    return choice in ("amino", "both"), choice in ("codon", "both")


# Which amino-acid trees codon models are fitted on (--codon_amino_trees; owner decision
# 2026-10-02: the best by BIC by default, every one on request).
CODON_AMINO_TREES = [
    ("best", "The best amino-acid tree (lowest BIC), per tree source"),
    ("all", "Every amino-acid tree (each amino-acid model's)"),
]


def codon_tree_params(search: bool, amino: bool, which: str | None) -> dict[str, Any]:
    """Job.params for the codon combinations' trees: their own IQ-TREE search and/or fits
    on the amino-acid trees ('none': only fits to the user's own trees)."""
    chosen = [key for key, on in (("search", search), ("amino", amino)) if on]
    params: dict[str, Any] = {"codon_trees": chosen or ["none"]}
    if amino:
        params["codon_amino_trees"] = which if which in dict(CODON_AMINO_TREES) else "best"
    return params


def code_label(code: codons.GeneticCode) -> str:
    return f"{code.ncbi}. {code.name} ({' and '.join(code.tools)})"


# Tree searches; keep in step with parseTools in pipeline/main.nf (tested). Labels: the
# New Job page's checkboxes.
TREE_BUILDERS = [
    ("iqtree", "IQ-TREE 3 (maximum-likelihood tree search)"),
    ("raxml", "RAxML-NG 2 (maximum-likelihood tree search; a second, independent heuristic)"),
]
# GeneRax's reconciled tree: a tree source of its own, offered in the
# reconciliation block rather than beside IQ-TREE and RAxML-NG (owner, 2026-10-17).
RECONCILED_TREE = "generax"
POLYTOMIES = [
    (GUIDED, "From the job's gene trees (recommended)"),
    (ARBITRARY, "At random, from the job's seed"),
]
# Gap placement in ancestors (src/phylobot/gaps.py); keep in step with main.nf (tested).
GAP_METHODS = [
    ("parsimony", "Parsimony (shared presence) · fast"),
    ("fitch", "Parsimony (Fitch) · fast; uses the outgroup root if given"),
    ("ml", "Maximum likelihood (with gap probabilities) · one extra IQ-TREE run per combination"),
    ("prank", "PRANK's own ancestors · alternative reconstruction · needs an outgroup"),
    (
        "grasp",
        "GRASP's own ancestors (indel-aware) · alternative reconstruction · needs an outgroup",
    ),
]
NEED_OUTGROUP = {"prank", "grasp"}  # keep in step with phylobot.gaps.NEED_OUTGROUP (tested)
MAX_PROVIDED = 5  # your own alignments, and trees, each; as max_provided in nextflow.config
STAGES = ("alignments", "trees", "ancestors")  # stop_after values, in pipeline order

# The New Job defaults (owner decision): enough alignments,
# models and gap versions that an ancestral library has something to compare. PRANK is
# left out for its run time; ModelFinder's pick may duplicate a fixed model (the results
# table merges the two rows).
DEFAULT_ALIGNERS = ["mafft", "muscle", "clustalo"]
DEFAULT_MODELS = ["LG+G4", "WAG+G4"]
DEFAULT_GAP_METHODS = ["parsimony", "fitch"]


class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class MultipleFileField(forms.FileField):
    """Several files from one input (Django's documented pattern)."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        kwargs.setdefault("widget", MultipleFileInput())
        super().__init__(*args, **kwargs)

    def clean(self, data: Any, initial: Any = None) -> list[Any]:
        single = super().clean
        if isinstance(data, (list, tuple)):
            return [single(item, initial) for item in data if item]
        return [single(data, initial)] if data else []


def chosen_alignments(owner: Any, keys: Sequence[str]) -> list[Alignment]:
    """The owner's alignments among `keys` (IDs), in the order given."""
    wanted = [k for k in dict.fromkeys(keys) if is_uuid(k)][: MAX_PROVIDED * 4]
    if owner is None or not wanted:
        return []
    found = {str(a.id): a for a in Alignment.objects.for_user(owner).select_related("job")
             .filter(id__in=wanted)}  # fmt: skip
    return [found[k] for k in wanted if k in found]


def is_uuid(text: str) -> bool:
    try:
        uuid.UUID(text)
    except ValueError:
        return False
    return True


STORE = {"hidden": True, "data-file-store": ""}  # a file input the page fills (no picker)


def used_files(files: list[Any], data: Any, kind: str) -> list[Any]:
    """The uploads a file list (`_file_list.html`) leaves ticked: a file the page has
    listed (`<kind>_listed`) but not ticked (`<kind>_use`) is left out; a file not yet
    listed (just chosen) is used. One file per name (the page keeps names unique)."""
    listed, use = set(data.getlist(f"{kind}_listed")), set(data.getlist(f"{kind}_use"))
    out, seen = [], set()
    for upload in files:
        name = upload.name or ""
        if name in seen or (name in listed and name not in use):
            continue
        seen.add(name)
        out.append(upload)
    return out


def label_for(filename: str) -> str:
    """A default label for an uploaded file: its name without the extension. VALIDATE
    checks labels (length, characters) like any user text."""
    return Path(filename).stem[:100] or "uploaded"


class SubmitForm(forms.Form):
    name = forms.CharField(max_length=100, required=False, label="Job name (optional)")
    sequences_text = forms.CharField(
        required=False,
        strip=False,  # keep the FASTA exactly as pasted
        label="Sequences (FASTA: protein or coding DNA)",
        widget=forms.Textarea(attrs={"rows": 12, "spellcheck": "false"}),
    )
    # The file inputs hold every file chosen so far (submit.html's script adds to them).
    sequences_file = MultipleFileField(
        required=False, label="…or upload FASTA files:", widget=MultipleFileInput(STORE),
    )  # fmt: skip
    # A name in more than one sequence file with different sequences (phylobot.inputs).
    duplicate_names = forms.ChoiceField(
        choices=[(inputs.RENAME, "Rename"), (inputs.FIRST, "First")], required=False,
    )  # fmt: skip
    seqtype = forms.ChoiceField(
        choices=SEQTYPES, initial="auto", required=False, widget=forms.RadioSelect,
        label="Sequence type",
    )  # fmt: skip
    genetic_code = forms.TypedChoiceField(
        choices=[(c.ncbi, code_label(c)) for c in codons.GENETIC_CODES],
        coerce=int,
        initial=1,
        required=False,
        label="Genetic code",
    )
    analyses = forms.ChoiceField(
        choices=ANALYSES, initial="both", required=False, widget=forms.RadioSelect,
        label="Analyse coding DNA with",
    )  # fmt: skip
    # Codon models (coding DNA), mirroring the protein options: ModelFinder among the
    # ticked candidates (all ticked: IQ-TREE's default codon set) and/or a full search
    # under each ticked codon model (none ticked).
    codon_modelfinder = forms.BooleanField(required=False, initial=True)
    codon_candidates = forms.MultipleChoiceField(
        choices=[(m, m) for m in catalog.all_codon_models()],
        required=False,
        initial=catalog.all_codon_models(),
    )
    codon_fixed = forms.BooleanField(required=False, initial=False)
    # The codon combinations' trees: a search under each codon model (slow) and/or a fit
    # of each codon model on the amino-acid trees (cheap; no codon search).
    codon_search = forms.BooleanField(required=False, initial=True)
    codon_on_amino = forms.BooleanField(required=False, initial=True)  # owner decision 2026-10-02
    codon_amino_trees = forms.ChoiceField(
        choices=CODON_AMINO_TREES, initial="best", required=False, widget=forms.RadioSelect
    )
    # PAML codeml (coding DNA, with codon models; owner decisions 2026-10-01): M0, M3,
    # M7 and M8 ticked, M1a and M2a offered for the positive-selection test; each
    # alignment's best codon combination by default.
    codeml = forms.BooleanField(required=False, initial=False, label="PAML codeml")
    codeml_models = forms.MultipleChoiceField(
        choices=[(m.key, m.label) for m in paml.SITE_MODELS],
        initial=paml.DEFAULT_MODELS,
        required=False,
        widget=forms.CheckboxSelectMultiple,
    )
    codeml_trees = forms.ChoiceField(
        choices=[
            ("best", "Each alignment's best codon combination (by BIC)"),
            ("all", "Every codon combination (slower)"),
        ],
        initial="best",
        required=False,
        widget=forms.RadioSelect,
    )
    codon_models = forms.MultipleChoiceField(
        choices=[(m, m) for m in catalog.all_codon_models()], required=False, initial=[]
    )
    aligners = forms.MultipleChoiceField(
        choices=ALIGNERS,
        initial=DEFAULT_ALIGNERS,
        required=False,
        widget=forms.CheckboxSelectMultiple,
    )
    alignment_files = MultipleFileField(
        required=False, label="…upload alignment (FASTA or PHYLIP):",
        widget=MultipleFileInput(STORE),
    )  # fmt: skip
    # Your own alignments with different sequences: use the ones they all share.
    shared_names = forms.BooleanField(required=False)
    # A name used twice in an uploaded alignment with different rows (phylobot.inputs).
    alignment_duplicates = forms.ChoiceField(
        choices=[(inputs.RENAME, "Rename"), (inputs.FIRST, "First")], required=False,
    )  # fmt: skip
    # From My alignments (phylobot.alignments): the choices are the posted (or initial)
    # IDs that name the user's own alignments; others are refused as invalid choices.
    saved_alignments = forms.MultipleChoiceField(required=False)
    run_trees = forms.BooleanField(required=False, initial=True, label="Infer trees")
    tree_builders = forms.MultipleChoiceField(
        choices=TREE_BUILDERS,
        initial=["iqtree"],
        required=False,
        widget=forms.CheckboxSelectMultiple,
    )
    tree_files = MultipleFileField(required=False, label="My own tree(s)")
    run_ancestors = forms.BooleanField(
        required=False, initial=True, label="Reconstruct ancestral sequences"
    )
    modelfinder = forms.BooleanField(required=False, initial=True)
    fixed_models = forms.BooleanField(required=False, initial=True)
    # Two grids: ModelFinder's candidates (all ticked: IQ-TREE's default set) and the
    # models that each get a full tree search (DEFAULT_MODELS ticked).
    candidates = forms.MultipleChoiceField(
        choices=[(m, m) for m in catalog.all_models()],
        required=False,
        initial=catalog.all_models(),
    )
    models = forms.MultipleChoiceField(
        choices=[(m, m) for m in catalog.all_fixed_models()],
        required=False,
        initial=DEFAULT_MODELS,
    )
    outgroup = forms.MultipleChoiceField(required=False, widget=forms.CheckboxSelectMultiple)
    # Reconciliation with the species tree (owner decisions): on
    # by default, any user can untick it; the species go to the Open Tree of Life only
    # once its consent box is ticked (never by default), and an uploaded species tree
    # needs no consent. The table's cells (species_for/_query/_taxid) are read by
    # species_web.edits, not as fields, since there is one per sequence.
    reconcile = forms.BooleanField(
        required=False, initial=True, label="Reconcile with the species tree"
    )
    species_consent = forms.BooleanField(
        required=False,
        initial=False,
        label="Send the species names to the Open Tree of Life",
        help_text=(
            "Agree to send the species named in your sequences' headers (and any taxids), "
            "not the sequences or their names, to the Open Tree of Life "
            "(api.opentreeoflife.org), which returns the species tree. Its terms apply. "
            "Not needed if you upload your own species tree."
        ),
    )
    accession_consent = forms.BooleanField(
        required=False,
        initial=False,
        label="Look up accessions at NCBI and UniProt",
        help_text=(
            "Agree to send the accession numbers in your sequence names (e.g. XP_012345.1, "
            "P32581), not the sequences, to NCBI (eutils.ncbi.nlm.nih.gov) and UniProt "
            "(rest.uniprot.org), which return each one's species. Their terms apply."
        ),
    )
    # A species list: pasted, or uploaded (e.g. the table's own download, edited).
    species_paste = forms.CharField(
        required=False,
        max_length=species_web.MAX_LIST_BYTES,
        widget=forms.Textarea(attrs={"rows": 4, "spellcheck": "false"}),
        label="Paste a list",
    )
    species_table_file = forms.FileField(required=False, label="…or upload one")
    reconciled_tree = forms.BooleanField(
        required=False, initial=False, label="Also search a reconciled tree (GeneRax)"
    )
    polytomies = forms.ChoiceField(
        choices=POLYTOMIES, initial=GUIDED, required=False, widget=forms.RadioSelect,
        label="Resolve the species tree's unresolved nodes",
    )  # fmt: skip
    species_tree_file = forms.FileField(required=False, label="My own species tree")
    # Taxids carried from a homolog search ("name<TAB>taxid" lines; species_web.hints).
    species_hints = forms.CharField(required=False, widget=forms.HiddenInput, max_length=2_000_000)
    gap_methods = forms.MultipleChoiceField(
        choices=GAP_METHODS,
        required=False,
        # Both parsimony rules (owner decision 2026-10-12), with or without an outgroup.
        # The pipeline's own default for an unset list is still gaps.default_methods.
        initial=DEFAULT_GAP_METHODS,
        widget=forms.CheckboxSelectMultiple,
    )
    # Tier-0 property posteriors. Ticked by default (owner decision) because it costs
    # milliseconds per ancestor; the pipeline's own default is none, so the form passes the
    # keys explicitly and golden runs stay untouched.
    properties = forms.BooleanField(
        required=False, initial=True, label="Measure properties over each ancestor's posterior"
    )
    # The neural predictors, opt-in (owner, 2026-10-21): they cost minutes to hours where the
    # three above cost seconds, so "Measure properties" alone never runs them. Measured
    # alongside the defaults, over the same draws, when ticked.
    extra_properties = forms.MultipleChoiceField(
        choices=[(p.key, p.label) for p in properties_mod.REGISTRY.values() if not p.default],
        required=False,
        initial=[],
        widget=forms.CheckboxSelectMultiple,
    )
    # How many sequences each property is measured on. More draws tighten the interval
    # and cost proportionally; blank means the default, as the replicate counts do.
    property_draws = forms.IntegerField(
        min_value=properties_mod.MIN_DRAWS,
        max_value=properties_mod.MAX_DRAWS,
        initial=properties_mod.DEFAULT_DRAWS,
        required=False,
    )
    # Branch support (owner decision 2026-09-30: SH-aLRT on, UFBoot off). Each measure has
    # its own replicate count; IQ-TREE refuses fewer than 1000. Counts are not required
    # (blank = 1000): the page disables a count while its box is unticked.
    sh_alrt = forms.BooleanField(required=False, initial=True, label="SH-aLRT")
    alrt = forms.IntegerField(min_value=1000, max_value=10000, initial=1000, required=False)
    ufboot = forms.BooleanField(required=False, initial=False, label="UFBoot")
    bootstrap = forms.IntegerField(min_value=1000, max_value=10000, initial=1000, required=False)
    seed = forms.IntegerField(min_value=1, max_value=2**31 - 1, initial=1)
    # "Duplicate and run again": the job this one copies (the view checks it is the
    # user's, and reuses its uploaded alignments and trees unless new ones are uploaded).
    from_job = forms.UUIDField(required=False, widget=forms.HiddenInput)
    # "Extend" (phylobot.extend): the finished job this one builds on (the view checks it
    # is the user's and can be extended; its sequences and alignments come from it).
    extend_job = forms.UUIDField(required=False, widget=forms.HiddenInput)
    # "Use in a new job" on a simulation's replicate: "<simulation job id>:<seed>". The
    # new job is marked as run on simulated data (phylobot.simulations).
    simulated_from = forms.CharField(required=False, max_length=80, widget=forms.HiddenInput)
    # The fingerprint of the ambiguous codons the user confirmed in the New Job page's
    # dialog (validate.ambiguity_report); coding DNA with ambiguous codons needs it.
    ambiguity_ack = forms.CharField(required=False, max_length=32, widget=forms.HiddenInput)

    def __init__(
        self,
        *args: Any,
        limits: Limits,
        reuse: dict[str, list[tuple[str, bytes]]] | None = None,
        owner: Any = None,
        extension: extend.Inheritance | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.limits = limits
        self.extension = extension  # the job this one extends, and what it inherits
        self.owner = owner
        self.simulation: tuple[Job, int] | None = None  # (simulation, replicate seed)
        self.ambiguity: dict[str, Any] | None = None  # validation.json's ambiguous_codons
        self.left_out: list[str] = []  # sequences not in every alignment (shared_names)
        if self.is_bound:
            asked = self.data.getlist("saved_alignments")  # type: ignore[attr-defined]
        else:
            asked = list(self.initial.get("saved_alignments") or [])
        self.saved = {str(a.id): a for a in chosen_alignments(owner, asked)}
        choice = self.fields["saved_alignments"]
        assert isinstance(choice, forms.MultipleChoiceField)
        choice.choices = [(key, a.title) for key, a in self.saved.items()]
        self.reuse = reuse or {}  # "alignment"/"tree": (label, content) of the copied job
        self.fasta = b""
        self.sequence_type = "protein"  # the sequence type in effect, set by clean()
        self.alignments: list[tuple[str, bytes]] = []  # (label, content) of uploads
        # Per alignment, an inherited aligner's key or extend.PROVIDED (extensions only).
        self.alignment_keys: list[str] = []
        self.additions = extend.Additions()  # an extension's, to every alignment
        self.trees: list[tuple[str, bytes]] = []
        self.species: species_web.Lookup | None = None  # set by clean() if reconciling
        if self.is_bound:
            # The outgroup's choices are the submitted names; validate() then checks that
            # each one is among the sequences.
            submitted = self.data.getlist("outgroup")  # type: ignore[attr-defined]
            outgroup = self.fields["outgroup"]
            assert isinstance(outgroup, forms.MultipleChoiceField)
            outgroup.choices = [(name, name) for name in submitted]

    def clean_simulated_from(self) -> str:
        value = str(self.cleaned_data.get("simulated_from") or "")
        if not value:
            return ""
        job_id, _, seed = value.partition(":")
        job = None
        if is_uuid(job_id) and seed.isdigit() and self.owner is not None:
            job = Job.objects.for_user(self.owner).filter(id=job_id).first()
        if job is None or not job.params.get("simulate"):
            raise forms.ValidationError("That simulation is not one of yours.")
        self.simulation = (job, int(seed))
        return value

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        if self.extension is not None:
            self.inherit(data)
        else:
            self.read_inputs(data)
        self.sequence_type = self.effective_seqtype(data)
        coding = self.sequence_type == "codon"
        if not coding and CODON_ALIGNERS & set(data.get("aligners") or []):
            raise forms.ValidationError(
                "PRANK's codon model needs coding DNA: untick it, or choose coding DNA."
            )
        amino, codon = analysis_flags(self.sequence_type, data.get("analyses") or "both")
        if not amino:  # codon models only: the amino-acid model options are ignored
            data["modelfinder"] = data["fixed_models"] = False
        if not codon:
            data["codon_modelfinder"] = data["codon_fixed"] = False
        if not data.get("codon_fixed"):
            data["codon_models"] = []
        if not data.get("codon_modelfinder"):
            data["codon_candidates"] = []
        if not codon or self.stop_after(data) != "ancestors":  # codeml needs codon models
            data["codeml"] = False
        if not data.get("codeml"):
            data["codeml_models"] = []
        stage = self.stop_after(data)
        if stage == "alignments" and self.trees:
            raise forms.ValidationError("Your own trees need “Infer trees” ticked.")
        if stage != "alignments":
            if not data.get("tree_builders") and not self.trees:
                raise forms.ValidationError("Tick IQ-TREE or RAxML-NG, or upload your own tree.")
            if amino and not data.get("modelfinder") and not data.get("fixed_models"):
                raise forms.ValidationError(
                    "Choose ModelFinder, a full search under each ticked model, or both."
                )
            if codon and not data.get("codon_modelfinder") and not data.get("codon_fixed"):
                raise forms.ValidationError(
                    "Choose codon ModelFinder, a full search under each ticked codon model, "
                    "or both."
                )
            if codon and data.get("codon_fixed") and not data.get("codon_models"):
                raise forms.ValidationError(
                    "Tick the codon models to search in full, or untick that option."
                )
            if data.get("codeml") and not data.get("codeml_models"):
                raise forms.ValidationError("Tick at least one codeml model, or untick codeml.")
            code = codons.CODES.get(data.get("genetic_code") or 1, codons.STANDARD)
            if data.get("codeml") and code.paml_icode is None:
                raise forms.ValidationError(
                    f"PAML doesn't have genetic code {code.ncbi} ({code.name}): untick codeml "
                    "or choose a code PAML offers."
                )
            builders = set(data["tree_builders"])
            search = data.get("codon_search") and "iqtree" in builders
            on_amino = data.get("codon_on_amino") and builders & {"iqtree", "raxml"}
            if codon and data.get("codon_on_amino") and not amino:
                raise forms.ValidationError(
                    "Fitting codon models on the amino-acid trees needs amino-acid models: "
                    "analyse with both, or untick that option."
                )
            if codon and not search and not on_amino and not self.trees:
                raise forms.ValidationError(
                    "Codon models need a tree: IQ-TREE's search under each codon model, a fit "
                    "on the amino-acid trees, or your own tree (RAxML-NG has no codon models)."
                )
            if data.get("fixed_models") and not data.get("models"):
                raise forms.ValidationError(
                    "Tick the models to search in full, or untick “Full ML tree search”."
                )
            missing = catalog.raxml_missing(self.used_models(data))
            mixtures = [m for m in missing if catalog.is_mixture(m)]
            if "raxml" in (data.get("tree_builders") or []) and mixtures:
                raise forms.ValidationError(
                    f"RAxML-NG has no site-profile mixtures ({', '.join(mixtures)}): untick "
                    "those models, or untick RAxML-NG."
                )
            if "raxml" in (data.get("tree_builders") or []) and missing:
                matrices = listing(sorted({catalog.split(m)[0] for m in missing}))
                raise forms.ValidationError(
                    f"RAxML-NG does not have the {matrices} matrix ({', '.join(missing)}): "
                    f"untick {matrices} in the model grids, or untick RAxML-NG."
                )
        if stage == "ancestors":
            needing = sorted(NEED_OUTGROUP & set(data.get("gap_methods") or []))
            if needing and not data.get("outgroup"):
                raise forms.ValidationError(
                    f"Gap placement by {' and '.join(m.upper() for m in needing)} "
                    "needs an outgroup."
                )
        if stage == "alignments" or not data.get("reconcile"):
            data["reconciled_tree"] = False
        if data.get("reconciled_tree"):
            if not amino:
                raise forms.ValidationError(
                    "GeneRax's reconciled tree needs amino-acid models (GeneRax has no codon "
                    "models): analyse with both, or untick it."
                )
        raxml = stage != "alignments" and "raxml" in (data.get("tree_builders") or [])
        self.check_sequences(
            self.used_models(data, "models"),
            self.used_models(data, "candidates"),
            data.get("outgroup") or [],
            raxml,
            (data.get("codon_models") or []) if stage != "alignments" else [],
            (data.get("codon_candidates") or []) if stage != "alignments" else [],
            data.get("codeml_models") or [],
        )
        confirmed = self.extension is not None  # confirmed when the parent was submitted
        if (
            self.ambiguity
            and not confirmed
            and data.get("ambiguity_ack") != self.ambiguity["fingerprint"]
        ):
            raise forms.ValidationError(
                "Your coding DNA has codons with ambiguous bases: review them (under "
                "Sequences) and confirm how they will be treated before submitting."
            )
        if stage != "alignments" and data.get("reconcile"):
            self.species = self.look_up_species(data)
        return data

    def read_inputs(self, data: dict[str, Any]) -> None:
        """A job's own inputs: sequences, alignments and trees, as posted."""
        uploaded = used_files(data.get("alignment_files") or [], self.data, "alignment")
        self.alignments = self.read_uploads(uploaded, "alignment")
        self.alignments = self.dedupe(self.alignments, data.get("alignment_duplicates") or "")
        self.alignments += self.read_saved(data.get("saved_alignments") or [], self.alignments)
        self.trees = self.read_uploads(data.get("tree_files") or [], "tree")
        # A duplicated job keeps its own uploads unless the user uploads replacements.
        self.alignments = self.alignments or self.reuse.get("alignment", [])
        self.trees = self.trees or self.reuse.get("tree", [])
        text = data.get("sequences_text") or ""
        files = used_files(data.get("sequences_file") or [], self.data, "sequences")
        if text.strip() and files:
            raise forms.ValidationError("Paste sequences or upload files (not both).")
        if not text.strip() and not files and not self.alignments:
            raise forms.ValidationError(
                "Paste or upload sequences, or upload or choose your own alignment."
            )
        if files:
            self.fasta = self.merge_files(files, data.get("duplicate_names") or "")
        elif text.strip():
            self.check_size(len(text.encode()))
            self.fasta = self.merge_texts(
                [(inputs.PASTED, text.encode())], data.get("duplicate_names") or ""
            )
        self.left_out = self.use_shared(data) if data.get("shared_names") else []
        if not data.get("aligners") and not self.alignments:
            raise forms.ValidationError("Tick an aligner, or upload or choose your own alignment.")

    def inherit(self, data: dict[str, Any]) -> None:
        """An extension's inputs: the parent's sequences and alignments (each under its
        aligner's key), sequence type, genetic code and seed, whatever was posted. Only
        additions are the user's: more aligners, alignments of their own, later steps."""
        assert self.extension is not None
        parent = self.extension.parent.params
        uploaded = used_files(data.get("alignment_files") or [], self.data, "alignment")
        added = self.read_uploads(uploaded, "alignment")
        added = self.dedupe(added, data.get("alignment_duplicates") or "")
        added += self.read_saved(data.get("saved_alignments") or [], added)
        inherited = self.extension.alignments
        self.fasta = self.extension.fasta
        self.alignments = [(a.label, a.content) for a in inherited] + added
        self.alignment_keys = [a.key for a in inherited] + [extend.PROVIDED] * len(added)
        if self.alignment_keys.count(extend.PROVIDED) > MAX_PROVIDED:
            raise forms.ValidationError(f"At most {MAX_PROVIDED} of your own alignments.")
        self.trees = self.read_uploads(data.get("tree_files") or [], "tree")
        data["seqtype"] = parent.get("seqtype") or "protein"
        data["genetic_code"] = parent.get("genetic_code") or 1
        data["seed"] = parent.get("seed", 1)  # so it gives what one job would have
        data["shared_names"] = False  # the parent's sequences stand
        data["aligners"] = [
            a for a in data.get("aligners") or [] if a not in self.extension.aligners
        ]
        if self.extension.mode == extend.RESULTS:
            self.inherit_settings(data, added)
        elif self.stop_after(data) == "alignments" and not data["aligners"] and not added:
            raise forms.ValidationError(
                "Nothing to add: tick another aligner, add an alignment of your own, or "
                "switch on trees."
            )

    def inherit_settings(self, data: dict[str, Any], added: list[tuple[str, bytes]]) -> None:
        """An extension of a job with trees (results mode): every later setting is the
        parent's, whatever was posted, so new alignments run exactly as its own did; its
        own trees run on them too. What it may add: aligners, and for every alignment,
        the inherited ones too, tree sources, fixed models and ModelFinder
        (extend.Additions; what the parent has is read, not run again)."""
        assert self.extension is not None
        if added:
            raise forms.ValidationError(
                "Only aligners can be added to a job that has trees: to bring your own "
                "alignment, duplicate the job instead."
            )
        parent = self.extension.parent
        adds = extend.additions(
            parent, data.get("tree_builders") or [], data.get("models") or [],
            bool(data.get("modelfinder")),
        )  # fmt: skip
        if not data["aligners"] and not adds:
            raise forms.ValidationError(
                "Nothing to add: tick an aligner the job has not run, or a tree source or "
                "model to add to every alignment."
            )
        self.additions = adds
        self.alignments, self.alignment_keys = [], []
        self.trees = list(self.extension.trees)
        initial, selected = initial_from_job(parent.params, "", "")
        for section in ("trees", "ancestors"):
            for name in SECTIONS[section]:
                if name in initial:
                    data[name] = initial[name]
        for grid, ticked in selected.items():
            data[grid] = sorted(ticked)
        data["analyses"] = initial["analyses"]
        data["outgroup"] = parent.params.get("outgroup") or []
        data["reconcile"] = data["reconciled_tree"] = False  # refused if the parent did
        if adds:  # checked as one job with every setting would be (RAxML-NG's matrices...)
            data["tree_builders"] = list(data["tree_builders"]) + list(adds.builders)
            fixed = data["models"] if parent.params.get("fixed_models") else []
            data["models"] = list(fixed) + list(adds.models)
            data["fixed_models"] = bool(data["models"])
            data["modelfinder"] = bool(data["modelfinder"] or adds.modelfinder)
            if adds.modelfinder:
                data["candidates"] = sorted(catalog.all_models())
            if "raxml" in adds.builders:
                # New alignments' ModelFinder leaves out what RAxML-NG lacks, as a job
                # with it from the start would; an inherited pick must not need it.
                data["candidates"] = [
                    m for m in data["candidates"] if not catalog.raxml_missing([m])
                ]
                reused = self.inherited_params().get("modelfinder_reused") or []
                if bad := extend.unusable_picks(self.extension, reused):
                    raise forms.ValidationError(
                        "RAxML-NG can't be added: ModelFinder picked "
                        + ", ".join(f"{model} for {aligner}" for aligner, model in bad)
                        + ", a matrix RAxML-NG doesn't have."
                    )
            self.alignments = extend.alignment_contents(self.extension)
            self.alignment_keys = [a.key for a in self.extension.alignments]

    def look_up_species(self, data: dict[str, Any]) -> species_web.Lookup:
        """The sequences' species and species tree, from the table as the user left it.
        Never an error: what is missing makes the job run without reconciliation."""
        upload = data.get("species_tree_file")
        text = ""
        if upload:
            self.check_size(upload.size)
            text = upload.read().decode("utf-8", errors="replace")
        listed = data.get("species_table_file")
        listed_text = ""
        if listed:
            self.check_size(listed.size)
            listed_text = listed.read().decode("utf-8", errors="replace")
        fasta = self.fasta or (ungapped(self.alignments[0][1]) if self.alignments else b"")
        gathered = species_web.gather(fasta, self.data, listed_text)
        return species_web.look_up(gathered.rows, bool(data.get("species_consent")), text)

    def analysis_options(self) -> list[dict[str, Any]]:
        """The coding-DNA analysis choices, each with whether it can be picked yet."""
        chosen = self["analyses"].value() or "both"
        return [
            {"value": value, "label": label, "checked": value == chosen,
             "available": value in AVAILABLE_ANALYSES}
            for value, label in ANALYSES
        ]  # fmt: skip

    def effective_seqtype(self, data: dict[str, Any]) -> str:
        """The chosen sequence type, or with "Detect" what the sequences look like (the
        sequences, or else the first uploaded alignment's rows)."""
        chosen = data.get("seqtype") or "auto"
        if chosen != "auto":
            return str(chosen)
        fasta = self.fasta or (self.alignments[0][1] if self.alignments else b"")
        return detect_seqtype(fasta.decode("utf-8", errors="replace"))

    @staticmethod
    def used_models(data: dict[str, Any], grid: str = "") -> list[str]:
        """The models a job uses from one grid ("models", "candidates") or both (""): a
        grid counts only while its option is ticked."""
        used: list[str] = []
        if grid in ("", "candidates") and data.get("modelfinder"):
            used += data.get("candidates") or []
        if grid in ("", "models") and data.get("fixed_models"):
            used += data.get("models") or []
        return list(dict.fromkeys(used))

    @staticmethod
    def stop_after(data: dict[str, Any]) -> str:
        """The last stage to run, from the two switches (ancestors need trees)."""
        if not data.get("run_trees"):
            return "alignments"
        return "ancestors" if data.get("run_ancestors") else "trees"

    def check_size(self, size: int) -> None:
        if size > self.limits.max_upload_bytes:
            raise forms.ValidationError(
                f"An input is {size:,} bytes; your limit is {self.limits.max_upload_bytes:,}."
            )

    def merge_files(self, files: list[Any], resolution: str) -> bytes:
        """Several sequence files as one FASTA (phylobot.inputs.merge_sequences)."""
        if len(files) > inputs.MAX_SEQUENCE_FILES:
            raise forms.ValidationError(f"At most {inputs.MAX_SEQUENCE_FILES} sequence files.")
        for upload in files:
            self.check_size(upload.size)
        merged = self.merge_texts([(f.name or "", f.read()) for f in files], resolution)
        self.check_size(len(merged))
        return merged

    def merge_texts(self, contents: list[tuple[str, bytes]], resolution: str) -> bytes:
        """Sequence files or pasted text with repeated names resolved (one copy of a
        sequence given twice; another sequence under a used name as the user chose)."""
        merged = inputs.merge_sequences(contents, resolution)
        if merged.error:
            raise forms.ValidationError(f"Sequences: {merged.error}.")
        if merged.unresolved:
            names = len({c.name for c in merged.clashes})
            raise forms.ValidationError(
                f"{names} sequence name{'s are' if names > 1 else ' is'} given to different "
                "sequences: choose what to do with them under Sequences."
            )
        return merged.fasta

    def dedupe(
        self, alignments: list[tuple[str, bytes]], resolution: str
    ) -> list[tuple[str, bytes]]:
        """Uploaded alignments with repeated names resolved (inputs.dedupe_alignment)."""
        out = []
        for label, content in alignments:
            merged = inputs.dedupe_alignment(label, content, resolution)
            if merged.error:  # VALIDATE says what is wrong with it
                out.append((label, content))
                continue
            if merged.unresolved:
                raise forms.ValidationError(
                    f"Your alignment {label} uses a sequence name more than once with "
                    "different rows: choose what to do with them under Alignments."
                )
            out.append((label, merged.fasta))
        return out

    def use_shared(self, data: dict[str, Any]) -> list[str]:
        """Restrict the sequences and every alignment to the names they all have (the
        user's choice when their alignments differ); returns the names left out."""
        sources = [content for _, content in self.alignments]
        names = [set(sequence_names(c)) for c in sources]
        if self.fasta:
            names.append(set(sequence_names(self.fasta)))
        keep = inputs.shared_names(names)
        if len(keep) < MIN_SEQUENCES:
            raise forms.ValidationError(
                f"Your alignments share {len(keep)} sequence{'s' * (len(keep) != 1)}; a job "
                f"needs at least {MIN_SEQUENCES}."
            )
        try:
            self.alignments = [(n, inputs.restrict(c, keep, True)) for n, c in self.alignments]
            if self.fasta:
                self.fasta = inputs.restrict(self.fasta, keep, False)
        except (InputError, UnicodeDecodeError) as exc:
            raise forms.ValidationError(f"Input error: {exc}") from None
        return sorted(set().union(*names) - keep)

    def read_uploads(self, files: list[Any], kind: str) -> list[tuple[str, bytes]]:
        if len(files) > MAX_PROVIDED:
            raise forms.ValidationError(f"At most {MAX_PROVIDED} of your own {kind}s.")
        uploads = []
        for upload in files:
            self.check_size(upload.size)
            content = upload.read()
            if kind == "alignment":  # PHYLIP and "?" read as the pipeline reads FASTA
                try:
                    content = alignments.normalise_upload(content)
                except alignments.AlignmentError as exc:
                    raise forms.ValidationError(f"{upload.name}: {exc}.") from None
            uploads.append((label_for(upload.name or ""), content))
        return uploads

    def read_saved(
        self, keys: list[str], uploads: list[tuple[str, bytes]]
    ) -> list[tuple[str, bytes]]:
        """Alignments chosen from My alignments, labelled by their names (made unique)."""
        if len(uploads) + len(keys) > MAX_PROVIDED:
            raise forms.ValidationError(f"At most {MAX_PROVIDED} of your own alignments.")
        taken = {label.casefold() for label, _ in uploads}
        out = []
        for key in keys:
            record = self.saved[key]
            content = storage.read_bytes(record.storage_key)
            self.check_size(len(content))
            label, n = record.title[:100], 2
            while label.casefold() in taken:
                label, n = f"{record.title[:94]} ({n})", n + 1
            taken.add(label.casefold())
            out.append((label, content))
        return out

    def check_sequences(
        self,
        models: list[str],
        candidates: list[str],
        outgroup: list[str],
        raxml: bool = False,
        codon_models: list[str] | None = None,
        codon_candidates: list[str] | None = None,
        codeml_models: list[str] | None = None,
    ) -> None:
        """Run the pipeline's own validation now, so errors show on the form."""
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            if self.fasta:
                (workdir / "input.fasta").write_bytes(self.fasta)
            for index, (_, content) in enumerate(self.alignments, 1):
                (workdir / f"alignment{index}.fasta").write_bytes(content)
            for index, (_, content) in enumerate(self.trees, 1):
                (workdir / f"tree{index}.nwk").write_bytes(content)
            manifest = {
                "sequences": bool(self.fasta),
                "alignments": [
                    {"filename": label, "label": label, **inherited_key(self.alignment_keys, i)}
                    for i, (label, _) in enumerate(self.alignments)
                ],
                "trees": [{"filename": label, "label": label} for label, _ in self.trees],
                "models": models,
                "candidates": candidates,
                "raxml": raxml,
                "outgroup": outgroup,
                "seqtype": self.sequence_type,
                "genetic_code": self.cleaned_data.get("genetic_code") or 1,
                "codon_models": codon_models or [],
                "codon_candidates": codon_candidates or [],
                "codeml_models": codeml_models or [],
                "max_sequences": self.limits.max_sequences,
                "max_length": self.limits.max_length,
                "max_provided": MAX_PROVIDED,
            }
            (workdir / "manifest.json").write_text(json.dumps(manifest))
            try:
                validate(workdir / "manifest.json", workdir)
            except InputError as exc:
                raise forms.ValidationError(f"Input error: {exc}") from exc
            report = json.loads((workdir / "validation.json").read_text())
            self.ambiguity = report.get("ambiguous_codons")

    def size(self) -> workload.Size:
        """The input's size for the work estimate: the sequences, or else the first
        uploaded alignment (gaps are not counted)."""
        fasta = self.fasta or (self.alignments[0][1] if self.alignments else b"")
        return workload.size_of(fasta, self.sequence_type == "codon") or workload.IME2_SIZE

    def params(self) -> dict[str, Any]:
        """Job.params: options only (the worker allowlists them; files go by fixed names)."""
        if self.extension is not None and self.extension.mode == extend.RESULTS:
            return self.inherited_params()
        data = self.cleaned_data
        stage = self.stop_after(data)
        params: dict[str, Any] = {
            "stop_after": stage,
            "aligners": data["aligners"] or ["none"],
            "outgroup": data["outgroup"],
            "seed": data["seed"],
            "seqtype": self.sequence_type,
        }
        if self.left_out:  # display only: the sequences not every alignment had
            params["left_out"] = self.left_out[:200]
            params["left_out_count"] = len(self.left_out)
        if self.sequence_type == "codon":
            params["genetic_code"] = data.get("genetic_code") or 1
        if self.simulation is not None:
            job, seed = self.simulation
            params["simulated_from"] = {"job": str(job.id), "seed": seed, "name": job.name}
        if self.simulation is not None or any(
            self.saved[k].simulated for k in data.get("saved_alignments") or [] if k in self.saved
        ):
            params["simulated"] = True  # never passed to Nextflow: not on PARAM_ALLOWLIST
        if self.alignments:
            params["alignment_labels"] = [label for label, _ in self.alignments]
        if self.extension is not None:
            parent = self.extension.parent
            params["alignment_keys"] = self.alignment_keys
            params["extends"] = {  # display only
                "job": str(parent.id),
                "name": parent.name,
                "alignments": len(self.extension.alignments),
            }
            for key in ("left_out", "left_out_count", "simulated", "simulated_from"):
                if key in parent.params:
                    params[key] = parent.params[key]
        # for My Jobs (not a pipeline option: not on the executors' allowlist)
        counted = workload.size_of(
            self.fasta or (ungapped(self.alignments[0][1]) if self.alignments else b"")
        )
        if counted is not None:
            params["sequence_count"] = counted.sequences
        if stage != "alignments":
            params |= {
                "tree_builders": data["tree_builders"] or ["none"],
                "modelfinder": data["modelfinder"],
                "fixed_models": data["fixed_models"],
                "models": self.used_models(data, "models") or None,
                "candidates": self.used_models(data, "candidates") or None,
                "alrt": (data["alrt"] or 1000) if data["sh_alrt"] else 0,
                "bootstrap": (data["bootstrap"] or 1000) if data["ufboot"] else 0,
            }
            if self.trees:
                params["tree_labels"] = [label for label, _ in self.trees]
            if data.get("reconcile"):
                params |= self.species_params(data)
            if data.get("codon_modelfinder"):
                params["codon_modelfinder"] = True
                params["codon_candidates"] = data.get("codon_candidates") or None
            if data.get("codon_models"):
                params["codon_models"] = data["codon_models"]
            if data.get("codon_modelfinder") or data.get("codon_models"):
                params |= codon_tree_params(
                    bool(data.get("codon_search")),
                    bool(data.get("codon_on_amino")),
                    data.get("codon_amino_trees"),
                )
        if stage == "ancestors":
            params["gap_methods"] = data["gap_methods"]
            if data.get("properties"):
                extra = data.get("extra_properties") or []
                params["properties"] = [
                    k for k in properties_mod.REGISTRY if k in properties_mod.DEFAULTS or k in extra
                ]
                params["property_draws"] = (
                    data.get("property_draws") or properties_mod.DEFAULT_DRAWS
                )
            if data.get("codeml_models"):
                params["codeml_models"] = data["codeml_models"]
                params["codeml_trees"] = data.get("codeml_trees") or "best"
        return params

    def inherited_params(self) -> dict[str, Any]:
        """An extension of a job with trees (extend.extended_params)."""
        assert self.extension is not None
        return extend.extended_params(
            self.extension,
            self.cleaned_data["aligners"],
            self.additions,
            [label for label, _ in self.trees],
        )

    def species_params(self, data: dict[str, Any]) -> dict[str, Any]:
        """Reconciliation's params: `reconcile` and how polytomies are resolved when the
        species are ready (plus GeneRax's tree source if ticked), else `reconcile_skipped`
        with the reason the job page gives. `species` describes the species tree."""
        found = self.species
        if found is None or not found.ready:
            return {"reconcile_skipped": found.skipped if found else "no species were read"}
        params: dict[str, Any] = {
            "reconcile": True,
            "polytomies": data.get("polytomies") or GUIDED,
            "species": {
                "source": found.provenance.get("source", ""),
                "synth_id": found.provenance.get("synth_id", ""),
                "species": found.species,
                "polytomies": found.polytomies,
            },
        }
        if data.get("reconciled_tree"):
            builders = [t for t in data["tree_builders"] or [] if t != "none"]
            params["tree_builders"] = [*builders, RECONCILED_TREE]
        return params

    def input_files(self) -> dict[str, bytes]:
        """Job inputs by their fixed names (the executor stages these)."""
        files = {"sequences.fasta": self.fasta} if self.fasta else {}
        if self.extension is not None and self.additions:  # the parent's ModelFinder picks
            reused = self.inherited_params().get("modelfinder_reused") or []
            files |= extend.picks(self.extension, reused)
        files |= {f"alignment{i}.fasta": c for i, (_, c) in enumerate(self.alignments, 1)}
        if self.stop_after(self.cleaned_data) != "alignments":
            files |= {f"tree{i}.nwk": c for i, (_, c) in enumerate(self.trees, 1)}
            found = self.species
            if self.cleaned_data.get("reconcile") and found is not None and found.ready:
                assert found.tree is not None  # ready means a tree
                files |= {
                    "species.tsv": species_table.write(found.rows).encode(),
                    "species_tree.nwk": (found.tree.newick() + "\n").encode(),
                    "species_provenance.json": (
                        json.dumps(found.provenance, indent=2) + "\n"
                    ).encode(),
                }
        return files


@dataclass
class InputCheck:
    """The New Job page's live check of the sequences: their type and size, or the
    first problem the pipeline's validation would refuse them for."""

    seqtype: str  # in effect: chosen, or detected
    detected: bool  # the type was detected, not chosen
    error: str = ""
    sequences: int = 0
    min_length: int = 0
    max_length: int = 0
    warnings: list[str] = field(default_factory=list)
    ambiguity: dict[str, Any] | None = None  # validation.json's ambiguous_codons

    @property
    def unit(self) -> str:
        return "codons" if self.seqtype == "codon" else "residues"

    @property
    def ambiguous_rows(self) -> list[dict[str, Any]]:
        """The ambiguous codons as the confirmation dialog lists them."""
        if not self.ambiguity:
            return []
        return [
            row | {"readings": readings(row["aminos"]), "as": amino_name(row["translation"])}
            for row in self.ambiguity["listed"]
        ]

    @property
    def ambiguous_unlisted(self) -> int:
        return self.ambiguity["codons"] - len(self.ambiguity["listed"]) if self.ambiguity else 0


AMINO_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "stop",
    "X": "X (unknown)",
}  # fmt: skip


def amino_name(letter: str) -> str:
    return AMINO_NAMES.get(letter, letter)


def readings(aminos: str) -> str:
    """What a codon's readings encode, e.g. "Asn or Thr"; every amino acid for NNN."""
    acids = aminos.replace("*", "")
    if len(acids) > 4:
        return f"any of {len(acids)} amino acids" + (" or a stop" if "*" in aminos else "")
    names = [amino_name(a) for a in acids] + ["stop"] * ("*" in aminos)
    return ", ".join(names[:-1]) + f" or {names[-1]}" if len(names) > 1 else names[0]


def check_input(fasta: bytes, seqtype: str, genetic_code: int, limits: Limits) -> InputCheck:
    """Run the pipeline's validation on the sequences alone (phylobot.pipeline.validate,
    the one implementation of every rule), as the user types."""
    text = fasta.decode("utf-8", errors="replace")
    chosen = seqtype if seqtype in ("protein", "codon") else "auto"
    effective = detect_seqtype(text) if chosen == "auto" else chosen
    result = InputCheck(effective, chosen == "auto")
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        (workdir / "input.fasta").write_bytes(fasta)
        manifest = {
            "sequences": True,
            "seqtype": effective,
            "genetic_code": genetic_code,
            "max_sequences": limits.max_sequences,
            "max_length": limits.max_length,
            "max_provided": MAX_PROVIDED,
        }
        (workdir / "manifest.json").write_text(json.dumps(manifest))
        try:
            validate(workdir / "manifest.json", workdir)
        except InputError as exc:
            result.error = str(exc)
            return result
        report = json.loads((workdir / "validation.json").read_text())
    result.sequences = report["sequences"]
    result.min_length, result.max_length = report["min_length"], report["max_length"]
    # The ambiguous codons have their own line and dialog (_sequence_check.html).
    result.warnings = [w for w in report["warnings"] if not w.startswith(AMBIGUITY_WARNING)]
    result.ambiguity = report.get("ambiguous_codons")
    return result


def initial_from_job(
    params: dict[str, Any], name: str, fasta: str
) -> tuple[dict[str, Any], dict[str, set[str]]]:
    """The New Job form's values and model-grid ticks for a copy of a job ("Duplicate and
    run again"): everything as the job was submitted, to change and submit."""
    stage = params.get("stop_after", "ancestors")
    coding = params.get("seqtype") == "codon"
    amino = bool(params.get("modelfinder") or params.get("fixed_models"))
    codon = bool(params.get("codon_modelfinder") or params.get("codon_models"))
    analyses = "both" if amino and codon else ("codon" if codon else "amino")
    initial: dict[str, Any] = {
        "name": name,
        "sequences_text": fasta,
        "seqtype": params.get("seqtype", "auto") if coding else "auto",
        "genetic_code": params.get("genetic_code", 1),
        "analyses": analyses if coding else "both",
        "aligners": [a for a in params.get("aligners") or [] if a != "none"],
        "run_trees": stage != "alignments",
        "run_ancestors": stage == "ancestors",
        "tree_builders": [
            t for t in params.get("tree_builders") or [] if t not in ("none", RECONCILED_TREE)
        ],
        "modelfinder": bool(params.get("modelfinder")),
        "fixed_models": bool(params.get("fixed_models")),
        "codon_modelfinder": bool(params.get("codon_modelfinder")),
        "codon_fixed": bool(params.get("codon_models")),
        "codon_search": "search" in (params.get("codon_trees") or ["search"]),
        "codon_on_amino": "amino" in (params.get("codon_trees") or []),
        "codon_amino_trees": params.get("codon_amino_trees", "best"),
        "outgroup": params.get("outgroup") or [],
        # The species are looked up again, so the consent box starts unticked as ever.
        "reconcile": bool(params.get("reconcile") or params.get("reconcile_skipped")),
        "reconciled_tree": RECONCILED_TREE in (params.get("tree_builders") or []),
        "polytomies": params.get("polytomies", GUIDED),
        "gap_methods": params.get("gap_methods") or [],
        "properties": bool(params.get("properties")),
        "extra_properties": [
            k for k in params.get("properties") or [] if k not in properties_mod.DEFAULTS
        ],
        "property_draws": params.get("property_draws") or properties_mod.DEFAULT_DRAWS,
        "sh_alrt": bool(params.get("alrt")),
        "alrt": params.get("alrt") or 1000,
        "ufboot": bool(params.get("bootstrap")),
        "bootstrap": params.get("bootstrap") or 1000,
        "seed": params.get("seed", 1),
        "codeml": bool(params.get("codeml_models")),
        "codeml_models": params.get("codeml_models") or paml.DEFAULT_MODELS,
        "codeml_trees": params.get("codeml_trees", "best"),
    }
    if (of := params.get("simulated_from")) and of.get("job"):  # not on an unlinked import
        initial["simulated_from"] = f"{of['job']}:{of['seed']}"
    # Unset candidates meant IQ-TREE's default set: the whole grid, as on a fresh form.
    selected = {
        "candidates": set(params.get("candidates") or catalog.all_models()),
        "models": set(params.get("models") or []),
        "codon_candidates": set(params.get("codon_candidates") or catalog.all_codon_models()),
        "codon_models": set(params.get("codon_models") or []),
    }
    return initial, selected


# Each New Job step's settings (not its data: sequences, uploads, outgroup).
SECTIONS = {
    "sequences": ["seqtype", "genetic_code", "analyses"],
    "alignments": ["aligners"],
    "trees": [
        "run_trees", "tree_builders", "modelfinder", "candidates", "fixed_models",
        "models", "codon_modelfinder", "codon_candidates", "codon_fixed",
        "codon_models", "codon_search", "codon_on_amino", "codon_amino_trees",
        "sh_alrt", "alrt", "ufboot", "bootstrap", "reconcile", "reconciled_tree",
        "polytomies",
    ],
    "ancestors": [
        "run_ancestors", "gap_methods", "properties", "extra_properties",
        "property_draws", "codeml", "codeml_models", "codeml_trees",
    ],
}  # fmt: skip


def section_defaults() -> dict[str, dict[str, Any]]:
    """Each New Job step's defaults, for its "Use defaults" button: the fields' initial
    values (and the model grids' default ticks), never the sequences, uploads or outgroup,
    which belong to the data rather than the settings."""
    fields = SubmitForm.base_fields
    return {
        section: {name: fields[name].initial for name in names}
        for section, names in SECTIONS.items()
    }


def grid(*, mixtures: bool = False) -> list[dict[str, Any]]:
    """Model grid rows for the template: groups of matrices, each with one cell per variant.
    `mixtures` adds the site-profile columns: the full-search grid, not ModelFinder's."""
    variants = catalog.VARIANTS + (catalog.MIXTURE_VARIANTS if mixtures else ())
    groups: dict[str, list[dict[str, Any]]] = {}
    for matrix in catalog.MATRICES:
        cells = [{"model": matrix.name + v.suffix, "variant": v,
                  "mixture": v in catalog.MIXTURE_VARIANTS} for v in variants]  # fmt: skip
        groups.setdefault(matrix.group, []).append({"matrix": matrix.name, "cells": cells})
    return [{"group": group, "rows": rows} for group, rows in groups.items()]


def codon_grid() -> list[dict[str, Any]]:
    """The codon model grid (as `grid`): a cell's model is None where IQ-TREE refuses
    that combination (MG and MGK with empirical codon frequencies)."""
    groups: dict[str, list[dict[str, Any]]] = {}
    for matrix in catalog.CODON_MATRICES:
        cells = [
            {"model": matrix.name + v.suffix if v.freq in matrix.freqs else None, "variant": v}
            for v in catalog.CODON_VARIANTS
        ]
        groups.setdefault(matrix.group, []).append(
            {"matrix": matrix.name, "title": matrix.label, "cells": cells}
        )
    return [{"group": group, "rows": rows} for group, rows in groups.items()]


def ungapped(alignment: bytes) -> bytes:
    """An alignment's sequences with gaps removed (names kept)."""
    lines = alignment.decode(errors="replace").splitlines()
    return "\n".join(ln if ln.startswith(">") else re.sub(r"[-.]", "", ln) for ln in lines).encode()


def sequence_names(fasta: bytes) -> list[str]:
    """Sequence names in a FASTA (for the outgroup checklist); [] if it can't be read."""
    from phylobot.pipeline.validate import parse_fasta

    try:
        return [name for name, _, _ in parse_fasta(fasta.decode("utf-8"), "sequences")]
    except InputError, UnicodeDecodeError:
        return []


def listing(items: list[str]) -> str:
    """'a', 'a and b', 'a, b and c'."""
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def inherited_key(keys: list[str], index: int) -> dict[str, str]:
    """A VALIDATE manifest entry's `key`: an inherited alignment's aligner (extensions)."""
    key = keys[index] if index < len(keys) else ""
    return {"key": key} if key and key != extend.PROVIDED else {}


def your(count: int, kind: str) -> str:
    return f"your {kind}" if count == 1 else f"your {count} {kind}s"


def describe(params: dict[str, Any], alignment_files: int, tree_files: int) -> list[str]:
    """What a job will do, one phrase per stage, for the New Job page's summary."""
    if params.get("simulate"):  # a simulation (phylobot.simulations)
        from phylobot.simulations import describe as describe_simulation

        return describe_simulation(params)
    if params.get("esr"):  # an ESR job (phylobot.esr_views)
        from phylobot.esr_views import describe as describe_esr

        return describe_esr(params)
    if params.get("properties_only"):  # a follow-up job from the ancestor viewer
        of = params.get("properties_of") or {}
        keys = [properties_mod.REGISTRY[k].label for k in of.get("properties", [])
                if k in properties_mod.REGISTRY]  # fmt: skip
        count = of.get("ancestors") or 1
        return [
            f"measure {listing(keys) or 'properties'} over the posterior of "
            f"{count} ancestor{'s' * (count != 1)} already reconstructed "
            f"({of.get('draws', 1000)} draws each)"
        ]
    if params.get("uncertainty_only"):  # a follow-up job from the Uncertainty page
        of = params.get("uncertainty_of") or {}
        count = of.get("combinations") or 1
        return [
            f"split the variation of the ancestors of {count} combination"
            f"{'s' * (count != 1)} into the posterior's own and each choice's"
        ]
    if params.get("branch_site"):  # a follow-up job from the tree viewer
        of = params.get("branch_site_of") or {}
        return [
            "test for positive selection with codeml's branch-site test (Model A against its "
            f"null), foreground {branch_site.describe(of)}, on the tree of "
            f"{of.get('label', '?')} and the codon alignment"
        ]
    labels = dict(ALIGNERS) | {
        "muscle": "MUSCLE", "prank": "PRANK", "prank_codon": "PRANK (codon model)",
    }  # fmt: skip
    aligners = [labels[a] for a in params.get("aligners") or [] if a in labels]
    steps = []
    if params.get("seqtype") == "codon":
        code = codons.CODES.get(int(params.get("genetic_code") or 1), codons.STANDARD)
        name = "the standard code" if code.ncbi == 1 else f"code {code.ncbi} ({code.name})"
        steps.append(f"check and translate the coding DNA ({name})")
    if count := params.get("left_out_count"):
        shown = ", ".join(params.get("left_out") or [])
        more = count - len(params.get("left_out") or [])
        steps.append(f"use only the sequences every alignment has ({count} left out: "
                     f"{shown}{f' and {more} more' if more > 0 else ''})")  # fmt: skip
    also = ""
    if extended := params.get("extends"):  # its first alignments are another job's
        n = extended.get("alignments") or 0
        if not extended.get("results"):  # inputs of this job, among its alignment_labels
            alignment_files -= n
        name = f"“{extended['name']}”" if extended.get("name") else "the job it extends"
        if extended.get("results"):  # read in place (phylobot.joined)
            steps.append(f"keep the results of the {n} alignment{'s' * (n != 1)} of {name}")
            if skipped := params.get("skip_combinations"):  # additions to every alignment
                reused = len(params.get("modelfinder_reused") or [])
                steps.append(
                    f"add to {'them' if n != 1 else 'it'} only what is new ({len(skipped)} "
                    f"combination{'s' * (len(skipped) != 1)} read, not rerun"
                    + (f"; {reused} ModelFinder pick{'s' * (reused != 1)} reused" if reused else "")
                    + ")"
                )
        else:
            steps.append(f"use the {n} alignment{'s' * (n != 1)} of {name}")
        also = "also "
    sources = aligners + ([your(alignment_files, "alignment")] if alignment_files else [])
    if aligners:
        codon = ", each also as a codon alignment" if params.get("seqtype") == "codon" else ""
        steps.append(f"{also}align with {listing(sources)}{codon}")
    elif sources:
        steps.append(f"{also}use {sources[0]}")
    stage = params.get("stop_after", "ancestors")
    if stage == "alignments":
        return [*steps, "stop (alignments only)"]
    models = []
    if params.get("modelfinder"):
        models.append("best model by ModelFinder")
    if params.get("fixed_models"):
        fixed = params.get("models") or []
        n = len(fixed)
        models.append("1 ticked model" if n == 1 else f"each of {n} ticked models")
    if params.get("codon_modelfinder"):
        models.append("best codon model by ModelFinder")
    codon_models = list(params.get("codon_models") or [])
    if len(codon_models) == 1:
        models.append(f"codon model {codon_models[0]}")
    elif codon_models:
        models.append(f"each of {len(codon_models)} ticked codon models")
    how = f" ({listing(models)})" if models else ""
    trees = []
    builders = [name for key, name in (("iqtree", "IQ-TREE"), ("raxml", "RAxML-NG"))
                if key in (params.get("tree_builders") or [])]  # fmt: skip
    if builders:
        trees.append(f"infer trees with {listing(builders)}{how}")
    codon_trees = params.get("codon_trees") or ["search"]
    has_codon = params.get("codon_modelfinder") or codon_models
    if has_codon and builders and "amino" in codon_trees:
        which = (
            "every amino-acid tree"
            if params.get("codon_amino_trees") == "all"
            else "the best amino-acid tree"
        )
        searched = "" if "search" in codon_trees else " (no codon tree search)"
        trees.append(f"fit the codon models on {which}{searched}")
    if tree_files:
        trees.append(f"fit {your(tree_files, 'tree')}{'' if trees else how} "
                     "(topology kept, branch lengths re-estimated)")  # fmt: skip
    supports = [name for name, key in (("SH-aLRT", "alrt"), ("UFBoot", "bootstrap"))
                if params.get(key)]  # fmt: skip
    if trees:
        support = f" with {listing(supports)} support" if supports else ", no branch support"
        steps.append(" and ".join(trees) + support)
    if skipped := params.get("reconcile_skipped"):
        steps.append(f"skip reconciliation with the species tree: {skipped}")
    if params.get("reconcile"):
        how = "from the gene trees" if params.get("polytomies") != "arbitrary" else "at random"
        steps.append("reconcile every tree with the species tree (GeneRax: duplications "
                     f"and losses; its unresolved nodes resolved {how})")  # fmt: skip
        if "generax" in (params.get("tree_builders") or []):
            steps.append("search a reconciled tree with GeneRax for each amino-acid model "
                         "and test every tree source's topology (AU test)")  # fmt: skip
    if stage == "trees":
        return [*steps, "stop (no ancestral reconstruction)"]
    names = dict(GAP_METHODS)
    methods = [names[m].split(" ·")[0] for m in params.get("gap_methods") or [] if m in names]
    gaps_text = f" with gap placement by {listing(methods)}" if methods else " (no gap placement)"
    steps.append(f"reconstruct ancestral sequences{gaps_text}")
    if params.get("properties"):
        draws = params.get("property_draws") or properties_mod.DEFAULT_DRAWS
        extra = [
            properties_mod.REGISTRY[k].label.lower().replace("esm-2", "ESM-2")
            for k in params.get("properties") or []
            if k in properties_mod.REGISTRY and k not in properties_mod.DEFAULTS
        ]
        also = f", including {listing(extra)}" if extra else ""
        steps.append(
            f"measure sequence properties of every ancestor over its posterior{also} "
            f"({draws} draws each)"
        )
    if models_run := [m for m in params.get("codeml_models") or [] if m in paml.MODELS]:
        where = (
            "every codon combination"
            if params.get("codeml_trees") == "all"
            else "each alignment's best codon combination"
        )
        steps.append(f"run PAML codeml ({', '.join(models_run)}) on {where}")
    return steps


class HomologFilterForm(forms.Form):
    """The homolog filters (phylobot.homologs.filters). Every value is bounded here and
    again by filters.Settings; the defaults are suggestions, shown with the histograms
    so the user can move a threshold where their hits say it belongs."""

    EVALUES = [(v, f"{v:g}") for v in (1e-3, 1e-5, 1e-10, 1e-20, 1e-50)]

    max_evalue = forms.TypedChoiceField(
        label="E-value at most", choices=EVALUES, coerce=float, initial=1e-5
    )
    min_coverage = forms.FloatField(
        label="Query coverage at least (%)",
        min_value=0,
        max_value=100,
        initial=round(filters.Settings.min_coverage),
        help_text=(
            "Suggested 70%. Lower it for proteins whose ends are poorly conserved: then "
            "true homologs align only over the conserved domain (Ime2: about 55%)."
        ),
    )
    min_identity = forms.FloatField(
        label="Identity from (%)",
        min_value=0,
        max_value=100,
        initial=round(filters.Settings.min_identity),
    )
    max_identity = forms.FloatField(
        label="to (%)", min_value=0, max_value=100, initial=round(filters.Settings.max_identity)
    )

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        if not self.errors:
            try:
                data["settings"] = filters.Settings(
                    max_evalue=data["max_evalue"],
                    min_coverage=data["min_coverage"],
                    min_identity=data["min_identity"],
                    max_identity=data["max_identity"],
                )
            except ValueError as exc:
                raise forms.ValidationError(str(exc)) from exc
        return data


def search_choices(mode: str) -> list[tuple[str, str]]:
    """'service:database' choices; services limited to local mode are left out otherwise."""
    return [
        (f"{key}:{database}", f"{source.label} · {label}")
        for key, source in sources.SOURCES.items()
        if mode == "local" or not source.local_only
        for database, label in source.databases.items()
    ]


class HomologSearchForm(forms.Form):
    """Start a remote homolog search from one sequence. The query is checked by
    remote.clean_query before anything is sent; only its residues leave the machine."""

    name = forms.CharField(label="Name (optional)", max_length=100, required=False)
    query = forms.CharField(
        label="Your sequence",
        widget=forms.Textarea(attrs={"rows": 6, "cols": 80, "spellcheck": "false"}),
        max_length=20_000,
    )
    database = forms.ChoiceField(label="Search", widget=forms.RadioSelect)
    consent = forms.BooleanField(
        label="Consent to send my sequence to the search service",
        help_text=(
            "Agree to send the residues of your sequence (not its name) to the service "
            "you choose: the search runs at EBI or NCBI, not on this machine, and the "
            "service's terms apply. Nothing else is sent."
        ),
        error_messages={
            "required": (
                "Tick this box to search: the search runs at EBI or NCBI, so your "
                "sequence must be sent there. Nothing has been sent."
            )
        },
    )

    def __init__(
        self,
        *args: Any,
        mode: str = "local",
        typical: dict[str, float] | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.typical = typical or {}  # "service:database" -> seconds, from past searches
        field = self.fields["database"]
        assert isinstance(field, forms.ChoiceField)
        field.choices = search_choices(mode)
        default = sources.SOURCES[sources.DEFAULT_SOURCE]
        field.initial = f"{sources.DEFAULT_SOURCE}:{default.default}"

    def database_options(self) -> list[dict[str, Any]]:
        """Each database choice with what the page says beside it: speed, timing, notes,
        and how long recent searches here took (if any finished)."""
        return [
            {
                "radio": radio,
                "about": sources.ABOUT[radio.data["value"]],
                "typical_seconds": self.typical.get(radio.data["value"]),
            }
            for radio in self["database"]
        ]

    def clean_query(self) -> str:
        try:
            return remote.clean_query(self.cleaned_data["query"]).residues
        except remote.SearchError as exc:
            raise forms.ValidationError(str(exc)) from exc


class SequenceLookupForm(forms.Form):
    """Find a protein's sequence by name and species (upstream of homolog search). Only the
    name and species are sent, to UniProt; they are cleaned by lookup.clean_term before the
    query is built (phylobot.homologs.lookup)."""

    protein = forms.CharField(
        label="Protein or gene name",
        max_length=120,
        widget=forms.TextInput(
            attrs={
                "placeholder": "Cytochrome c oxidase subunit 1",
                "autocomplete": "off",
                "list": "gene-list",
                # Type-ahead (sequence_lookup.js): gene symbols from UniProt, in the species
                # below if one is entered.
                "data-suggest": "gene",
                "data-suggest-list": "gene-list",
                "data-suggest-species": "id_species",
            }  # fmt: skip
        ),
    )
    species = forms.CharField(
        label="Species (optional)",
        max_length=120,
        required=False,
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g. Homo sapiens — or leave blank",
                "autocomplete": "off",
                "list": "species-list",
                # Type-ahead: species whose Latin or common name matches.
                "data-suggest": "species",
                "data-suggest-list": "species-list",
            }  # fmt: skip
        ),
    )
    reviewed_only = forms.BooleanField(
        label="Reviewed entries only (Swiss-Prot)", required=False, initial=True
    )
    consent = forms.BooleanField(
        label="Consent to send search to UniProt",
        # A distinct id: this form and HomologSearchForm both have a "consent" field on the
        # same page, so the default id_consent would collide (one label toggling both boxes).
        widget=forms.CheckboxInput(attrs={"id": "lookup_consent"}),
        error_messages={
            "required": (
                "Tick this box to look up the sequence: the name and species are sent to "
                "UniProt. Your sequences are not."
            )
        },
    )


class SiteSettingsForm(forms.Form):
    """The Settings page (phylobot.site_settings): one field per setting, built from its
    registry. Blank = use the environment variable or built-in default."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        for s in site_settings.SETTINGS:
            field: forms.Field
            if s.field == "search_email":
                field = forms.EmailField(required=False, max_length=254)
            elif s.switch:
                field = forms.TypedChoiceField(
                    required=False,
                    choices=[("", "(default)"), ("1", "On"), ("0", "Off")],
                    coerce=lambda v: v == "1",
                    empty_value=None,
                )
            elif s.real:
                field = forms.FloatField(
                    required=False,
                    min_value=s.low,
                    max_value=s.high,
                    widget=forms.NumberInput(attrs={"step": "any"}),
                )
            else:
                field = forms.IntegerField(required=False, min_value=s.low, max_value=s.high)
            field.label = s.label
            field.help_text = s.help
            self.fields[s.field] = field
        # The browser's own limit, too: no more than the worker detected.
        if (limits := machine.detected_limits()) is not None:
            self.fields["machine_cpus"].widget.attrs["max"] = limits[0]
            self.fields["machine_memory_gb"].widget.attrs["max"] = limits[1]

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        errors = machine.check_overrides(
            data.get("machine_cpus"),
            data.get("machine_memory_gb"),
            data.get("reserved_memory_gb"),
            data.get("task_memory_gb"),
        )
        for name, message in errors.items():
            self.add_error(name, message)
        return data


class CalibrationForm(forms.Form):
    """A user's own calibration coefficients (My calibration): blank = the site's value.
    Bounded as the Settings page bounds the site's."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        bound = site_settings.CALIBRATION_BOUND
        for term in site_settings.CALIBRATION_TERMS:
            self.fields[term] = forms.FloatField(
                required=False,
                min_value=-bound,
                max_value=bound,
                widget=forms.NumberInput(attrs={"step": "any"}),
            )


class AlignmentUploadForm(forms.Form):
    """An alignment for My alignments (phylobot.alignments): a FASTA or PHYLIP file or
    pasted text (PHYLIP becomes FASTA, alignments.phylip_to_fasta), checked by
    alignments.parse as the New Job page checks its inputs."""

    name = forms.CharField(max_length=100, required=False, label="Name (optional)")
    file = forms.FileField(required=False, label="FASTA or PHYLIP file")
    text = forms.CharField(
        required=False,
        strip=False,
        label="…or paste it (FASTA or PHYLIP)",
        widget=forms.Textarea(attrs={"rows": 6, "spellcheck": "false"}),
    )
    seqtype = forms.ChoiceField(
        choices=[("auto", "Detect"), ("protein", "Protein"), ("codon", "Coding DNA (codons)")],
        initial="auto",
        required=False,
        widget=forms.RadioSelect,
        label="Sequence type",
    )
    genetic_code = forms.TypedChoiceField(
        choices=[(c.ncbi, code_label(c)) for c in codons.GENETIC_CODES],
        coerce=int,
        initial=1,
        required=False,
        label="Genetic code (coding DNA)",
    )
    # A name used twice with different rows (phylobot.inputs): rename or keep the first.
    duplicates = forms.ChoiceField(
        choices=[(inputs.RENAME, "Rename"), (inputs.FIRST, "First")], required=False,
    )  # fmt: skip

    def __init__(self, *args: Any, limits: Limits, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.limits = limits
        self.rows: list[Any] = []  # alignments.Row, set by clean()
        self.sequence_type = "protein"
        self.repeats: inputs.Merged | None = None  # repeated names, as resolved

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        upload, text = data.get("file"), data.get("text") or ""
        if upload and text.strip():
            raise forms.ValidationError("Choose a file or paste an alignment (not both).")
        if not upload and not text.strip():
            raise forms.ValidationError("Choose a file or paste an alignment.")
        size = upload.size if upload else len(text.encode())
        if size > self.limits.max_upload_bytes:
            raise forms.ValidationError(
                f"The alignment is {size:,} bytes; your limit is {self.limits.max_upload_bytes:,}."
            )
        raw = upload.read() if upload else text.encode()
        try:
            content = raw.decode("utf-8")
        except UnicodeDecodeError:
            raise forms.ValidationError("The file is not UTF-8 text.") from None
        source = upload.name if upload else "the pasted alignment"
        try:
            raw = alignments.normalise_upload(raw)
        except alignments.AlignmentError as exc:
            raise forms.ValidationError(f"Not saved: {exc}.") from None
        content = raw.decode("utf-8")
        self.repeats = inputs.dedupe_alignment(source, raw, data.get("duplicates") or "")
        if not self.repeats.error:  # else alignments.parse says what is wrong
            if self.repeats.unresolved:
                names = len({c.name for c in self.repeats.clashes})
                raise forms.ValidationError(
                    f"Not saved: {names} sequence name{'s are' if names > 1 else ' is'} used "
                    "for different rows: choose what to do with them."
                )
            content = self.repeats.fasta.decode("utf-8")
        chosen = data.get("seqtype") or "auto"
        self.sequence_type = alignments.detect_seqtype(content) if chosen == "auto" else chosen
        code = codons.CODES.get(data.get("genetic_code") or 1, codons.STANDARD)
        try:
            self.rows = alignments.parse(content, self.sequence_type, code, self.limits)
        except alignments.AlignmentError as exc:
            raise forms.ValidationError(f"Not saved: {exc}.") from None
        data["name"] = alignments.clean_name(data.get("name") or "") or (
            label_for(upload.name or "") if upload else "Pasted alignment"
        )
        data["genetic_code"] = code.ncbi
        return data
