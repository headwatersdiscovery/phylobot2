"""The Simulations page: a form that queues a simulation job, the
user's simulations beneath it, and the simulation's section of the job page.

The form's six steps mirror New Job's: name and seed; a tree (one of the user's jobs'
trees, or a random one with the user's parameters); the root sequence (an ancestor of
that job, a random one of a given length, or pasted); the model (a catalog matrix with
its rate variation); indels; and replicates. Everything the user chooses ends up on
AliSim's command line, so the form bounds every number and the spec is checked again by
`phylobot.simulate.parse` before the job is queued.
"""

import json
from dataclasses import dataclass
from typing import Any

from django import forms
from django.contrib.auth.models import User
from django.db import transaction
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse
from django.views.decorators.http import require_GET

from phylobot import ancestors, job_types, jobs, simulations, site_settings, storage, workload
from phylobot.forms import is_uuid
from phylobot.models import Job, Limits, limits_for
from phylobot.simulate import SpecError
from phylobot.simulate.spec import SHAPES

TREE_SOURCES = [("job", "A tree from one of my jobs"), ("random", "A random tree")]
ROOT_SOURCES = [
    ("ancestor", "An ancestor from that job"),
    ("random", "A random sequence, of a length I choose"),
    ("paste", "A sequence I paste"),
]


def current_user(request: HttpRequest) -> User:
    """The signed-in user (local mode: the local user); middleware guarantees one."""
    assert isinstance(request.user, User)
    return request.user


class SimulationForm(forms.Form):
    name = forms.CharField(max_length=100, required=False, label="Simulation name (optional)")
    seed = forms.IntegerField(min_value=1, max_value=2**31 - 1_001, initial=1)

    tree_source = forms.ChoiceField(choices=TREE_SOURCES, initial="job", widget=forms.RadioSelect)
    job = forms.UUIDField(required=False)
    combination = forms.CharField(required=False, max_length=200)
    shape = forms.ChoiceField(choices=[(k, v) for k, v in SHAPES.items()], initial="yh")
    taxa = forms.IntegerField(min_value=3, initial=20, required=False, label="Sequences")
    length_min = forms.FloatField(min_value=0.0001, max_value=10, initial=0.02, required=False)
    length_mean = forms.FloatField(min_value=0.0001, max_value=10, initial=0.1, required=False)
    length_max = forms.FloatField(min_value=0.0001, max_value=10, initial=0.3, required=False)
    birth = forms.FloatField(min_value=0.0001, max_value=100, initial=0.1, required=False)
    death = forms.FloatField(min_value=0, max_value=100, initial=0.05, required=False)

    root_source = forms.ChoiceField(
        choices=ROOT_SOURCES, initial="ancestor", widget=forms.RadioSelect
    )
    node = forms.CharField(required=False, max_length=100)
    root_length = forms.IntegerField(min_value=10, initial=300, required=False)
    root_text = forms.CharField(required=False, widget=forms.Textarea(attrs={"rows": 4}))

    matrix = forms.ChoiceField(choices=[(m, m) for m in simulations.MATRICES], initial="LG")
    rate = forms.ChoiceField(choices=list(simulations.RATES.items()), initial="G4")
    alpha = forms.FloatField(min_value=0.01, max_value=100, initial=0.5, required=False,
                             label="Gamma shape (alpha)")  # fmt: skip
    invariant = forms.FloatField(min_value=0, max_value=0.99, initial=0.1, required=False,
                                 label="Proportion of invariant sites")  # fmt: skip

    indels = forms.BooleanField(required=False, initial=False, label="Simulate indels")
    insertion_rate = forms.FloatField(min_value=0.0001, max_value=1, initial=0.03, required=False)
    deletion_rate = forms.FloatField(min_value=0.0001, max_value=1, initial=0.03, required=False)
    indel_mean = forms.FloatField(min_value=1, max_value=50, initial=3, required=False,
                                  label="Mean indel length")  # fmt: skip

    replicates = forms.IntegerField(min_value=1, initial=1)

    def __init__(self, *args: Any, owner: User, limits: Limits, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.owner, self.limits = owner, limits
        self.max_replicates = site_settings.simulation_max_replicates()
        replicates = self.fields["replicates"]
        assert isinstance(replicates, forms.IntegerField)
        replicates.max_value = self.max_replicates
        replicates.validators.append(_at_most(self.max_replicates))
        taxa = self.fields["taxa"]
        assert isinstance(taxa, forms.IntegerField)
        taxa.validators.append(_at_most(limits.max_sequences))
        length = self.fields["root_length"]
        assert isinstance(length, forms.IntegerField)
        length.validators.append(_at_most(limits.max_length))
        self.files_out: dict[str, bytes] = {}  # input file name -> content
        self.spec: dict[str, Any] = {}
        self.source: dict[str, Any] = {}  # what the simulation was built on (display only)
        self.size = workload.Size(0, 0)

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        if self.errors:
            return data
        tree, taxa = self._tree(data)
        root, length = self._root(data)
        rate = data["rate"]
        alpha = _given(data, "alpha", 0.5)
        invariant = _given(data, "invariant", 0.1)
        model = simulations.model_string(data["matrix"], rate, alpha, invariant)
        indels = None
        if data.get("indels"):
            indels = (
                _given(data, "insertion_rate", 0.03),
                _given(data, "deletion_rate", 0.03),
                _given(data, "indel_mean", 3.0),
            )
        replicates = data["replicates"]
        written = simulations.residues_written(taxa, length, replicates)
        cap = site_settings.simulation_max_residues()
        if written > cap:
            raise forms.ValidationError(
                f"That would write {written:,} residues ({replicates} replicate(s) x "
                f"{2 * taxa - 1} rows, ancestors included, x {length} sites); the limit is "
                f"{cap:,}. Use fewer replicates, sequences or sites."
            )
        try:
            self.spec = simulations.build_spec(
                description=data.get("name") or "",
                tree=tree,
                root=root,
                model=model,
                indels=indels,
                seed=data["seed"],
                replicates=replicates,
            )
        except (SpecError, ValueError) as error:
            raise forms.ValidationError(str(error)) from error
        self.files_out["spec.json"] = json.dumps(self.spec, indent=2).encode()
        self.size = workload.Size(taxa, length)
        return data

    def _tree(self, data: dict[str, Any]) -> tuple[dict[str, Any], int]:
        if data["tree_source"] == "random":
            taxa = data.get("taxa")
            if not taxa:
                raise forms.ValidationError("Say how many sequences the random tree has.")
            low, mean, high = (_given(data, f"length_{k}", 0.1) for k in ("min", "mean", "max"))
            if not low <= mean <= high or low == high:
                raise forms.ValidationError(
                    "Branch lengths: minimum <= mean <= maximum, the minimum below the maximum."
                )
            tree: dict[str, Any] = {"random": data["shape"], "taxa": taxa,
                                    "lengths": [low, mean, high]}  # fmt: skip
            if data["shape"] == "bd":
                birth, death = _given(data, "birth", 0.1), data.get("death") or 0
                if death >= birth:
                    raise forms.ValidationError("The death rate must be below the birth rate.")
                tree |= {"birth": birth, "death": death}
            if data["root_source"] == "ancestor":
                raise forms.ValidationError(
                    "A random tree has no ancestors to start from: choose a random or pasted "
                    "root sequence."
                )
            self.source = {"kind": "random"}
            return tree, taxa
        job = self._job(data)
        choice = next(
            (c for c in simulations.tree_choices(job) if c.key == data.get("combination")), None
        )
        if choice is None:
            raise forms.ValidationError("Choose one of that job's trees.")
        newick = storage.read_bytes(f"{job.storage_prefix}results/{choice.path}").decode()
        renamed, names = simulations.rename_tree(newick)
        if len(names) > self.limits.max_sequences:
            raise forms.ValidationError(
                f"That tree has {len(names)} sequences; your limit is {self.limits.max_sequences}."
            )
        self.files_out["sim_tree.nwk"] = renamed.encode()
        self.files_out["names.tsv"] = (
            "simulated\toriginal\n" + "".join(f"{s}\t{o}\n" for s, o in names.items())
        ).encode()
        self.source = {"kind": "job", "job": str(job.id), "combination": choice.key,
                       "label": f"{job.name or 'job'} · {choice.label}"}  # fmt: skip
        return {}, len(names)

    def _job(self, data: dict[str, Any]) -> Job:
        job = None
        if data.get("job"):
            job = Job.objects.for_user(self.owner).filter(id=data["job"]).first()
        if job is None or job.status != Job.Status.SUCCEEDED:
            raise forms.ValidationError("Choose one of your finished jobs.")
        return job

    def _root(self, data: dict[str, Any]) -> tuple[dict[str, Any], int]:
        source = data["root_source"]
        if source == "random":
            length = data.get("root_length")
            if not length:
                raise forms.ValidationError("Say how long the random root sequence is.")
            return {"length": length}, length
        if source == "paste":
            try:
                sequence = simulations.clean_root(_fasta_body(data.get("root_text") or ""))
            except ValueError as error:
                raise forms.ValidationError(str(error)) from error
        else:
            job = self._job(data)
            versions = simulations.ancestor_versions(job, data.get("combination") or "")
            version = simulations.default_version(job, versions)
            node = data.get("node") or ""
            if node not in versions.get(version, {}):
                raise forms.ValidationError("Choose one of the tree's ancestors as the root.")
            sequence = simulations.clean_root(versions[version][node])
            self.source |= {"node": node, "version": version}
        if len(sequence) > self.limits.max_length:
            raise forms.ValidationError(
                f"The root sequence has {len(sequence)} residues; your limit is "
                f"{self.limits.max_length}."
            )
        self.files_out["sim_root.fasta"] = f">root\n{sequence}\n".encode()
        return {"name": "root"}, len(sequence)

    def params(self) -> dict[str, Any]:
        """Job.params: the spec, and what it was built on (display only)."""
        return {"simulate": True, "spec": self.spec, "simulation_of": self.source}


def _at_most(limit: int) -> Any:
    def check(value: int | None) -> None:
        if value is not None and value > limit:
            raise forms.ValidationError(f"At most {limit}.")

    return check


def _given(data: dict[str, Any], key: str, default: float) -> float:
    value = data.get(key)
    return float(default if value is None else value)


def _fasta_body(text: str) -> str:
    """A pasted sequence, with or without a FASTA header (one sequence only)."""
    lines = [line.strip() for line in text.strip().splitlines()]
    if sum(line.startswith(">") for line in lines) > 1:
        raise ValueError("paste one sequence")
    return "".join(line for line in lines if not line.startswith(">"))


@dataclass
class NodeOption:
    name: str
    label: str


def node_options(job: Job, combination: str) -> list[NodeOption]:
    """A combination's ancestors, root first, for the root-sequence menu."""
    choice = next((c for c in simulations.tree_choices(job) if c.key == combination), None)
    if choice is None:
        return []
    versions = simulations.ancestor_versions(job, combination)
    version = simulations.default_version(job, versions)
    names = list(versions.get(version, {}))
    newick = storage.read_bytes(f"{job.storage_prefix}results/{choice.path}").decode()
    return [NodeOption(c.name, c.label) for c in ancestors.node_choices(newick, names, {})]


def simulations_page(request: HttpRequest) -> HttpResponse:
    owner = current_user(request)
    limits = limits_for(owner)
    if request.method == "POST":
        form = SimulationForm(request.POST, owner=owner, limits=limits)
        if form.is_valid():
            with transaction.atomic():
                params = form.params()
                job = Job.objects.create(
                    owner=owner,
                    name=form.cleaned_data["name"],
                    params=params,
                    work_units=workload.estimate(params, form.size).units,
                    input_size=workload.size_record(form.size),
                )
                for filename, content in form.files_out.items():
                    storage.put_bytes(f"{job.storage_prefix}inputs/{filename}", content)
                transaction.on_commit(lambda: jobs.enqueue(job))
            return redirect("job_detail", job_id=job.id)
    else:
        initial: dict[str, Any] = {}
        if (job_id := request.GET.get("job")) and is_uuid(job_id):
            initial["job"] = job_id
            initial["combination"] = request.GET.get("combination", "")
        form = SimulationForm(initial=initial, owner=owner, limits=limits)
    chosen = _chosen_job(owner, form)
    combination = str(form["combination"].value() or "")
    choices = simulations.tree_choices(chosen) if chosen else []
    if choices and combination not in {c.key for c in choices}:
        combination = choices[0].key
    context = {
        "form": form,
        "has_jobs": simulations.jobs_with_trees(owner).exists(),
        "chosen_job": chosen,
        "tree_choices": choices,
        "combination": combination,
        "nodes": node_options(chosen, combination) if chosen and combination else [],
        "node": str(form["node"].value() or ""),
        "shapes": simulations.SHAPE_HELP,
        "max_replicates": form.max_replicates,
        "limits": limits,
        "simulations": Job.objects.for_user(owner).filter(params__simulate=True)[:100],
    }
    return render(request, "phylobot/simulations.html", context)


def _chosen_job(owner: User, form: SimulationForm) -> Job | None:
    value = str(form["job"].value() or "")
    if not is_uuid(value):
        return None
    job = Job.objects.for_user(owner).filter(id=value, status=Job.Status.SUCCEEDED).first()
    return job


@require_GET
def tree_options(request: HttpRequest) -> HttpResponse:
    """The tree and ancestor menus for the chosen job (HTMX, as the job or tree changes)."""
    owner = current_user(request)
    job_id = request.GET.get("job", "")
    job = None
    if is_uuid(job_id):
        job = Job.objects.for_user(owner).filter(id=job_id, status=Job.Status.SUCCEEDED).first()
    choices = simulations.tree_choices(job) if job else []
    combination = request.GET.get("combination", "")
    if combination not in {c.key for c in choices}:
        combination = choices[0].key if choices else ""
    nodes = node_options(job, combination) if job and combination else []
    context = {"tree_choices": choices, "combination": combination, "nodes": nodes,
               "node": request.GET.get("node", ""), "chosen_job": job}  # fmt: skip
    return render(request, "phylobot/_simulation_tree_options.html", context)


PICKER_PAGE = 30


@require_GET
def job_picker(request: HttpRequest) -> HttpResponse:
    """The job picker's rows (HTMX): finished jobs with trees, newest first, filtered as
    the user types (`q`), a page at a time as the list scrolls (`page`)."""
    text = request.GET.get("q", "").strip()[:100]
    page = request.GET.get("page", "0")
    number = int(page) if page.isdigit() and len(page) < 6 else 0
    start = number * PICKER_PAGE
    found = simulations.jobs_with_trees(current_user(request), text).prefetch_related("tags")
    rows = list(found[start : start + PICKER_PAGE + 1])
    context = {
        "rows": [(job, job_types.types_of(job.params)) for job in rows[:PICKER_PAGE]],
        "next_page": number + 1 if len(rows) > PICKER_PAGE else None,
        "q": text,
        "first": number == 0,
    }
    return render(request, "phylobot/_simulation_job_rows.html", context)


def job_prefill(request: HttpRequest) -> dict[str, Any] | None:
    """For the New Job page opened with ?simulation=<job>&replicate=<seed>: that
    replicate's unaligned sequences as step 2 (submitted through the normal form)."""
    job_id, seed = request.GET.get("simulation", ""), request.GET.get("replicate", "")
    if not job_id:
        return None
    if not is_uuid(job_id) or not seed.isdigit():
        raise Http404
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if not simulations.is_simulation(job) or job.status != Job.Status.SUCCEEDED:
        raise Http404
    try:
        fasta = simulations.replicate_input(job, int(seed))
    except Exception as error:
        raise Http404 from error
    return {"job": job, "seed": int(seed), "fasta": fasta,
            "back_url": reverse("job_detail", args=[job.id])}  # fmt: skip


def context(job: Job, keys: set[str]) -> dict[str, Any]:
    """The job page's simulation section: the replicates and their files."""
    replicates = simulations.replicates(job, keys) if job.status == Job.Status.SUCCEEDED else []
    names = simulations.original_names(job)
    return {
        "simulation": {
            "source": job.params.get("simulation_of") or {},
            "spec": job.params.get("spec") or {},
            "replicates": replicates,
            "names": sorted(names.items()),
            "new_job": reverse("submit"),
        }
    }
