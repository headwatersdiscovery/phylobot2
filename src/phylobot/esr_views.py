"""The ESR page: a form that queues an extant-sequence reconstruction job, the
user's ESR jobs beneath it, and the ESR section of the job page.

An ESR job holds sequences out of an alignment, reconstructs each as if it were an
ancestor (`phylobot.esr`) and scores the posterior against the real sequence. It runs on
one alignment, tree and model: one of the user's finished combinations (its published
alignment.fasta, tree.nwk and model.txt), or an alignment and tree the user uploads with
a catalog model, whose branch lengths are then fitted first.

Steps, as the Simulations page's: name and seed; the alignment, tree and model; which
sequences to hold out; the levels (E0, E1). The uploads are hostile: they are parsed
here and checked by the same `phylobot.pipeline.esr.check` the pipeline runs again, and
no name in them ever reaches a command line (the pipeline renames sequences t1, t2, ...).
"""

import csv
import io
import json
import random
import re
from dataclasses import dataclass, field
from typing import Any

from django import forms
from django.contrib.auth.models import User
from django.db import transaction
from django.http import Http404, HttpRequest, HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.views.decorators.http import require_GET

from phylobot import (
    bundle,
    catalog,
    esr,
    esr_lengths,
    jobs,
    joined,
    manifest,
    simulations,
    storage,
    workload,
)
from phylobot.executors import ESR_DIR
from phylobot.forms import is_uuid
from phylobot.models import Job, Limits, limits_for
from phylobot.pipeline import esr as esr_step
from phylobot.trees import load_tree

SOURCES = [
    ("job", "A combination from one of my jobs"),
    ("upload", "An alignment and tree I upload"),
]
TIP_CHOICES = [
    ("all", "Every sequence"),
    ("sample", "A random sample (drawn with the seed)"),
    ("list", "These sequences"),
]
LEVEL_HELP = {
    "e0": "Sennett & Theobald's sequence-wise ESR: the tree, branch lengths and model as "
    "fitted with the sequence in. Quick (about a second a sequence).",
    "e1": "The model and branch lengths refitted without the sequence, which goes back on "
    "the same branch: the sequence no longer informs its own reconstruction. A model fit "
    "per sequence (seconds to a minute each, depending on the model).",
}
COMBINATION_FILES = ("alignment.fasta", "tree.nwk", "model.txt")
# The job page's "?" beside E0 and E1 (HTML: _help.html shows it unescaped).
LEVELS_POPOVER = (
    "<span><b>E0</b> (Sennett &amp; Theobald's sequence-wise ESR): the tree, branch "
    "lengths and model as they were fitted with the sequence in. The sequence still "
    "informs its own reconstruction, through its branch length and the model's "
    "parameters.</span>"
    "<span><b>E1</b>: the model and branch lengths refitted without the sequence, which "
    "then goes back on the same branch. Closer to a real ancestor, whose sequence no "
    "data informs directly.</span>"
    "<span>If E1 is less accurate than E0, the sequence was helping to reconstruct "
    "itself.</span>"
)


def parse_fasta(text: str) -> dict[str, str]:
    """{name: sequence}, each name the whole header line (as `seqio.read_fasta`: the
    tree names sequences the same way); ValueError on a name given twice."""
    records: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            if name in records:
                raise ValueError(f"the name {name!r} is given twice")
            records[name] = ""
        elif name is not None:
            records[name] += "".join(line.split())
    return records


def current_user(request: HttpRequest) -> User:
    assert isinstance(request.user, User)
    return request.user


@dataclass(frozen=True)
class CombinationChoice:
    key: str
    label: str
    model: str  # IQ-TREE's model string

    @property
    def text(self) -> str:
        """The menu's text: the label, with the fitted model when the label does not name
        it (ModelFinder's pick)."""
        return self.label if self.model in self.label else f"{self.label} · {self.model}"


def combination_choices(job: Job) -> list[CombinationChoice]:
    """A finished job's combinations ESR can run on: amino-acid, with a model whose
    parameters can be fixed, and the three files published (an extension's through
    `joined`)."""
    files = joined.results(job)
    if "summary.json" not in files.keys:
        return []
    out = []
    for combo in files.summary.get("combinations", []):
        key = combo.get("key", "")
        model = str((combo.get("fit") or {}).get("model") or "")
        if combo.get("data", "amino") != "amino" or not esr.supported(model):
            continue
        if all(f"{key}/{name}" in files.keys for name in COMBINATION_FILES):
            out.append(CombinationChoice(key, combo.get("label", key), model))
    return out


MATRICES = [(m.name, m.name) for m in catalog.MATRICES]
VARIANTS = [(v.suffix, v.label) for v in catalog.VARIANTS]


class ESRForm(forms.Form):
    name = forms.CharField(max_length=100, required=False, label="Name (optional)")
    seed = forms.IntegerField(min_value=1, max_value=2**31 - 1, initial=1)

    source = forms.ChoiceField(choices=SOURCES, initial="job", widget=forms.RadioSelect)
    job = forms.UUIDField(required=False)
    combination = forms.CharField(required=False, max_length=200)
    alignment_file = forms.FileField(required=False, label="Alignment (FASTA)")
    tree_file = forms.FileField(required=False, label="Tree (Newick)")
    matrix = forms.ChoiceField(choices=MATRICES, initial="LG")
    variant = forms.ChoiceField(choices=VARIANTS, initial="+G4", required=False)

    tips = forms.ChoiceField(choices=TIP_CHOICES, initial="all", widget=forms.RadioSelect)
    sample = forms.IntegerField(min_value=1, initial=20, required=False, label="How many")
    tip_names = forms.CharField(required=False, widget=forms.Textarea(attrs={"rows": 4}),
                                label="Sequence names, one per line")  # fmt: skip

    e0 = forms.BooleanField(required=False, initial=True, label="E0")
    e1 = forms.BooleanField(required=False, initial=True, label="E1")

    def __init__(self, *args: Any, owner: User, limits: Limits, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.owner, self.limits = owner, limits
        self.files_out: dict[str, bytes] = {}
        self.plan: dict[str, Any] = {}
        self.simulated = False
        self.size = workload.Size(0, 0)

    def clean(self) -> dict[str, Any]:
        data = super().clean() or {}
        if self.errors:
            return data
        levels = [lvl for lvl in esr.LEVELS if data.get(lvl)]
        if not levels:
            raise forms.ValidationError("Choose at least one level: E0, E1 or both.")
        if data["source"] == "job":
            alignment, newick, model, source = self._from_job(data)
            fit_lengths = False
        else:
            alignment, newick, model = self._uploaded(data)
            source, fit_lengths = {"kind": "upload"}, True
        try:
            esr_step.check(alignment, newick, model)
        except ValueError as error:
            raise forms.ValidationError(f"These inputs can't be used: {error}.") from error
        if "e1" in levels and not esr_step.binary(newick):
            raise forms.ValidationError(
                "E1 prunes each sequence from the tree, so it needs a binary tree; this one has "
                "polytomies. Untick E1, or resolve them first."
            )
        tips = self._tips(data, list(alignment))
        settings = {"levels": levels, "tips": tips, "seed": data["seed"],
                    "fit_lengths": fit_lengths}  # fmt: skip
        records = "".join(f">{name}\n{seq}\n" for name, seq in alignment.items())
        self.files_out = {
            "alignment.fasta": records.encode(),
            "tree.nwk": (newick.strip() + "\n").encode(),
            "model.txt": f"{model}\n".encode(),
            "settings.json": (json.dumps(settings, indent=2) + "\n").encode(),
        }
        length = len(next(iter(alignment.values())))
        self.size = workload.Size(len(alignment), length)
        self.plan = source | {
            "model": model, "levels": levels, "tips": len(tips), "tip_choice": data["tips"],
            "fit_lengths": fit_lengths, "sequences": len(alignment), "columns": length,
        }  # fmt: skip
        return data

    def _from_job(self, data: dict[str, Any]) -> tuple[dict[str, str], str, str, dict[str, Any]]:
        job = None
        if data.get("job"):
            job = Job.objects.for_user(self.owner).filter(id=data["job"]).first()
        if job is None or job.status != Job.Status.SUCCEEDED:
            raise forms.ValidationError("Choose one of your finished jobs.")
        choice = next((c for c in combination_choices(job) if c.key == data.get("combination")),
                      None)  # fmt: skip
        if choice is None:
            raise forms.ValidationError("Choose one of that job's combinations.")
        files = joined.results(job)
        read = files.read_many([f"{choice.key}/{name}" for name in COMBINATION_FILES])
        alignment = parse_fasta(read[f"{choice.key}/alignment.fasta"])
        model = read[f"{choice.key}/model.txt"].strip()
        self.simulated = bool(job.params.get("simulated") or job.params.get("simulate"))
        source = {"kind": "job", "job": str(job.id), "combination": choice.key,
                  "label": f"{job.name or 'job'} · {choice.label}"}  # fmt: skip
        return alignment, read[f"{choice.key}/tree.nwk"], model, source

    def _uploaded(self, data: dict[str, Any]) -> tuple[dict[str, str], str, str]:
        uploads = {}
        for key, what in (("alignment_file", "an alignment"), ("tree_file", "a tree")):
            upload = data.get(key)
            if upload is None:
                raise forms.ValidationError(f"Upload {what}.")
            if upload.size > self.limits.max_upload_bytes:
                raise forms.ValidationError(
                    f"{upload.name} is {upload.size:,} bytes; your limit is "
                    f"{self.limits.max_upload_bytes:,}."
                )
            try:
                uploads[key] = upload.read().decode()
            except UnicodeDecodeError as error:
                raise forms.ValidationError(f"{upload.name} is not a text file.") from error
        try:
            alignment = {
                n: seq.upper() for n, seq in parse_fasta(uploads["alignment_file"]).items()
            }
        except ValueError as error:
            raise forms.ValidationError(f"The alignment can't be used: {error}.") from error
        if not alignment:
            raise forms.ValidationError("The alignment has no sequences (FASTA: >name lines).")
        if "" in alignment:
            raise forms.ValidationError("Every sequence in the alignment needs a name.")
        if len(alignment) > self.limits.max_sequences:
            raise forms.ValidationError(
                f"The alignment has {len(alignment)} sequences; your limit is "
                f"{self.limits.max_sequences}."
            )
        if max(len(s) for s in alignment.values()) > self.limits.max_length:
            raise forms.ValidationError(
                f"Alignments may be at most {self.limits.max_length} columns."
            )
        try:
            tree = load_tree(uploads["tree_file"])
        except Exception as error:
            raise forms.ValidationError(
                "The tree is not a Newick tree that can be read."
            ) from error
        newick = esr.plain(tree.as_string(schema="newick", suppress_rooting=True))
        model = data["matrix"] + (data.get("variant") or "")
        if not catalog.is_model(model) or not esr.supported(model):
            raise forms.ValidationError("Choose a model from the menus.")
        return alignment, newick, model

    def _tips(self, data: dict[str, Any], names: list[str]) -> list[str]:
        match data["tips"]:
            case "all":
                return names
            case "sample":
                n = data.get("sample")
                if not n:
                    raise forms.ValidationError("Say how many sequences to sample.")
                return sorted(random.Random(data["seed"]).sample(names, min(n, len(names))),
                              key=names.index)  # fmt: skip
        wanted = [line.strip() for line in (data.get("tip_names") or "").splitlines()]
        wanted = list(dict.fromkeys(w for w in wanted if w))
        if not wanted:
            raise forms.ValidationError("Name at least one sequence to hold out.")
        if unknown := [w for w in wanted if w not in names]:
            shown = ", ".join(unknown[:5]) + ("…" if len(unknown) > 5 else "")
            raise forms.ValidationError(f"Not in the alignment: {shown}.")
        return wanted

    def params(self) -> dict[str, Any]:
        """Job.params: `esr` and the seed, and what the run was built on (display, the
        step graph's tip count and the work estimate)."""
        params: dict[str, Any] = {"esr": True, "seed": self.cleaned_data["seed"],
                                  "esr_of": self.plan}  # fmt: skip
        if self.simulated:
            params["simulated"] = True  # simulated data never passes as real
        return params


def esr_page(request: HttpRequest) -> HttpResponse:
    owner = current_user(request)
    limits = limits_for(owner)
    if request.method == "POST":
        form = ESRForm(request.POST, request.FILES, owner=owner, limits=limits)
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
                    storage.put_bytes(f"{job.storage_prefix}inputs/{ESR_DIR}/{filename}", content)
                transaction.on_commit(lambda: jobs.enqueue(job))
            return redirect("job_detail", job_id=job.id)
    else:
        initial: dict[str, Any] = {}
        if (job_id := request.GET.get("job")) and is_uuid(job_id):
            initial["job"] = job_id
            initial["combination"] = request.GET.get("combination", "")
        form = ESRForm(initial=initial, owner=owner, limits=limits)
    chosen = _chosen_job(owner, str(form["job"].value() or ""))
    choices = combination_choices(chosen) if chosen else []
    combination = str(form["combination"].value() or "")
    if choices and combination not in {c.key for c in choices}:
        combination = choices[0].key
    context = {
        "form": form,
        "has_jobs": simulations.jobs_with_trees(owner).exists(),
        "chosen_job": chosen,
        "combination_choices": choices,
        "combination": combination,
        "level_help": LEVEL_HELP,
        "limits": limits,
        "esr_jobs": Job.objects.for_user(owner).filter(params__esr=True)[:100],
    }
    return render(request, "phylobot/esr.html", context)


def _chosen_job(owner: User, value: str) -> Job | None:
    if not is_uuid(value):
        return None
    return Job.objects.for_user(owner).filter(id=value, status=Job.Status.SUCCEEDED).first()


@require_GET
def combination_options(request: HttpRequest) -> HttpResponse:
    """The combination menu for the chosen job (HTMX, as the job changes)."""
    job = _chosen_job(current_user(request), request.GET.get("job", ""))
    choices = combination_choices(job) if job else []
    combination = request.GET.get("combination", "")
    if combination not in {c.key for c in choices}:
        combination = choices[0].key if choices else ""
    context = {"combination_choices": choices, "combination": combination, "chosen_job": job}
    return render(request, "phylobot/_esr_combination_options.html", context)


# --- the job page --------------------------------------------------------------------

SERIES = {"e0": 1, "e1": 2}  # level -> --viz slot
# The slots' light values, as SVG attributes so a downloaded figure keeps them (plots.js
# saves the SVG without the page's CSS); dark.css remaps them on the page.
COLOURS = {1: "#2a78d6", 2: "#eb6834"}
# The reliability chart, in SVG user units.
CHART_W, CHART_H = 400, 320
LEFT, RIGHT, TOP, BOTTOM = 52, 16, 14, 44


@dataclass
class Point:
    x: float
    y: float
    pp: float
    accuracy: float
    sites: int
    low: float
    high: float


@dataclass
class Series:
    level: str
    slot: int
    points: list[Point]
    legend_y: float = 0.0  # its key, drawn inside the figure so downloads keep it

    @property
    def colour(self) -> str:
        return COLOURS[self.slot]

    @property
    def path(self) -> str:
        return " ".join(f"{'M' if i == 0 else 'L'}{p.x:.1f},{p.y:.1f}"
                        for i, p in enumerate(self.points))  # fmt: skip


@dataclass
class Chart:
    width: int = CHART_W
    height: int = CHART_H
    left: int = LEFT
    top: int = TOP
    right_edge: int = CHART_W - RIGHT
    bottom_edge: int = CHART_H - BOTTOM
    series: list[Series] = field(default_factory=list)
    ticks: list[dict[str, Any]] = field(default_factory=list)

    @property
    def legend_x(self) -> int:
        return self.left + 12

    @property
    def x_middle(self) -> float:
        return (self.left + self.right_edge) / 2

    @property
    def y_middle(self) -> float:
        return (self.top + self.bottom_edge) / 2

    @property
    def diagonal(self) -> str:
        return f"M{self.left},{self.bottom_edge} L{self.right_edge},{self.top}"


def chart(levels: dict[str, Any]) -> Chart:
    """Observed accuracy against the ML residue's posterior, a point per non-empty bin
    (at the bin's mean posterior), one series per level; the diagonal is perfect
    calibration."""
    out = Chart()
    w, h = out.right_edge - out.left, out.bottom_edge - out.top

    def x(v: float) -> float:
        return out.left + v * w

    def y(v: float) -> float:
        return out.bottom_edge - v * h

    out.ticks = [
        {"value": f"{t / 10:.1f}", "x": x(t / 10), "y": y(t / 10)} for t in range(0, 11, 2)
    ]
    for level, slot in SERIES.items():
        if level not in levels:
            continue
        points = [
            Point(
                x(b["pp"]),
                y(b["accuracy"]),
                b["pp"],
                b["accuracy"],
                b["sites"],
                b["low"],
                b["high"],
            )  # fmt: skip
            for b in levels[level]["bins"]
            if b["sites"]
        ]
        out.series.append(Series(level, slot, points, out.top + 12 + 20 * len(out.series)))
    return out


@dataclass
class TipRow:
    tip: str
    pendant: str
    sites: str
    cells: list[tuple[str, str]]  # (accuracy, mean posterior) per level shown


def tip_rows(text: str, levels: list[str]) -> list[TipRow]:
    """esr_tips.tsv as one row per held-out sequence, its levels side by side."""
    records: dict[str, dict[str, dict[str, str]]] = {}
    lines = text.splitlines()
    header = lines[0].split("\t") if lines else []
    for line in lines[1:]:
        record = dict(zip(header, line.split("\t"), strict=False))
        records.setdefault(record["tip"], {})[record["level"]] = record
    rows = []
    for tip, by_level in records.items():
        first = next(iter(by_level.values()))
        cells = [(by_level[lvl]["accuracy"], by_level[lvl]["mean_pp"]) if lvl in by_level
                 else ("–", "–") for lvl in levels]  # fmt: skip
        rows.append(TipRow(tip, first["pendant"], first["sites"], cells))
    return rows


TIPS_KEY = "esr/esr_tips.tsv"


def tips_export(request: HttpRequest, job_id: str, fmt: str) -> HttpResponse:
    """The job page's table of held-out sequences as CSV or Excel: a row per sequence,
    each level's figures side by side, as esr_tips.tsv has them (not the page's text)."""
    job = get_object_or_404(Job.objects.for_user(request.user), id=job_id)
    if fmt not in ("csv", "xlsx") or not job.params.get("esr"):
        raise Http404
    if job.status != Job.Status.SUCCEEDED or TIPS_KEY not in manifest.paths(job, "esr/"):
        raise Http404
    text = storage.read_bytes(f"{job.storage_prefix}results/{TIPS_KEY}").decode()
    records = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    levels = [lvl for lvl in esr.LEVELS if any(r["level"] == lvl for r in records)]
    columns = ["Sequence", "Pendant branch length", "Sites"]
    for level in levels:
        columns += [
            f"{level.upper()} {name}" for name in ("accuracy", "mean posterior", "log score")
        ]
    by_tip: dict[str, dict[str, dict[str, str]]] = {}
    for record in records:
        by_tip.setdefault(record["tip"], {})[record["level"]] = record
    rows = []
    for tip, found in by_tip.items():
        first = next(iter(found.values()))
        row = [tip, first["pendant"], first["sites"]]
        for level in levels:
            r = found.get(level, {})
            row += [r.get("accuracy", ""), r.get("mean_pp", ""), r.get("log_score", "")]
        rows.append(row)
    stem = re.sub(r"[^A-Za-z0-9._+-]", "_", f"{bundle.folder_name(job)}_esr_sequences")
    if fmt == "csv":
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
        response = HttpResponse(out.getvalue(), content_type="text/csv; charset=utf-8")
    else:
        from phylobot import sites_excel  # openpyxl: only needed here

        plan = job.params.get("esr_of") or {}
        about = [("Job", job.name or str(job.id)), ("Job ID", str(job.id)),
                 ("Table", "Each sequence held out (ESR)"),
                 ("On", str(plan.get("label") or "an uploaded alignment and tree")),
                 ("Model", str(plan.get("model", "")))]  # fmt: skip
        response = HttpResponse(
            sites_excel.results_workbook(columns, rows, about, name="Sequences"),
            content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        )
    response["Content-Disposition"] = f'attachment; filename="{stem}.{fmt}"'
    return response


def context(job: Job, keys: set[str]) -> dict[str, Any]:
    """The job page's ESR section."""
    plan = job.params.get("esr_of") or {}
    out: dict[str, Any] = {
        "of": plan,
        "levels": [],
        "chart": None,
        "tips": [],
        "fit": None,
        "lengths": None,
        "level_help": LEVELS_POPOVER,
    }
    summary_key, tips_key = "esr/esr_summary.json", TIPS_KEY
    if job.status == Job.Status.SUCCEEDED and summary_key in keys:
        prefix = f"{job.storage_prefix}results/"
        summary = json.loads(storage.read_bytes(prefix + summary_key))
        levels = summary.get("levels") or {}
        out["levels"] = [
            {
                "key": key,
                "label": key.upper(),
                "slot": SERIES[key],
                **levels[key],
                "gap": levels[key]["mean_pp"] - levels[key]["accuracy"],
            }  # fmt: skip
            for key in esr.LEVELS
            if key in levels
        ]
        out["chart"] = chart(levels)
        out["fit"] = summary.get("fit")
        if tips_key in keys:
            shown = [level["key"] for level in out["levels"]]
            tips_text = storage.read_bytes(prefix + tips_key).decode()
            out["tips"] = tip_rows(tips_text, shown)
            out["lengths"] = esr_lengths.for_job(prefix, keys, tips_text, shown)
            if out["lengths"]:
                for level in out["levels"]:
                    level |= out["lengths"]["headline"].get(level["key"], {})
        out["files"] = [k for k in ("esr/esr_tips.tsv", "esr/esr_sites.tsv", "esr/tree.nwk",
                                    summary_key) if k in keys]  # fmt: skip
    return {"esr_run": out}


def describe(params: dict[str, Any]) -> list[str]:
    """The ESR job in words, for the job page's plan line."""
    plan = params.get("esr_of") or {}
    n = plan.get("tips") or 1
    levels = " and ".join(lvl.upper() for lvl in plan.get("levels") or ["e0"])
    if plan.get("kind") == "job":
        on = f"the alignment, tree and model of {plan.get('label', 'a job')}"
    else:
        on = f"the uploaded alignment and tree under {plan.get('model', '?')}"
    steps = [f"fit {on}"]
    if plan.get("fit_lengths"):
        steps[0] += " (the tree's branch lengths first)"
    steps.append(f"hold out {n} sequence{'s' * (n != 1)} and reconstruct each as an "
                 f"ancestor ({levels})")  # fmt: skip
    steps.append("score the posteriors against the real sequences")
    return steps
