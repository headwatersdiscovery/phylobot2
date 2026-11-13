import re
import uuid
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Self

from django.conf import settings
from django.contrib.auth.models import AbstractBaseUser, AnonymousUser
from django.db import models
from django.db.models.functions import Lower


class UserLimits(models.Model):
    """Per-user overrides of the site defaults (site_settings.default_limits).
    Blank fields fall back to the default. Edited in the Django admin."""

    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="limits"
    )
    max_sequences = models.PositiveIntegerField(null=True, blank=True)
    max_length = models.PositiveIntegerField(null=True, blank=True, help_text="residues")
    max_upload_bytes = models.PositiveIntegerField(null=True, blank=True)
    max_concurrent_jobs = models.PositiveIntegerField(null=True, blank=True)

    class Meta:
        verbose_name_plural = "user limits"

    def __str__(self) -> str:
        return f"Limits for {self.user}"


class UserCalibration(models.Model):
    """A user's own calibration coefficients (phylobot.calibration), over the site's
    defaults (SiteSettings.calibration_*): a blank field uses the site's value. Edited on
    My calibration."""

    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="calibration"
    )
    intercept = models.FloatField(null=True, blank=True)
    branch_length = models.FloatField(null=True, blank=True)
    tip_distance = models.FloatField(null=True, blank=True)
    site_rate = models.FloatField(null=True, blank=True)
    entropy = models.FloatField(null=True, blank=True)
    gap_fraction = models.FloatField(null=True, blank=True)
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self) -> str:
        return f"Calibration for {self.user}"


@dataclass(frozen=True)
class Limits:
    max_sequences: int
    max_length: int
    max_upload_bytes: int
    max_concurrent_jobs: int


def limits_for(user: AbstractBaseUser | AnonymousUser) -> Limits:
    from phylobot import site_settings  # it imports this module

    values: dict[str, Any] = site_settings.default_limits()
    if not user.is_authenticated:
        return Limits(**values)
    overrides = UserLimits.objects.filter(user_id=user.pk).first()
    if overrides is not None:
        for field in values:
            if (value := getattr(overrides, field)) is not None:
                values[field] = value
    return Limits(**values)


TAG_PATTERN = re.compile(r"^[a-z0-9][a-z0-9_-]{0,39}$")


def parse_tags(text: str) -> list[str]:
    """'#Ime2, review  #draft' -> ['ime2', 'review', 'draft']. Raises ValueError naming
    any invalid tag. Tags are lowercase letters, digits, '-' and '_' (max 40)."""
    tags: list[str] = []
    for raw in re.split(r"[\s,]+", text.strip()):
        tag = raw.lstrip("#").lower()
        if not tag:
            continue
        if not TAG_PATTERN.match(tag):
            raise ValueError(
                f"#{tag} is not a valid tag: use letters, digits, '-' or '_' (up to 40)"
            )
        if tag not in tags:
            tags.append(tag)
    return tags


class Tag(models.Model):
    """A user's hashtag for organising jobs. Names are stored normalised (see parse_tags)."""

    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="tags"
    )
    name = models.CharField(max_length=40)

    class Meta:
        ordering = ["name"]
        constraints = [models.UniqueConstraint(fields=["owner", "name"], name="unique_tag")]

    def __str__(self) -> str:
        return f"#{self.name}"


class OwnedQuerySet[M: models.Model](models.QuerySet[M]):
    """Records a user owns (jobs, homolog searches). Deleting one keeps its record (the
    paper trail for usage and billing) with `deleted_at` set."""

    def for_user(self, user: AbstractBaseUser | AnonymousUser) -> Self:
        """The only way views should look up these records: a user sees only their own,
        and never deleted ones."""
        if not user.is_authenticated:
            return self.none()
        return self.filter(owner_id=user.pk, deleted_at__isnull=True)

    def deleted_for(self, user: AbstractBaseUser | AnonymousUser) -> Self:
        """A user's deleted records (their files are gone)."""
        if not user.is_authenticated:
            return self.none()
        return self.filter(owner_id=user.pk, deleted_at__isnull=False)


class JobQuerySet(OwnedQuerySet["Job"]):
    pass


class Job(models.Model):
    class Status(models.TextChoices):
        QUEUED = "queued"
        RUNNING = "running"
        SUCCEEDED = "succeeded"
        FAILED = "failed"
        CANCELLED = "cancelled"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="jobs"
    )
    name = models.CharField(max_length=100, blank=True)
    status = models.CharField(max_length=16, choices=Status.choices, default=Status.QUEUED)
    # Pipeline parameters (tool keys, labels, limits) as passed to Nextflow's params file.
    params = models.JSONField(default=dict)
    created_at = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    finished_at = models.DateTimeField(null=True, blank=True)
    heartbeat_at = models.DateTimeField(null=True, blank=True)
    worker = models.CharField(max_length=100, blank=True)
    progress = models.JSONField(default=dict, blank=True)
    # Shown to the user: validation messages only, never internal errors.
    error = models.TextField(blank=True)
    # Compute used: sum over pipeline steps of runtime x CPUs, over every attempt (a
    # resumed attempt counts only the steps it ran, not those it reused).
    cpu_seconds = models.FloatField(null=True, blank=True)
    # The basis for billing: each attempt's completed steps in full plus its failed steps
    # at the rate in force when it ended (JobAttempt.billable_cpu_seconds, summed).
    billable_cpu_seconds = models.FloatField(null=True, blank=True)
    # The current attempt (JobAttempt.number). Resuming a failed job starts the next one
    # in the same work directory; every write from a runner names its attempt, so a
    # superseded runner's late writes change nothing (phylobot.jobs).
    attempt = models.PositiveIntegerField(default=1)
    # Whether it ran from its own folder, with Nextflow's resume cache there: a resume
    # reuses finished steps only then (jobs from before resuming existed start over).
    resume_cache = models.BooleanField(default=False)
    # The kept local folder of a failed job (cache and every step's files): its size,
    # measured by the worker, and when the worker removed it after the Settings page's
    # retention period (a resume then starts from the beginning). Both clear on resume.
    work_bytes = models.BigIntegerField(null=True, blank=True)
    work_expired_at = models.DateTimeField(null=True, blank=True)
    # The attempt whose folder is archived in object storage (phylobot.work_archive), so a
    # resume survives losing the worker's disk, and the archive's size.
    work_archive_attempt = models.PositiveIntegerField(null=True, blank=True)
    work_archive_bytes = models.BigIntegerField(null=True, blank=True)
    # Work estimated at submission (phylobot.workload; null for older jobs): shown beside
    # cpu_seconds, and the ratio of the two calibrates time estimates for this server.
    work_units = models.FloatField(null=True, blank=True)
    # The input size that estimate used ({"sequences": n, "length": mean residues}), so
    # each stage's estimate can be set beside its CPU time (Work records); {} for older jobs.
    input_size = models.JSONField(default=dict, blank=True)
    # Bytes under storage_prefix, measured once the job ends (jobs.measure_storage; its
    # files no longer change); null until measured. Kept after the files are purged, as
    # the record of what the job held until `deleted_at` (My usage's storage history), so
    # what a job occupies *now* is this only while `purged_at` is null. Stored so pages
    # never list storage per request.
    storage_bytes = models.BigIntegerField(null=True, blank=True)
    # Deleting a job removes its files but keeps this record (the paper trail).
    deleted_at = models.DateTimeField(null=True, blank=True)
    purged_at = models.DateTimeField(null=True, blank=True, help_text="files removed")
    # When its ResultFile rows were last written in full (phylobot.manifest); null while
    # it runs, while its results go up, and for jobs from before the manifest existed.
    results_indexed_at = models.DateTimeField(null=True, blank=True)
    tags = models.ManyToManyField(Tag, blank=True, related_name="jobs")
    # "Duplicate and run again" (the New Job page, ?from=<job>) starts a new job from this
    # one's inputs and settings, leaving this one as it is; the new job points back here.
    rerun_of = models.ForeignKey(
        "self", null=True, blank=True, on_delete=models.SET_NULL, related_name="reruns"
    )
    # "Extend" (the New Job page, ?extend=<job>): a new job that builds on this one's
    # results (its alignments, kept under their aligners' keys) and only adds to them,
    # leaving this one as it is. A finished job's results never change (cited, exported).
    extends = models.ForeignKey(
        "self", null=True, blank=True, on_delete=models.SET_NULL, related_name="extensions"
    )
    # Imported from another PhyloBot's ZIP (phylobot.job_import): when, and what the ZIP
    # said about where it came from (its id there, PhyloBot version, export time, and the
    # CPU it used there). Its results were not computed here and a ZIP can be edited, so
    # it is badged wherever it appears, never resumed, and its CPU is never ours to count:
    # cpu_seconds and billable_cpu_seconds stay null.
    imported_at = models.DateTimeField(null=True, blank=True)
    imported_from = models.JSONField(default=dict, blank=True)

    objects = JobQuerySet.as_manager()

    class Meta:
        ordering = ["-created_at"]
        indexes = [models.Index(fields=["status", "created_at"])]

    def __str__(self) -> str:
        return f"{self.name or 'Job'} {self.id} ({self.status})"

    @property
    def storage_prefix(self) -> str:
        return f"jobs/{self.id}/"

    @property
    def exported_at(self) -> datetime | None:
        """When an imported job's ZIP was made, as the ZIP says (shown on the job page)."""
        try:
            return datetime.fromisoformat(self.imported_from.get("exported_at") or "")
        except AttributeError, TypeError, ValueError:
            return None

    def add_tags(self, names: list[str]) -> None:
        tags = [Tag.objects.get_or_create(owner_id=self.owner_id, name=n)[0] for n in names]
        self.tags.add(*tags)

    def remove_tags(self, names: list[str]) -> None:
        self.tags.remove(*self.tags.filter(name__in=names))


class ResultFile(models.Model):
    """One file under a job's results/ in object storage (phylobot.manifest): what pages
    read instead of listing storage, which takes 32 s on a job of 524,069 files."""

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="result_files")
    path = models.TextField()  # relative to results/
    size = models.BigIntegerField()

    class Meta:
        constraints = [models.UniqueConstraint(fields=["job", "path"], name="unique_result_file")]

    def __str__(self) -> str:
        return f"{self.job_id} {self.path}"


class AncestorName(models.Model):
    """A user's name for an ancestor (phylobot.ancestor_names): it belongs to the clade
    (`ancestor_map.tsv`'s content-derived ID), so it shows on every combination of the
    job that has that ancestor. User text: never in pipeline files or paths, applied
    only where the web tier displays and exports."""

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="ancestor_names")
    # The map's matching key: the clade ID with an outgroup, else the parts ID.
    clade = models.CharField(max_length=16)
    # The parts ID of the node it was named on: where another combination's node has the
    # same clade but other branching inside, the name is marked as approximate.
    parts = models.CharField(max_length=16)
    name = models.CharField(max_length=40)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["job", "clade"], name="unique_ancestor_clade"),
            models.UniqueConstraint(Lower("name"), "job", name="unique_ancestor_name"),
        ]

    def __str__(self) -> str:
        return f"{self.name} ({self.job_id} {self.clade})"


class AncestorNote(models.Model):
    """A user's description of an ancestor (phylobot.ancestor_names): free text, keyed by
    the clade as AncestorName is, so every combination with the same descendants shares
    it. Apart from the name so a description needs no name, and an empty name can never
    reach a label or a FASTA id."""

    LIMIT = 20_000

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="ancestor_notes")
    clade = models.CharField(max_length=16)
    parts = models.CharField(max_length=16)
    text = models.TextField()
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["job", "clade"], name="unique_ancestor_note_clade"),
        ]

    def __str__(self) -> str:
        return f"note on {self.job_id} {self.clade}"


class CladeCollapse(models.Model):
    """A clade the user folded into a triangle in the tree viewer (phylobot.clade_collapse).
    Like AncestorName it belongs to the clade (`ancestor_map.tsv`'s content-derived ID),
    not one tree's node label, so the fold, its label and its triangle style follow the
    clade onto every combination of the job that has it, and persist across page loads.
    Its presence means collapsed; removing the row expands the clade again. User text and
    colours are applied only where the web tier draws the tree, never in pipeline files."""

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="clade_collapses")
    # The map's matching key: the clade ID with an outgroup, else the parts ID.
    clade = models.CharField(max_length=16)
    # The parts ID of the node it was folded on (as AncestorName): a combination whose node
    # has the same clade but other branching inside is folded all the same.
    parts = models.CharField(max_length=16)
    label = models.CharField(max_length=60, blank=True, default="")  # drawn beside the triangle
    fill = models.CharField(max_length=16, default="#cfd8e6")  # hex or "none"
    border_colour = models.CharField(max_length=16, default="#2457a6")
    border_width = models.FloatField(default=1.6)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(fields=["job", "clade"], name="unique_clade_collapse"),
        ]

    def __str__(self) -> str:
        return f"collapse {self.clade} ({self.job_id})"


class AlignmentTracks(models.Model):
    """The tracks a user last chose over one of a job's alignments in the alignment viewer
    (alignment_views.track_list): adding, removing or reordering a track writes the list,
    and the viewer opens that alignment with it next time. An empty list is a choice too
    (no tracks). A display preference only: it is not part of a job's ZIP."""

    owner = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)
    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="alignment_tracks")
    aligner = models.CharField(max_length=64)  # the alignment's key in the job
    tracks = models.JSONField(default=list)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["owner", "job", "aligner"], name="unique_alignment_tracks"
            ),
        ]

    def __str__(self) -> str:
        return f"tracks of {self.aligner} ({self.job_id})"


class SavedPath(models.Model):
    """A mutational path the user pinned (phylobot.saved_paths), to come back to without
    tracing it again: the combination it was traced in, its two ends as that tree names
    them, and the view it was saved with (the path page's own query, re-checked against
    the job every time it is opened). A finished job's results never change, so the names
    keep meaning the same nodes. `name` is user text, shown only on pages."""

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="saved_paths")
    name = models.CharField(max_length=80)
    combination = models.CharField(max_length=255)
    start = models.CharField(max_length=255)  # the ancestor it starts from
    end = models.CharField(max_length=255)  # an ancestor or an extant sequence below it
    view = models.JSONField(default=dict)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ["name", "id"]
        constraints = [
            models.UniqueConstraint(
                fields=["job", "combination", "start", "end"], name="unique_saved_path"
            ),
        ]

    def __str__(self) -> str:
        return f"{self.name} ({self.job_id} {self.start} to {self.end})"


class Structure(models.Model):
    """A protein structure attached to a job (phylobot.structures): a PDB
    entry, an AlphaFold model or the user's own file, to see the job's ancestors on. The
    file (PDB format or mmCIF) is in object storage at `storage_key`; it is never part of a job's
    results, and nothing in the pipeline reads it. Kept on `ancestor_names.holder(job)`,
    as ancestor names are, so a job and its extensions share them. `label` and `title`
    are text from the user or from a database: shown on pages only."""

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="structures")
    label = models.CharField(max_length=120)
    source = models.CharField(max_length=16)  # structures.remote.SOURCES
    accession = models.CharField(max_length=32, blank=True)  # PDB ID or UniProt accession
    title = models.CharField(max_length=300, blank=True)
    chains = models.JSONField(default=list)  # [{"id": "A", "residues": 450}]
    size = models.PositiveIntegerField(default=0)  # bytes
    format = models.CharField(max_length=3, default="pdb")  # structures.pdb.FORMATS
    # A predicted model (AlphaFold's, or an upload the user says is one): its B-factor
    # column holds pLDDT, and is coloured in pLDDT's bands.
    predicted = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["created_at", "id"]

    def __str__(self) -> str:
        return f"{self.label} ({self.job_id})"

    @property
    def storage_key(self) -> str:
        return f"{self.job.storage_prefix}structures/{self.id}.{self.format}"


class JobAttempt(models.Model):
    """One run of a job: the first, and each resume (phylobot.jobs.resume). The Job row
    shows the current attempt; these keep every attempt's times, outcome and compute,
    which are the record for billing."""

    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="attempts")
    number = models.PositiveIntegerField()
    worker = models.CharField(max_length=100, blank=True)
    started_at = models.DateTimeField(null=True, blank=True)
    finished_at = models.DateTimeField(null=True, blank=True)
    status = models.CharField(max_length=16, choices=Job.Status.choices, blank=True)
    error = models.TextField(blank=True)
    # From the attempt's trace: steps it ran to completion, steps that failed or were
    # aborted, and how many it reused from earlier attempts (charged there, not here).
    completed_cpu_seconds = models.FloatField(null=True, blank=True)
    failed_cpu_seconds = models.FloatField(null=True, blank=True)
    tasks_run = models.PositiveIntegerField(null=True, blank=True)
    tasks_reused = models.PositiveIntegerField(null=True, blank=True)
    # The admin's rate for failed steps when the attempt's usage was recorded (percent of
    # the normal rate), fixed then so a later change never re-prices past work.
    failed_rate_percent = models.PositiveIntegerField(null=True, blank=True)
    billable_cpu_seconds = models.FloatField(null=True, blank=True)
    # The share of the machine Nextflow was given (phylobot.machine), set at launch.
    cpus = models.PositiveIntegerField(null=True, blank=True)
    memory_bytes = models.BigIntegerField(null=True, blank=True)

    class Meta:
        ordering = ["job", "number"]
        constraints = [models.UniqueConstraint(fields=["job", "number"], name="unique_job_attempt")]

    def __str__(self) -> str:
        return f"{self.job_id} attempt {self.number} ({self.status or 'running'})"


class HomologSearch(models.Model):
    """A homolog search at a remote service (phylobot.homologs), followed step by step
    by phylobot.searches (run_searches); not a pipeline Job: it uses no compute here and
    mostly waits. Metadata only: the query and the service's raw result live in object
    storage (storage_prefix). Like jobs, deleting one keeps the record."""

    class Status(models.TextChoices):
        QUEUED = "queued"
        RUNNING = "running"  # submitted; waiting on the service
        SUCCEEDED = "succeeded"
        FAILED = "failed"
        CANCELLED = "cancelled"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="homolog_searches"
    )
    name = models.CharField(max_length=100, blank=True)
    status = models.CharField(max_length=16, choices=Status.choices, default=Status.QUEUED)
    service = models.CharField(max_length=8)  # phylobot.homologs.sources key
    database = models.CharField(max_length=40)  # the service's database name
    query_length = models.PositiveIntegerField()
    remote_id = models.CharField(max_length=100, blank=True, help_text="the service's job ID")
    created_at = models.DateTimeField(auto_now_add=True)
    submitted_at = models.DateTimeField(null=True, blank=True)
    finished_at = models.DateTimeField(null=True, blank=True)
    # When the runner next looks at it; claiming pushes it forward (a lease), so a runner
    # that dies leaves the search to be picked up again.
    next_check_at = models.DateTimeField(null=True, blank=True)
    check_failures = models.PositiveSmallIntegerField(default=0)
    last_checked_at = models.DateTimeField(null=True, blank=True)  # shown while it waits
    # Shown to the user.
    error = models.TextField(blank=True)
    # From the finished search.
    program = models.CharField(max_length=40, blank=True)
    database_release = models.CharField(max_length=40, blank=True)
    database_sequences = models.PositiveBigIntegerField(null=True, blank=True)
    hit_count = models.PositiveIntegerField(null=True, blank=True)
    truncated = models.BooleanField(null=True, blank=True, help_text="hit the hit limit")
    # A duplicate (searches.duplicate) copies a finished search's files; nothing is sent
    # to the service again, so copies count in neither the daily limit nor the timings.
    copied_from = models.ForeignKey(
        "self", null=True, blank=True, on_delete=models.SET_NULL, related_name="copies"
    )
    copied_from_example = models.CharField(max_length=40, blank=True)
    deleted_at = models.DateTimeField(null=True, blank=True)
    purged_at = models.DateTimeField(null=True, blank=True, help_text="files removed")

    objects = OwnedQuerySet["HomologSearch"].as_manager()

    class Meta:
        ordering = ["-created_at"]
        indexes = [models.Index(fields=["status", "next_check_at"])]

    def __str__(self) -> str:
        return f"{self.name or 'Search'} {self.id} ({self.status})"

    @property
    def is_copy(self) -> bool:
        return self.copied_from_id is not None or bool(self.copied_from_example)

    @property
    def storage_prefix(self) -> str:
        return f"searches/{self.id}/"

    @property
    def query_key(self) -> str:
        return f"{self.storage_prefix}query.fasta"

    @property
    def result_key(self) -> str:
        return f"{self.storage_prefix}result.xml"

    @property
    def sequences_key(self) -> str:
        """The hits' full-length sequences (BLAST reports only the aligned stretch)."""
        return f"{self.storage_prefix}sequences.fasta"


class Library(models.Model):
    """An ancestral library design over one of the user's jobs (phylobot.library): the
    targets collected from the viewers, and later the axes, budget and strategy. Not a
    Job: it uses no compute, only Python over files the pipeline already wrote (the
    reasoning that made HomologSearch not a Job either). Members are not stored: the
    job's files never change and the generator is deterministic, so `spec` reproduces
    them. Editable for as long as the user likes; each export leaves a LibraryExport
    receipt instead. Like jobs, deleting one keeps the record."""

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="libraries"
    )
    job = models.ForeignKey(Job, on_delete=models.CASCADE, related_name="libraries")
    name = models.CharField(max_length=100)
    # The first segment of every member name (phylobot.library.member_name).
    code = models.CharField(max_length=6)
    # The user's free text about the library (a description, notes, anything pasted).
    NOTES_LIMIT = 20_000
    notes = models.TextField(blank=True, default="")
    spec = models.JSONField(default=dict)
    spec_version = models.PositiveSmallIntegerField(default=1)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    deleted_at = models.DateTimeField(null=True, blank=True)

    objects = OwnedQuerySet["Library"].as_manager()

    class Meta:
        ordering = ["-updated_at"]
        verbose_name_plural = "libraries"
        indexes = [models.Index(fields=["owner", "job", "updated_at"])]

    def __str__(self) -> str:
        return f"{self.name} ({self.code}) on {self.job_id}"

    @property
    def targets(self) -> list[dict[str, Any]]:
        from phylobot import library

        return library.spec_targets(self.spec)

    @property
    def storage_prefix(self) -> str:
        """Exports live under the library's own prefix: the web tier never writes into a
        job's results."""
        return f"libraries/{self.id}/"


class LibraryExport(models.Model):
    """The receipt of one download of a library: what was exported, from which spec,
    with which PhyloBot. Never edited. Its files (FASTA, manifest, names.tsv, as
    downloaded) live at `storage_prefix`; its spec can be restored into the library."""

    library = models.ForeignKey(Library, on_delete=models.CASCADE, related_name="exports")
    number = models.PositiveIntegerField()  # 1, 2, ... per library
    created_at = models.DateTimeField(auto_now_add=True)
    spec = models.JSONField()  # a copy, as it was
    spec_version = models.PositiveSmallIntegerField()
    member_count = models.PositiveIntegerField()
    phylobot_version = models.CharField(max_length=40)

    class Meta:
        ordering = ["library", "-number"]
        constraints = [
            models.UniqueConstraint(fields=["library", "number"], name="unique_export_number")
        ]

    def __str__(self) -> str:
        return f"{self.library.name} export {self.number}"

    @property
    def storage_prefix(self) -> str:
        return f"{self.library.storage_prefix}exports/{self.number}/"


class AlignmentQuerySet(OwnedQuerySet["Alignment"]):
    def for_user(self, user: AbstractBaseUser | AnonymousUser) -> Self:
        """As for other records, and a job's alignments go with the job when it is
        deleted (its files are gone)."""
        return (
            super()
            .for_user(user)
            .filter(models.Q(job__isnull=True) | models.Q(job__deleted_at__isnull=True))
        )


class Alignment(models.Model):
    """An alignment in the user's workspace ("My alignments", phylobot.alignments):
    uploaded, saved from the alignment editor, or published by one of their jobs.

    The editor's "Save a copy" makes a new alignment (`derived_from`, with the edits that
    made it); "Save (overwrite)" replaces a saved alignment's rows and bumps `version`
    (jobs that used it keep their own copy of their inputs). Saved alignments live under
    their own prefix; a job's are registered in place (`job`, `path`), are read-only, and
    leave with the job. Like jobs, deleting one keeps the record."""

    class Origin(models.TextChoices):
        UPLOAD = "upload"
        EDIT = "edit"
        JOB = "job"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    owner = models.ForeignKey(
        settings.AUTH_USER_MODEL, on_delete=models.CASCADE, related_name="alignments"
    )
    # A job's alignment shows as "<job name> · <name>", the name being the aligner's label.
    name = models.CharField(max_length=100)
    origin = models.CharField(max_length=8, choices=Origin.choices)
    seqtype = models.CharField(max_length=8, default="protein")  # "protein" or "codon"
    genetic_code = models.PositiveSmallIntegerField(default=1)  # NCBI number (codons)
    sequences = models.PositiveIntegerField()
    columns = models.PositiveIntegerField(help_text="sites: codons for a codon alignment")
    job = models.ForeignKey(
        Job, null=True, blank=True, on_delete=models.CASCADE, related_name="alignments"
    )
    job_key = models.CharField(max_length=40, blank=True)  # the job's alignment key (mafft)
    path = models.CharField(max_length=200, blank=True)  # under the job's results/
    filename = models.CharField(max_length=200, blank=True)  # an upload's file name
    derived_from = models.ForeignKey(
        "self", null=True, blank=True, on_delete=models.SET_NULL, related_name="derivatives"
    )
    # How this came from `derived_from` (or, for an overwritten upload, from the uploaded
    # file): the columns kept, sequences removed, row order (alignments.provenance and
    # alignments.compose), and the parent's version then.
    edits = models.JSONField(default=dict, blank=True)
    # 1, then +1 at each overwrite: caches, undo histories and copies' provenance name it.
    version = models.PositiveIntegerField(default=1)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    # Opened in the editor, downloaded or used in a job (alignments.touch): the New Job
    # page's picker lists the most recently used first.
    last_used_at = models.DateTimeField(null=True, blank=True)
    deleted_at = models.DateTimeField(null=True, blank=True)
    # Simulated data (a simulation's replicate, or anything saved from one): permanent,
    # never cleared, and badged wherever the alignment appears.
    simulated = models.BooleanField(default=False)

    objects = AlignmentQuerySet.as_manager()

    class Meta:
        ordering = ["-created_at"]
        indexes = [models.Index(fields=["owner", "created_at"])]
        constraints = [
            models.UniqueConstraint(
                fields=["job", "job_key"],
                condition=models.Q(job__isnull=False),
                name="unique_job_alignment",
            )
        ]

    def __str__(self) -> str:
        return self.title

    @property
    def from_job(self) -> bool:
        return self.job_id is not None

    @property
    def title(self) -> str:
        if self.job is not None:
            return f"{self.job.name or 'Job'} · {self.name}"
        return self.name

    @property
    def storage_prefix(self) -> str:
        return f"alignments/{self.id}/"

    @property
    def storage_key(self) -> str:
        if self.job is not None:
            return f"{self.job.storage_prefix}results/{self.path}"
        return f"{self.storage_prefix}alignment.fasta"


class WorkerMachine(models.Model):
    """What a worker host found when it last looked (phylobot.machine.detect): shown on
    the Settings page, which cannot ask Docker itself. One row per host."""

    host = models.CharField(max_length=255, unique=True)
    cpus = models.PositiveIntegerField()  # usable: isolated ones (or the cpuset) applied
    isolated = models.CharField(max_length=255, blank=True)
    memory_bytes = models.BigIntegerField()
    note = models.CharField(max_length=255, blank=True)
    detected_at = models.DateTimeField()

    def __str__(self) -> str:
        return f"{self.host}: {self.cpus} CPUs, {self.memory_bytes / 1024**3:.1f} GB"


class SiteSettings(models.Model):
    """Site-wide settings edited on the Settings page (phylobot.site_settings): one row.
    A blank field means "use the environment variable or built-in default"; a value
    here overrides it, and takes effect without restarting anything."""

    search_email = models.EmailField(blank=True)
    searches_per_day = models.PositiveIntegerField(null=True, blank=True)
    searches_in_flight_per_user = models.PositiveIntegerField(null=True, blank=True)
    search_time_limit_minutes = models.PositiveIntegerField(null=True, blank=True)
    max_sequences = models.PositiveIntegerField(null=True, blank=True)
    max_length = models.PositiveIntegerField(null=True, blank=True)
    max_upload_megabytes = models.PositiveIntegerField(null=True, blank=True)
    max_concurrent_jobs = models.PositiveIntegerField(null=True, blank=True)
    job_time_limit_hours = models.PositiveIntegerField(null=True, blank=True)
    worker_concurrency = models.PositiveIntegerField(null=True, blank=True)
    # The machine as PhyloBot treats it, over what the worker detects (phylobot.machine).
    machine_cpus = models.PositiveIntegerField(null=True, blank=True)
    machine_memory_gb = models.FloatField(null=True, blank=True)
    reserved_memory_gb = models.FloatField(null=True, blank=True)
    task_memory_gb = models.FloatField(null=True, blank=True)
    failed_task_rate_percent = models.PositiveIntegerField(null=True, blank=True)
    keep_failed_work_days = models.PositiveIntegerField(null=True, blank=True)
    archive_failed_work = models.BooleanField(null=True, blank=True)
    simulation_max_replicates = models.PositiveIntegerField(null=True, blank=True)
    simulation_max_megaresidues = models.PositiveIntegerField(null=True, blank=True)
    import_max_megabytes = models.PositiveIntegerField(null=True, blank=True)
    import_max_unpacked_megabytes = models.PositiveIntegerField(null=True, blank=True)
    import_max_files = models.PositiveIntegerField(null=True, blank=True)
    import_max_record_megabytes = models.PositiveIntegerField(null=True, blank=True)
    import_max_libraries = models.PositiveIntegerField(null=True, blank=True)
    import_max_exports = models.PositiveIntegerField(null=True, blank=True)
    # Calibration coefficients (phylobot.calibration): ln T = intercept + sum(beta x).
    calibration_intercept = models.FloatField(null=True, blank=True)
    calibration_branch_length = models.FloatField(null=True, blank=True)
    calibration_tip_distance = models.FloatField(null=True, blank=True)
    calibration_site_rate = models.FloatField(null=True, blank=True)
    calibration_entropy = models.FloatField(null=True, blank=True)
    calibration_gap_fraction = models.FloatField(null=True, blank=True)
    updated_at = models.DateTimeField(auto_now=True)
    updated_by = models.ForeignKey(
        settings.AUTH_USER_MODEL, null=True, blank=True, on_delete=models.SET_NULL, related_name="+"
    )

    class Meta:
        verbose_name_plural = "site settings"

    def __str__(self) -> str:
        return "Site settings"
