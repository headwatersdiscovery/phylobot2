"""Submit a job from the command line: python manage.py submit_job SEQUENCES.fasta [options]

Creates a Job owned by the given user (default: the local user in local mode), stores
the input in object storage and queues it. Prints the job ID. Options map to pipeline
parameters; the worker passes only allowlisted ones to Nextflow.
"""

from pathlib import Path
from typing import Any

from django.conf import settings
from django.contrib.auth import get_user_model
from django.contrib.auth.models import User
from django.core.management.base import BaseCommand, CommandError, CommandParser
from django.db import transaction

from phylobot import catalog, gaps, jobs, storage, workload
from phylobot.models import Job
from phylobot.modes import local_user


class Command(BaseCommand):
    help = "Queue a PhyloBot job for a protein FASTA file."

    def add_arguments(self, parser: CommandParser) -> None:
        parser.add_argument("fasta", type=Path)
        parser.add_argument("--name", default="")
        parser.add_argument("--user", help="owner username (default: local user, local mode)")
        parser.add_argument("--aligners", default="mafft")
        parser.add_argument("--tree-builders", default="iqtree")
        parser.add_argument("--models", default="", help="comma-separated catalog models")
        parser.add_argument("--fixed-models", action="store_true")
        parser.add_argument("--no-modelfinder", action="store_true")
        parser.add_argument(
            "--gap-methods", help="comma-separated; default fitch with an outgroup, else parsimony"
        )
        parser.add_argument("--outgroup", default="", help="comma-separated sequence names")
        parser.add_argument(
            "--stop-after", choices=["alignments", "trees", "ancestors"], default="ancestors"
        )
        parser.add_argument(
            "--species-map", type=Path, help="species.tsv, for reconciliation with a species tree"
        )
        parser.add_argument(
            "--species-tree", type=Path, help="Newick species tree (with --species-map)"
        )
        parser.add_argument("--polytomies", choices=["guided", "arbitrary"], default="guided")

    def handle(self, *args: Any, **options: Any) -> None:
        fasta: Path = options["fasta"]
        if not fasta.is_file():
            raise CommandError(f"no such file: {fasta}")
        owner: User
        if options["user"]:
            owner = get_user_model().objects.get(username=options["user"])
        elif settings.PHYLOBOT_MODE == "local":
            owner = local_user()
        else:
            raise CommandError("--user is required in cloud mode")
        models = [m for m in options["models"].split(",") if m]
        unknown = [m for m in models if not catalog.is_model(m)]
        if unknown:
            raise CommandError(f"unknown models: {', '.join(unknown)}")
        outgroup = [n for n in options["outgroup"].split(",") if n]
        species = [options["species_map"], options["species_tree"]]
        if any(species) and not all(species):
            raise CommandError("reconciliation needs both --species-map and --species-tree")
        for path in species:
            if path is not None and not path.is_file():
                raise CommandError(f"no such file: {path}")
        params = {
            "stop_after": options["stop_after"],
            "aligners": options["aligners"],
            "tree_builders": options["tree_builders"],
            "modelfinder": not options["no_modelfinder"],
            "fixed_models": options["fixed_models"],
            "models": models or None,
            "outgroup": outgroup,
            "gap_methods": (
                gaps.default_methods(bool(outgroup))
                if options["gap_methods"] is None
                else [m for m in options["gap_methods"].split(",") if m]
            ),
        }
        if all(species):  # checked by the pipeline's VALIDATE, like every input
            params |= {"reconcile": True, "polytomies": options["polytomies"]}
        size = workload.size_of(fasta.read_bytes()) or workload.IME2_SIZE
        with transaction.atomic():
            job = Job.objects.create(
                owner=owner,
                name=options["name"],
                params=params,
                work_units=workload.estimate(params, size).units,
                input_size=workload.size_record(size),
            )
            storage.put_bytes(f"{job.storage_prefix}inputs/sequences.fasta", fasta.read_bytes())
            if all(species):
                inputs = f"{job.storage_prefix}inputs/"
                storage.put_bytes(inputs + "species.tsv", options["species_map"].read_bytes())
                storage.put_bytes(inputs + "species_tree.nwk", options["species_tree"].read_bytes())
            transaction.on_commit(lambda: jobs.enqueue(job))
        self.stdout.write(str(job.id))
