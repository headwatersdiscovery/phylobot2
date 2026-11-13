"""Tier-0 property posteriors for one combination's ancestors.

Reads a finished combination's published files and, for every ancestral node and gap
version, draws from the per-column posteriors and writes each property's values per draw
plus a summary. Nothing here needs a structure, a network or a GPU, which is why tier 0
can run inline on every job.

    phylobot-properties <combination dir> <out dir> --properties charge,moment
                        [--draws 1000] [--seed 1] [--gap-versions fitch,parsimony]
                        [--nodes Node1,Node4]

`--nodes` measures only those ancestral nodes. Unset means every node in the staged
files, which is what the pipeline's inline run wants; a measurement asked for after the
job has finished names the nodes it was asked for, so a scope of
one ancestor costs one ancestor.

Output, under `<out>/` (`phylobot.property_store`): draws.parquet (every draw of every
node and gap version, a column per property), summaries.json (each node's n, mean, sd,
interval and its width, f(ML), the shift, cv) and index.json.

Per-draw values are kept rather than only the summary, so a later viewer can pivot
across nodes, combinations and properties without recomputing, and
so a figure can be rebuilt from what was published.
"""

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from phylobot import ancestors as ancestors_mod
from phylobot import library, properties, property_store
from phylobot.seqio import read_fasta


def node_distributions(combination: Path, node: str) -> list[list[tuple[str, float]]]:
    """`node`'s per-column distributions, read through the table's byte index where it
    has one: reading the whole table once per node would be the dominant cost.
    """
    whole = combination / "ancestors_posteriors.tsv"
    if not whole.exists():
        return []
    index = combination / ancestors_mod.index_name(whole.name)
    if index.exists():
        found = ancestors_mod.read_index(index.read_text())
        text = found.text(node, ancestors_mod.file_range(whole))
        return ancestors_mod.distributions(text, node) if text else []
    return ancestors_mod.distributions(whole.read_text(), node)


def versions_present(combination: Path, chosen: list[str]) -> list[str]:
    """The gap versions to measure: those the user asked for that this combination has.

    PRANK's and GRASP's reconstructions carry no posteriors, so no distribution can be
    drawn from them and they are skipped here, as the library's targets skip them.
    """
    out = []
    for version in [ancestors_mod.ORIGINAL, *chosen]:
        if version in ("prank", "grasp") or version not in ancestors_mod.VERSIONS:
            continue
        if (combination / ancestors_mod.VERSIONS[version][1]).exists() and version not in out:
            out.append(version)
    return out


def chosen_properties(keys: list[str]) -> list[properties.Property]:
    unknown = [k for k in keys if k not in properties.REGISTRY]
    if unknown:
        known = ", ".join(sorted(properties.REGISTRY))
        raise SystemExit(f"unknown properties: {', '.join(unknown)} (known: {known})")
    return [properties.REGISTRY[k] for k in keys]


def measure_node(
    target: library.Target,
    props: list[properties.Property],
    draws: int,
    seed: str,
) -> property_store.Measured:
    """One node and gap version: every property's draws and their summary."""
    values = properties.draw_values(target, props, draws, seed)
    summaries = {}
    for prop in props:
        rows = values[prop.key]
        summary = properties.summarise(rows, ml=prop.fn(target.template))
        summaries[prop.key] = {
            "label": prop.label,
            "units": prop.units,
            **({"method": prop.method} if prop.method else {}),
            **{
                field: round(number, 6) if isinstance(number, float) else number
                for field, number in asdict(summary).items()
            },
            "width": round(summary.width, 6),
            "shift": round(summary.shift, 6),
            "standardised_shift": round(summary.standardised_shift, 6),
            "cv": round(summary.cv, 6),
        }
    version = target.source.gap or ancestors_mod.ORIGINAL
    summary_of = {
        "node": target.source.node,
        "gaps": version,
        "draws": draws,
        "seed": seed,
        # Site-independent draws: stated here as well as in the UI, since a figure made
        # from these files must carry the caveat.
        "sampling": "per-column marginal posteriors, drawn independently",
        "weighted": False,
        "properties": summaries,
    }
    # As written to text before (6 significant digits), so stored values did not change.
    kept = {key: [float(f"{v:.6g}") for v in rows] for key, rows in values.items()}
    return property_store.Measured(version, target.source.node, kept, summary_of)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("combination", type=Path, help="a published <aligner>/<tree>/<model>")
    parser.add_argument("out", type=Path, help="where to write properties/")
    parser.add_argument("--properties", required=True, help="comma-separated property keys")
    parser.add_argument("--draws", type=int, default=properties.DEFAULT_DRAWS)
    parser.add_argument("--seed", default="1")
    parser.add_argument("--gap-versions", default="", help="comma-separated gap methods run")
    parser.add_argument("--nodes", default="", help="comma-separated nodes; unset means all")
    args = parser.parse_args(argv)

    keys = [k for k in args.properties.split(",") if k]
    if not keys:
        raise SystemExit("no properties asked for")
    props = chosen_properties(keys)
    if args.draws < 1:
        raise SystemExit(f"--draws must be at least 1, got {args.draws}")

    wanted = {n for n in args.nodes.split(",") if n}
    chosen_versions = [v for v in args.gap_versions.split(",") if v]
    versions = versions_present(args.combination, chosen_versions)
    written = 0
    measured: set[str] = set()
    results: list[property_store.Measured] = []
    for version in versions:
        rows = read_fasta(args.combination / ancestors_mod.VERSIONS[version][1])
        for node, template in sorted(rows.items()):
            if wanted and node not in wanted:
                continue
            dists = node_distributions(args.combination, node)
            if not dists or len(dists) != len(template):
                continue  # an extant row, or a reconstruction without posteriors
            source = library.Source(
                combination=0,
                node=node,
                gap="" if version == ancestors_mod.ORIGINAL else version,
            )
            target = library.Target(source=source, template=template, dists=dists)
            results.append(measure_node(target, props, args.draws, f"{args.seed}:{node}:{version}"))
            written += 1
            measured.add(node)
    missing = sorted(wanted - measured)
    if missing:
        # Not an error: the rest were measured, and a run that asked for one ancestor
        # too many is better finished than failed. It is said out loud, though, because
        # a node this job cannot see means the caller's plan and the staged files
        # disagree.
        print(f"asked for nodes not found here: {', '.join(missing)}", file=sys.stderr)
    if not written:
        # Not an error: a trees-only job or a combination without posteriors has nothing
        # to measure, and the workflow should not fail for it.
        print("no ancestral nodes with posteriors found", file=sys.stderr)
    args.out.mkdir(parents=True, exist_ok=True)
    if results:
        property_store.write(args.out, results, [p.key for p in props], args.draws)
    (args.out / property_store.INDEX).write_text(
        json.dumps(
            {
                "properties": [p.key for p in props],
                "gap_versions": versions,
                "nodes": written // max(len(versions), 1),
                "asked": sorted(wanted),
                "draws": args.draws,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
