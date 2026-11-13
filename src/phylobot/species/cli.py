"""Species for reconciliation, by hand (pipeline-only users) or in a pipeline step.

    phylobot-species resolve FASTA OUTDIR [--table species.tsv] [--tree species.nwk]
                             [--cache DIR]
    phylobot-species polytomies SPECIES_TREE SPECIES_TSV OUT_TREE OUT_REPORT
                             --gene-trees TREE... [--outgroup FILE] [--seed N]
                             [--method guided|arbitrary]

`resolve` reads each sequence's species from its FASTA header (or from `--table`, an
edited species.tsv) and resolves them against the Open Tree of Life, which receives the
species names and taxids only, never sequences; with `--tree` it uses the given species
tree instead and contacts nothing. It writes OUTDIR/species.tsv, species_tree.nwk and
provenance.json, and exits 1 if any sequence still needs a species: edit the `query`
(a name, or "ott<id>") or `taxid` of those rows in species.tsv and run it again with
`--table`.

`polytomies` resolves the species tree's polytomies, guided by gene trees (whose leaves
are the sequences of species.tsv) unless `--method arbitrary`.
"""

import argparse
import json
import sys
from pathlib import Path

from phylobot.pipeline.validate import InputError, parse_fasta
from phylobot.species import polytomies, resolve, speciestree, table
from phylobot.species.opentree import Cache, Client, OpenTreeError
from phylobot.species.polytomies import PolytomyError
from phylobot.species.speciestree import SpeciesTreeError
from phylobot.species.table import OK, TableError

ERRORS = (InputError, TableError, SpeciesTreeError, OpenTreeError, PolytomyError, ValueError)


def run_resolve(args: argparse.Namespace) -> int:
    records = parse_fasta(args.fasta.read_text(), str(args.fasta))
    if args.table:
        rows = table.read(args.table.read_text())
        listed = {row.sequence for row in rows}
        detected = resolve.detect([(n, d) for n, d, _ in records if n not in listed])
        names = {n for n, _, _ in records}
        rows = [row for row in rows if row.sequence in names] + detected
    else:
        rows = resolve.detect([(name, description) for name, description, _ in records])
    if args.tree:
        result = resolve.with_upload(rows, args.tree.read_text())
    else:
        client = Client(cache=Cache(args.cache) if args.cache else None)
        result = resolve.with_open_tree(rows, client)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "species.tsv").write_text(table.write(result.rows))
    (args.out / "provenance.json").write_text(json.dumps(result.provenance, indent=2) + "\n")
    tree_path = args.out / "species_tree.nwk"
    if result.tree is not None:
        tree_path.write_text(result.tree.newick() + "\n")
    elif tree_path.exists():
        tree_path.unlink()
    pending = [row for row in result.rows if row.status != OK]
    species = len(result.tree.leaves()) if result.tree else 0
    print(f"{len(result.rows) - len(pending)} of {len(result.rows)} sequences resolved "
          f"to {species} species; {len(result.tree.polytomies()) if result.tree else 0} "
          f"polytomies in the species tree")  # fmt: skip
    for row in pending:
        print(f"  {row.sequence}: {row.status}: {row.query or '-'}: {row.note}", file=sys.stderr)
    if pending or result.tree is None:
        print(f"Edit {args.out / 'species.tsv'} and run again with --table.", file=sys.stderr)
        return 1
    return 0


def run_polytomies(args: argparse.Namespace) -> int:
    rows = table.read(args.table.read_text())
    lacking = [row.sequence for row in rows if not row.species]
    if lacking:
        raise TableError(f"sequences without a species: {', '.join(lacking[:10])}")
    species_of = {row.sequence: row.species for row in rows}
    tree = speciestree.parse(args.species_tree.read_text())
    outgroup = set()
    if args.outgroup:
        outgroup = {line.strip() for line in args.outgroup.read_text().splitlines() if line.strip()}
    resolved, report = polytomies.resolve(
        tree,
        [path.read_text() for path in args.gene_trees],
        species_of,
        seed=args.seed,
        method=args.method,
        outgroup=outgroup,
    )
    args.out_tree.write_text(resolved.newick() + "\n")
    args.out_report.write_text(json.dumps(report, indent=2) + "\n")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    steps = parser.add_subparsers(dest="step", required=True)
    r = steps.add_parser("resolve")
    r.add_argument("fasta", type=Path)
    r.add_argument("out", type=Path)
    r.add_argument("--table", type=Path, help="an edited species.tsv")
    r.add_argument("--tree", type=Path, help="your own species tree (Newick)")
    r.add_argument("--cache", type=Path, help="directory for Open Tree answers")
    p = steps.add_parser("polytomies")
    p.add_argument("species_tree", type=Path)
    p.add_argument("table", type=Path)
    p.add_argument("out_tree", type=Path)
    p.add_argument("out_report", type=Path)
    p.add_argument("--gene-trees", type=Path, nargs="*", default=[])
    p.add_argument("--outgroup", type=Path, help="outgroup sequence names, one per line")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--method", choices=polytomies.METHODS, default=polytomies.GUIDED)
    args = parser.parse_args()
    try:
        code = run_resolve(args) if args.step == "resolve" else run_polytomies(args)
    except ERRORS as exc:
        print(f"phylobot-species: {exc}", file=sys.stderr)
        code = 1
    sys.exit(code)


if __name__ == "__main__":
    main()
