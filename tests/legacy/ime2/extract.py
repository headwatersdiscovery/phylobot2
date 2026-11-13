"""Extract the 2014 PhyloBot Ime2 run from its SQLite database into plain files.

One-time provenance script; the extracted files are committed. Re-run with:

    uv run python tests/legacy/ime2/extract.py \
        ../phylobot-django/phylobot/static/media/anclibs/asr.db

Control1.IME2 (an 8-residue test sequence, YQLIRAPE) is dropped from sequences,
alignments and ancestors, but kept in the trees, which are copied verbatim.
"""

import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
CONTROL = "Control1.IME2"
MODEL_NAMES = {"PROTGAMMALG": "LG", "PROTGAMMAWAG": "WAG"}


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def main(db_path: str) -> None:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    aligners = dict(con.execute("SELECT id, name FROM AlignmentMethods"))
    models = {k: MODEL_NAMES[v] for k, v in con.execute("SELECT modelid, name FROM PhyloModels")}

    seqs = con.execute(
        "SELECT t.fullname, o.sequence FROM Taxa t JOIN OriginalSequences o ON o.taxonid = t.id "
        "WHERE t.fullname != ? ORDER BY t.id",
        (CONTROL,),
    ).fetchall()
    write_fasta(OUT / "sequences.fasta", seqs)

    for al_id, al_name in aligners.items():
        rows = con.execute(
            "SELECT t.fullname, s.alsequence FROM AlignedSequences s "
            "JOIN Taxa t ON t.id = s.taxonid WHERE s.almethod = ? AND t.fullname != ? "
            "ORDER BY t.id",
            (al_id, CONTROL),
        ).fetchall()
        write_fasta(OUT / "alignments" / f"{al_name}.fasta", rows)

    trees = con.execute(
        "SELECT u.almethod, u.phylomodelid, u.newick, c.newick FROM UnsupportedMlPhylogenies u "
        "JOIN AncestralCladogram c ON c.unsupportedmltreeid = u.id"
    ).fetchall()
    for al_id, model_id, ml_newick, cladogram in trees:
        stem = f"{aligners[al_id]}_{models[model_id]}"
        (OUT / "trees").mkdir(exist_ok=True)
        (OUT / "trees" / f"{stem}.ml.nwk").write_text(ml_newick.strip() + "\n")
        (OUT / "trees" / f"{stem}.cladogram.nwk").write_text(cladogram.strip() + "\n")

    # Most probable state per site for each named ancestor.
    named = con.execute(
        "SELECT aa.alias, a.id, a.name, a.almethod, a.phylomodel FROM AncestorsAlias aa "
        "JOIN Ancestors a ON a.id = aa.ancid ORDER BY aa.alias, a.almethod, a.phylomodel"
    ).fetchall()
    fasta: list[tuple[str, str]] = []
    tsv = ["alias\tnode\taligner\tmodel\tsite\tstate\tpp"]
    for alias, anc_id, node, al_id, model_id in named:
        best: dict[int, tuple[str, float]] = {}
        for site, state, pp in con.execute(
            "SELECT site, state, pp FROM AncestralStates WHERE ancid = ?", (anc_id,)
        ):
            if site not in best or pp > best[site][1]:
                best[site] = (state, pp)
        label = f"{alias}|{node}|{aligners[al_id]}|{models[model_id]}"
        fasta.append((label, "".join(best[s][0] for s in sorted(best))))
        for site in sorted(best):
            state, pp = best[site]
            tsv.append(
                f"{alias}\t{node}\t{aligners[al_id]}\t{models[model_id]}\t{site}\t{state}\t{pp:.6g}"
            )
    write_fasta(OUT / "ancestors.fasta", fasta)
    (OUT / "ancestors_ml_pp.tsv").write_text("\n".join(tsv) + "\n")

    groups: dict[str, list[str]] = defaultdict(list)
    for group, taxon in con.execute(
        "SELECT g.name, t.fullname FROM GroupsTaxa gt JOIN TaxaGroups g ON g.id = gt.groupid "
        "JOIN Taxa t ON t.id = gt.taxonid"
    ):
        groups[group].append(taxon)
    settings = con.execute("SELECT keyword, value FROM Settings").fetchall()
    lines = ["# Original run settings (keyword\tvalue)"] + [f"{k}\t{v}" for k, v in settings]
    lines += ["", "# Ancestor definitions (group\ttaxa)"]
    lines += [f"{g}\t{', '.join(sorted(t))}" for g, t in sorted(groups.items())]
    (OUT / "run_settings.tsv").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
