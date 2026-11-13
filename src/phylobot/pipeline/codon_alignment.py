"""Both forms of one alignment of coding sequences: protein and codon.

    phylobot-codon-alignment KEY FORM ALIGNMENT CDS GENETIC_CODE OUTDIR

FORM is the alignment's form as the aligner wrote it. 'protein' (MAFFT, MUSCLE, Clustal
Omega, PRANK on the translations): each row is back-translated with its sequence's own
codons. 'codon' (PRANK -codon): each row is checked to hold its coding sequence and in
whole codons, then translated; codons with an ambiguity code that PRANK wrote as NNN
get the user's own codon back (`codons.restore_codons`). Either way OUTDIR gets
alignment.fasta (protein) and alignment.cds.fasta (codons), and any other change an
aligner made to a sequence fails the step: PRANK, for one, exits 0 after silently
replacing a stop codon with NNN.
KEY (the aligner's key) only labels messages.
"""

import sys
from pathlib import Path

from phylobot import codons
from phylobot.seqio import read_fasta, write_fasta


class AlignmentError(Exception):
    pass


def codon_alignment(
    key: str, form: str, alignment: Path, cds_path: Path, code: codons.GeneticCode, outdir: Path
) -> None:
    if not alignment.exists() or not alignment.stat().st_size:
        raise AlignmentError(f"{key} wrote no alignment")
    rows = read_fasta(alignment)
    cds = read_fasta(cds_path)
    if set(rows) != set(cds):
        raise AlignmentError(f"{key}: the alignment's sequences are not the input's")
    protein, codon = {}, {}
    for name in cds:  # input order
        row = rows[name].upper()
        try:
            if form == "protein":
                codon[name] = codons.backtranslate(row, cds[name], code)
                protein[name] = row
            else:
                row = codons.restore_codons(row, cds[name])
                protein[name] = codons.translate_row(row, code)
                codon[name] = row
        except codons.CodonError as exc:
            raise AlignmentError(f"{key}: sequence {name}: {exc}") from None
    if len({len(row) for row in codon.values()}) != 1:
        raise AlignmentError(f"{key}: rows have different lengths")
    outdir.mkdir(parents=True, exist_ok=True)
    write_fasta(outdir / "alignment.fasta", protein)
    write_fasta(outdir / "alignment.cds.fasta", codon)


def main() -> None:
    if len(sys.argv) != 7 or sys.argv[2] not in ("protein", "codon"):
        sys.exit(__doc__)
    key, form, alignment, cds, code, outdir = sys.argv[1:]
    try:
        codon_alignment(
            key, form, Path(alignment), Path(cds), codons.genetic_code(code), Path(outdir)
        )
    except (AlignmentError, codons.CodonError) as exc:
        print(f"Alignment error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
