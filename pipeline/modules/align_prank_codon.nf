// Codon aligner (coding DNA only):
//   in:  tuple val(meta), path(coding sequences with safe IDs)
//   out: tuple val(meta + [aligner: 'prank_codon']), path(codon alignment)
// PRANK with its codon model, aligning the DNA directly rather than the translations.
// PRANK exits 0 even when it fails (a length that is not whole codons: a message and no
// alignment), so the script checks that an alignment was written; CODON_ALIGNMENT then
// checks every row against its input (PRANK silently replaces stop codons with NNN).
process ALIGN_PRANK_CODON {
    tag 'prank_codon'
    container params.containers.prank

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta + [aligner: 'prank_codon']), path('alignment.cds.fasta'), emit: alignment
    path 'versions.tsv', emit: versions

    script:
    """
    prank -d=${fasta} -o=prank -f=fasta -codon -seed=${params.seed}
    test -s prank.best.fas || { echo 'PRANK -codon wrote no alignment' >&2; exit 1; }
    mv prank.best.fas alignment.cds.fasta
    printf 'prank\\t%s\\n' "\$(prank -help 2>&1 | sed -n 's/^prank v\\.\\([0-9]*\\).*/\\1/p' | head -1)" > versions.tsv
    """
}
