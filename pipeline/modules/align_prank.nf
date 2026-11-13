// Aligner interface:
//   in:  tuple val(meta), path(unaligned FASTA with safe IDs)
//   out: tuple val(meta + [aligner: KEY]), path(aligned FASTA)
// PRANK is phylogeny-aware: it keeps independent insertions in separate columns.
// It seeds from the clock unless told otherwise, so the job's seed is passed for
// reproducible output. Its -version flag checks online for updates, so the version is
// read from the help header instead. Slow: expect minutes, not seconds.
process ALIGN_PRANK {
    tag 'prank'
    container params.containers.prank

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta + [aligner: 'prank']), path('alignment.fasta'), emit: alignment
    path 'versions.tsv', emit: versions

    script:
    """
    prank -d=${fasta} -o=prank -f=fasta -seed=${params.seed}
    mv prank.best.fas alignment.fasta
    printf 'prank\\t%s\\n' "\$(prank -help 2>&1 | sed -n 's/^prank v\\.\\([0-9]*\\).*/\\1/p' | head -1)" > versions.tsv
    """
}
