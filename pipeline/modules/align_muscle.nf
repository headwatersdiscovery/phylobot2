// Aligner interface:
//   in:  tuple val(meta), path(unaligned FASTA with safe IDs)
//   out: tuple val(meta + [aligner: KEY]), path(aligned FASTA)
// MUSCLE 5 reorders rows by its guide tree; downstream steps match rows by name.
process ALIGN_MUSCLE {
    tag 'muscle'
    container params.containers.muscle

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta + [aligner: 'muscle']), path('alignment.fasta'), emit: alignment
    path 'versions.tsv', emit: versions

    script:
    """
    muscle -align ${fasta} -output alignment.fasta -threads ${task.cpus}
    printf 'muscle\\t%s\\n' "\$(muscle -version 2>/dev/null | head -1 | cut -d' ' -f2 | sed 's/\\.linux64//')" > versions.tsv
    """
}
