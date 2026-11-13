// Aligner interface:
//   in:  tuple val(meta), path(unaligned FASTA with safe IDs)
//   out: tuple val(meta + [aligner: KEY]), path(aligned FASTA)
process ALIGN_MAFFT {
    tag 'mafft'
    container params.containers.mafft

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta + [aligner: 'mafft']), path('alignment.fasta'), emit: alignment
    path 'versions.tsv', emit: versions

    script:
    """
    mafft --auto --quiet --thread ${task.cpus} ${fasta} > alignment.fasta
    printf 'mafft\\t%s\\n' "\$(mafft --version 2>&1 | cut -d' ' -f1)" > versions.tsv
    """
}
