// Aligner interface:
//   in:  tuple val(meta), path(unaligned FASTA with safe IDs)
//   out: tuple val(meta + [aligner: KEY]), path(aligned FASTA)
process ALIGN_CLUSTALO {
    tag 'clustalo'
    container params.containers.clustalo

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta + [aligner: 'clustalo']), path('alignment.fasta'), emit: alignment
    path 'versions.tsv', emit: versions

    script:
    """
    clustalo -i ${fasta} -o alignment.fasta --outfmt=fasta --threads=${task.cpus} --force
    printf 'clustalo\\t%s\\n' "\$(clustalo --version)" > versions.tsv
    """
}
