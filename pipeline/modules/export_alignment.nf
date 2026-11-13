// Publishes one alignment, with the original sequence names, to <outdir>/<aligner>/
// alignment.fasta, and for coding DNA its codon form, alignment.cds.fasta ([] if none).
// Runs for every job; a job that stops after alignments ends here.
process EXPORT_ALIGNMENT {
    tag "${meta.aligner}"
    container params.containers.tools
    publishDir path: { "${params.outdir}/${meta.aligner}" }, mode: 'copy',
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(alignment), path(codon, stageAs: 'codon/*')
    path names

    output:
    tuple val(meta), path('out/alignment.fasta'), emit: alignment
    tuple val(meta), path('out/alignment.cds.fasta'), emit: codon, optional: true

    script:
    """
    phylobot-export-alignment ${names} ${alignment} out ${codon ?: ''}
    """
}
