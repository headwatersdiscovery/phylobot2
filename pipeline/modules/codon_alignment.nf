// Both forms of one alignment of coding sequences (coding DNA only): protein
// (alignment.fasta, for every amino-acid step) and codon (alignment.cds.fasta). `form` is
// what the aligner wrote: 'protein' rows are back-translated with each sequence's own
// codons, 'codon' rows (PRANK -codon) translated. Fails if an aligner changed any
// sequence (src/phylobot/pipeline/codon_alignment.py).
process CODON_ALIGNMENT {
    tag "${meta.aligner}"
    container params.containers.tools

    input:
    tuple val(meta), path(alignment, stageAs: 'aligned.fasta'), val(form)
    path cds

    output:
    tuple val(meta), path('out/alignment.fasta'),     emit: protein
    tuple val(meta), path('out/alignment.cds.fasta'), emit: codon

    script:
    """
    phylobot-codon-alignment ${meta.aligner} ${form} aligned.fasta ${cds} ${params.genetic_code} out
    """
}
