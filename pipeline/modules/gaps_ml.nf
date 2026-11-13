// ML gap placement, step 1 of 2 (EXPORT applies the result): reconstruct the gap (0) /
// residue (1) pattern of the alignment on the combination's own tree. Binary GTR2 with
// branch lengths rescaled by one fitted factor (-blscale): gaps change far more slowly
// than residues, and one rate is better determined than every branch length.
//   in:  tuple val(meta), path(alignment), path(asr tree)
//   out: tuple val(meta), path(gap .state), path(gap tree)
process GAPS_ML {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.iqtree

    input:
    tuple val(meta), path(alignment), path(tree)

    output:
    tuple val(meta), path('gaps.state'), path('gaps.treefile'), emit: gaps

    script:
    """
    sed '/^>/!{s/[^-]/1/g; s/-/0/g;}' ${alignment} > gaps.fasta
    iqtree3 -s gaps.fasta -st BIN -te ${tree} -m GTR2 -blscale --ancestral \\
        --seed ${params.seed} -T ${task.cpus} --prefix gaps --quiet
    """
}
