// PRANK gap placement: PRANK's own ancestors on the combination's alignment, kept as is
// (-keep), and tree. An alternative reconstruction (PRANK's residues, not IQ-TREE's);
// EXPORT renames its nodes to ours. PRANK roots the tree by splitting a branch; that
// extra node matches none of ours and is dropped.
//   in:  tuple val(meta), path(alignment), path(asr tree)
//   out: tuple val(meta), path(ancestors), path(node-labelled tree)
process GAPS_PRANK {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.prank

    input:
    tuple val(meta), path(alignment), path(tree)

    output:
    tuple val(meta), path('prank.anc.fas'), path('prank.anc.dnd'), emit: ancestors

    script:
    """
    sed -E 's/\\)Node[0-9]+/)/g' ${tree} > guide.nwk
    prank -d=${alignment} -t=guide.nwk -keep -showanc -f=fasta -o=prank -seed=${params.seed}
    """
}
