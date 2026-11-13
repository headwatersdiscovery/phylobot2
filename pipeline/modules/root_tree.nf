// Roots a combination's tree on the user's outgroup (only run when one was given).
// Records whether the outgroup is monophyletic; a scattered outgroup is a warning on
// the job page, not an error (users often pick a few representatives).
//   in:  tuple val(meta), path(asr tree); path(outgroup.txt)
//   out: tuple val(meta), path(rooted tree), path(rooting.json)
process ROOT_TREE {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.tools

    input:
    tuple val(meta), path(tree)
    path outgroup

    output:
    tuple val(meta), path('rooted.nwk'), path('rooting.json'), emit: rooted

    script:
    """
    phylobot-root ${tree} ${outgroup} rooted.nwk rooting.json
    """
}
