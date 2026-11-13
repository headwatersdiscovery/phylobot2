// A combination's tree without ancestors (jobs that stop after trees), published to
// <outdir>/<aligner>/<tree>/<model>/: alignment.fasta, tree.nwk (internal nodes named
// Node1...), supports.tsv, model.txt, fit.tsv, modelfinder.tsv, and with an outgroup
// tree.rooted.nwk and rooting.json.
process EXPORT_TREE {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.tools
    publishDir path: { "${params.outdir}/${meta.aligner}/${meta.tree}/${meta.model}" },
        mode: 'copy',
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(alignment), path(tree), path(model), path(report), path(ranking, stageAs: 'modelfinder/*')
    path names
    path outgroup  // outgroup.txt from VALIDATE (safe IDs; empty = none)

    output:
    tuple val(meta), path('out/*'), emit: results
    tuple val(meta), path('out/fit.tsv'), emit: fit

    script:
    // Support measures in the tree's node labels, in IQ-TREE's order: SH-aLRT on every
    // tree, UFBoot only from IQ-TREE's own search (FIT_IQTREE cannot compute it).
    def supports = ((params.alrt as Integer) > 0 ? ['sh_alrt'] : []) +
        (meta.tree == 'iqtree' && (params.bootstrap as Integer) > 0 ? ['ufboot'] : [])
    def rankingFile = ranking ? "--modelfinder-report ${ranking}" : ''
    """
    phylobot-export-tree ${names} ${alignment} ${tree} ${model} ${report} ${outgroup} out \
        --supports '${supports.join(',')}' ${rankingFile}
    """
}
