// A combination's reconciliation on its published tree, with the sequences' and
// species' own names (src/phylobot/pipeline/reconcile_export.py): reconciliation.tsv,
// .json and .nwk, published beside tree.nwk.
//   in:  tuple val(meta), path(RecPhyloXML), path(GeneRax's stats), path(generax_model.txt),
//        path(the tree as published: the reconstruction's, or a trees-only job's search
//        tree, whose nodes are then named as EXPORT_TREE names them); names.tsv;
//        species_nodes.tsv
process RECONCILE_EXPORT {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.tools
    publishDir path: { "${params.outdir}/${meta.aligner}/${meta.tree}/${meta.model}" },
        mode: 'copy',
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(xml), path(stats), path(model), path(tree)
    path names
    path species_nodes
    val name_nodes  // true when `tree` is a search tree (trees-only jobs)

    output:
    tuple val(meta), path('out/*'), emit: results

    script:
    def naming = name_nodes ? '--name-nodes' : ''
    """
    phylobot-reconcile-export ${names} ${species_nodes} ${tree} ${xml} ${stats} ${model} out ${naming}
    """
}
