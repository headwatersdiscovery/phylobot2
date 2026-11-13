// The job's species tree, made ready for reconciliation: once per job, after every tree
// search, its polytomies are resolved (guided by every combination's gene tree unless
// --polytomies arbitrary; src/phylobot/species/polytomies.py) and GeneRax's species map
// written. Publishes <outdir>/species/: species_tree.nwk (display names),
// species_nodes.tsv and polytomies.json.
//   in:  VALIDATE's species.safe.tsv, species_tree.safe.nwk, species_names.tsv,
//        outgroup.txt; every gene tree (safe IDs)
//   out: the resolved species tree (safe IDs, every node labelled), mapping.link,
//        species_nodes.tsv
process SPECIES_TREE {
    container params.containers.tools
    publishDir path: "${params.outdir}/species", mode: 'copy',
        saveAs: { filename -> filename in ['out/species_tree.nwk', 'out/species_nodes.tsv', 'out/polytomies.json'] ? filename.substring(4) : null }

    input:
    path mapping
    path species_tree
    path species_names
    path outgroup
    path gene_trees, stageAs: 'genes/tree*.nwk'

    output:
    path 'out/species_tree.resolved.nwk', emit: tree
    path 'out/mapping.link',              emit: mapping
    path 'out/species_nodes.tsv',         emit: nodes
    path 'out/*',                         emit: results

    script:
    // params.polytomies is checked in main.nf (guided or arbitrary).
    """
    phylobot-species-tree ${mapping} ${species_tree} ${species_names} ${outgroup} out \\
        --seed ${params.seed} --method ${params.polytomies} --gene-trees genes/*
    """
}
