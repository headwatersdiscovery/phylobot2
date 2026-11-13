// RAxML-NG tree search: a second, independent ML search heuristic. Only the topology is
// kept: FIT_IQTREE then fits the same model and branch lengths on it, as for a provided
// tree, so every tree source is compared (and its ancestors reconstructed) under one
// likelihood implementation. IQ-TREE's model strings mean the same in RAxML-NG 2.0.3
// (checked: same log-likelihoods on a fixed tree), except JTTDCMut, which it lacks;
// VALIDATE refuses that combination and leaves JTTDCMut out of ModelFinder's candidates.
//   in:  tuple val(meta), path(alignment), path(picked model)   -- as TREE_IQTREE
//   out: tuple val(meta), path(best tree)   -- same meta as the input (main.nf rejoins it
//        with its alignment and model for FIT_IQTREE)
process TREE_RAXML {
    tag "${meta.aligner}/${meta.model}"
    container params.containers.raxml

    input:
    tuple val(meta), path(alignment), path(picked, stageAs: 'picked/*')

    output:
    tuple val(meta), path('raxml.nwk'), emit: tree
    path 'versions.tsv', emit: versions

    script:
    // Fixed: a catalog model (validated). ModelFinder: MODELFINDER's pick (checked there).
    def model = meta.model == 'modelfinder' ? "\"\$(cat ${picked})\"" : "'${meta.model}'"
    """
    raxml-ng --search --msa ${alignment} --model ${model} --seed ${params.seed} \\
        --threads ${task.cpus} --prefix search > raxml.out
    cp search.raxml.bestTree raxml.nwk
    printf 'raxml-ng\\t%s\\n' "\$(raxml-ng --version | sed -n 's/^RAxML-NG v\\. \\([^ ]*\\).*/\\1/p' | head -1)" > versions.tsv
    """
}
