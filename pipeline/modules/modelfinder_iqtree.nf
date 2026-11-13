// ModelFinder, once per alignment: picks the model (BIC) that every tree source then uses
// for that alignment's 'modelfinder' combinations (TREE_IQTREE, FIT_IQTREE). Scored on
// ModelFinder's own initial tree, as inside a full IQ-TREE run.
//   in:  tuple val(meta), path(alignment)   -- meta: [aligner]
//        path(modelfinder_args)             -- candidate flags (empty = IQ-TREE's default)
//   out: tuple val(meta), path(model.txt: the pick), path(report: the candidate ranking)
process MODELFINDER {
    tag "${meta.aligner}"
    container params.containers.iqtree

    input:
    tuple val(meta), path(alignment)
    path modelfinder_args

    output:
    tuple val(meta), path('model.txt'), path('modelfinder.iqtree'), emit: model
    path 'versions.tsv', emit: versions

    script:
    // The pick goes into other steps' command lines (quoted); it must look like a model.
    // Codon ModelFinder (MODELFINDER_CODON; meta.codon) runs on the codon alignment.
    def seqtype = meta.codon ? "-st CODON${params.genetic_code}" : ''
    """
    iqtree3 -s ${alignment} ${seqtype} -m MF \$(cat ${modelfinder_args}) --seed ${params.seed} -T ${task.cpus} --prefix modelfinder --quiet
    sed -n 's/^Best-fit model according to BIC: //p' modelfinder.iqtree > model.txt
    grep -Eqx '[A-Za-z0-9.+]+' model.txt
    printf 'iqtree\\t%s\\n' "\$(iqtree3 --version | head -1 | cut -d' ' -f3)" > versions.tsv
    """
}
