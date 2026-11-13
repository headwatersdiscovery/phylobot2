// Tree interface:
//   in:  tuple val(meta), path(alignment), path(picked model)
//        -- meta.model: a catalog model or 'modelfinder'; picked model: MODELFINDER's
//           model.txt for 'modelfinder' combinations, empty ([]) for fixed models
//   out: tuple val(meta + [tree: KEY]), path(tree, nodes labelled with branch support),
//        path(model.txt), path(IQ-TREE report: final fit statistics)
include { mixtureFlag } from './mixtures'

process TREE_IQTREE {
    tag "${meta.aligner}/${meta.model}"
    container params.containers.iqtree

    input:
    tuple val(meta), path(alignment), path(picked, stageAs: 'picked/*')

    output:
    tuple val(meta + [tree: 'iqtree']), path('tree.treefile'), path('model.txt'), path('tree.iqtree'), emit: tree
    path 'versions.tsv', emit: versions

    script:
    // as Integer: a command-line value is text. main.nf checks both ranges.
    def bootstrap = (params.bootstrap as Integer) > 0 ? "-B ${params.bootstrap as Integer}" : ''
    def alrt = (params.alrt as Integer) > 0 ? "-alrt ${params.alrt as Integer}" : ''
    // Fixed: a catalog model (validated). ModelFinder: MODELFINDER's pick (checked there).
    def model = meta.model.startsWith('modelfinder') ? "\"\$(cat ${picked})\"" : "'${meta.model}'"
    // Codon models (meta.codon) run on the codon alignment, in the job's genetic code.
    def seqtype = meta.codon ? "-st CODON${params.genetic_code}" : ''
    def weights = mixtureFlag(meta)  // appended: off, the command is unchanged
    """
    iqtree3 -s ${alignment} ${seqtype} -m ${model}${weights} ${bootstrap} ${alrt} --seed ${params.seed} -T ${task.cpus} --prefix tree --quiet
    echo ${model} > model.txt
    test -s model.txt
    printf 'iqtree\\t%s\\n' "\$(iqtree3 --version | head -1 | cut -d' ' -f3)" > versions.tsv
    """
}
