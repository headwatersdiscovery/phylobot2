// Tree interface, for a given topology (provided, RAxML-NG's, or for a codon model an
// amino-acid tree): fits the model and branch lengths on the alignment without searching
// for a new tree. Internal node labels (another analysis's branch supports) are removed
// first: IQ-TREE would otherwise keep them in its output as if they were this fit's.
// Branch support: SH-aLRT only (IQ-TREE refuses UFBoot on a fixed topology).
// A provided tree marked 'keep' (--tree_branch_lengths) keeps its branch lengths (-blfix).
//   in:  tuple val(meta), path(alignment), path(provided tree), path(picked model)
//        -- meta.model as in TREE_IQTREE; picked model as in TREE_IQTREE
//   out: tuple val(meta), path(fitted tree), path(model.txt), path(IQ-TREE report)
include { keepsLengths } from './branch_lengths'
include { mixtureFlag } from './mixtures'

process FIT_IQTREE {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.iqtree

    input:
    tuple val(meta), path(alignment), path(tree), path(picked, stageAs: 'picked/*')

    output:
    tuple val(meta), path('fit.treefile'), path('model.txt'), path('fit.iqtree'), emit: tree
    path 'versions.tsv', emit: versions

    script:
    // Fixed: a catalog model (validated). ModelFinder: MODELFINDER's pick (checked there).
    def model = meta.model.startsWith('modelfinder') ? "\"\$(cat ${picked})\"" : "'${meta.model}'"
    def alrt = (params.alrt as Integer) > 0 ? "-alrt ${params.alrt as Integer}" : ''  // as Integer: see TREE_IQTREE
    def seqtype = meta.codon ? "-st CODON${params.genetic_code}" : ''  // as in TREE_IQTREE
    def blfix = keepsLengths(meta) ? ' -blfix' : ''  // appended: off, the command is unchanged
    def weights = mixtureFlag(meta)  // as blfix
    """
    sed -E 's/\\)[^:,;()]+/)/g' ${tree} > topology.nwk
    iqtree3 -s ${alignment} ${seqtype} -te topology.nwk${blfix} -m ${model}${weights} ${alrt} --seed ${params.seed} -T ${task.cpus} --prefix fit --quiet
    echo ${model} > model.txt
    test -s model.txt
    printf 'iqtree\\t%s\\n' "\$(iqtree3 --version | head -1 | cut -d' ' -f3)" > versions.tsv
    """
}
