// Ancestral interface:
//   in:  tuple val(meta), path(alignment), path(tree), path(model.txt)
//   out: tuple val(meta), path(node-labelled tree), path(IQ-TREE .state); and (sites)
//        per-column rates (-wsr; header only without rate heterogeneity) and
//        log-likelihoods (-wsl) on the final tree, at no measurable cost
// (converted to PhyloBot's format by EXPORT). A provided tree marked 'keep'
// (--tree_branch_lengths) keeps its branch lengths (-blfix): without it, --ancestral
// re-estimates them.
include { keepsLengths } from './branch_lengths'
include { mixtureFlag } from './mixtures'

process ANCESTRAL_IQTREE {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.iqtree

    input:
    tuple val(meta), path(alignment), path(tree), path(model)

    output:
    tuple val(meta), path('asr.treefile'), path('asr.state'), emit: ancestors
    tuple val(meta), path('asr.rate'), path('asr.sitelh'), emit: sites
    path 'versions.tsv', emit: versions

    script:
    def seqtype = meta.codon ? "-st CODON${params.genetic_code}" : ''  // codon posteriors
    def blfix = keepsLengths(meta) ? ' -blfix' : ''  // as in FIT_IQTREE
    def weights = mixtureFlag(meta)  // as blfix; --ancestral refits the model it is given
    """
    iqtree3 -s ${alignment} ${seqtype} -te ${tree}${blfix} -m "\$(cat ${model})"${weights} --ancestral -wsr -wsl \\
        --seed ${params.seed} -T ${task.cpus} --prefix asr --quiet
    printf 'iqtree\\t%s\\n' "\$(iqtree3 --version | head -1 | cut -d' ' -f3)" > versions.tsv
    """
}
