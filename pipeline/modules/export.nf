// Converts one combination's output to PhyloBot's results files under
// <outdir>/<aligner>/<tree>/<model>/, restoring the original sequence names, and writes
// ancestors.gaps-<method>.fasta for each chosen gap placement method (masking only;
// ancestors.fasta stays unmasked). Inputs for ml, prank and grasp come from GAPS_ML,
// GAPS_PRANK and GAPS_GRASP; for methods not chosen they are empty.
process EXPORT {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.tools
    publishDir path: { "${params.outdir}/${meta.aligner}/${meta.tree}/${meta.model}" },
        mode: 'copy',
        saveAs: { filename -> filename.tokenize('/').last() }

    input:
    tuple val(meta), path(asr_tree), path(state), path(alignment), path(support_tree), path(model), path(report), path(ranking, stageAs: 'modelfinder/*'), path(gap_state, stageAs: 'gaps/*'), path(gap_tree, stageAs: 'gaps/*'), path(prank_anc, stageAs: 'prank/*'), path(prank_tree, stageAs: 'prank/*'), path(grasp_anc, stageAs: 'grasp/*'), path(grasp_tree, stageAs: 'grasp/*'), path(rooted_tree, stageAs: 'rooted/*'), path(rooting, stageAs: 'rooted/*'), path(codon_alignment, stageAs: 'codon/*'), path(site_rates, stageAs: 'sites/*'), path(site_lh, stageAs: 'sites/*')
    path names
    val gap_methods  // comma-separated, validated by the workflow ('' for none)

    output:
    tuple val(meta), path('out/*'), emit: results
    tuple val(meta), path('out/fit.tsv'), emit: fit

    script:
    // Support measures in the tree's node labels, in IQ-TREE's order: SH-aLRT on every
    // tree, UFBoot only from IQ-TREE's own search (FIT_IQTREE cannot compute it).
    def supports = ((params.alrt as Integer) > 0 ? ['sh_alrt'] : []) +
        (meta.tree == 'iqtree' && (params.bootstrap as Integer) > 0 ? ['ufboot'] : [])
    def rankingFile = ranking ? "--modelfinder-report ${ranking}" : ''
    def mlFiles = gap_state ? "--gap-state ${gap_state} --gap-tree ${gap_tree}" : ''
    def prankFiles = prank_anc ? "--prank-ancestors ${prank_anc} --prank-tree ${prank_tree}" : ''
    def graspFiles = grasp_anc ? "--grasp-ancestors ${grasp_anc} --grasp-tree ${grasp_tree}" : ''
    def rootedFiles = rooted_tree ? "--rooted-tree ${rooted_tree} --rooting ${rooting}" : ''
    // A codon model: STATEFILE holds codon posteriors (summed into amino acids).
    def codonFiles = codon_alignment ? "--codon-alignment ${codon_alignment} --genetic-code ${params.genetic_code}" : ''
    // Written to out/ so no output can overwrite a staged input of the same name.
    """
    phylobot-export ${names} ${alignment} ${asr_tree} ${support_tree} ${state} ${model} ${report} out \\
        --gap-methods '${gap_methods}' --supports '${supports.join(',')}' ${rankingFile} ${mlFiles} ${prankFiles} ${graspFiles} ${rootedFiles} ${codonFiles} \\
        --site-rates ${site_rates} --site-lh ${site_lh}
    """
}
