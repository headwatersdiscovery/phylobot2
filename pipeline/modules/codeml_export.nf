// codeml's results for one codon combination: the best start per model, ancestors under
// the combination's node names, likelihood-ratio tests and BEB sites
// (src/phylobot/pipeline/codeml_export.py), published to <combination>/codeml/; with the
// free-ratio model also its per-branch dN/dS (FR/branches.tsv, original sequence names).
process CODEML_EXPORT {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.tools
    publishDir path: { "${params.outdir}/${meta.aligner}/${meta.tree}/${meta.model}/codeml" },
        mode: 'copy', saveAs: { filename -> filename.startsWith('out/') ? filename.substring(4) : null }

    input:
    tuple val(meta), path(mlcs, stageAs: 'runs/*'), path(rsts, stageAs: 'runs/*'), path(asr_tree)
    path names

    output:
    tuple val(meta), path('out/**'), emit: results

    script:
    """
    phylobot-codeml-export ${asr_tree} ${params.genetic_code} runs out ${names}
    """
}
