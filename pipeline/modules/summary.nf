// Job-level summary.json (alignments, combinations, primary, fit statistics, the last
// stage run) and model_comparison.tsv (BIC ranking within each alignment).
process SUMMARY {
    container params.containers.tools
    publishDir params.outdir, mode: 'copy'

    input:
    path fits         // one row per combination: aligner, tree, model, fit.tsv values
    val primary       // preferred "aligner/tree/model" key
    path sources      // sources.tsv: labels for provided inputs
    val alignments    // the job's alignment keys, comma-separated
    val stop_after    // alignments, trees or ancestors

    output:
    path 'summary.json'
    path 'model_comparison.tsv'

    script:
    """
    phylobot-summary ${sources} ${fits} '${primary}' '${alignments}' ${stop_after}
    """
}
