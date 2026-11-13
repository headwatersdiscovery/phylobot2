// Tier-0 property posteriors for one combination's ancestors.
// Draws from each node's per-column posteriors and writes every property's values per
// draw plus a summary, under <outdir>/<aligner>/<tree>/<model>/properties/: draws.parquet,
// summaries.json and index.json (phylobot.property_store).
// No structure, no network, no GPU. The first properties are pure Python (phylobot-tools);
// the neural predictors need phylobot-ml, so the container follows
// the keys asked for (`ml_keys`, from main.nf's mlPropertyKeys).
// Only runs when params.properties names at least one property.
// `nodes` restricts it to named ancestral nodes (empty = all): the inline run measures
// every node, while a run asked for afterwards names its scope's
// nodes, so measuring one ancestor costs one ancestor. The names are checked against a
// safe pattern by the workflow before they reach this command line.
//   in:  tuple val(meta), path(the combination's published files), val(nodes)
//   out: tuple val(meta), path(properties/)
process PROPERTIES {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container { keys.tokenize(',').any { it in ml_keys } ? params.containers.ml : params.containers.tools }
    // The predictors need more than the 2 GB default: torch's own footprint, plus a
    // model, plus activations quadratic in sequence length. Measured 2026-10-21 on
    // ~1050-column ancestors: 1.05 GB for metapredict over 1000 draws, 2.15 GB for ESM-2
    // before its batch was cut to one window per pass.
    memory { keys.tokenize(',').any { it in ml_keys } ? '4 GB' : '2 GB' }
    publishDir path: { "${params.outdir}/${meta.aligner}/${meta.tree}/${meta.model}" },
        mode: 'copy'

    input:
    tuple val(meta), path(results, stageAs: 'combination/*'), val(nodes)
    val keys          // comma-separated property keys
    val gap_versions  // comma-separated gap methods the job ran ('' for none)
    val ml_keys       // the keys that need phylobot-ml

    output:
    tuple val(meta), path('properties'), emit: properties

    script:
    """
    phylobot-properties combination properties \\
        --properties '${keys}' \\
        --draws ${params.property_draws} \\
        --seed ${params.seed} \\
        --gap-versions '${gap_versions}' \\
        --nodes '${nodes}'
    """
}
