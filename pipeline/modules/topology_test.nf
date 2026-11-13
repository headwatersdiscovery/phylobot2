// Is each tree source's topology acceptable to the sequences alone? One IQ-TREE
// topology test per alignment and model over every tree source's tree (-z, -n 0,
// 10,000 RELL replicates, AU test), under the combination's model. GeneRax's reconciled
// tree trades sequence likelihood for reconciliation likelihood, so it always loses on
// BIC; the AU p says whether it is significantly worse. Run only
// when GeneRax is a tree source. Publishes <aligner>/topology_tests/<model>.tsv.
//   in:  tuple val(meta: aligner, model), path(alignment), path(model.txt),
//        val(tree keys), path(trees) -- keys and trees in the same order
//   out: tuple val(meta), path(topology_test.tsv): tree, logL, deltaL, bp_rell, p_kh,
//        p_sh, c_elw, p_au (one row per tree source)
include { mixtureFlag } from './mixtures'

process TOPOLOGY_TEST {
    tag "${meta.aligner}/${meta.model}"
    container params.containers.iqtree
    publishDir path: { "${params.outdir}/${meta.aligner}/topology_tests" }, mode: 'copy',
        saveAs: { filename -> filename == 'topology_test.tsv' ? "${meta.model}.tsv" : null }

    input:
    tuple val(meta), path(alignment), path(model), val(keys), path(trees, stageAs: 'trees/tree*.nwk')

    output:
    tuple val(meta), path('topology_test.tsv'), emit: test
    path 'versions.tsv', emit: versions

    script:
    // Tree keys are short stable keys (iqtree, raxml, generax, provided1...), never user text.
    def order = keys.join(' ')
    def staged = trees instanceof List ? trees.join(' ') : trees.toString()
    def weights = mixtureFlag(meta)  // appended: off, the command is unchanged
    """
    for t in ${staged}; do sed -E 's/\\)[^:,;()]+/)/g' "\$t"; done > candidates.nwk
    iqtree3 -s ${alignment} -m "\$(cat ${model})"${weights} -z candidates.nwk -n 0 -zb 10000 -au \\
        --seed ${params.seed} -T ${task.cpus} --prefix au --quiet
    printf '%s\\n' ${order} > keys.txt
    printf 'tree\\tlogL\\tdeltaL\\tbp_rell\\tp_kh\\tp_sh\\tc_elw\\tp_au\\n' > topology_test.tsv
    awk 'NR == FNR { key[NR] = \$1; next }
         /^Tree +logL/ { table = 1; next }
         table && /^ *[0-9]+ +-?[0-9]/ { print key[\$1] "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$6 "\\t" \$8 "\\t" \$10 "\\t" \$12 }
         table && /^\$/ { table = 0 }' keys.txt au.iqtree >> topology_test.tsv
    test "\$(wc -l < topology_test.tsv)" -eq \$((\$(wc -l < keys.txt) + 1))
    printf 'iqtree\\t%s\\n' "\$(iqtree3 --version | head -1 | cut -d' ' -f3)" > versions.tsv
    """
}
