// ESR, step 2 of 3: one held-out tip (an ID from ESR_FIT, never a sequence name),
// reconstructed at every level the settings ask for.
//   in:  tuple path(fit/), val(tip)
//   out: path(tip_<id>/) with esr_tips.tsv and esr_sites.tsv
process ESR_TIP {
    tag "${tip}"
    container params.containers.tools

    input:
    tuple path(fit, stageAs: 'fit'), val(tip)

    output:
    path "tip_${tip}", emit: result

    script:
    """
    phylobot-esr tip fit ${tip} tip_${tip}
    """
}
