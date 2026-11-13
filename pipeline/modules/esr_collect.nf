// ESR, step 3 of 3: every tip's rows with the real names put back, and the job page's
// figures, published to <outdir>/esr/.
//   in:  path(fit/); path(tip_*/) collected
//   out: path(esr/)
process ESR_COLLECT {
    tag 'collect'
    container params.containers.tools
    publishDir path: params.outdir, mode: 'copy'

    input:
    path fit, stageAs: 'fit'
    path tips, stageAs: 'tips/*'

    output:
    path 'esr/**', emit: results

    script:
    """
    phylobot-esr collect fit esr tips/*
    """
}
