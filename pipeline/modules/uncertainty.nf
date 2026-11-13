// Where a finished job's ancestral uncertainty comes from (src/phylobot/uncertainty_run.py):
// per ancestor and site, the posterior's own variance and each choice's Shapley share of
// the rest, published to <outdir>/uncertainty/.
process UNCERTAINTY {
    tag 'uncertainty'
    container params.containers.tools
    memory '2 GB'
    publishDir path: params.outdir, mode: 'copy'

    input:
    path published, stageAs: 'published'
    path settings, stageAs: 'settings.json'

    output:
    path 'uncertainty/**', emit: results

    script:
    """
    phylobot-uncertainty published settings.json .
    """
}
