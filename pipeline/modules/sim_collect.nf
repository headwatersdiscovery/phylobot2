// Simulation, step 3 of 3: reads AliSim's outputs into one replicate's truth bundle,
// published to <outdir>/replicates/<seed>/ (the golden layout; src/phylobot/simulate).
//   in:  tuple val(seed), path(run/), path(version.txt); path(spec.json)
//   out: path(<seed>/)
process SIM_COLLECT {
    tag "seed ${seed}"
    container params.containers.tools
    publishDir path: "${params.outdir}/replicates", mode: 'copy'

    input:
    tuple val(seed), path(run), path(version)
    path spec, stageAs: 'spec.json'

    output:
    path "${seed}", emit: replicate

    script:
    """
    phylobot-simulate collect ${spec} ${run} ${seed} --version ${version} --image ${params.containers.iqtree}
    """
}
