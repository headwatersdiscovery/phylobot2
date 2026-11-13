// Simulation, step 1 of 3 (--simulate; src/phylobot/simulate): checks the spec and
// writes AliSim's inputs and planned commands for one seed. The tree and root sequence
// come as staged files ([] for a random tree or root), never as paths in the spec.
//   in:  val(seed); path(spec.json); path(tree) or []; path(root) or []
//   out: tuple val(seed), path(prep/)
process SIM_PREPARE {
    tag "seed ${seed}"
    container params.containers.tools

    input:
    val seed
    path spec, stageAs: 'spec.json'
    path tree, stageAs: 'input_tree.nwk'
    path root, stageAs: 'input_root.fasta'

    output:
    tuple val(seed), path('prep'), emit: prepared

    script:
    """
    phylobot-simulate prepare ${spec} ${seed} prep ${tree ? "--tree ${tree}" : ''} ${root ? "--root ${root}" : ''}
    """
}
