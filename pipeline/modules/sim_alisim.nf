// Simulation, step 2 of 3: runs the planned AliSim commands in order, in the pinned
// IQ-TREE container. Each command<i>.args holds NUL-separated arguments written by
// SIM_PREPARE, so no argument is ever parsed by the shell.
//   in:  tuple val(seed), path(prep/)
//   out: tuple val(seed), path(run/) (prep/ plus AliSim's outputs), path(version.txt)
process SIM_ALISIM {
    tag "seed ${seed}"
    container params.containers.iqtree

    input:
    tuple val(seed), path(prep)

    output:
    tuple val(seed), path('run'), path('version.txt'), emit: simulated

    script:
    """
    mkdir run
    cp -R ${prep}/. run/
    cd run
    n=\$(ls command*.args | wc -l)
    i=0
    while [ \$i -lt \$n ]; do
        xargs -0 -x -a command\$i.args iqtree3 > alisim\$i.log 2>&1
        i=\$((i + 1))
    done
    cd ..
    iqtree3 --version | head -1 > version.txt
    """
}
