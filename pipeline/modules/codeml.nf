// PAML codeml, one site model from one starting omega, on a codon combination's
// alignment and fixed tree (the IQ-TREE tree, internal node names removed; its branch
// lengths are codeml's starting values). The control file comes from VALIDATE
// (phylobot.codeml.write_controls): models are defined in one place.
//   in:  tuple val(meta), path(codon alignment), path(asr tree), val(model), val(start), path(ctl)
//   out: tuple val(meta), path(<model>.start<n>.mlc), path(<model>.start<n>.rst); the
//        control file as <model>.start<n>.ctl
process CODEML {
    tag "${meta.aligner}/${meta.tree}/${meta.model} ${model}.${start}"
    container params.containers.paml

    input:
    tuple val(meta), path(alignment, stageAs: 'input.fasta'), path(tree, stageAs: 'input.nwk'), val(model), val(start), path(ctl, stageAs: 'control.ctl')

    output:
    tuple val(meta), path("${model}.${start}.mlc"), path("${model}.${start}.rst"), emit: run
    tuple val(meta), path("${model}.${start}.ctl"), emit: control  // as run, for publishing
    path 'versions.tsv', emit: versions

    script:
    """
    cp ${alignment} aln.fasta
    # codeml's tree file: "<sequences> <trees>", then the tree (without the header it
    # analyses the tree, then exits 1 looking for another one).
    { echo "\$(grep -c '>' aln.fasta) 1"; sed -E 's/\\)Node[0-9]+/)/g' ${tree}; } > tree.nwk
    cp ${ctl} codeml.ctl
    cp ${ctl} ${model}.${start}.ctl
    codeml codeml.ctl > codeml.log
    # (the branch-site models run without ancestors: an empty rst)
    test -s mlc && test -e rst
    mv mlc ${model}.${start}.mlc
    mv rst ${model}.${start}.rst
    printf 'paml\\t%s\\n' "\$(sed -n 's/^CODONML (in paml version \\([^,]*\\).*/\\1/p' ${model}.${start}.mlc | head -1)" > versions.tsv
    """
}
