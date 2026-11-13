// Validates every user input and renames to safe IDs. User-controlled file names and
// labels arrive only inside manifest.json; files are staged under fixed names.
process VALIDATE {
    container params.containers.tools

    input:
    path manifest
    path sequences,  stageAs: 'input.fasta'
    path alignments, stageAs: 'alignment?.fasta'
    path trees,      stageAs: 'tree?.nwk'
    // Reconciliation only (else empty): the species table and species tree.
    path species_map,  stageAs: 'species.tsv'
    path species_tree, stageAs: 'species_tree.nwk'

    output:
    path 'sanitized.fasta',               emit: fasta
    path 'names.tsv',                     emit: names
    path 'sources.tsv',                   emit: sources
    path 'modelfinder_args.txt',          emit: modelfinder_args
    path 'codon_modelfinder_args.txt',    emit: codon_modelfinder_args
    path 'outgroup.txt',                  emit: outgroup
    path 'provided/alignments/*.fasta',   emit: alignments, optional: true
    path 'provided/trees/*.nwk',          emit: trees,      optional: true
    // Coding DNA only (--seqtype codon): the coding sequences, and provided codon alignments.
    path 'sanitized.cds.fasta',           emit: cds,        optional: true
    path 'provided/codon_alignments/*.fasta', emit: codon_alignments, optional: true
    path 'codeml/*/*.ctl',                emit: codeml_controls, optional: true
    // Reconciliation only: species renamed sp001..., tree nodes spn001..., and their names.
    path 'species.safe.tsv',              emit: species_map,   optional: true
    path 'species_tree.safe.nwk',         emit: species_tree,  optional: true
    path 'species_names.tsv',             emit: species_names, optional: true
    path 'validation.json',               emit: report
    path 'versions.tsv',                  emit: versions

    script:
    """
    phylobot-validate ${manifest}
    printf 'phylobot\\t%s\\n' "\$(python3 -c 'import importlib.metadata as m; print(m.version("phylobot"))')" > versions.tsv
    printf 'python\\t%s\\n' "\$(python3 --version | cut -d' ' -f2)" >> versions.tsv
    """
}
