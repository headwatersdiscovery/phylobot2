// GRASP gap placement: GRASP's joint reconstruction with indel inference (its default,
// BEP) on the combination's alignment and tree. An alternative reconstruction; EXPORT
// renames its nodes (N0...) to ours. GRASP offers JTT, Dayhoff, LG and WAG among our
// matrices; other models fall back to LG (recorded in grasp_model.txt).
//   in:  tuple val(meta), path(alignment), path(asr tree), path(model.txt)
//   out: tuple val(meta), path(ancestors), path(node-labelled tree)
process GAPS_GRASP {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.grasp

    input:
    tuple val(meta), path(alignment), path(tree), path(model)

    output:
    tuple val(meta), path('grasp/aln_ancestors.fa'), path('grasp/aln_ancestors.nwk'), emit: ancestors
    path 'grasp_model.txt', emit: model

    script:
    """
    sed -E 's/\\)Node[0-9]+/)/g' ${tree} > guide.nwk
    cp ${alignment} aln.fasta
    matrix=\$(cut -d+ -f1 ${model} | tr '[:lower:]' '[:upper:]')
    case "\$matrix" in
        LG) grasp_model=LG ;; WAG) grasp_model=WAG ;; JTT) grasp_model=JTT ;;
        DAYHOFF) grasp_model=Dayhoff ;; *) grasp_model=LG ;;
    esac
    echo "\$grasp_model (combination model: \$(cat ${model}))" > grasp_model.txt
    grasp -a aln.fasta -n guide.nwk -s \$grasp_model -j --seed ${params.seed} \\
        -o grasp --save-as FASTA TREE
    """
}
