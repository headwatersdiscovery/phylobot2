// GeneRax reconciles one combination's gene tree with the species tree, keeping its
// topology (--strategy EVAL; UndatedDL): it roots the tree, infers duplication and
// loss rates, and labels every node with its event. Events on a fixed topology do not
// depend on the substitution model, so GeneRax gets the combination's matrix +G if it
// has it, else LG+G (LG+G for codon models: the protein alignment is used), recorded in
// generax_model.txt. Internal node labels (supports) are stripped from the tree.
//   in:  tuple val(meta), path(protein alignment), path(tree), path(model.txt);
//        the resolved species tree; mapping.link (both from SPECIES_TREE)
//   out: tuple val(meta), path(RecPhyloXML), path(GeneRax's stats), path(generax_model.txt)
process RECONCILE {
    tag "${meta.aligner}/${meta.tree}/${meta.model}"
    container params.containers.generax

    input:
    tuple val(meta), path(alignment), path(tree), path(model)
    path species_tree
    path mapping

    output:
    tuple val(meta), path('reconciliation.xml'), path('generax_stats.txt'), path('generax_model.txt'), emit: reconciliation
    path 'versions.tsv', emit: versions

    script:
    def codon = meta.codon ? 'true' : 'false'
    """
    sed -E 's/\\)[^:,;()]+/)/g' ${tree} > gene.nwk
    matrix=\$(cut -d+ -f1 ${model})
    case "\$matrix" in
        LG|WAG|JTT|Dayhoff|DCMut|VT|PMB|Blosum62|mtREV|mtART|mtZOA|mtMAM|cpREV|rtREV|HIVb|HIVw|FLU) ;;
        JTTDCMut) matrix=JTT ;;
        *) matrix=LG ;;
    esac
    if [ ${codon} = true ]; then matrix=LG; fi
    echo "\$matrix+G" > generax_model.txt
    printf '[FAMILIES]\\n- family\\nstarting_gene_tree = gene.nwk\\nalignment = %s\\nmapping = %s\\nsubst_model = %s\\n' \\
        ${alignment} ${mapping} "\$(cat generax_model.txt)" > families.txt
    generax -f families.txt -s ${species_tree} -p gx --rec-model UndatedDL --strategy EVAL \\
        --seed ${params.seed} > generax.out 2>&1
    cp gx/reconciliations/family_reconciliated.xml reconciliation.xml
    cp gx/results/family/stats.txt generax_stats.txt
    printf 'generax\\t%s\\n' "\$(grep -m1 -o 'GeneRax [0-9.]*' generax.out | cut -d' ' -f2)" > versions.tsv
    """
}
