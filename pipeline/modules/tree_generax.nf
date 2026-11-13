// GeneRax tree search (--strategy SPR; UndatedDL): the gene tree that best fits the
// sequences and the species tree jointly, starting from the alignment's own ML tree
// (IQ-TREE's, else RAxML-NG's, else the first provided tree; main.nf chooses). Only the
// topology is kept: FIT_GENERAX (FIT_IQTREE) then fits the combination's own model on it,
// as for RAxML-NG. GeneRax's substitution model only steers the search: the combination's
// matrix +G if GeneRax has it, else LG+G, recorded in generax_model.txt (as in RECONCILE).
// Amino-acid combinations only.
//   in:  tuple val(meta), path(alignment), path(starting tree), path(model.txt);
//        the resolved species tree; mapping.link (both from SPECIES_TREE)
//   out: tuple val(meta), path(GeneRax's tree), path(generax_model.txt)
process TREE_GENERAX {
    tag "${meta.aligner}/${meta.model}"
    container params.containers.generax

    input:
    tuple val(meta), path(alignment), path(tree), path(model)
    path species_tree
    path mapping

    output:
    tuple val(meta), path('generax.nwk'), path('generax_model.txt'), emit: tree
    path 'versions.tsv', emit: versions

    script:
    """
    sed -E 's/\\)[^:,;()]+/)/g' ${tree} > start.nwk
    matrix=\$(cut -d+ -f1 ${model})
    case "\$matrix" in
        LG|WAG|JTT|Dayhoff|DCMut|VT|PMB|Blosum62|mtREV|mtART|mtZOA|mtMAM|cpREV|rtREV|HIVb|HIVw|FLU) ;;
        JTTDCMut) matrix=JTT ;;
        *) matrix=LG ;;
    esac
    echo "\$matrix+G" > generax_model.txt
    printf '[FAMILIES]\\n- family\\nstarting_gene_tree = start.nwk\\nalignment = %s\\nmapping = %s\\nsubst_model = %s\\n' \\
        ${alignment} ${mapping} "\$(cat generax_model.txt)" > families.txt
    generax -f families.txt -s ${species_tree} -p gx --rec-model UndatedDL --strategy SPR \\
        --seed ${params.seed} > generax.out 2>&1
    cp gx/results/family/geneTree.newick generax.nwk
    printf 'generax\\t%s\\n' "\$(grep -m1 -o 'GeneRax [0-9.]*' generax.out | cut -d' ' -f2)" > versions.tsv
    """
}
