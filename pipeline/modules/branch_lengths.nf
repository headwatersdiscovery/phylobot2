// Whether a combination's tree keeps the branch lengths it was given: a provided tree
// (provided1, provided2, ...) marked 'keep' in --tree_branch_lengths (checked in main.nf).
// FIT_IQTREE and ANCESTRAL_IQTREE then pass -blfix; both are needed, as IQ-TREE's
// --ancestral re-estimates the lengths of the tree it is given.
def keepsLengths(meta) {
    if (!(meta.tree ==~ /provided[0-9]+/)) {
        return false
    }
    def values = params.tree_branch_lengths instanceof List
        ? params.tree_branch_lengths
        : (params.tree_branch_lengths ?: '').toString().tokenize(',')
    def index = (meta.tree.substring('provided'.length()) as Integer) - 1
    return index < values.size() && values[index].toString().trim() == 'keep'
}
