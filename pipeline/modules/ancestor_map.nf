// Job-level ancestor_map.tsv: every combination's ancestors by content-derived IDs
// (phylobot.analogs), so the viewers match the same ancestor across combinations
// without reading every tree. One line per combination in `trees`: key, "rooted" or
// "unrooted", and its published tree's Newick.
process ANCESTOR_MAP {
    container params.containers.tools
    publishDir params.outdir, mode: 'copy'

    input:
    path trees

    output:
    path 'ancestor_map.tsv'

    script:
    """
    phylobot-ancestor-map ${trees} ancestor_map.tsv
    """
}
