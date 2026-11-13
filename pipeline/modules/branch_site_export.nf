// The branch-site test's results (src/phylobot/pipeline/branch_site_export.py): the best
// start of Model A and of its null, the likelihood-ratio test, BEB sites, the marked
// tree under the original names and codeml's own files (codeml/: rerun as they ran),
// published to <outdir>/branch_site/ (summary.json to <outdir>/).
process BRANCH_SITE_EXPORT {
    tag "${meta.aligner}/${meta.tree}"
    container params.containers.tools
    publishDir path: params.outdir, mode: 'copy', saveAs: { filename ->
        filename == 'out/summary.json' ? 'summary.json'
            : filename.startsWith('out/') ? 'branch_site/' + filename.substring(4) : null
    }

    input:
    tuple val(meta), path(mlcs, stageAs: 'runs/*'), path(rsts, stageAs: 'runs/*'), path(controls, stageAs: 'runs/*'), path(alignment, stageAs: 'codons.fasta'), path(tree)
    path names

    output:
    tuple val(meta), path('out/**'), emit: results

    script:
    """
    phylobot-branch-site-export ${tree} runs out ${names} ${alignment}
    """
}
