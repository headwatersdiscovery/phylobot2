// ESR, step 1 of 3 (--esr; src/phylobot/pipeline/esr.py): checks the staged
// alignment, tree, model and settings, renames the sequences t1, t2, ..., fits the model
// on the tree (its branch lengths first, for an uploaded tree) and lists the tips to
// hold out.
//   in:  path(staged/)
//   out: path(fit/), the fit; path(fit/tips.txt), one tip ID per line
process ESR_FIT {
    tag 'fit'
    container params.containers.tools

    input:
    path staged, stageAs: 'staged'

    output:
    path 'fit', emit: fit
    path 'fit/tips.txt', emit: tips

    script:
    """
    phylobot-esr fit staged fit
    """
}
