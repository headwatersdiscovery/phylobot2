"""Species of each sequence, and the species tree they sit on.

A gene tree is reconciled against a species tree, which needs two things the sequences
do not carry on their own: which species each sequence comes from (`table`, written as
species.tsv), and a rooted tree of those species (`speciestree`). The species are read
from FASTA headers (`headers`) or given by the user, resolved against the Open Tree of
Life (`opentree`, `resolve`), or taken from a species tree the user uploads.

Open Tree's synthesis leaves many nodes unresolved, and GeneRax refuses polytomies, so
`polytomies` resolves them, guided by the job's own gene trees by default (owner
decision 2026-10-16; the user can choose a seeded arbitrary resolution instead, or
upload a resolved tree). That step needs gene trees, so it runs after the tree searches;
everything else runs before the job (in the web tier, or by hand with
`phylobot-species`), so the pipeline never calls a remote service.
"""
