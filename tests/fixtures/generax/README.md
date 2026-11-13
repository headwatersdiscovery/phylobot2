GeneRax 2.1.3 output (bioconda build hf316886_2, `--rec-model UndatedDL --strategy EVAL
--seed 1`, 2026-10-16) on the Ime2 sequences: `discordant_gene.nwk` is a deliberately
wrong gene tree (each Saccharomycetaceae gene beside a CTG-clade gene) reconciled with
`species.nwk` (the Ime2 species tree from the Open Tree of Life, polytomies resolved from
the real gene trees, internal nodes labelled spn001...). GeneRax infers 3 duplications
and 12 losses. `discordant.xml` is its RecPhyloXML, `stats.txt` its per-family stats.
