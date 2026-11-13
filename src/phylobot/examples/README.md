# Bundled example sequence sets

These FASTA files back the New Job page's **Load examples** drawer
(`phylobot.example_sets`, `views.example`). Each is an unaligned sequence set loaded
into the Sequences box; the `.subset.fasta` files are evenly-strided samples for a quick
run. The registry (`example_sets.EXAMPLES`) holds the label, sequence count, sequence
type and citation for each — keep it in step with the files here.

| Key | File(s) | n (full / subset) | Type | Source |
|---|---|---|---|---|
| `ime2` | `ime2.fasta` | 7 | protein | Bundled smoke-test set (2014 PhyloBot run); Hanson-Smith & Johnson 2016, PLoS Comput. Biol. 12:e1004976 |
| `gfp` | `gfp.fasta`, `gfp.subset.fasta` | 109 / 20 | coding DNA | Alieva et al. 2008, PLoS ONE 3:e2680 (degapped) |
| `gba1` | `gba1.fasta`, `gba1.subset.fasta` | 131 / 20 | protein | Curated for PhyloBot 2 |
| `maguk` | `maguk.fasta`, `maguk.subset.fasta` | 135 / 20 | protein | Anderson et al. 2016, eLife 5:e10147 (degapped; names cleaned and de-duplicated) |
| `pink1` | `pink1.fasta`, `pink1.subset.fasta` | 105 / 20 | protein | Curated for PhyloBot 2 |
| `rubisco` | `rubisco.fasta`, `rubisco.subset.fasta` | 81 / 20 | protein | Kacar et al. 2017, Geobiology 15:628–640 |
| `snca` | `snca.fasta`, `snca.subset.fasta` | 115 / 20 | protein | Curated for PhyloBot 2 |
| `src_abl` | `src_abl.fasta`, `src_abl.subset.fasta` | 76 / 20 | protein | Sennett & Theobald 2024, J. Mol. Evol. 92:181–206 |
| `sr238` | `sr238.fasta`, `sr238.subset.fasta` | 238 / 20 | protein | Eick et al. 2012, PLoS Genet. 8:e1003072; Bridgham, Carroll & Thornton 2006, Science 312:97–101 (SR238 steroid-receptor set) |
| `vatpase` | `vatpase.fasta`, `vatpase.subset.fasta` | 140 / 20 | protein | Finnigan et al. 2012, Nature 481:360–364 (V-ATPase ring subunits 3, 11, 16; degapped) |
| `lysin` | `lysin.fasta`, `lysin.subset.fasta` | 25 / 12 | coding DNA | Lee, Ota & Vacquier 1995, Mol. Biol. Evol. 12:231–238 |

## Subsets

The `.subset.fasta` files were made by evenly striding across the full set (always keeping
the first and last sequence) to roughly 20 sequences (12 for lysin). They are a quick-start
convenience, not a curated representative sample; the full set is authoritative.

Aligned sources (gfp, maguk, vatpase) were degapped per sequence. gfp's coding-DNA frames
were checked after degapping (all in-frame, no internal stop codons).
