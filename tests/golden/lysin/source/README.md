# Abalone sperm lysin: source data for the codon golden scenarios

Real data and fitted parameters behind `lysin/codon-easy`, `codon-selection` and
`codon-indels` (owner sign-off 2026-10-01). Lysin is a small, classic dataset for codon
models: 25 *Haliotis* species, 135 codons, with dN/dS varying strongly among sites.

## Files

| File | What |
|---|---|
| `lysin.fasta` | The real codon alignment (25 x 405 nt; one gap codon, in *H. rufescens*). |
| `lysin.nwk` | The published topology (Lee et al. 1995), no branch lengths. |
| `gy_g4.nwk` | `lysin.nwk` with branch lengths from IQ-TREE GY+F3X4+G4 (tree length 8.94). |
| `m3.nwk` | `lysin.nwk` with branch lengths from codeml M3 (tree length 8.94). |
| `../root.fasta` | Root of every simulation: *H. rufescens* lysin (the 1LIS crystal structure) without its gap codon, 134 codons. |

## Source

Copied from PAML's example data (`examples/lysin/lysin.nuc`, `lysin.trees`;
github.com/abacus-gene/paml, 2026-10-01), as used by Yang, Swanson & Vacquier (2000).
Names were cleaned (`__1_H._rufescens` -> `Haliotis_rufescens`, `_24_H._t.coccinea` ->
`Haliotis_tuberculata_coccinea`) and the `#1` branch mark removed from the tree.

- Lee YH, Ota T, Vacquier VD. *Positive selection is a general phenomenon in the
  evolution of abalone sperm lysin.* Mol Biol Evol 12(2): 231-238 (1995).
- Yang Z, Swanson WJ, Vacquier VD. *Maximum-likelihood analysis of molecular adaptation
  in abalone sperm lysin reveals variable selective pressures among lineages and sites.*
  Mol Biol Evol 17(10): 1446-1455 (2000).

**Licence:** the sequences are public GenBank data; the PAML repository is GPL-3.0 and
PhyloBot 2 is BSD-3. The owner decided (2026-10-01) that this copy is fine to bundle,
including as the New Job page's codon example.

## Fits (all on the published topology)

IQ-TREE 3.1.3 (`quay.io/biocontainers/iqtree:3.1.3--h8471819_0`):

    iqtree3 -s lysin.fasta -st CODON1 -te lysin.nwk -m GY+F3X4      # M0 check
    iqtree3 -s lysin.fasta -st CODON1 -te lysin.nwk -m GY+F3X4+G4   # -> gy_g4.nwk

PAML 4.10.10 (`quay.io/biocontainers/paml:4.10.10--h87e0c26_2`), codeml with
`seqtype = 1`, `CodonFreq = 2` (F3X4), `model = 0`, `NSsites = 0 3`, `ncatG = 3`,
`icode = 0`, `cleandata = 0`; M3's tree with branch lengths -> `m3.nwk`.

| Model | lnL | Parameters | Tree length |
|---|---|---|---|
| GY+F3X4 (IQ-TREE) | -4682.373 | omega 0.9303, kappa 1.5607 | 7.81 |
| M0 (codeml) | -4682.373 | kappa 1.5613 (the same model: the programs agree) | 7.81 |
| GY+F3X4+G4 (IQ-TREE) | -4490.576 | omega 1.0172, kappa 1.6247, alpha 1.0997 | 8.94 |
| M3 (codeml) | -4464.307 | p 0.32902 / 0.40231 / 0.26868, omega 0.08519 / 0.91120 / 3.06543, kappa 1.57610 | 8.94 |

F3X4 frequencies (observed; A C G T): position 1 0.34149 0.22514 0.26352 0.16984;
position 2 0.30553 0.16440 0.22393 0.30613; position 3 0.21154 0.24902 0.25960 0.27984.

## Checks on the simulations

- Mean pairwise amino-acid identity of the leaves: real 0.63; `codon-selection` 0.55
  (0.52-0.60 over 10 replicates); `codon-easy` 0.52; a single-dN/dS simulation 0.43.
- GY+F3X4 fitted to `codon-selection` replicates 302/303 gives omega 1.04 / 1.00 and
  tree length 8.3 / 8.4 (real data: 0.93, 7.8).
- In replicate 302 the three site classes (M3 rates 0.26 / 0.80 / 2.20) differ between
  leaves at 24% / 64% / 85% of codon pairs, 32% / 87% / 94% of those non-synonymous.
- `codon-indels`: 35% gaps, leaves 111-152 codons; no stop codons and no partial-codon
  gaps in any replicate (`generate.py` refuses both).
