# Backdrop structure images

Twenty protein structure images from the **RCSB PDB** (<https://www.rcsb.org>), used by
the landing and About pages' backdrop (`backdrop.js`, style "PDB structures"). Refresh or
change the set with:

    uv run python scripts/fetch-structures.py

Each file is `<pdb id>_assembly-1.jpeg` as RCSB serves it, downloaded unmodified;
`structures.json` lists each entry's PDB ID and title. Five of the twenty are proteins
whose sequences were themselves reconstructed from a phylogeny -- what PhyloBot is for --
and the rest are classics picked for legible, open silhouettes.

**Why they are committed rather than fetched at run time.** PhyloBot runs on a laptop
with no internet as a first-class mode, and pinning every input is a project rule. Linking
to `cdn.rcsb.org` would also send every visitor's address to a third party just to draw a
background.

**Licence.** RCSB PDB content is in the public domain under
[CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/); see RCSB's
[usage policies](https://www.rcsb.org/pages/policies). The images are redistributed here
unmodified. The browser keys out their white background and tints what is left to suit the
skin, so nothing is baked into the files.

Cite the PDB itself as:

> Berman HM, Westbrook J, Feng Z, Gilliland G, Bhat TN, Weissig H, Shindyalov IN,
> Bourne PE. The Protein Data Bank. *Nucleic Acids Research* 28(1): 235-242 (2000).
> [doi:10.1093/nar/28.1.235](https://doi.org/10.1093/nar/28.1.235)
