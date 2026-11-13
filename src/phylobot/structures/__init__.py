"""Ancestors on protein structures.

A job can have structures attached (`phylobot.models.Structure`): PDB entries, AlphaFold
models or the user's own files. The structure page (`phylobot.structure_views`) draws one
in the browser (3Dmol.js) coloured by an ancestor: its residues, their posterior
probabilities, where it differs from the structure's own protein, or what changed between
two ancestors.

Nothing here changes a structure's coordinates. A structure's residues are tied to
alignment columns (`mapping`), and the column is what an ancestor's state belongs to.

- `pdb`: the residues of a PDB-format file (hostile input: a file the user uploads).
- `align`: local pairwise alignment (numpy).
- `mapping`: structure residue -> extant sequence -> alignment column.
- `colours`: the colour schemes, computed on the server so the legend and the picture agree.
- `remote`: RCSB, AlphaFold DB and UniProt's cross-references.
"""
