# Bundled species data

- `uniprot_species.tsv.gz`: UniProt's controlled vocabulary of species (speclist.txt),
  code, kind, NCBI taxid, official name, common names and synonyms, read by
  `phylobot.species.known` so that UniProt entry names (IME2_YEAST), binomials and common
  names ("zebrafish") give their species offline. Copyright the UniProt
  Consortium, distributed under CC BY 4.0 (https://www.uniprot.org/terms); the release
  is on the file's first line. Refresh with `uv run python scripts/fetch-uniprot-species.py`
  (byte-identical for the same release).
