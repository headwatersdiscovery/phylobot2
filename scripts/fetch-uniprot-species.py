#!/usr/bin/env python3
"""Download UniProt's species codes and write the compact table PhyloBot bundles.

UniProt entry names end in a species code (IME2_YEAST, Q6C4M5_YARLI): `speclist.txt`,
UniProt's "controlled vocabulary of species", maps each code to an NCBI taxid and an
official name. PhyloBot reads species from sequence names offline with it
(`phylobot.species.codes`), so nothing is sent anywhere. UniProt distributes it under
CC BY 4.0 (https://www.uniprot.org/terms); the release is recorded in the table's
first line.

    uv run python scripts/fetch-uniprot-species.py

Writes src/phylobot/species/bundled/uniprot_species.tsv.gz: one line per code,
`code<TAB>kind<TAB>taxid<TAB>name<TAB>common names<TAB>synonyms` (several separated by
"|"), kind "species" for real organism codes and "group" for UniProt's virtual codes
(e.g. 9SACH, Saccharomycetaceae), which stand for a taxon above species.
"""

import gzip
import re
import sys
import urllib.request
from pathlib import Path

URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt"
OUT = Path(__file__).resolve().parent.parent / "src/phylobot/species/bundled/uniprot_species.tsv.gz"
ENTRY = re.compile(r"^([A-Z0-9]{1,5})\s+[ABEVO]\s+(\d{1,7}): N=(.+)$")
MORE = re.compile(r"^\s{10,}([CS])=(.+)$")  # a common name or synonym of the entry above


def main() -> None:
    text = urllib.request.urlopen(URL, timeout=120).read().decode()
    release = re.search(r"^Release:\s+(.+)$", text, re.M)
    if release is None or "CC BY 4.0" not in text:
        sys.exit("speclist.txt changed format (no release line, or a different licence)")
    kind = "species"
    entries: list[list[str]] = []
    for line in text.splitlines():
        if line.startswith('(2) "Virtual" codes'):
            kind = "group"
        if match := ENTRY.match(line):
            code, taxid, name = match.groups()
            entries.append([code, kind, taxid, name.strip(), "", ""])
        elif entries and (more := MORE.match(line)):
            column = 4 if more[1] == "C" else 5
            value = more[2].strip().replace("|", "/").replace("\t", " ")
            entries[-1][column] = "|".join(filter(None, [entries[-1][column], value]))
    rows = ["\t".join(entry) for entry in entries]
    if len(rows) < 20_000:
        sys.exit(f"only {len(rows)} codes read; refusing to write")
    header = f"# UniProt speclist.txt, release {release[1].strip()}; CC BY 4.0, UniProt Consortium"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0: the same release always gives the same bytes
    OUT.write_bytes(gzip.compress(("\n".join([header, *rows]) + "\n").encode(), mtime=0))
    print(f"{len(rows)} codes, release {release[1].strip()} -> {OUT}")


if __name__ == "__main__":
    main()
