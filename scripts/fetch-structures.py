#!/usr/bin/env python3
"""Download the backdrop's protein structure images from the RCSB PDB.

The landing and About pages rain these structures behind the content (``backdrop.js``).
They come from the RCSB PDB, whose content
is in the public domain (CC0 1.0), and are committed unmodified so the site works
offline and every run shows the same twenty structures: the browser keys out the white
background and tints what is left to suit the skin, so nothing has to be baked in.

To change the set, edit ENTRIES below and re-run; that is the whole protocol. Images for
entries you removed are deleted, and ``structures.json`` -- which is what the page reads --
is rewritten, so nothing else has to be touched.

    uv run python scripts/fetch-structures.py
"""

import json
import pathlib
import sys
import urllib.error
import urllib.request

# Twenty entries: five proteins whose sequences were themselves reconstructed from a
# phylogeny -- what PhyloBot is for -- and fifteen classics chosen for legible, open
# silhouettes (big assemblies turn into blobs at backdrop size).
ENTRIES = [
    "9XYT",
    "6Z1M",
    "3VVT",
    "5HIF",
    "6NEE",  # reconstructed ancestral proteins
    "1UBQ",
    "1CRN",
    "1MBN",
    "2LYZ",
    "1TIM",
    "1EMA",
    "2POR",
    "1BL8",
    "1PGB",
    "1ENH",
    "1SHG",
    "1BDD",
    "5PTI",
    "1HRC",
    "1ZAA",
]

IMAGE = "https://cdn.rcsb.org/images/structures/{lower}_assembly-1.jpeg"
ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{id}"
OUT = pathlib.Path(__file__).resolve().parent.parent / "src/phylobot/static/phylobot/backdrop"


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "phylobot2-backdrop"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    structures = []
    for entry in ENTRIES:
        lower = entry.lower()
        try:
            meta = json.loads(fetch(ENTRY.format(id=entry)))
            image = fetch(IMAGE.format(lower=lower))
        except (urllib.error.URLError, TimeoutError) as error:
            print(f"{entry}: {error}", file=sys.stderr)
            return 1
        (OUT / f"{lower}.jpeg").write_bytes(image)
        title = meta["struct"]["title"].strip()
        structures.append({"id": entry, "title": title, "file": f"{lower}.jpeg"})
        print(f"{entry}  {len(image) // 1024:>3} KB  {title[:60]}")

    (OUT / "structures.json").write_text(json.dumps(structures, indent=2) + "\n")

    keep = {entry["file"] for entry in structures}
    for stale in sorted(set(OUT.glob("*.jpeg")) - {OUT / name for name in keep}):
        stale.unlink()
        print(f"removed {stale.name} (no longer in ENTRIES)")

    print(f"\n{len(structures)} structures -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
