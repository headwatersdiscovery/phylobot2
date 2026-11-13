"""From sequences to species.tsv and a species tree.

`detect` fills each row's query from its header. `with_open_tree` resolves the queries
against the Open Tree of Life and fetches the induced species tree; `with_upload` maps
them onto the user's own species tree instead. Either way the rows say, one by one,
what was found and what still needs the user (`table.Row.status`), and only rows that
are "ok" are in the tree.

The user corrects a row by editing its `query` (a name, or "ott<id>" to pin an Open
Tree taxon) or `taxid`, and resolves again; `species`, `ott_id`, `status` and `note`
are always recomputed.
"""

import re
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from typing import Any

from phylobot.species import headers
from phylobot.species.opentree import Client, Match, Taxon
from phylobot.species.speciestree import Node, from_open_tree, from_upload, key
from phylobot.species.table import AMBIGUOUS, NOT_IN_TREE, OK, UNRESOLVED, Row

PINNED = re.compile(r"^ott(\d{1,10})$")
MAX_TAXID_LOOKUPS = 500  # taxids not found by name are looked up one request each
MAX_LISTED = 5
# Open Tree ranks at or below species: a taxon of one of these ranks that contains
# another of the job's taxa absorbs it (strains and subspecies join their species).
SPECIES_OR_BELOW = {
    "species", "subspecies", "varietas", "variety", "subvariety", "forma", "form",
    "strain", "isolate", "no rank - terminal",
}  # fmt: skip
CITATIONS = (
    "Open Tree of Life synthetic tree (release named in synth_id), https://tree.opentreeoflife.org",
    "Open Tree Taxonomy (version named in taxonomy_version)",
    "Hinchliff CE et al. (2015) Synthesis of phylogeny and taxonomy into a comprehensive "
    "tree of life. PNAS 112: 12764-12769. doi:10.1073/pnas.1423041112",
    "Redelings BD, Holder MT (2017) A supertree pipeline for summarizing phylogenetic and "
    "taxonomic information for millions of species. PeerJ 5: e3058. doi:10.7717/peerj.3058",
)


@dataclass
class Resolution:
    rows: list[Row]
    tree: Node | None  # over the species of the "ok" rows; None if fewer than two
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return self.tree is not None and all(row.status == OK for row in self.rows)


def detect(records: list[tuple[str, str]]) -> list[Row]:
    """One row per (sequence name, description), its query read from the header."""
    rows = []
    for name, description in records:
        found = headers.detect(name, description)
        if found is None:
            rows.append(Row(sequence=name, note="no species found in the header"))
        else:
            rows.append(
                Row(sequence=name, query=found.name, taxid=found.taxid, source=found.source,
                    note=found.note)
            )  # fmt: skip
    return rows


def two_words(query: str) -> str | None:
    words = query.split()
    return " ".join(words[:2]) if len(words) > 2 else None


def describe(t: Taxon) -> str:
    return f"{t.name} (ott{t.ott_id}, {t.rank or 'no rank'})"


def left_out(t: Taxon) -> tuple[Taxon, str, str]:
    """A taxon the synthetic tree leaves out: kept (its ott id goes in the row, so the
    page can offer its parent), but not "ok"."""
    return t, NOT_IN_TREE, f"{describe(t)} is left out of the synthetic tree"


def listing(items: list[str]) -> str:
    shown = items[:MAX_LISTED]
    more = len(items) - len(shown)
    return "; ".join(shown) + (f"; and {more} more" if more > 0 else "")


def choose(matches: list[Match], query: str = "") -> tuple[Taxon | None, str, str]:
    """(taxon, status, note) for one query's exact matches."""
    if not matches:
        return None, UNRESOLVED, "no taxon of that name in the Open Tree taxonomy"
    in_synthesis = [m for m in matches if m.taxon.in_synthesis]
    if not in_synthesis:
        flags = sorted({f for m in matches for f in m.taxon.flags})
        return (
            matches[0].taxon,
            NOT_IN_TREE,
            (
                f"{describe(matches[0].taxon)} is left out of the synthetic tree"
                + (f" ({', '.join(flags)})" if flags else "")
            ),
        )
    if len({m.taxon.ott_id for m in in_synthesis}) > 1:
        accepted = [m for m in in_synthesis if not m.synonym]
        if len({m.taxon.ott_id for m in accepted}) != 1:
            return (
                None,
                AMBIGUOUS,
                "several taxa have this name: "
                + listing([describe(m.taxon) for m in in_synthesis])
                + "; give a taxid or 'ott<id>' as the query",
            )
        in_synthesis = accepted
    match = in_synthesis[0]
    note = f"{query!r} is now {match.taxon.name!r}" if match.synonym and query else ""
    return match.taxon, OK, note


def with_open_tree(rows: list[Row], client: Client) -> Resolution:
    """Resolve every row against the Open Tree of Life and fetch the species tree."""
    synthesis = client.synthesis()
    queries = {row.query for row in rows if row.query and not PINNED.match(row.query)}
    queries |= {w for row in rows if row.query and (w := two_words(row.query))}
    matched = client.match_names(sorted(queries)) if queries else {}

    found: dict[str, tuple[Taxon | None, str, str]] = {}  # sequence -> (taxon, status, note)
    taxid_lookups: dict[int, Taxon | None] = {}
    for row in rows:
        if pinned := PINNED.match(row.query):
            t = client.taxon_by_ott(int(pinned[1]))
            if t is None:
                found[row.sequence] = (None, UNRESOLVED, f"no Open Tree taxon {row.query}")
            elif not t.in_synthesis:
                found[row.sequence] = left_out(t)
            else:
                found[row.sequence] = (t, OK, "")
            continue
        exact = matched.get(row.query, [])
        taxon, status, note = choose(exact, row.query)
        if status == UNRESOLVED and (shorter := two_words(row.query)):
            taxon, status, note = choose(matched.get(shorter, []), shorter)
            if status == OK:
                note = f"matched as {shorter!r}" + (f"; {note}" if note else "")
        if row.taxid is not None:
            # A taxid outranks the name: take the name's match only if it is that taxid.
            by_taxid = [m for m in exact if row.taxid in m.taxon.ncbi]
            if by_taxid:
                taxon, status, note = choose(by_taxid, row.query)
            else:
                if row.taxid not in taxid_lookups:
                    if len(taxid_lookups) >= MAX_TAXID_LOOKUPS:
                        note = "too many taxids to look up one by one"
                        found[row.sequence] = (None, UNRESOLVED, note)
                        continue
                    taxid_lookups[row.taxid] = client.taxon_from_ncbi(row.taxid)
                t = taxid_lookups[row.taxid]
                if t is None:
                    if status != OK:
                        note = f"NCBI taxid {row.taxid} is not in the Open Tree taxonomy"
                elif not t.in_synthesis:
                    taxon, status, note = left_out(t)
                else:
                    differs = taxon.name if taxon is not None and taxon.ott_id != t.ott_id else ""
                    taxon, status = t, OK
                    note = f"from taxid {row.taxid}" + (
                        f"; the name {row.query!r} would give {differs}" if differs else ""
                    )
        if not row.query and row.taxid is None:
            note = row.note or "no species given"
        found[row.sequence] = (taxon, status, note)

    unresolved = sorted(
        {row.query for row in rows if row.query and found[row.sequence][1] == UNRESOLVED}
    )
    if unresolved:
        suggestions = client.match_names(unresolved, approximate=True)
        for row in rows:
            taxon, status, note = found[row.sequence]
            if status == UNRESOLVED and (near := suggestions.get(row.query)):
                close = list(dict.fromkeys(m.taxon.name for m in near if m.taxon.in_synthesis))
                if close:
                    found[row.sequence] = (taxon, status, f"{note}; did you mean {listing(close)}?")

    # Species names: Open Tree's names, made unique where two taxa share one.
    taxa = {t.ott_id: t for t, status, _ in found.values() if t is not None and status == OK}
    names: dict[str, list[int]] = {}
    for t in taxa.values():
        names.setdefault(t.name, []).append(t.ott_id)
    species_of = {
        t.ott_id: t.name if len(names[t.name]) == 1 else f"{t.name} ott{t.ott_id}"
        for t in taxa.values()
    }

    tree: Node | None = None
    provenance: dict[str, Any] = {
        "source": "opentree",
        **synthesis,
        "fetched_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "citations": list(CITATIONS),
    }
    if len(taxa) >= 2:
        induced = client.induced_subtree(sorted(taxa))
        tree, nested = from_open_tree(induced, species_of)
        # A strain or subspecies inside a species of the job joins that species; a
        # species inside a genus of the job leaves the genus's sequences to the user.
        within = {
            inner: outer
            for outer, inners in nested.items()
            if taxa[outer].rank in SPECIES_OR_BELOW
            for inner in inners
        }
        if within:
            outermost: dict[int, int] = {}
            for inner, outer in within.items():
                seen = {inner}
                while outer in within and outer not in seen:  # two taxa on one node loop
                    seen.add(outer)
                    outer = within[outer]
                if outer not in seen:
                    outermost[inner] = outer
            within = outermost
            tree, nested = from_open_tree(
                induced, {i: n for i, n in species_of.items() if i not in within}
            )
            for inner, outer in within.items():
                species_of[inner] = species_of[outer]
        provenance |= {
            "requested_ott_ids": sorted(taxa),
            "pruned_ott_ids": list(induced.pruned),
            "broken": {f"ott{k}": v for k, v in induced.broken.items()},
            "supporting_studies": list(induced.studies),
        }
        for sequence, (t, status, note) in list(found.items()):
            if t is None or status != OK:
                continue
            if t.ott_id in induced.pruned:
                note = f"{describe(t)} is pruned from the synthetic tree"
                found[sequence] = (t, NOT_IN_TREE, note)
            elif t.ott_id in nested:
                inside = [species_of[i] for i in nested[t.ott_id]]
                found[sequence] = (
                    None,
                    AMBIGUOUS,
                    (
                        f"{describe(t)} contains other species of this job "
                        f"({listing(sorted(set(inside)))}); give this sequence's species"
                    ),
                )
            elif t.ott_id in within:
                extra = f"placed with {species_of[t.ott_id]}, which contains it"
                found[sequence] = (t, OK, f"{note}; {extra}" if note else extra)
            elif t.ott_id in induced.broken:
                extra = "not monophyletic in the synthetic tree; placed at its stand-in"
                found[sequence] = (t, OK, f"{note}; {extra}" if note else extra)
        if len(tree.leaves()) < 2:
            tree = None

    out = []
    for row in rows:
        t, status, note = found[row.sequence]
        out.append(
            replace(
                row,
                ott_id=t.ott_id if t is not None else None,
                species=species_of[t.ott_id] if t is not None and status == OK else "",
                status=status,
                note=note,
            )
        )
    return Resolution(out, tree, provenance)


def with_upload(rows: list[Row], newick: str) -> Resolution:
    """Map each row's query onto the leaves of the user's species tree (names compared
    with underscores as spaces)."""
    tree = from_upload(newick, {row.query for row in rows if row.query})
    leaves = {key(leaf): leaf for leaf in tree.leaves()}
    out = []
    for row in rows:
        if not row.query:
            out.append(replace(row, ott_id=None, species="", status=UNRESOLVED,
                               note=row.note or "no species given"))  # fmt: skip
        elif key(row.query) not in leaves:
            out.append(replace(row, ott_id=None, species="", status=NOT_IN_TREE,
                               note="not a leaf of the uploaded species tree"))  # fmt: skip
        else:
            species = leaves[key(row.query)]
            out.append(replace(row, ott_id=None, species=species, status=OK, note=""))
    provenance = {
        "source": "upload",
        "fetched_at": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    return Resolution(out, tree, provenance)
