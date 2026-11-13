"""Search hits: one record per database entry, whatever the search's output format.

`parse_ebi_xml` reads EBI's XML for a BLAST search (EBIApplicationResult);
`parse_ncbi_xml` reads NCBI's BLAST XML (BlastOutput). The XML comes from a remote
service, so it is size-capped; ElementTree never fetches external entities or DTDs,
and the bundled expat refuses entity-expansion bombs.
"""

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass

# UniRef cluster names: "<name> n=<members> Tax=<common taxon> TaxID=<id> RepID=<entry>".
UNIREF_DESCRIPTION = re.compile(
    r"^(?P<name>.*?)\s+n=(?P<members>\d+)\s+Tax=(?P<taxon>.*?)\s+TaxID=(?P<taxid>\d+)"
    r"\s+RepID=(?P<repid>\S+)\s*$"
)
# UniProtKB (EBI's Swiss-Prot hits): "V-type proton ATPase subunit c OS=Saccharomyces
# cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=VMA3 PE=1 SV=1".
UNIPROT_DESCRIPTION = re.compile(r"^(?P<name>.*?)\s+OS=(?P<taxon>.*?)\s+OX=(?P<taxid>\d+)(\s|$)")


@dataclass(frozen=True)
class Hit:
    id: str  # the database's ID, e.g. UniRef90_P32581
    name: str  # protein name
    taxon: str  # for UniRef, the common taxon of the cluster's members
    taxid: int | None
    members: int | None  # sequences in the cluster (UniRef only)
    length: int
    evalue: float  # best HSP
    bits: float  # best HSP
    identity: float  # % identical over the best HSP's alignment
    coverage: float  # % of the query covered by all HSPs together


@dataclass(frozen=True)
class SearchResult:
    program: str  # e.g. "BLASTP 2.16.0+"
    database: str  # e.g. "uniref90"
    database_created: str  # the release's date, as EBI states it; NCBI does not say
    database_sequences: int
    query_length: int
    max_hits: int  # the search's hit limit
    hits: list[Hit]

    @property
    def truncated(self) -> bool:
        """True when the search hit its limit, so weaker homologs were never reported."""
        return len(self.hits) >= self.max_hits


def parse_description(description: str) -> tuple[str, str, int | None, int | None]:
    """(name, taxon, taxid, members) from a UniRef or UniProtKB description (members only
    for UniRef); other descriptions are kept whole as the name."""
    match = UNIREF_DESCRIPTION.match(description.strip())
    if not match:
        entry = UNIPROT_DESCRIPTION.match(description.strip())
        if entry:
            return entry["name"], entry["taxon"], int(entry["taxid"]), None
        return description.strip(), "", None, None
    return (
        match["name"],
        match["taxon"],
        int(match["taxid"]),
        int(match["members"]),
    )


def coverage(query_length: int, intervals: list[tuple[int, int]]) -> float:
    """% of query positions inside any of the (start, end) intervals (1-based, inclusive;
    either order, as BLAST writes reverse ranges for some programs)."""
    if query_length <= 0:
        return 0.0
    covered: set[int] = set()
    for start, end in intervals:
        low, high = sorted((start, end))
        covered.update(range(max(low, 1), min(high, query_length) + 1))
    return 100.0 * len(covered) / query_length


MAX_XML_BYTES = 100_000_000
NS = {"e": "http://www.ebi.ac.uk/schema"}


class ParseError(Exception):
    pass


def _find(element: ET.Element, path: str) -> ET.Element:
    found = element.find(path, NS)
    if found is None:
        raise ParseError(f"BLAST XML lacks {path}")
    return found


def _number(element: ET.Element, path: str) -> float:
    text = _find(element, path).text
    try:
        return float(text or "")
    except ValueError as exc:
        raise ParseError(f"BLAST XML: {path} is not a number: {text!r}") from exc


def _parse(content: bytes) -> ET.Element:
    if len(content) > MAX_XML_BYTES:
        raise ParseError("BLAST XML is too large")
    try:
        return ET.fromstring(content)
    except ET.ParseError as exc:
        raise ParseError(f"BLAST XML is malformed: {exc}") from exc


def parse_ebi_xml(content: bytes) -> SearchResult:
    """Read EBI's XML result for a blastp search against a single database."""
    root = _parse(content)
    header = _find(root, "e:Header")
    params = _find(header, "e:parameters")
    databases = _find(params, "e:databases")
    database = _find(databases, "e:database")
    query = _find(params, "e:sequences/e:sequence")
    query_length = int(query.get("length", "0"))
    if query_length <= 0:
        raise ParseError("BLAST XML: the query has no length")

    hits = []
    for element in root.iterfind("e:SequenceSimilaritySearchResult/e:hits/e:hit", NS):
        hsps = element.findall("e:alignments/e:alignment", NS)
        if not hsps:
            continue
        best = max(hsps, key=lambda hsp: _number(hsp, "e:bits"))
        intervals = []
        for hsp in hsps:
            span = _find(hsp, "e:querySeq")
            intervals.append((int(span.get("start", "0")), int(span.get("end", "0"))))
        name, taxon, taxid, members = parse_description(element.get("description", ""))
        hits.append(
            Hit(
                id=element.get("id", ""),
                name=name,
                taxon=taxon,
                taxid=taxid,
                members=members,
                length=int(element.get("length", "0")),
                evalue=_number(best, "e:expectation"),
                bits=_number(best, "e:bits"),
                identity=_number(best, "e:identity"),
                coverage=coverage(query_length, intervals),
            )
        )
    program = _find(header, "e:program")
    return SearchResult(
        program=program.get("version", program.get("name", "")),
        database=database.get("name", ""),
        database_created=database.get("created", ""),
        database_sequences=int(databases.get("sequences", "0")),
        query_length=query_length,
        max_hits=int(_find(params, "e:alignments").text or "0"),
        hits=hits,
    )


# NCBI's descriptions: "<name> [<organism>]", several joined by " >" when identical
# sequences (or, in ClusteredNR, a cluster's members) share an entry. The first is used.
NCBI_ORGANISM = re.compile(r"^(?P<name>.*?)\s*\[(?P<organism>[^\[\]]+)\]\s*$")


def _text(element: ET.Element, path: str) -> str:
    found = element.find(path)
    if found is None or found.text is None:
        raise ParseError(f"BLAST XML lacks {path}")
    return found.text


def _int(element: ET.Element, path: str) -> int:
    text = _text(element, path)
    try:
        return int(text)
    except ValueError as exc:
        raise ParseError(f"BLAST XML: {path} is not a whole number: {text!r}") from exc


def _float(element: ET.Element, path: str) -> float:
    text = _text(element, path)
    try:
        return float(text)
    except ValueError as exc:
        raise ParseError(f"BLAST XML: {path} is not a number: {text!r}") from exc


def parse_ncbi_xml(content: bytes, max_hits: int) -> SearchResult:
    """Read NCBI's BLAST XML (one query). `max_hits` is the limit the search asked for;
    the XML does not record it."""
    root = _parse(content)
    if root.tag != "BlastOutput":
        raise ParseError("BLAST XML lacks BlastOutput")
    query_length = _int(root, "BlastOutput_query-len")
    if query_length <= 0:
        raise ParseError("BLAST XML: the query has no length")
    iteration = root.find("BlastOutput_iterations/Iteration")
    if iteration is None:
        raise ParseError("BLAST XML lacks BlastOutput_iterations/Iteration")
    hits = []
    for element in iteration.iterfind("Iteration_hits/Hit"):
        hsps = element.findall("Hit_hsps/Hsp")
        if not hsps:
            continue
        best = max(hsps, key=lambda hsp: _float(hsp, "Hsp_bit-score"))
        aligned = _int(best, "Hsp_align-len")
        intervals = [(_int(h, "Hsp_query-from"), _int(h, "Hsp_query-to")) for h in hsps]
        description = _text(element, "Hit_def").split(" >")[0].strip()
        match = NCBI_ORGANISM.match(description)
        hits.append(
            Hit(
                id=_text(element, "Hit_accession"),
                name=match["name"] if match else description,
                taxon=match["organism"] if match else "",
                taxid=None,
                members=None,
                length=_int(element, "Hit_len"),
                evalue=_float(best, "Hsp_evalue"),
                bits=_float(best, "Hsp_bit-score"),
                identity=100.0 * _int(best, "Hsp_identity") / aligned if aligned else 0.0,
                coverage=coverage(query_length, intervals),
            )
        )
    stat = iteration.find("Iteration_stat/Statistics")
    return SearchResult(
        program=_text(root, "BlastOutput_version"),
        database=_text(root, "BlastOutput_db"),
        database_created="",
        database_sequences=_int(stat, "Statistics_db-num") if stat is not None else 0,
        query_length=query_length,
        max_hits=max_hits,
        hits=hits,
    )
