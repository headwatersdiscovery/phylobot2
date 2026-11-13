"""Ancestors on structures, the pure part (phylobot.structures): reading PDB and mmCIF
files as hostile text, aligning a chain with the job's sequences, tying residues to
alignment columns, laying one structure on another, and the colour schemes. No network,
no database."""

import json
import math

import numpy as np
import pytest

from phylobot.homologs import remote as homolog_remote
from phylobot.homologs.remote import SearchError
from phylobot.structures import align, colours, mapping, pdb, remote, superpose

ONE_TO_THREE = {one: three for three, one in pdb.THREE_TO_ONE.items() if three != "MSE"} | {
    "M": "MET", "C": "CYS", "K": "LYS", "S": "SER", "T": "THR", "Y": "TYR", "P": "PRO",
}  # fmt: skip
Point = tuple[float, float, float]


def helix(index: int) -> Point:
    """Where the C-alpha of residue `index` of an ideal helix is (so that structures
    made here have a shape, and one can be laid on another)."""
    turn = math.radians(100.0) * index
    return (2.3 * math.cos(turn), 2.3 * math.sin(turn), 1.5 * index)


def moved(point: Point) -> Point:
    """`point` turned a quarter about x and shifted: another crystal's frame."""
    x, y, z = point
    return (x + 10.0, -z + 20.0, y - 5.0)


def atom(serial: int, name: str, chain: str, number: int, b: float = 20.0, **extra: str) -> str:
    """One C-alpha ATOM record, in the PDB format's fixed columns."""
    record = extra.get("record", "ATOM  ")
    atom_name = extra.get("atom", " CA ")
    icode = extra.get("icode", " ")
    x, y, z = helix(serial)
    return (
        f"{record}{serial:5d} {atom_name} {name:>3} {chain}{number:4d}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{b:6.2f}           C"
    )


def pdb_text(chains: dict[str, str], start: int = 1, b: float = 20.0) -> bytes:
    lines = ["TITLE     A TEST STRUCTURE"]
    serial = -1
    for chain, sequence in chains.items():
        for i, aa in enumerate(sequence):
            serial += 1
            lines.append(atom(serial, ONE_TO_THREE[aa], chain, start + i, b))
    return ("\n".join(lines) + "\nEND\n").encode()


def cif_text(sequence: str, chain: str = "A", b: float = 90.0, title: str = "A model") -> bytes:
    """A predicted model as AlphaFold 3 writes one: mmCIF, an N and a C-alpha per residue,
    pLDDT in the B-factor column, in a frame of its own (`moved`)."""
    names = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
        "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
        "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_asym_id",
        "pdbx_PDB_model_num",
    ]  # fmt: skip
    lines = ["data_model", "#", f"_struct.title '{title}'", "_ma_qa_metric.id 1", "#", "loop_"]
    lines += [f"_atom_site.{name}" for name in names]
    serial = 0
    for i, aa in enumerate(sequence):
        x, y, z = moved(helix(i))
        for element, name in (("N", "N"), ("C", "CA")):
            serial += 1
            lines.append(
                f"ATOM {serial} {element} {name} . {ONE_TO_THREE[aa]} {chain} 1 {i + 1} ? "
                f"{x:.3f} {y:.3f} {z:.3f} 1.00 {b:.2f} {i + 1} {chain} 1"
            )
    return ("\n".join(lines) + "\n#\n").encode()


def test_pdb_reads_chains_by_their_c_alphas() -> None:
    lines = [
        "TITLE     KINASE DOMAIN",
        *(atom(i, "ALA", "A", i) for i in range(1, 6)),
        atom(6, "MSE", "A", 6, record="HETATM"),  # selenomethionine: a residue
        atom(7, "GLY", "A", 6),  # the same residue again (alternate location): once
        atom(8, "GLY", "A", 7, icode="B"),
        atom(9, "XYZ", "A", 8),  # unknown residue name
        atom(10, " CA", "A", 900, record="HETATM", atom="CA  "),  # a calcium ion
        atom(11, "ALA", "A", 9, atom=" CB "),  # not a C-alpha
        "ATOM  bad line",
        *(atom(20 + i, "LYS", "B", i) for i in range(1, 4)),  # too short to be a chain
        "ENDMDL",
        *(atom(40 + i, "ALA", "C", i) for i in range(1, 9)),  # a second model: ignored
    ]
    parsed = pdb.parse("\n".join(lines).encode())
    assert parsed.title == "KINASE DOMAIN" and [c.id for c in parsed.chains] == ["A"]
    assert parsed.format == "pdb"
    chain = parsed.chains[0]
    assert chain.sequence == "AAAAAMGX"
    assert pdb.key("A", chain.residues[6]) == "A:7B" and chain.residues[0].b == 20.0
    assert chain.residues[0].z == 1.5  # the C-alpha's position is kept


def test_mmcif_is_read_as_3dmol_names_its_residues() -> None:
    data = cif_text("TAYIAKQRQ", b=91.5, title="Model of Node2")
    parsed = pdb.parse(data)
    assert parsed.format == "cif" and parsed.title == "Model of Node2"
    chain = parsed.chains[0]
    assert chain.sequence == "TAYIAKQRQ" and chain.residues[0].b == 91.5
    assert [pdb.key("A", r) for r in chain.residues[:2]] == ["A:1", "A:2"]
    assert (chain.residues[1].x, chain.residues[1].y) == pytest.approx(
        moved(helix(1))[:2], abs=1e-3
    )
    assert pdb.looks_predicted(data) and not pdb.looks_predicted(pdb_text({"A": "TAYIAKQ"}))
    # a quoted atom name is read; a chain of one residue is no chain; the first model only
    rows = (
        "ATOM 97 C 'CA' . ALA A 1 10 ? 1 2 3 1.00 5 10 A 1\n"
        "ATOM 98 C CA . ALA X 1 9 ? 1 2 3 1.00 5 99 X 1\n"
        "ATOM 99 C CA . ALA A 1 50 ? 1 2 3 1.00 5 50 A 2\n"
        "ATOM 100 C CA . ALA A 1 51 ? 1 2 3 1.00 5 51 A 1\n"
    )
    again = pdb.parse(
        (data.decode()[:-2] + rows + "#\nloop_\n_atom_site_anisotrop.id\n1\n").encode()
    )
    assert [c.id for c in again.chains] == ["A"]
    assert [r.number for r in again.chains[0].residues][-2:] == [9, 10]


@pytest.mark.parametrize(
    ("data", "why"),
    [
        (b"<!doctype html><title>403</title>403 Forbidden", "no protein chain"),
        (b"data_2SRC\nloop_\n_atom_site.id\n1\n", "no protein chain"),
        (b"ATOM\0\0", "not a text file"),
        (b"", "no protein chain"),
    ],
)
def test_what_is_not_a_structure_is_refused(data: bytes, why: str) -> None:
    with pytest.raises(pdb.StructureError, match=why):
        pdb.parse(data)


def test_local_alignment_finds_the_shared_stretch() -> None:
    assert (align.BLOSUM62 == align.BLOSUM62.T).all() and align.BLOSUM62[17, 17] == 11  # W
    found = align.local("MKTAYIAKQRQISFVKSHFSRQ", "GGGGTAYIAKQRSFVKSHFSRQPPP")
    # a's QI has no partner in b: one gap, not a run of mismatches
    assert found.pairs[0] == (2, 4) and found.pairs[-1] == (21, 21) and found.identity == 1.0
    assert (10, 12) not in found.pairs and (12, 12) in found.pairs
    assert align.local("", "MKT").pairs == () and align.local("WWWW", "GGGG").score == 0


def test_a_chain_is_tied_to_columns_through_its_closest_sequence() -> None:
    alignment = {
        "far": "MSTNPKPQRK--TKRNTNRRPQDVKF",
        "near": "MKTAYIAK--QRQISFVKSHFSRQ--",
        "nearer": "MKTAYLAK--QRQISFVRSHFSRQ--",
    }
    parsed = pdb.parse(pdb_text({"A": "TAYIAKQRQISFVKSH", "B": "GGGGGGGGGG"}, start=10))
    maps = mapping.chain_maps(parsed, mapping.ungapped(alignment))
    assert list(maps) == ["A"]  # chain B is like nothing in the job
    found = maps["A"]
    assert found.reference == "near" and found.identity == 1.0 and found.mapped == 16
    cols = mapping.columns(found, alignment["near"])
    assert cols[:7] == [2, 3, 4, 5, 6, 7, 10]  # the alignment's gap columns are stepped over

    row = "MKTAYLAKGGQ-QISFVKSHFSRQ--"  # an ancestor: L at column 6, residues at 9-10, gap at 12
    other = "MKTAYIAK--QRQISFIKSHFSRQ--"
    sites = mapping.sites(parsed.chains[0], cols, row, [0.5] * len(row), other)
    assert [s.key for s in sites[:2]] == ["A:10", "A:11"]
    by_key = {s.key: s for s in sites}
    assert by_key["A:13"].differs == "sub" and by_key["A:13"].change == "sub"  # I in the structure
    assert by_key["A:17"].state == "-" and by_key["A:17"].differs == "absent"
    assert by_key["A:17"].change == "gained" and by_key["A:17"].pp is None
    assert by_key["A:22"].change == "sub" and by_key["A:22"].differs == ""
    between = mapping.between(sites, row, other)
    assert [(b.after, b.before, b.count, b.other) for b in between] == [("A:15", "A:16", 2, 0)]


def test_one_structure_is_laid_on_another() -> None:
    # a small globule: C-alphas on a 4 x 4 x 3 grid, 3.8 apart
    fixed = [(3.8 * (i % 4), 3.8 * (i // 4 % 4), 3.8 * (i // 16)) for i in range(48)]
    moving = [moved(p) for p in fixed]
    fit = superpose.fit(moving, fixed)
    assert fit is not None and fit.rmsd < 1e-6 and (fit.kept, fit.pairs) == (48, 48)
    rotation, translation = np.array(fit.rotation), np.array(fit.translation)
    assert np.allclose(np.array(moving) @ rotation.T + translation, np.array(fixed), atol=1e-6)
    assert np.isclose(np.linalg.det(rotation), 1.0) and len(fit.matrix) == 12
    # a part that has moved (a domain, a loop) is left out of the fit rather than spoiling it
    bent = moving[:40] + [(x + 14.0, y, z) for x, y, z in moving[40:]]
    fit = superpose.fit(bent, fixed)
    assert fit is not None and fit.kept == 40 and fit.rmsd < 1e-6
    assert superpose.fit(moving[:3], fixed[:3]) is None  # too few pairs to mean anything


def test_colour_schemes() -> None:
    parsed = pdb.parse(pdb_text({"A": "TAYIAKQ"}, b=95.0))
    row, other = "TAYL-KQ", "TAYLAK-"
    cols: list[int | None] = [0, 1, 2, 3, 4, 5, None]
    sites = mapping.sites(parsed.chains[0], cols, row, [0.0, 1, 1, 1, 1, 1, 1], other)
    assert colours.available(compared=False, has_pp=False) == [
        "differs", "chemistry", "hydropathy", "confidence",
    ]  # fmt: skip
    assert "steps" in colours.available(True, True, has_steps=True)
    assert set(colours.ABOUT) == set(colours.SCHEMES)  # the page's "?" explains each
    pp = colours.paint("pp", sites, "coolwarm", False)
    assert pp.colours[0] == "#d73027" and pp.colours[1] == "#4575b4"
    assert pp.colours[4] == colours.ABSENT and pp.colours[6] == colours.UNMAPPED
    differs = colours.paint("differs", sites, "coolwarm", False).colours
    assert differs[:5] == [colours.SAME] * 3 + [colours.SUBSTITUTED, colours.ABSENT]
    changes = colours.paint("changes", sites, "coolwarm", False).colours
    assert changes[3:6] == [colours.SAME, colours.GAINED, colours.SAME]
    # a panel that paints the second sequence holds the two the other way round: the
    # change keeps its name
    swapped = mapping.sites(parsed.chains[0], cols, other, None, row)
    assert colours.paint("changes", swapped, "coolwarm", False, reverse=True).colours == changes
    # an extant sequence has no probabilities: it is not drawn as "no residue"
    observed = colours.paint("pp", swapped, "coolwarm", False)
    assert observed.colours[0] == colours.SAME and "Observed" in observed.legend[-3][0]
    # a path of three branches: each changed column in its last branch's colour
    steps = colours.Steps(["n1 → n2", "n2 → n3", "n3 → x"], {1: 0, 3: 2})
    stepped = colours.paint("steps", sites, "coolwarm", False, steps=steps)
    assert stepped.colours[0] == colours.SAME and stepped.colours[6] == colours.UNMAPPED
    assert stepped.colours[1] == colours.step_colour(0, 3) != stepped.colours[3]
    assert [label for label, _ in stepped.legend[:3]] == steps.branches
    # a predicted model's B-factors are pLDDT, in AlphaFold's bands
    assert colours.paint("confidence", sites, "coolwarm", True).colours[0] == "#0053d6"
    plain = colours.paint("confidence", sites, "coolwarm", False)
    assert "B-factor" in plain.legend[0][0]
    for scheme in colours.SCHEMES:  # every colour is one structure.js can use as it is
        painted = colours.paint(scheme, sites, "coolwarm", False, steps=steps)
        assert all(len(c) == 7 and c.startswith("#") for c in painted.colours)
        assert len(painted.colours) == len(sites) and painted.legend


def test_identifiers_are_checked_before_they_reach_a_url() -> None:
    assert remote.clean_id(" 2src ") == (remote.RCSB, "2SRC")
    assert remote.clean_id("p12931") == (remote.ALPHAFOLD, "P12931")
    assert remote.clean_id("AF-A0A024R161-F1") == (remote.ALPHAFOLD, "A0A024R161")
    for bad in ("", "2SRC/../x", "https://evil.example/2SRC.pdb", "SRC", "12345", "2SR?"):
        with pytest.raises(SearchError):
            remote.clean_id(bad)
    with pytest.raises(SearchError):
        remote.fetch("upload", "2SRC")
    # several at once, each checked, duplicates once
    assert remote.clean_ids("2src, 1OPL;P12931\n2SRC", 8) == [
        (remote.RCSB, "2SRC"), (remote.RCSB, "1OPL"), (remote.ALPHAFOLD, "P12931"),
    ]  # fmt: skip
    for bad, most in (("2SRC, ../x", 8), ("2SRC 1OPL 1FMK", 2), (" , ", 8)):
        with pytest.raises(SearchError):
            remote.clean_ids(bad, most)
    names = ["Homo_sapiens_P00519", "Monosiga_ovata_BAF02917", "sp|P12931|SRC_HUMAN", "AcaSchAR"]
    assert remote.accessions_in(names) == {
        "Homo_sapiens_P00519": "P00519",
        "sp|P12931|SRC_HUMAN": "P12931",
    }


def test_an_entry_too_large_for_pdb_format_comes_as_mmcif(monkeypatch: pytest.MonkeyPatch) -> None:
    asked: list[str] = []

    def answer(url: str, data: object = None) -> bytes:
        asked.append(url)
        if url.endswith(".pdb"):
            raise SearchError("the search service answered 404: not found")
        return cif_text("TAYIAKQRQ") if url.endswith(".cif") else b"{}"

    monkeypatch.setattr(homolog_remote, "request", answer)
    fetched = remote.fetch_rcsb("4v6x")
    assert pdb.parse(fetched.data).format == "cif"
    assert asked[:2] == [
        "https://files.rcsb.org/download/4V6X.pdb",
        "https://files.rcsb.org/download/4V6X.cif",
    ]


def test_alphafold_files_come_only_from_alphafold(monkeypatch: pytest.MonkeyPatch) -> None:
    asked: list[str] = []

    def answer(url: str, data: object = None) -> bytes:
        asked.append(url)
        return json.dumps([{"pdbUrl": "https://evil.example/files/x.pdb"}]).encode()

    monkeypatch.setattr(homolog_remote, "request", answer)
    with pytest.raises(SearchError, match="no model"):
        remote.fetch_alphafold("P12931")
    assert asked == ["https://alphafold.ebi.ac.uk/api/prediction/P12931"]  # nothing else fetched


def test_uniprot_structures_are_parsed_best_resolution_first() -> None:
    payload = {
        "results": [
            {"primaryAccession": "Q00001", "entryType": "UniProtKB unreviewed (TrEMBL)",
             "organism": {"scientificName": "Mus musculus"}, "sequence": {"length": 50}},
            {"primaryAccession": "P12931", "uniProtkbId": "SRC_HUMAN",
             "entryType": "UniProtKB reviewed (Swiss-Prot)",
             "organism": {"scientificName": "Homo sapiens"}, "sequence": {"length": 536},
             "genes": [{"geneName": {"value": "SRC"}}],
             "uniProtKBCrossReferences": [
                 {"database": "PDB", "id": "1FMK", "properties": [
                     {"key": "Method", "value": "X-ray"}, {"key": "Resolution", "value": "1.50 A"},
                     {"key": "Chains", "value": "A=86-536"}]},
                 {"database": "PDB", "id": "1HCS", "properties": [
                     {"key": "Method", "value": "NMR"}, {"key": "Resolution", "value": "-"}]},
                 {"database": "PDB", "id": "<script>", "properties": []},
                 {"database": "AlphaFoldDB", "id": "P12931"}]},
            {"primaryAccession": "not an accession"},
        ]
    }  # fmt: skip
    src, other = remote.parse(json.dumps(payload).encode())
    assert src.accession == "P12931" and src.alphafold and src.reviewed and src.gene == "SRC"
    assert [(e.pdb_id, e.resolution) for e in src.entries] == [("1FMK", "1.50 A"), ("1HCS", "")]
    assert other.accession == "Q00001" and not other.entries and not other.alphafold
    with pytest.raises(SearchError):
        remote.parse(b"<html>")
