"""The protein chains of a structure file, PDB format or mmCIF: each residue's number,
name, B-factor and C-alpha position.

The file is hostile input (the user may upload it, and a fetched one is still somebody
else's text): PDB is read by fixed columns and mmCIF by its `_atom_site` table alone,
nothing in it is trusted to be well formed, and anything odd is skipped or refused with a
`StructureError` that names the reason. Only what the page needs is kept: residues that
have a C-alpha atom, in file order, of the first model.

mmCIF is what AlphaFold 3 and AlphaFold Server write (they give no PDB file), and the only
format RCSB has for its largest entries. A residue is named as 3Dmol.js names it, or the
page's colours would land on the wrong residues: the author's chain and number
(`auth_asym_id`, `auth_seq_id`), and no insertion code, which 3Dmol's mmCIF reader drops.
"""

import re
from dataclasses import dataclass
from typing import Any

MAX_BYTES = 25 * 1024 * 1024
MAX_RESIDUES = 20000
MAX_CHAINS = 200
PDB, CIF = "pdb", "cif"
FORMATS = (PDB, CIF)
MIN_RESIDUES = 5  # shorter chains are peptides and tags, not something to map

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common modified residues, read as the residue they were made from.
    "MSE": "M", "SEC": "C", "PYL": "K", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C",
    "HYP": "P", "MLY": "K", "CME": "C", "KCX": "K", "LLP": "K",
}  # fmt: skip


class StructureError(Exception):
    pass


@dataclass(frozen=True)
class Residue:
    number: int
    icode: str  # insertion code, "" for none
    name: str  # three-letter name as the file gives it
    aa: str  # one letter; X for a residue this does not know
    b: float  # the C-alpha's B-factor (pLDDT in a predicted model)
    x: float = 0.0  # the C-alpha's position
    y: float = 0.0
    z: float = 0.0

    @property
    def label(self) -> str:
        return f"{self.number}{self.icode}"


@dataclass(frozen=True)
class Chain:
    id: str
    residues: tuple[Residue, ...]

    @property
    def sequence(self) -> str:
        return "".join(r.aa for r in self.residues)


@dataclass(frozen=True)
class Parsed:
    title: str
    chains: tuple[Chain, ...]
    format: str = PDB

    def chain(self, chain_id: str) -> Chain | None:
        return next((c for c in self.chains if c.id == chain_id), None)


def key(chain_id: str, residue: Residue) -> str:
    """How the page and structure.js name one residue ("A:295", "A:52B")."""
    return f"{chain_id}:{residue.label}"


def _text(data: bytes) -> str:
    if len(data) > MAX_BYTES:
        raise StructureError(f"the file is over {MAX_BYTES // (1024 * 1024)} MB")
    if b"\0" in data:
        raise StructureError("this is not a text file in PDB format")
    return data.decode("ascii", errors="replace")


def is_cif(text: str) -> bool:
    return text.lstrip().startswith("data_") or "\n_atom_site." in text[:2_000_000]


def looks_predicted(data: bytes) -> bool:
    """Whether the file says it is a predicted model (its B-factors are then pLDDT):
    ModelCIF's quality table, or the name of a predictor near the top."""
    head = data[:200_000].lower()
    return b"_ma_qa_metric" in data or any(
        name in head for name in (b"alphafold", b"colabfold", b"esmfold", b"openfold")
    )


class _Reader:
    """Collects residues as either format's reader finds C-alpha atoms."""

    def __init__(self, icodes: bool) -> None:
        self.icodes = icodes
        self.found: dict[str, list[Residue]] = {}
        self.seen: set[tuple[str, int, str]] = set()
        self.total = 0

    def add(self, chain_id: str, number: int, icode: str, name: str, b: float, xyz: Any) -> None:
        if not self.icodes:
            icode = ""
        if not (chain_id.isalnum() or chain_id == "_") or len(chain_id) > 4:
            return
        if icode and not icode.isalnum():
            return
        if (chain_id, number, icode) in self.seen:
            return  # an alternate location of a residue already read
        self.seen.add((chain_id, number, icode))
        if chain_id not in self.found and len(self.found) >= MAX_CHAINS:
            return
        self.total += 1
        if self.total > MAX_RESIDUES:
            raise StructureError(f"the structure has over {MAX_RESIDUES} residues")
        residue = Residue(number, icode, name[:3], THREE_TO_ONE.get(name, "X"), b, *xyz)
        self.found.setdefault(chain_id, []).append(residue)

    def chains(self) -> tuple[Chain, ...]:
        return tuple(
            Chain(chain_id, tuple(residues))
            for chain_id, residues in self.found.items()
            if len(residues) >= MIN_RESIDUES
        )


def _number(text: str) -> float:
    try:
        value = float(text)
    except ValueError:
        return 0.0
    return value if -1e6 < value < 1e6 else 0.0  # also turns nan and inf away


def _parse_pdb(text: str) -> Parsed:
    title: list[str] = []
    reader = _Reader(icodes=True)
    for line in text.splitlines():
        record = line[:6]
        if record == "TITLE ":
            title.append(line[10:80].strip())
        elif record == "ENDMDL":
            break  # the first model only (NMR ensembles)
        elif record in ("ATOM  ", "HETATM") and line[12:16].strip() == "CA":
            name = line[17:20].strip().upper()
            if record == "HETATM" and name not in THREE_TO_ONE:
                continue  # a calcium ion, not a C-alpha
            try:
                number = int(line[22:26])
            except ValueError:
                continue
            xyz = (_number(line[30:38]), _number(line[38:46]), _number(line[46:54]))
            reader.add(
                line[21:22].strip() or "_", number, line[26:27].strip(), name,
                _number(line[60:66]), xyz,
            )  # fmt: skip
    return Parsed(" ".join(part for part in title if part)[:300], reader.chains(), PDB)


# A CIF value: quoted (a quote ends only before white space), or a run of non-blanks.
_TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'|"(?:[^"]|"(?=\S))*"|\S+""")


def _tokens(line: str) -> list[str]:
    return [
        t[1:-1] if len(t) >= 2 and t[0] in "'\"" and t[-1] == t[0] else t
        for t in _TOKEN.findall(line)
    ]


def _parse_cif(text: str) -> Parsed:
    lines = text.splitlines()
    title = ""
    names: list[str] = []
    start = len(lines)
    for i, line in enumerate(lines):
        if line.startswith("_struct.title"):
            rest = line[len("_struct.title") :].strip()
            if rest:
                title = " ".join(_tokens(rest))
            elif i + 1 < len(lines) and lines[i + 1].startswith(";"):
                title = lines[i + 1][1:].strip()
        elif (
            line.strip() == "loop_"
            and lines[i + 1 : i + 2]
            and lines[i + 1].startswith("_atom_site.")
        ):
            j = i + 1
            while j < len(lines) and lines[j].startswith("_atom_site."):
                names.append(lines[j].split()[0][len("_atom_site.") :])
                j += 1
            start = j
            break
    column = {name: i for i, name in enumerate(names)}

    def first(*wanted: str) -> int | None:
        return next((column[name] for name in wanted if name in column), None)

    atom_at = first("label_atom_id", "auth_atom_id")
    comp_at = first("auth_comp_id", "label_comp_id")
    chain_at = first("auth_asym_id", "label_asym_id")
    number_at = first("auth_seq_id", "label_seq_id")
    xyz_at = [column.get(f"Cartn_{axis}") for axis in "xyz"]
    if None in (atom_at, comp_at, chain_at, number_at, *xyz_at):
        return Parsed(title[:300], (), CIF)
    assert atom_at is not None and comp_at is not None and chain_at is not None
    assert number_at is not None
    group_at, b_at, model_at = (
        column.get("group_PDB"),
        column.get("B_iso_or_equiv"),
        column.get("pdbx_PDB_model_num"),
    )
    reader = _Reader(icodes=False)
    model = None
    for line in lines[start:]:
        if not line or line[0] in "#_" or line.startswith(("loop_", "data_")):
            break
        if "CA" not in line:
            continue  # only C-alphas are wanted: most rows are never split
        row = _tokens(line)
        if len(row) != len(names) or row[atom_at] != "CA":
            continue
        if model_at is not None:
            if model is None:
                model = row[model_at]
            elif row[model_at] != model:
                break  # the first model only
        name = row[comp_at].upper()
        if group_at is not None and row[group_at] == "HETATM" and name not in THREE_TO_ONE:
            continue  # a calcium ion
        try:
            number = int(row[number_at])
        except ValueError:
            continue
        xyz = tuple(_number(row[at]) for at in xyz_at if at is not None)
        b = _number(row[b_at]) if b_at is not None else 0.0
        reader.add(row[chain_at], number, "", name, b, xyz)
    return Parsed(title[:300], reader.chains(), CIF)


def parse(data: bytes) -> Parsed:
    """The file's protein chains; StructureError if it has none or is neither format."""
    text = _text(data)
    parsed = _parse_cif(text) if is_cif(text) else _parse_pdb(text)
    if not parsed.chains:
        raise StructureError("no protein chain was found (is this a PDB or mmCIF file?)")
    return parsed
