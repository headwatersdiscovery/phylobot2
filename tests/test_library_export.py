"""The documents a library hands to the bench (phylobot.library_export, step L5b).

Pure text in, pure text out: a handful of members is enough to check every file.
"""

import pytest

from phylobot import library, library_export


def target(template: str, *dists: dict[str, float]) -> library.Target:
    return library.Target(
        source=library.Source(1, "Node3", "parsimony", weight=0.8),
        template=template,
        dists=[sorted(d.items(), key=lambda kv: -kv[1]) for d in dists],
    )


ONE = target("MK-Q", {"M": 0.9, "L": 0.1}, {"K": 0.6, "R": 0.4}, {"-": 1.0}, {"Q": 0.55, "E": 0.45})


def rows_for(axes: library.Axes) -> list[library_export.Row]:
    members = library.reduce(library.generate([ONE], axes, seed=1))
    source = library_export.Source(
        combination="mafft/iqtree/LG",
        label="MAFFT × IQ-TREE × LG",
        node="Node3",
        version="parsimony",
        version_label="Parsimony (shared presence)",
        short_label="Parsimony (shared)",
        weight=0.8,
    )
    return [
        library_export.Row(
            name=member.name("IME2"),
            member=member,
            source=source,
            template=ONE.template,
            selected=index < 2,
            mean_pp=member.mean_posterior(ONE),
        )
        for index, member in enumerate(members)
    ]


def test_fasta_holds_the_protein_that_would_be_ordered() -> None:
    rows = rows_for(library.Axes(threshold=0.2))
    text = library_export.fasta(rows)
    lines = text.splitlines()
    assert lines[0].startswith(">IME2-C1-N3-ML-P Most likely of Node3 · MAFFT")
    # Ungapped: the alignment's gap column is not part of the protein.
    assert lines[1] == "MKQ"
    assert text.endswith("\n") and "-" not in lines[1]
    assert library_export.fasta([]) == ""
    # Long sequences wrap, as every other FASTA the pipeline writes does.
    assert library_export.wrapped("A" * 130).count("\n") == 2


def test_manifest_names_the_sites_in_both_numberings() -> None:
    """A single-site alternate is named by its alignment column; the manifest adds the
    residue number, so the protein's own numbering is there too (naming decision)."""
    rows = rows_for(library.Axes(threshold=0.2, max_substitutions=1))
    text = library_export.manifest(rows)
    header, *lines = text.splitlines()
    assert header.split("\t") == list(library_export.MANIFEST_COLUMNS)
    assert len(lines) == len(rows)
    swapped = [line for line in lines if "K2R" in line.split("\t")[1]]
    assert swapped, "the K->R alternate should be a member"
    cells = dict(zip(library_export.MANIFEST_COLUMNS, swapped[0].split("\t"), strict=True))
    assert cells["kind"] == "single"
    assert cells["sites"] == "K2R (residue 2)"  # column 2 is the second residue here
    assert cells["node"] == "Node3" and cells["combination"] == "mafft/iqtree/LG"
    assert cells["changes_from_ml"] == "1"
    # Selected is spelled out rather than left as a bare boolean.
    assert {line.split("\t")[3] for line in lines} <= {"yes", "no"}


def test_a_gap_column_is_said_to_be_a_gap_rather_than_given_a_number() -> None:
    assert library.residue_number("MK-Q", 3) == 0  # the gap itself
    assert library.residue_number("MK-Q", 4) == 3  # Q is the third residue
    assert library.residue_number("MK-Q", 99) == 0  # beyond the row
    assert library.parse_swaps("K300RL22M") == [("K", 300, "R"), ("L", 22, "M")]
    assert library.parse_swaps("ALT") == []


def test_names_table_maps_each_name_to_what_it_is() -> None:
    rows = rows_for(library.Axes(threshold=0.2))
    lines = library_export.names(rows).splitlines()
    assert lines[0] == "name\tdescription"
    assert all(len(line.split("\t")) == 2 for line in lines[1:])
    assert any("AltAll of Node3" in line for line in lines[1:])


def methods_for(design: library.Design) -> str:
    rows = rows_for(design.axes)
    return library_export.paragraph(
        library_export.Methods(
            library_name="Library 2",
            code="IME2",
            job_name="Ime2 test",
            design=design,
            sources=[rows[0].source],
            members=len(rows),
            selected=2,
            ambiguous=2,
            coverage="3 of 4 states appear in a member; 1 does not.",
            seed=7,
            version="0.1.0",
        )
    )


def test_methods_paragraph_states_the_settings_and_the_caveats() -> None:
    """The caveats are what make the library defensible: draws come from the marginal
    per-site posteriors, and single-site alternates say nothing about epistasis."""
    text = methods_for(
        library.Design(
            axes=library.Axes(threshold=0.25, samples=10, max_substitutions=1),
            budget=24,
            strategy="spanning",
        )
    )
    assert "Library 2 (IME2)" in text and "Ime2 test" in text
    assert "0.25" in text and "10 sequence(s) drawn" in text
    assert "maximally dissimilar subset" in text  # the strategy, in words
    assert "budget of 24" in text
    assert "3 of 4 states appear in a member" in text  # the coverage sentence, verbatim
    assert "marginal per-site posterior probabilities" in text
    assert "no joint distribution over sites is available" in text
    assert "do not address epistasis" in text
    assert "Random seed 7; PhyloBot 0.1.0" in text
    assert "Eick GN" in text and "Hanson-Smith V" in text  # the two citations
    # Prose uses the short reconstruction name, not the menu's explanatory one.
    assert "Parsimony (shared)" in text and "shared presence" not in text


def test_methods_paragraph_leaves_out_axes_that_were_not_used() -> None:
    text = methods_for(
        library.Design(
            axes=library.Axes(alt_all=False, threshold=0.2), budget=4, strategy="supported"
        )
    )
    assert "AltAll" not in text
    assert "drawn per ancestor" not in text
    assert "single-site alternate(s) at those sites" not in text
    assert "the greatest number of reconstructions agreed on" in text


@pytest.mark.parametrize(
    ("name", "expected"),
    [("Library 2", "Library-2"), ("  ", "library"), ("a/b:c", "a-b-c"), ("x" * 90, "x" * 60)],
)
def test_download_names_are_safe(name: str, expected: str) -> None:
    assert library_export.safe_name(name) == expected


def test_documents_writes_every_file_an_export_keeps() -> None:
    rows = rows_for(library.Axes(threshold=0.2))
    design = library.Design(axes=library.Axes(threshold=0.2), budget=4)
    methods = library_export.Methods(
        library_name="L", code="C", job_name="J", design=design, sources=[rows[0].source],
        members=len(rows), selected=2, ambiguous=2, coverage="", seed=1, version="0.1.0",
    )  # fmt: skip
    written = library_export.documents(rows, methods, '{"seed": 1}')
    assert set(written) == set(library_export.FILES)
    assert all(text for text in written.values())
