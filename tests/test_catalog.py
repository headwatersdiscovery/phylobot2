import pytest

from phylobot import catalog


def test_catalog_size_and_names() -> None:
    models = catalog.all_models()
    assert len(models) == 28 * 9
    assert "LG+G4" in models and "Q.pfam+F+R4" in models and "WAG" in models
    assert len(set(models)) == len(models)


@pytest.mark.parametrize("model", ["LG+G4", "WAG", "mtREV+I+R4", "Q.yeast+F"])
def test_valid_models(model: str) -> None:
    assert catalog.is_model(model)


@pytest.mark.parametrize(
    "model", ["", "lg+G4", "LG+G", "LG+G4+I", "LG+R5", "NQ.pfam", "LG+G4;rm", "LG+G4 -redo"]
)
def test_invalid_models(model: str) -> None:
    assert not catalog.is_model(model)


def test_rectangular_selection_is_exact() -> None:
    ticked = ["LG+G4", "LG+F+G4", "WAG+G4", "WAG+F+G4"]
    candidates = catalog.modelfinder_candidates(ticked)
    assert candidates.args() == ["-mset", "LG,WAG", "-mrate", "G", "-mfreq", "FU,F"]
    assert sorted(candidates.models()) == sorted(ticked)
    assert candidates.count == 4


def test_non_rectangular_selection_widens() -> None:
    candidates = catalog.modelfinder_candidates(["LG+G4", "WAG+R4"])
    assert candidates.args() == ["-mset", "LG,WAG", "-mrate", "G,R", "-mfreq", "FU"]
    assert sorted(candidates.models()) == ["LG+G4", "LG+R4", "WAG+G4", "WAG+R4"]


def test_raxml_lacks_only_jttdcmut() -> None:
    assert catalog.raxml_missing(["LG+G4", "JTTDCMut+G4", "JTT+R4"]) == ["JTTDCMut+G4"]
    assert catalog.modelfinder_args([]) == []  # IQ-TREE's default set
    matrices = catalog.modelfinder_args([], raxml=True)
    assert matrices[0] == "-mset" and len(matrices) == 2  # default rate/frequency types
    assert "JTTDCMut" not in matrices[1].split(",") and "JTT" in matrices[1].split(",")
    assert len(matrices[1].split(",")) == len(catalog.MATRICES) - 1
    ticked = catalog.modelfinder_args(["LG+G4"], raxml=True)  # ticked: as without RAxML-NG
    assert ticked == catalog.modelfinder_args(["LG+G4"])


def test_candidates_follow_catalog_order() -> None:
    candidates = catalog.modelfinder_candidates(["WAG+R4", "LG+G4", "LG"])
    assert candidates.matrices == ["LG", "WAG"]
    assert candidates.rates == ["E", "G", "R"]


def test_no_models_rejected() -> None:
    with pytest.raises(ValueError):
        catalog.modelfinder_candidates([])
    with pytest.raises(ValueError, match="unknown model"):
        catalog.modelfinder_candidates(["LG+G4", "XX"])


def test_the_whole_grid_is_iqtrees_default_set() -> None:
    everything = catalog.all_models()
    assert catalog.is_default_set([]) and catalog.is_default_set(everything)
    assert catalog.modelfinder_args(everything) == []
    # The cross product is what counts: one cell less still covers every axis...
    assert catalog.is_default_set(everything[1:])
    # ...a whole matrix less does not.
    assert not catalog.is_default_set([m for m in everything if not m.startswith("LG")])
    without_raxml_missing = [m for m in everything if not m.startswith("JTTDCMut")]
    assert catalog.is_default_set(without_raxml_missing, raxml=True)
    assert not catalog.is_default_set(without_raxml_missing)
    assert catalog.modelfinder_args(without_raxml_missing, raxml=True) == catalog.modelfinder_args(
        [], raxml=True
    )


def test_every_rate_and_frequency_type_keeps_iqtrees_shortcuts() -> None:
    every = [m for m in catalog.all_models() if not m.startswith("mtREV")]
    candidates = catalog.modelfinder_candidates(every)
    assert candidates.shortcuts and candidates.args() == ["-mset", ",".join(candidates.matrices)]
    assert candidates.fitted == 27 * 22 * 2  # IQ-TREE: "test up to 1188 protein models"
    # One cell less (LG+F+R4 is the only F+R cell) still covers every type...
    assert catalog.modelfinder_candidates(every[1:]).shortcuts
    # ...one whole type less does not: then every flag is explicit, as before.
    no_i = [m for m in every if "+I" not in m]
    assert not catalog.modelfinder_candidates(no_i).shortcuts
    assert "-mrate" in catalog.modelfinder_candidates(no_i).args()


def test_codon_models_are_explicit_and_apart_from_protein_ones() -> None:
    models = catalog.all_codon_models()
    assert len(models) == len(set(models)) == 76
    assert catalog.is_codon_model("GY+F3X4+G4") and "GY+F3X4+G4" in catalog.DEFAULT_CODON_MODELS
    assert not catalog.is_codon_model("MG+F")  # IQ-TREE refuses MG with empirical codon freqs
    assert not catalog.is_codon_model("GY")  # frequencies always spelled out
    assert not catalog.is_model("GY+F3X4") and not catalog.is_codon_model("LG+G4")
    assert set(models).isdisjoint(catalog.all_models())


def test_codon_modelfinder_candidates() -> None:
    assert catalog.codon_modelfinder_args([]) == []  # IQ-TREE's default codon set
    assert catalog.codon_modelfinder_args(catalog.all_codon_models()) == []
    ticked = ["MG+F1X4+G4", "GY+F3X4", "ECMK07+F3X4+G4"]
    assert catalog.codon_modelfinder_args(ticked) == [
        "-mset",
        "MG,GY,ECMK07",
        "-mrate",
        "E,G",
        "-mfreq",
        "F3X4,F1X4",
    ]  # fmt: skip  -- the cross product, as IQ-TREE tests it (checked: 12 models)
    assert catalog.codon_candidate_count(ticked) == 12
    assert catalog.is_modelfinder("modelfinder_codon") and catalog.label("modelfinder_codon") == (
        "ModelFinder (codon)"
    )


MIXTURES = ["LG+C20+G4", "Q.pfam+C20+F+G4"]


@pytest.mark.parametrize("model", MIXTURES)
def test_mixtures_are_fixed_models_only(model: str) -> None:
    assert catalog.is_model(model) and catalog.is_mixture(model)
    assert model in catalog.all_fixed_models() and model not in catalog.all_models()
    assert catalog.raxml_missing([model]) == [model]
    assert catalog.candidate_models([model, "LG+G4"]) == ["LG+G4"]
    with pytest.raises(ValueError, match="not ModelFinder candidates"):
        catalog.modelfinder_candidates([model, "LG+G4"])


def test_mixtures_leave_modelfinders_grid_alone() -> None:
    assert len(catalog.all_fixed_models()) == 28 * (9 + 2)
    assert not any(catalog.is_mixture(m) for m in catalog.all_models())
    assert catalog.is_default_set(catalog.all_models())


@pytest.mark.parametrize("model", ["LG+C20", "LG+C10+G4", "LG+C60+G4", "LG+G4+C20", "LG+F+C20+G4"])
def test_only_the_offered_mixtures_are_models(model: str) -> None:
    assert not catalog.is_model(model) and not catalog.is_mixture(model)
