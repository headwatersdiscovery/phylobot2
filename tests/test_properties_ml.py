"""The neural property predictors (`phylobot.properties_ml`).

These need torch, metapredict and the ESM-2 weights, which only the phylobot-ml image
carries, so `make test` skips them and `make ml-test` runs them inside that image. The
window arithmetic and the sequence cleaning are pure and run everywhere.

Expectations are about *direction and order*, never exact floats:
neural scores drift across platforms, their orderings do not.
"""

import math
import random

import pytest

from phylobot import properties, properties_ml

# alpha-synuclein (UniProt P37840), disordered as a monomer; hen egg-white lysozyme
# (P00698, mature chain), a compact folded domain.
SYNUCLEIN = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAG"
    "SIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
)
LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSS"
    "DITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)


class TestCleaning:
    def test_gaps_and_nonstandard_letters_are_dropped(self) -> None:
        assert properties_ml.protein("m-KU*x") == "MK"


class TestWindows:
    @pytest.mark.parametrize("length", [1, 500, 1022, 1023, 1500, 2100, 5000])
    def test_kept_ranges_tile_the_sequence_once(self, length: int) -> None:
        spans = properties_ml.windows(length)
        covered = [i for _, _, a, b in spans for i in range(a, b)]
        assert covered == list(range(length))

    @pytest.mark.parametrize("length", [1023, 1500, 2100, 5000])
    def test_every_window_fits_the_model_and_holds_what_it_keeps(self, length: int) -> None:
        for start, end, keep_from, keep_to in properties_ml.windows(length):
            assert end - start <= properties_ml.ESM_WINDOW
            assert 0 <= start <= keep_from <= keep_to <= end <= length

    def test_a_short_sequence_is_one_window(self) -> None:
        assert properties_ml.windows(300) == [(0, 300, 0, 300)]


def mutate(sequence: str, share: float, seed: int) -> str:
    rng = random.Random(seed)
    letters = list(sequence)
    for i in rng.sample(range(len(letters)), int(len(letters) * share)):
        letters[i] = rng.choice(sorted(properties_ml.STANDARD))
    return "".join(letters)


class TestDisorder:
    def setup_method(self) -> None:
        pytest.importorskip("metapredict")

    def test_synuclein_is_disordered_and_lysozyme_is_not(self) -> None:
        synuclein, lysozyme = properties_ml.disorder_fraction([SYNUCLEIN, LYSOZYME])
        assert synuclein > 0.8
        assert lysozyme < 0.2

    def test_it_is_a_fraction(self) -> None:
        for value in properties_ml.disorder_fraction([SYNUCLEIN, LYSOZYME, "MK", ""]):
            assert 0.0 <= value <= 1.0

    def test_gaps_do_not_change_it(self) -> None:
        gapped = "--" + LYSOZYME[:40] + "---" + LYSOZYME[40:]
        assert properties_ml.disorder_fraction([gapped]) == properties_ml.disorder_fraction(
            [LYSOZYME]
        )

    def test_the_registry_entry_agrees_with_the_batch(self) -> None:
        prop = properties.REGISTRY["disorder"]
        assert prop.values([SYNUCLEIN, LYSOZYME, SYNUCLEIN]) == [
            prop.fn(SYNUCLEIN),
            prop.fn(LYSOZYME),
            prop.fn(SYNUCLEIN),
        ]


class TestEsm:
    def setup_method(self) -> None:
        pytest.importorskip("esm")

    def test_a_real_protein_beats_its_mutant_and_its_shuffle(self) -> None:
        """What makes it a plausibility score: a real protein scores
        above a degraded copy of itself, which held for 7 of 7 Ime2 sequences in the
        2026-10-21 probe.

        Deliberately *not* asserting mutant > shuffle: the probe found that order broke
        for 2 of 7 sequences, since a shuffle keeps composition and a point mutant keeps
        order. The score is a one-sided test, not a ruler of how degraded a sequence is.
        """
        shuffled = "".join(random.Random(1).sample(LYSOZYME, len(LYSOZYME)))
        real, mutant, shuffle = properties_ml.esm_log_likelihood(
            [LYSOZYME, mutate(LYSOZYME, 0.1, 2), shuffled]
        )
        assert real > mutant
        assert real > shuffle

    def test_it_is_a_log_probability_per_residue(self) -> None:
        for value in properties_ml.esm_log_likelihood([SYNUCLEIN, LYSOZYME]):
            assert -math.log(33) < value < 0.0  # better than uniform over the vocabulary

    def test_no_residues_is_nan(self) -> None:
        assert math.isnan(properties_ml.esm_log_likelihood(["---"])[0])

    def test_batching_does_not_change_a_score(self) -> None:
        alone = properties_ml.esm_log_likelihood([LYSOZYME])[0]
        together = properties_ml.esm_log_likelihood([SYNUCLEIN, LYSOZYME, "MKV"])[1]
        assert together == pytest.approx(alone, abs=1e-4)

    def test_a_protein_longer_than_the_context_is_scored(self) -> None:
        long = "".join(mutate(LYSOZYME, 0.3, seed) for seed in range(9))  # 1161 residues
        assert len(long) > properties_ml.ESM_WINDOW
        (value,) = properties_ml.esm_log_likelihood([long])
        assert -math.log(33) < value < 0.0
