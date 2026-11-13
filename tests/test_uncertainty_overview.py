"""The Uncertainty page's overview figures (phylobot.uncertainty_overview): the map places
reconstructions by the distances between their posterior means, the explainer's bar
splits the job's total, and the strip names the largest source per ancestor."""

from __future__ import annotations

from typing import Any

import numpy as np

from phylobot import uncertainty_overview as o
from phylobot.uncertainty_views import COMPONENTS


def summary_with(means: np.ndarray, own: float) -> dict[str, Any]:
    """Cells at these mean points; D_ij = own + |m_i - m_j|^2 (D_ii = own)."""
    diff = ((means[:, None, :] - means[None, :, :]) ** 2).sum(-1)
    labels = ["A", "A", "B", "B", "C"][: len(means)]
    return {
        "nodes": 3,
        "factors": [{"key": "alignment", "label": "Alignment", "varied": True,
                     "levels": [{"key": k, "label": k} for k in dict.fromkeys(labels)]}],
        "cells": [{"key": f"c{i}", "label": f"cell {i}", "levels": [lab]}
                  for i, lab in enumerate(labels)],
        "pairs": {"posterior": (own + diff).tolist()},
    }  # fmt: skip


PARTS = {"all": {"total": 30.0, "within": 10.0, "alignment": 20.0, "tree": 0.0,
                 "model": 0.0, "gaps": 0.0, "other": 0.0}}  # fmt: skip


def test_the_map_keeps_the_distances_between_means() -> None:
    means = np.array([[0.0, 0], [0.5, 0], [6, 0], [6, 0.5], [0, 8]])
    m = o.reconstruction_map(summary_with(means, 4.0), "posterior", PARTS, COMPONENTS)
    assert m is not None and m.by == "alignment" and m.shown == 100
    xy = np.array([[p.x, p.y] for p in m.points])
    drawn = np.linalg.norm(xy[:, None] - xy[None], axis=-1)
    true = np.linalg.norm(means[:, None] - means[None], axis=-1)
    ratio = drawn[np.triu_indices(5, 1)] / true[np.triu_indices(5, 1)]
    assert np.allclose(ratio, ratio[0], rtol=1e-6)  # one scale: an exact embedding
    # One halo per cluster, sqrt(own / 2) on the same scale; one label per level.
    assert len(m.halos) == 3 and np.isclose(m.halos[0][2], (4.0 / 2) ** 0.5 * ratio[0])
    assert [label.text for label in m.labels] == ["A", "B", "C"]
    assert [(level, n) for _, level, n in m.legend] == [("A", 2), ("B", 2), ("C", 1)]


def test_groups_that_sit_together_share_a_label() -> None:
    means = np.array([[0.0, 0], [0.1, 0], [0.05, 0.05], [0.1, 0.1], [9, 0]])
    m = o.reconstruction_map(summary_with(means, 0.5), "posterior", PARTS, COMPONENTS)
    assert m is not None
    assert sorted(label.text for label in m.labels) == ["A, B", "C"]


def test_hiding_a_level_lays_the_rest_out_again() -> None:
    means = np.array([[0.0, 0], [0.1, 0], [6, 0], [6, 0.1], [0, 80]])
    full = o.reconstruction_map(summary_with(means, 1.0), "posterior", PARTS, COMPONENTS)
    shown = {"alignment:A", "alignment:B"}
    part = o.reconstruction_map(summary_with(means, 1.0), "posterior", PARTS, COMPONENTS,
                                shown=shown)  # fmt: skip
    assert full is not None and part is not None
    assert part.hidden == 1 and len(part.points) == 4
    assert [t for _, t, ticked in part.choices[0].levels if not ticked] == ["C"]

    def width(m: o.ReconstructionMap) -> float:
        xs = [p.x for p in m.points if p.level != "C"]
        return max(xs) - min(xs)

    assert width(part) > 2 * width(full)  # A and B now fill the frame
    two = o.reconstruction_map(summary_with(means, 1.0), "posterior", PARTS, COMPONENTS,
                               shown={"alignment:C"})  # fmt: skip
    assert two is not None and two.note and not two.points


def test_nothing_placed_when_every_reconstruction_is_the_same() -> None:
    m = o.reconstruction_map(summary_with(np.zeros((4, 2)), 2.0), "posterior", PARTS, COMPONENTS)
    assert m is not None and not m.points and "all the same" in m.note


def test_the_explainer_splits_the_total() -> None:
    e = o.explainer(summary_with(np.zeros((5, 2)), 1.0), PARTS, COMPONENTS)
    assert e.cells == 5 and e.design == "3 alignments" and e.total == "30.0"
    assert sum(w for _, w, _, _ in e.bar) == o.BAR_WIDTH
    assert np.isclose(e.posterior_end, o.BAR_WIDTH / 3)
    assert e.choices_label == "alignment 20.0"


def test_pies_scale_by_area_and_split_like_the_bars() -> None:
    nodes = {
        "N1": {"all": {"total": 8.0}, "residue": {"within": 2.0, "alignment": 2.0},
               "presence": {"alignment": 4.0}},
        "N2": {"all": {"total": 2.0}, "residue": {"within": 2.0}, "presence": {}},
    }  # fmt: skip
    places = {"N1": (10.0, 20.0), "N2": (30.0, 40.0), "Leafless": (0.0, 0.0)}
    drawn = {p.node: p for p in o.pies(places, nodes, {"N1": "?node=N1"}, COMPONENTS, "pt")}
    assert set(drawn) == {"N1", "N2"}
    assert drawn["N1"].r == o.PIE_MAX and np.isclose(drawn["N2"].r, o.PIE_MAX / 2)  # area
    fills = [fill for _, fill, _ in drawn["N1"].slices]
    assert fills == ["#2a78d6", "#eb6834", "url(#pt-hatch-2)"]  # bar order, presence hatched
    assert drawn["N1"].href == "?node=N1"
    assert drawn["N2"].slices[0][0].count(" a ") == 2  # one slice: a whole circle


def test_an_arc_spans_its_fraction() -> None:
    assert o.arc(0, 0, 10, 0, 0.25) == "M 0.00 0.00 L 0.00 -10.00 A 10.00 10.00 0 0 1 10.00 0.00 Z"
    assert " 0 1 1 " in o.arc(0, 0, 10, 0, 0.75)  # the large arc


def test_dots_at_one_place_are_fanned_out() -> None:
    dots = [
        o.MapPoint(100.0, 100.0, "#2a78d6", "A", ""),
        o.MapPoint(100.5, 100.0, "#eb6834", "B", ""),
        o.MapPoint(200.0, 100.0, "#2a78d6", "A", ""),
    ]
    out = o.fanned(dots)
    assert out[2] == dots[2]  # alone: unmoved
    assert np.isclose(np.hypot(out[0].x - out[1].x, out[0].y - out[1].y), 2 * o.FAN_RING)


def test_nearest_tip_depths() -> None:
    depths = o.nearest_tip("((a:0.1,b:0.3)N2:0.2,(c:0.5,d:0.4)N3:0.05)N1;")
    assert np.isclose(depths["N2"], 0.1) and np.isclose(depths["N3"], 0.4)
    assert np.isclose(depths["N1"], 0.3)  # through N2 to a


def test_the_depth_scatter_and_its_trend() -> None:
    rows = [(f"N{i}", i / 10, 2 + i, {"all": {"total": 10.0, "within": float(i)}})
            for i in range(10)]  # fmt: skip
    s = o.depth_scatter(rows, "within", "tip", "share", COMPONENTS, {"N1": "?node=N1"})
    assert s is not None and len(s.dots) == 10 and s.trend.startswith("M ")
    assert s.dots[0][1] > s.dots[-1][1]  # deeper: more posterior, higher on the page
    assert "a median 10% of the total in the third of ancestors nearest" in s.note
    assert o.depth_scatter(rows, "nonsense", "tip", "share", COMPONENTS, {}) is None


def test_sliders_and_stretching() -> None:
    assert o.slider("9", o.TREE_SPACING) == o.TREE_SPACING[1]
    assert o.slider("nan", o.TREE_SPACING) == 1.0 and o.slider(None, o.PIE_SIZES) == 1.0
    from phylobot import treeview

    drawing = treeview.layout("((a:1,b:1)N2:1,c:1)N1;", "true")
    tall = o.stretch(treeview.layout("((a:1,b:1)N2:1,c:1)N1;", "true"), 2.0, 16, 44)
    gap = [lf.y for lf in drawing.leaves]
    wide = [lf.y for lf in tall.leaves]
    assert np.isclose(wide[1] - wide[0], 2 * (gap[1] - gap[0]))


def test_the_takeaways_say_what_it_means() -> None:
    parts = {"all": {"total": 30.0, "within": 10.0, "alignment": 20.0},
             "presence": {"alignment": 15.0}}  # fmt: skip
    said = o.takeaways({"nodes": 3, "columns": 100}, parts, {"N1": 5.0, "N2": 30.0, "N3": 55.0},
                       COMPONENTS, "posterior")  # fmt: skip
    text = " ".join(said)
    assert "on average they differ at 30.0 of its 100 alignment columns" in text
    assert "10.0 of those positions (33%)" in text and "20.0 (67%) come from your choices" in text
    assert "whether a residue is there at all" in text and "N3, at 55.0" in text


def test_the_glance_holds_the_takeaways_numbers() -> None:
    """The page's opening figures: the same numbers the sentences give, as tiles."""
    parts = {
        "all": {"total": 30.0, "within": 10.0, "alignment": 15.0, "tree": 1.0, "model": 3.0,
                "gaps": 1.0, "other": 0.0},
        "presence": {"total": 12.0, "alignment": 11.0},
        "residue": {"total": 18.0},
    }  # fmt: skip
    components = (
        ("within", "Posterior", 1), ("alignment", "Alignment", 2), ("tree", "Tree", 3),
        ("model", "Model", 4), ("gaps", "Gap placement", 5), ("other", "Other", 6),
    )  # fmt: skip
    totals = {"N1": 5.0, "N2": 30.0, "N3": 55.0}
    g = o.glance({"nodes": 3, "columns": 100}, parts, totals, components, "posterior")
    assert (g.total, g.columns, g.nodes, g.share) == ("30.0", 100, 3, "30%")
    tiles = {t.key: t for t in g.tiles}
    assert list(tiles) == ["posterior", "choices", "largest", "presence", "most"]
    assert (tiles["posterior"].value, tiles["posterior"].colour) == ("10.0", 1)
    assert tiles["choices"].value == "20.0" and "67%" in tiles["choices"].unit
    assert tiles["largest"].label == "Mostly the alignment" and tiles["largest"].colour == 2
    assert tiles["largest"].value == "15.0" and "75%" in tiles["largest"].unit
    assert tiles["presence"].value == "12.0"
    assert (tiles["most"].node, tiles["most"].value) == ("N3", "55.0")
    # Most probable sequences only: there is no posterior part to show.
    point = o.glance({"nodes": 3, "columns": 100}, parts, totals, components, "map")
    assert "posterior" not in {t.key for t in point.tiles}
    # Nothing differs: the number alone, and the spread over the ancestors.
    same = {"all": dict.fromkeys(parts["all"], 0.0)}
    nothing = o.glance({}, same, {"N1": 0.0}, components, "posterior")
    assert [t.key for t in nothing.tiles] == ["most"]
