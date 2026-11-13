"""The structure page's colour schemes: a colour for every residue of the structure, and
the legend that explains them. Computed on the server, so the picture, the legend and the
sequence strip under it can never disagree; structure.js only paints.

Categorical colours are the Okabe-Ito set (safe for colour-blind readers). The posterior
probability scheme uses the ancestor viewer's palettes, so a residue has the colour its
tile has there.
"""

from collections.abc import Sequence
from dataclasses import dataclass

from phylobot import ancestors
from phylobot.properties import KYTE_DOOLITTLE
from phylobot.structures.mapping import GAP, Site

UNMAPPED = "#f1f2f4"  # the structure's residue has no column in this alignment
ABSENT = "#5d6470"  # the ancestor has no residue at the column (a gap)
SAME = "#b4bac4"
SUBSTITUTED = "#d55e00"
LOST = "#0072b2"
GAINED = "#009e73"
INSERTION = "#cc79a7"  # the marker for residues the structure has no place for

CHEMISTRY = {
    "Acidic D E": ("DE", "#d55e00"),
    "Basic K R H": ("KRH", "#0072b2"),
    "Aliphatic A V L I M": ("AVLIM", "#e69f00"),
    "Aromatic F W Y": ("FWY", "#cc79a7"),
    "Polar S T N Q": ("STNQ", "#009e73"),
    "Special G P C": ("GPC", "#56b4e9"),
}
_CHEMISTRY = {aa: colour for letters, colour in CHEMISTRY.values() for aa in letters}
HYDROPHILIC, NEUTRAL, HYDROPHOBIC = (59, 111, 182), (244, 244, 244), (217, 98, 43)
# AlphaFold's own bands for pLDDT, the colours its database uses.
PLDDT = [
    (90.0, "Very high (> 90)", "#0053d6"),
    (70.0, "Confident (70 to 90)", "#65cbf3"),
    (50.0, "Low (50 to 70)", "#ffdb13"),
    (0.0, "Very low (< 50)", "#ff7d45"),
]
COLD, WARM = (49, 104, 176), (214, 69, 52)

SCHEMES = {
    "pp": "Posterior probability of the ancestral state",
    "differs": "Where the ancestor differs from this structure",
    "changes": "Changes between the two sequences",
    "steps": "When along the path each site changed",
    "chemistry": "Ancestral residue, by chemistry",
    "hydropathy": "Ancestral residue, by hydropathy",
    "confidence": "The structure's own confidence",
}
# What each scheme shows, for the "?" beside the page's Colour menu.
ABOUT = {
    "pp": "How sure the reconstruction is of each residue: the posterior probability of the "
    "ancestor's most probable state at the site, in the ancestor viewer's colours. Red sites "
    "are the ones to doubt.",
    "differs": "The ancestor against the protein the structure is of: where it has another "
    "residue, and where it has none. What you would change in that protein to make the "
    "ancestor.",
    "changes": "The sites that differ between the page's two sequences (the two ends of a "
    "branch or of a mutational path): substitutions, residues lost and residues gained.",
    "steps": "For a path of several branches: each changed site in the colour of the branch "
    "it last changed on, early branches dark, late ones bright.",
    "chemistry": "The ancestor's residue by chemical class: acidic, basic, aliphatic, "
    "aromatic, polar, and glycine, proline and cysteine.",
    "hydropathy": "The ancestor's residue on the Kyte-Doolittle scale, from hydrophilic "
    "(blue) to hydrophobic (orange).",
    "confidence": "Nothing of the ancestor: the structure's own B-factors, or for a "
    "predicted model its pLDDT, in AlphaFold's bands. Where a model is unsure, so is "
    "anything painted on it.",
}
DEFAULT_SCHEME = "pp"


@dataclass(frozen=True)
class Painted:
    colours: list[str]  # one per site
    legend: list[tuple[str, str]]  # (label, colour)
    note: str = ""


def hexed(rgb: tuple[int, int, int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def blend(a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> tuple[int, int, int]:
    t = min(max(t, 0.0), 1.0)
    return (
        round(a[0] + (b[0] - a[0]) * t),
        round(a[1] + (b[1] - a[1]) * t),
        round(a[2] + (b[2] - a[2]) * t),
    )


def hydropathy(aa: str) -> str:
    value = KYTE_DOOLITTLE.get(aa)
    if value is None:
        return UNMAPPED
    if value < 0:
        return hexed(blend(NEUTRAL, HYDROPHILIC, -value / 4.5))
    return hexed(blend(NEUTRAL, HYDROPHOBIC, value / 4.5))


@dataclass(frozen=True)
class Steps:
    """A path's branches, oldest first, and for each column that changes along it the
    branch on which it last did."""

    branches: list[str]  # "Node10 → Node11"
    of_column: dict[int, int]


# Viridis without its darkest end, which is lost on the dark skins' stage.
_RAMP = [(62, 74, 137), (49, 104, 142), (33, 145, 140), (53, 183, 121), (144, 215, 67),
         (253, 231, 37)]  # fmt: skip
MAX_STEP_LEGEND = 8


def step_colour(index: int, count: int) -> str:
    """The colour of a path's branch `index` of `count`: early dark blue, late yellow."""
    t = index / (count - 1) if count > 1 else 1.0
    at = t * (len(_RAMP) - 1)
    low = min(int(at), len(_RAMP) - 2)
    return hexed(blend(_RAMP[low], _RAMP[low + 1], at - low))


def available(compared: bool, has_pp: bool, has_steps: bool = False) -> list[str]:
    """The schemes this view can draw."""
    keys = []
    for key in SCHEMES:
        if (key == "pp" and not has_pp) or (key == "changes" and not compared):
            continue
        if key == "steps" and not has_steps:
            continue
        keys.append(key)
    return keys


def _of_state(sites: Sequence[Site], colour_of: dict[str, str] | None) -> list[str]:
    out = []
    for site in sites:
        if not site.mapped:
            out.append(UNMAPPED)
        elif site.state == GAP:
            out.append(ABSENT)
        elif colour_of is None:
            out.append(hydropathy(site.state))
        else:
            out.append(colour_of.get(site.state, UNMAPPED))
    return out


def paint(
    scheme: str,
    sites: Sequence[Site],
    palette: str,
    predicted: bool,
    reverse: bool = False,
    steps: Steps | None = None,
) -> Painted:
    """Every site's colour in `scheme`. `predicted`: the structure is a predicted model,
    whose B-factor column is pLDDT. `reverse`: the sites hold the page's *second* sequence
    as their state and the first as `other` (a panel that paints the second one), so a
    change is still named in the first-to-second direction. `steps`: the path between the
    two, for the "steps" scheme."""
    missing = [("The ancestor has no residue here", ABSENT), ("Not in the alignment", UNMAPPED)]
    if scheme == "pp":
        colours = []
        for site in sites:
            if not site.mapped:
                colours.append(UNMAPPED)
            elif site.state == GAP:
                colours.append(ABSENT)
            elif site.pp is None:  # an extant sequence: observed, not reconstructed
                colours.append(SAME)
            else:
                colours.append(hexed(ancestors.rgb(site.pp, palette)))
        ramp = [(f"{pp:g}", hexed(ancestors.rgb(pp, palette))) for pp in (0.2, 0.4, 0.6, 0.8, 1.0)]
        if any(s.mapped and s.state != GAP and s.pp is None for s in sites):
            ramp.append(("Observed, not reconstructed", SAME))
        return Painted(colours, [*ramp, *missing], "Posterior probability")
    if scheme == "chemistry":
        legend = [(label, colour) for label, (_, colour) in CHEMISTRY.items()]
        return Painted(_of_state(sites, _CHEMISTRY), [*legend, *missing])
    if scheme == "hydropathy":
        legend = [
            ("Hydrophilic (−4.5)", hexed(HYDROPHILIC)),
            ("0", hexed(NEUTRAL)),
            ("Hydrophobic (+4.5)", hexed(HYDROPHOBIC)),
        ]
        return Painted(_of_state(sites, None), [*legend, *missing])
    if scheme == "differs":
        by = {"": SAME, "sub": SUBSTITUTED, "absent": ABSENT}
        colours = [by[s.differs] if s.mapped else UNMAPPED for s in sites]
        legend = [
            ("Same residue as the structure", SAME),
            ("Another residue in the ancestor", SUBSTITUTED),
            *missing,
        ]
        return Painted(colours, legend)
    if scheme == "changes":
        by = {"": SAME, "sub": SUBSTITUTED, "lost": LOST, "gained": GAINED}
        if reverse:
            by |= {"lost": GAINED, "gained": LOST}
        colours = [by[s.change] if s.mapped else UNMAPPED for s in sites]
        legend = [
            ("Unchanged", SAME),
            ("Substitution", SUBSTITUTED),
            ("Residue lost (deletion)", LOST),
            ("Residue gained (insertion)", GAINED),
            ("Not in the alignment", UNMAPPED),
        ]
        return Painted(colours, legend)
    if scheme == "steps" and steps is not None:
        count = len(steps.branches)
        colours = []
        for site in sites:
            step = steps.of_column.get(site.column) if site.column is not None else None
            if not site.mapped:
                colours.append(UNMAPPED)
            else:
                colours.append(SAME if step is None else step_colour(step, count))
        shown = range(count)
        if count > MAX_STEP_LEGEND:  # a long path: its first, middle and last branches
            shown = sorted({0, count // 2, count - 1})  # type: ignore[assignment]
        legend = [(steps.branches[i], step_colour(i, count)) for i in shown]
        note = "The branch each site last changed on"
        return Painted(
            colours, [*legend, ("Unchanged", SAME), ("Not in the alignment", UNMAPPED)], note
        )
    # confidence: what the file's B-factor column holds
    values = [site.residue.b for site in sites]
    if predicted:
        colours = [next(c for least, _, c in PLDDT if b >= least or least == 0.0) for b in values]
        return Painted(colours, [(label, c) for _, label, c in PLDDT], "pLDDT, the model's own")
    low, high = (min(values), max(values)) if values else (0.0, 0.0)
    span = (high - low) or 1.0
    colours = [hexed(blend(COLD, WARM, (b - low) / span)) for b in values]
    legend = [(f"B-factor {low:.0f} (ordered)", hexed(COLD)), (f"{high:.0f} (mobile)", hexed(WARM))]
    return Painted(colours, legend, "B-factors of the C-alpha atoms")


def text_colour(colour: str) -> str:
    rgb = (int(colour[1:3], 16), int(colour[3:5], 16), int(colour[5:7], 16))
    return ancestors.text_colour(rgb)
