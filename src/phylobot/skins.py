"""Skins: the site's look, chosen per browser from the menus in the top bar.

Each skin is a stylesheet under ``static/phylobot/skins/`` whose rules are scoped to
``:root[data-skin="<key>"]``, so every skin's stylesheet is loaded on every page and
switching is instant (``skins.js`` sets the attribute and the cookie). "classic" (labelled
Skinless) has no stylesheet of its own: it is the base ``site.css`` alone, the pilot look.

The landing and About pages' backdrop (``backdrop.js``) is chosen on the Settings page,
under Display, and remembered in a cookie. Only Clean, Dark, Spectrum, Lagoon and Abyss
show one at all: on the others the choice is ignored and nothing is drawn, which the skin
menu's options carry to the browser so switching skins turns it off and on without a
reload.
"""

from typing import NamedTuple

from django.http import HttpRequest

COOKIE = "phylobot_skin"
BACKDROP_COOKIE = "phylobot_backdrop"
SPEED_COOKIE = "phylobot_backdrop_speed"
DENSITY_COOKIE = "phylobot_backdrop_density"


class Skin(NamedTuple):
    key: str
    label: str
    stylesheet: str | None  # path under static/, or None for the base styles only
    # The backdrop this skin shows. "none" means it never shows one, whatever the reader
    # picked: the quiet skins stay quiet.
    backdrop: str = "none"


class Backdrop(NamedTuple):
    key: str
    label: str


class Level(NamedTuple):
    """One step of a backdrop dial. ``value`` is the factor itself, so the cookie, the
    attribute ``backdrop.js`` reads and the menu's selected option are all one string."""

    value: str
    label: str


# In the order their stylesheets load, which the cascade depends on: Clean's layout
# before the skins built on it, Dark's before Matrix, Dim and Abyss. The menu has its
# own order (MENU).
SKINS = (
    Skin("clean", "Clean", "phylobot/skins/clean.css", "structures"),
    Skin("dark", "Dark", "phylobot/skins/dark.css", "structures"),  # Clean's layout; loads after it
    Skin("legible", "High legibility", "phylobot/skins/legible.css"),  # likewise; never a backdrop
    Skin("matrix", "Matrix", "phylobot/skins/matrix.css"),  # Dark's rules too; loads after
    Skin("dim", "Dark II", "phylobot/skins/dim.css"),  # likewise; quiet, so no backdrop
    Skin("abyss", "Abyss", "phylobot/skins/abyss.css", "structures"),  # likewise
    Skin("spectrum", "Spectrum", "phylobot/skins/spectrum.css", "structures"),  # Clean's layout
    Skin("lagoon", "Lagoon", "phylobot/skins/lagoon.css", "structures"),  # Clean's layout
    Skin("broadsheet", "Broadsheet", "phylobot/skins/broadsheet.css"),  # likewise; print is still
    Skin("retro95", "PhyloBot '95", "phylobot/skins/retro95.css"),  # Clean's layout
    Skin("sciencefair", "Science Fair", "phylobot/skins/sciencefair.css"),  # Clean's layout
    Skin("classic", "Skinless", None),  # the key predates the label: readers' cookies keep it
)
# The top bar's menu: the everyday skins first, then the playful ones, then none at all.
# (Science Fair is defined above but kept out of the menu.)
MENU_ORDER = (
    "clean",
    "dark",
    "dim",  # labelled "Dark II"
    "spectrum",
    "legible",
    "abyss",
    "broadsheet",
    "lagoon",
    "matrix",
    "retro95",
    "classic",
)
MENU = tuple(next(s for s in SKINS if s.key == key) for key in MENU_ORDER)
BACKDROPS = (
    Backdrop("structures", "PDB structures"),
    Backdrop("none", "None"),
)
# How fast the shapes fall and how many there are, as multipliers of the skin's own
# settings. The Display panel on the Settings page sets them, per browser.
SPEEDS = (
    Level("0.4", "Very slow"),
    Level("0.7", "Slow"),
    Level("1", "Normal"),
    Level("1.5", "Brisk"),
    Level("2.2", "Fast"),
)
DENSITIES = (
    Level("0.4", "Sparse"),
    Level("0.7", "Light"),
    Level("1", "Normal"),
    Level("1.5", "Busy"),
    Level("2.2", "Dense"),
)
NORMAL = "1"
DEFAULT = SKINS[0].key
_KEYS = {skin.key for skin in SKINS}
_BACKDROP_KEYS = {backdrop.key for backdrop in BACKDROPS}


# Top-bar section for each page (by URL name), so the bar can mark where the reader is.
_SECTIONS = {
    "submit": "new",
    "my_jobs": "jobs",
    "job_detail": "jobs",
    "tree_viewer": "jobs",
    "ancestor_viewer": "jobs",
    "my_libraries": "libraries",
    "library_page": "libraries",
    "my_alignments": "alignments",
    "alignment_editor": "alignments",
    "homologs": "homologs",
    "homolog_search": "homologs",
    "homolog_example": "homologs",
    "simulations": "simulations",
    "esr": "esr",
    "my_usage": "usage",
    "site_settings": "settings",
    "work_estimate": "settings",
    "work_records": "settings",
    "about": "about",
}

# The pages worth returning to, by URL name: what the home page's "Jump Back In" calls each
# kind (recent.js remembers them in the browser; phylobot.home_views). A new viewer belongs
# here, or that list never offers it.
PAGES = {
    "job_detail": "Job",
    "tree_viewer": "Tree",
    "ancestor_viewer": "Ancestors",
    "path_viewer": "Mutational Path",
    "structure_viewer": "Structure",
    "alignment_viewer": "Alignment",
    "job_properties": "Properties",
    "uncertainty": "Uncertainty",
    "library_page": "Library",
    "homolog_search": "Homolog Search",
    "alignment_editor": "Alignment Editor",
}


def current(request: HttpRequest) -> str:
    """The reader's skin: the cookie's value if it names a skin, else the default."""
    key = request.COOKIES.get(COOKIE, "")
    return key if key in _KEYS else DEFAULT


def backdrop_preference(request: HttpRequest) -> str:
    """What the reader picked under Settings -> Display, for the skins that show one."""
    key = request.COOKIES.get(BACKDROP_COOKIE, "")
    return key if key in _BACKDROP_KEYS else BACKDROPS[0].key


def current_backdrop(request: HttpRequest, skin: str) -> str:
    """The backdrop to draw: the reader's choice, or nothing on a skin that shows none."""
    if next(s.backdrop for s in SKINS if s.key == skin) == "none":
        return "none"
    return backdrop_preference(request)


def _level(request: HttpRequest, cookie: str, levels: tuple[Level, ...]) -> str:
    """One backdrop dial's factor: the cookie if it names a step, else Normal."""
    value = request.COOKIES.get(cookie, "")
    return value if any(level.value == value for level in levels) else NORMAL


def context(request: HttpRequest) -> dict[str, object]:
    """Template context processor: the skin and the top bar's current section."""
    match = getattr(request, "resolver_match", None)
    skin = current(request)
    url_name = match.url_name if match else ""
    return {
        "skin": skin,
        "skins": SKINS,
        "skin_menu": MENU,
        "backdrop": current_backdrop(request, skin),
        "backdrop_preference": backdrop_preference(request),
        "backdrops": BACKDROPS,
        "backdrop_speed": _level(request, SPEED_COOKIE, SPEEDS),
        "backdrop_speeds": SPEEDS,
        "backdrop_density": _level(request, DENSITY_COOKIE, DENSITIES),
        "backdrop_densities": DENSITIES,
        "nav_section": _SECTIONS.get(url_name, ""),
        # What "Jump Back In" on the home page calls this page ("" if it is not listed).
        "page_kind": PAGES.get(url_name or "", ""),
    }
