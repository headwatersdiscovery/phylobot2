"""Build phylobot.com: one static page, written to website/dist/ (``make website``).

The page is the app's title page and its About page in one, for people who have not
installed PhyloBot yet. Nothing here is a second copy of the app: the About text comes
from the app's own partials (``src/phylobot/templates/phylobot/about/``), the citations
and the comparison table from ``phylobot.about``, the skins from ``phylobot.skins``, and
the styles, fonts, logo and backdrop from the app's static folder, copied as they are. So
a rebuild picks up whatever changed there; ``tests/test_website.py`` builds the site and
checks that it holds together.

The output is plain files with relative links: any static host serves it, from the root
of a domain or from a folder. It makes no request to another site until a reader presses
play on the tutorial video.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import tomllib
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT / "src") not in sys.path:  # run straight from a checkout
    sys.path.insert(0, str(ROOT / "src"))

import django  # noqa: E402
from django.conf import settings  # noqa: E402

APP_TEMPLATES = ROOT / "src/phylobot/templates"
APP_STATIC = ROOT / "src/phylobot/static/phylobot"
DIST = HERE / "dist"

# What the page uses of the app's static folder, copied to dist/static/phylobot/ under the
# same names (the skins' stylesheets find the fonts by relative path).
APP_FILES = (
    "site.css",
    "backdrop.js",
    "copy.js",
    "matrix-rain.js",
    "phylobot-logo_new.png",
)
APP_FOLDERS = ("skins", "fonts")
OWN_FILES = ("website.css", "website.js")
SHOTS = HERE / "static" / "shots"

# The tools named under the pipeline, each of which must be one phylobot.about cites: the
# page may not claim a method the app does not list.
TOOLS = (
    "MAFFT", "MUSCLE", "Clustal Omega", "PRANK", "IQ-TREE", "RAxML-NG", "ModelFinder",
    "GeneRax", "GRASP", "Nextflow",
)  # fmt: skip

YOUTUBE_ID = re.compile(r"[A-Za-z0-9_-]{11}")
POSTER = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.(?:png|jpe?g|webp)")


class ConfigError(ValueError):
    """site.toml says something the page cannot be built from."""


def engine() -> Any:
    """Django's template engine alone, over this folder's templates and then the app's:
    no database, no apps, and none of the app's settings (an engine of its own, so the
    tests can build the site from inside the configured app)."""
    if not settings.configured:
        settings.configure()
        django.setup()
    from django.template import Engine

    return Engine(dirs=[str(HERE / "templates"), str(APP_TEMPLATES)])


def read_config(path: Path = HERE / "site.toml", youtube_id: str | None = None) -> dict[str, Any]:
    """site.toml, checked. ``youtube_id`` overrides the file's (to preview the layout)."""
    config = tomllib.loads(path.read_text())
    for key in ("url", "repository", "contact"):
        if not isinstance(config.get(key), str) or not config[key]:
            raise ConfigError(f"site.toml: {key} is missing")
    if not config["url"].startswith("https://") or not config["url"].endswith("/"):
        raise ConfigError("site.toml: url must start with https:// and end with /")
    tutorial = dict(config.get("tutorial") or {})
    video = tutorial.get("youtube_id", "") if youtube_id is None else youtube_id
    if video and not YOUTUBE_ID.fullmatch(video):
        raise ConfigError(
            "the YouTube ID is the eleven characters after 'v=' in the video's address, "
            f"not {video!r}"
        )
    poster = tutorial.get("poster", "")
    if poster and not (POSTER.fullmatch(poster) and (HERE / "static" / poster).is_file()):
        raise ConfigError(f"site.toml: poster {poster!r} is not an image in website/static/")
    config["shots"] = read_shots(config)
    config["tutorial"] = (
        {
            "youtube_id": video,
            "title": tutorial.get("title") or "Watch the Tutorial",
            "blurb": tutorial.get("blurb", ""),
            "poster": poster,
        }
        if video
        else None
    )
    return config


def read_shots(config: dict[str, Any]) -> list[dict[str, Any]]:
    """The screenshots of "A Look Inside" (site.toml's [[shots]]): each needs its light and
    dark picture in static/shots/. A page with none simply has no such section."""
    shots = []
    for shot in config.get("shots") or []:
        key = str(shot.get("key", ""))
        if not re.fullmatch(r"[a-z][a-z0-9-]*", key):
            raise ConfigError(f"site.toml: a screenshot's key must be a short name, not {key!r}")
        for theme in ("light", "dark"):
            if not (SHOTS / f"{key}-{theme}.webp").is_file():
                raise ConfigError(f"site.toml: static/shots/{key}-{theme}.webp is missing")
        for field in ("title", "text", "caption", "alt", "width", "height"):
            if not shot.get(field):
                raise ConfigError(f"site.toml: screenshot {key!r} has no {field}")
        shots.append(dict(shot))
    return shots


def structures_script() -> str:
    """The backdrop's structures as a script: the list, each image inline as a data:
    address (backdrop.js reads window.phylobotStructures). A page opened straight from disk
    can neither fetch the list nor read the pixels of an image file beside it, and the
    backdrop needs both; inline, it rains there too."""
    import base64

    folder = APP_STATIC / "backdrop"
    entries = []
    for entry in json.loads((folder / "structures.json").read_text()):
        data = base64.b64encode((folder / entry["file"]).read_bytes()).decode("ascii")
        entries.append(
            {"id": entry["id"], "title": entry["title"], "src": f"data:image/jpeg;base64,{data}"}
        )
    return "window.phylobotStructures = " + json.dumps(entries) + ";\n"


def context(config: dict[str, Any]) -> dict[str, Any]:
    from phylobot import about, skins

    cited = {method.name for _, group in about.METHODS for method in group}
    if missing := [tool for tool in TOOLS if tool not in cited]:
        raise ConfigError(f"phylobot.about cites no {', '.join(missing)}")

    menu = skins.MENU
    offered = {skin.key for skin in menu}
    default = config.get("default_skin") or skins.DEFAULT
    if default not in offered:
        raise ConfigError(
            f"site.toml: default_skin {default!r} is not in the skin menu "
            f"({', '.join(skin.key for skin in menu)})"
        )
    return {
        "site": config,
        "tutorial": config["tutorial"],
        "shots": config["shots"],
        "tools": TOOLS,
        "other_tools": len(about.COMPARISON_TOOLS) - 1,
        # The menu's skins only, in the order their stylesheets must load.
        "stylesheets": [
            skin.stylesheet for skin in skins.SKINS if skin.key in offered and skin.stylesheet
        ],
        "skin_menu": menu,
        "default_skin": default,
        "default_backdrop": next(skin.backdrop for skin in menu if skin.key == default),
        # For the script in <head>: each skin and whether it shows the backdrop. Our own
        # constants, so safe to write into the page as JSON.
        "skins_json": json.dumps({skin.key: skin.backdrop for skin in menu}),
        "reference": about.PHYLOBOT_2016,
        "methods": about.METHODS,
        "matrices": about.MATRICES,
        "variants": about.VARIANTS,
        "raxml_missing": about.RAXML_MISSING,
        "comparison_tools": about.COMPARISON_TOOLS,
        "comparison": about.COMPARISON,
        "comparison_legend": about.COMPARISON_LEGEND,
    }


def render(config: dict[str, Any]) -> str:
    from django.template import Context

    template = engine().get_template("index.html")
    return str(template.render(Context(context(config))))


def build(out: Path = DIST, youtube_id: str | None = None) -> Path:
    """Write the whole site to ``out`` (emptied first) and return it."""
    config = read_config(youtube_id=youtube_id)
    page = render(config)
    if out.exists():
        shutil.rmtree(out)
    app = out / "static" / "phylobot"
    app.mkdir(parents=True)
    for name in APP_FILES:
        shutil.copy2(APP_STATIC / name, app / name)
    for name in APP_FOLDERS:
        shutil.copytree(APP_STATIC / name, app / name)
    for name in OWN_FILES:
        shutil.copy2(HERE / "static" / name, out / "static" / name)
    (out / "static" / "structures.js").write_text(structures_script())
    if config["shots"]:
        (out / "static" / "shots").mkdir()
        for shot in config["shots"]:
            for theme in ("light", "dark"):
                name = f"{shot['key']}-{theme}.webp"
                shutil.copy2(SHOTS / name, out / "static" / "shots" / name)
    if config["tutorial"] and config["tutorial"]["poster"]:
        poster = config["tutorial"]["poster"]
        shutil.copy2(HERE / "static" / poster, out / "static" / poster)
    (out / "index.html").write_text(page)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument("--out", type=Path, default=DIST, help="where to write the site")
    parser.add_argument(
        "--youtube-id",
        default=None,
        help="build with this tutorial video, whatever site.toml says (to preview the layout)",
    )
    args = parser.parse_args()
    try:
        out = build(args.out, args.youtube_id)
    except ConfigError as error:
        raise SystemExit(str(error)) from error
    size = sum(path.stat().st_size for path in out.rglob("*") if path.is_file())
    print(f"Built {out} ({size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
