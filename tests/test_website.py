"""phylobot.com (website/): the static page builds from the app's own content and holds
together: every file it names is there, every in-page link has a target, and nothing is
fetched from another site."""

import html
import importlib.util
import json
import re
from pathlib import Path
from types import ModuleType

import pytest

from phylobot import about, skins

WEBSITE = Path(__file__).resolve().parent.parent / "website"
APP_STATIC = Path(__file__).resolve().parent.parent / "src/phylobot/static/phylobot"


@pytest.fixture(scope="module")
def builder() -> ModuleType:
    spec = importlib.util.spec_from_file_location("website_build", WEBSITE / "build.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def site(builder: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out: Path = builder.build(tmp_path_factory.mktemp("website") / "dist")
    return out


def attributes(page: str, tag: str, name: str) -> list[str]:
    return re.findall(rf'<{tag}\b[^>]*\s{name}="([^"]*)"', page)


def test_the_page_is_a_landing_page_built_from_the_apps_own_content(site: Path) -> None:
    page = (site / "index.html").read_text()
    assert page.lstrip().startswith("<!doctype html>")
    assert "{#" not in page and "#}" not in page and "{%" not in page and "{{" not in page
    # The title page, as the app draws it, without the app's bar.
    assert '<h1>PhyloBot <span class="two">2</span></h1>' in page
    assert "Every protein has a past." in page and "phylobot-logo_new.png" in page
    assert "<header" not in page and "My Jobs" not in page and "htmx" not in page
    # What a first visit needs: how to get it, and how to cite it.
    assert "make up" in page and "curl -L" in page and "make storage-init" not in page
    assert "git clone https://github.com/headwatersdiscovery/phylobot2.git" in page
    assert "releases/latest/download/phylobot.tar.gz" in page
    assert f"doi.org/{about.PHYLOBOT_2016.doi}" in page
    for heading in (
        "From Sequences to Ancestors", "A Look Inside", "Why Reconstruct Ancestors",
        "Confidence You Can Measure", "Quick Start", "How to Cite PhyloBot", "Included Methods",
    ):  # fmt: skip
        assert f">{heading}</h2>" in page, heading
    # The About page's long text and tables are all here, folded, from the same partials
    # and the same data: nothing on this page is a second copy.
    assert "Robustness across alignments, trees and models, by design." in page
    assert "Feature Comparison" in page and "FireProtASR" in page
    for _, group in about.METHODS:
        for method in group:
            assert f"<b>{method.name}</b>" in page, method.name
    for tool in about.COMPARISON_TOOLS:
        assert f">{tool}</th>" in page, tool
    # In the order a visitor asks: what is it, what does it look like, how do I get it.
    order = ['class="slogan"', 'id="how"', 'id="tour"', 'id="quick-start"', 'id="cite"']
    assert [page.index(mark) for mark in order] == sorted(page.index(mark) for mark in order)


def test_the_pipeline_demo_names_only_what_the_app_cites(builder: ModuleType, site: Path) -> None:
    page = (site / "index.html").read_text()
    cited = {method.name for _, group in about.METHODS for method in group}
    assert set(builder.TOOLS) <= cited
    for tool in builder.TOOLS:
        assert f"<span>{tool}</span>" in page, tool
    # The five steps are the page's own text, readable without the script that plays them.
    for step in ("Sequences", "Alignment", "Tree", "Ancestors", "Paths"):
        assert f"<b>{step}</b>" in page, step
    # The invented family is consistent: each change the path names is one the sequences
    # on either side of that branch really differ by.
    script = (site / "static/website.js").read_text()
    sequences = dict(re.findall(r'(n\d): \{ label: "[^"]+", x: \d+, seq: "([A-Z]+)"', script))
    tips = re.findall(r'\{ name: "seq \d+", seq: "([A-Z-]+)" \}', script)
    assert len(sequences) == 5 and len(tips) == 6
    path = re.findall(r'\["(n\d)", (?:"(n\d)"|(\d)), "([A-Z0-9 ]+)"\]', script)
    assert len(path) == 3
    for start, node, tip, changes in path:
        before, after = sequences[start], sequences[node] if node else tips[int(tip)]
        differ = {i for i in range(len(before)) if before[i] != after[i]}
        named = {int(change[1:-1]) - 1: change for change in changes.split()}
        assert differ == set(named), (start, node or tip)
        for column, change in named.items():
            assert (before[column], after[column]) == (change[0], change[-1]), change


def test_the_screenshots_come_in_both_themes(builder: ModuleType, site: Path) -> None:
    page = (site / "index.html").read_text()
    shots = builder.read_config()["shots"]
    assert len(shots) >= 3
    for shot in shots:
        for theme in ("light", "dark"):
            assert (site / f"static/shots/{shot['key']}-{theme}.webp").is_file(), shot["key"]
            assert f'src="static/shots/{shot["key"]}-{theme}.webp"' in page
        assert f"<h3>{shot['title']}</h3>" in page
        assert f"<figcaption>{html.escape(shot['caption'])}</figcaption>" in page
    # A screenshot site.toml names but the folder lacks stops the build.
    with pytest.raises(builder.ConfigError):
        builder.read_shots({"shots": [{"key": "missing", "title": "x", "text": "x", "alt": "x"}]})
    with pytest.raises(builder.ConfigError):
        builder.read_shots({"shots": [{"key": "../x"}]})


def test_every_skin_in_the_menu_comes_along(site: Path) -> None:
    page = (site / "index.html").read_text()
    assert f'data-skin="{skins.DEFAULT}"' in page and 'class="backdrop"' in page
    for skin in skins.MENU:
        assert f'<option value="{skin.key}" data-backdrop="{skin.backdrop}"' in page, skin.key
        assert f'"{skin.key}": "{skin.backdrop}"' in page, skin.key  # the script in <head>
        if skin.stylesheet:
            assert (site / "static" / skin.stylesheet).is_file(), skin.key
    # In the order the cascade depends on, after the base styles and before the page's own.
    sheets = attributes(page, "link", "href")
    sheets = [href for href in sheets if href.endswith(".css")]
    assert sheets[0] == "static/phylobot/site.css" and sheets[-1] == "static/website.css"
    expected = [s.stylesheet for s in skins.SKINS if s in skins.MENU and s.stylesheet]
    assert sheets[1:-1] == [f"static/{sheet}" for sheet in expected]
    # Fonts are found by the stylesheets' own relative paths.
    for sheet in (site / "static/phylobot/skins").glob("*.css"):
        for font in re.findall(r'url\("([^")]+\.woff2)"\)', sheet.read_text()):
            assert (sheet.parent / font).resolve().is_file(), (sheet.name, font)


def test_the_site_is_self_contained(site: Path) -> None:
    page = (site / "index.html").read_text()
    # Every file the page names is in the folder, by a relative path (any host, any folder).
    local = attributes(page, "script", "src") + attributes(page, "img", "src")
    local += attributes(page, "link", "href")
    assert local
    for path in local:
        if path.startswith("https://phylobot.com/"):
            continue  # the canonical link
        assert not re.match(r"[a-z]+:|/", path), path
        assert (site / path).is_file(), path
    # Nothing loads from another site: its only outside addresses are links to follow.
    for tag, name in (("script", "src"), ("img", "src"), ("iframe", "src")):
        assert not [u for u in attributes(page, tag, name) if "//" in u], tag
    stylesheets = re.findall(r'<link rel="stylesheet" href="([^"]*)"', page)
    assert stylesheets and not [u for u in stylesheets if "//" in u]
    # Every link within the page has somewhere to land.
    ids = set(re.findall(r'\sid="([^"]+)"', page))
    anchors = [href[1:] for href in attributes(page, "a", "href") if href.startswith("#")]
    assert anchors and set(anchors) <= ids, set(anchors) - ids
    assert not [href for href in attributes(page, "a", "href") if href.startswith("/")]
    # The backdrop's structures came too, inline: the list and each image as a data:
    # address, so the page rains even when it is opened straight from disk.
    script = (site / "static/structures.js").read_text()
    assert script.startswith("window.phylobotStructures = [")
    listed = json.loads((APP_STATIC / "backdrop/structures.json").read_text())
    assert script.count('"src": "data:image/jpeg;base64,') == len(listed) >= 8
    assert "window.phylobotStructures" in (site / "static/phylobot/backdrop.js").read_text()


def test_the_tutorial_video_has_a_place_once_it_exists(builder: ModuleType, tmp_path: Path) -> None:
    # site.toml as committed has no video, so the page says nothing of one.
    assert builder.read_config()["tutorial"] is None
    without = builder.render(builder.read_config())
    assert 'id="tutorial"' not in without and "youtube" not in without.lower()
    # With an ID: a section under the title, a button beside Quick Start, a contents entry.
    config = builder.read_config(youtube_id="abcDEF_-123")
    page = builder.render(config)
    assert (
        page.index('class="slogan"') < page.index('id="tutorial"') < page.index('id="quick-start"')
    )
    assert page.count('href="#tutorial"') == 2  # the title's button and the section links
    assert 'data-video="abcDEF_-123"' in page
    # Until the reader presses play it is a link, not a player: nothing comes from YouTube.
    assert 'href="https://www.youtube.com/watch?v=abcDEF_-123"' in page
    assert "<iframe" not in page and "youtube-nocookie.com" not in page
    assert "youtube-nocookie.com/embed/" in (WEBSITE / "static/website.js").read_text()
    # Only a video ID reaches the page.
    for bad in ("https://youtu.be/abcDEF_-123", "abc", '"><script>x</script>', "abcDEF_-1234"):
        with pytest.raises(builder.ConfigError):
            builder.read_config(youtube_id=bad)
    # A poster must be an image that is really in website/static/.
    toml = tmp_path / "site.toml"
    toml.write_text(
        (WEBSITE / "site.toml").read_text().replace('poster = ""', 'poster = "../build.py"')
    )
    with pytest.raises(builder.ConfigError):
        builder.read_config(toml, youtube_id="abcDEF_-123")
