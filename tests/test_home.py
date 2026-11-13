import json
import pathlib
from datetime import UTC, datetime, timedelta

import pytest
from django.contrib.auth.models import User
from django.test import Client
from django.utils import timezone

import phylobot
from phylobot import home_views, skins, urls
from phylobot.models import HomologSearch, Job


def local_user(client: Client) -> User:
    client.get("/")  # local mode makes its one user on the first request
    return User.objects.get(username="local")


@pytest.mark.django_db
def test_home_is_for_someone_who_has_phylobot_running(client: Client) -> None:
    response = client.get("/")
    assert response.status_code == 200
    page = response.content.decode()
    assert '<h1>PhyloBot <span class="two">2</span></h1>' in page
    assert "phylobot/phylobot-logo_new.png" in page and "Every protein has a past." in page
    # Ways in, not installation: that is phylobot.com's subject (website/).
    assert '<a class="button" href="/jobs/new/">Start a New Job</a>' in page
    for url in ("/homologs/", "/alignments/", "/simulations/", "/about/#reference"):
        assert f'href="{url}"' in page, url
    for gone in ("git clone", "make up", "curl -L", "Quick Start"):
        assert gone not in page, gone
    assert "htmx-2.0.11.min.js" in page
    # Nothing to list yet: the tagline stays, and the list says where the examples are.
    assert 'class="intro"' in page and "No jobs yet." in page
    assert "hx-trigger" not in page.split('id="home-work"')[1].split(">")[0]


@pytest.mark.django_db
def test_home_lists_the_latest_jobs_and_searches(client: Client) -> None:
    owner = local_user(client)
    now = timezone.now()
    done = Job.objects.create(
        owner=owner, name="Done <b>", status="succeeded", finished_at=now - timedelta(hours=3)
    )
    running = Job.objects.create(owner=owner, name="Running", status="running", started_at=now)
    HomologSearch.objects.create(
        owner=owner, name="A search", status="succeeded", service="ebi", database="uniref90",
        query_length=300, finished_at=now - timedelta(days=1, hours=1),
    )  # fmt: skip
    gone = Job.objects.create(owner=owner, name="Deleted", status="succeeded", deleted_at=now)
    theirs = Job.objects.create(owner=User.objects.create_user("other"), name="Not mine")
    page = client.get("/").content.decode()
    panel = page.split('id="home-work"')[1].split("</section>")[0]
    # Unfinished work first, then by when each last changed; names are escaped.
    assert panel.index("Running") < panel.index("Done &lt;b&gt;") < panel.index("A search")
    assert f'href="/jobs/{running.id}/"' in panel and f'href="/jobs/{done.id}/"' in panel
    assert "finished 3 h ago" in panel and "finished yesterday" in panel and "ETA:" in panel
    assert "Homolog search" in panel and 'href="/jobs/">All My Jobs</a>' in panel
    assert str(gone.id) not in page and str(theirs.id) not in page and "Not mine" not in page
    assert 'class="intro"' not in page  # a returning reader needs no pitch
    # While something is unfinished the list asks for itself again; then it stops.
    assert 'hx-get="/home/work/" hx-trigger="every 10s"' in panel
    again = client.get("/home/work/").content.decode()
    assert again.lstrip().startswith("<section") and "Running" in again
    running.status = "succeeded"
    running.finished_at = now
    running.save()
    assert "hx-trigger" not in client.get("/home/work/").content.decode()
    assert len(home_views.recent_work(owner, limit=2)) == 2


def test_how_long_ago_is_short() -> None:
    now = datetime(2026, 11, 11, 12, 0, tzinfo=UTC)
    for elapsed, text in (
        (timedelta(seconds=20), "just now"),
        (timedelta(minutes=5), "5 min ago"),
        (timedelta(hours=3, minutes=59), "3 h ago"),
        (timedelta(days=1, hours=2), "yesterday"),
        (timedelta(days=6), "6 days ago"),
        (timedelta(days=40), "2026-10-02"),
        (timedelta(seconds=-30), "just now"),  # a clock a little ahead
    ):
        assert home_views.ago(now - elapsed, now) == text, elapsed


@pytest.mark.django_db
def test_jump_back_in_is_remembered_by_the_browser_and_checked_by_the_server(
    client: Client,
) -> None:
    owner = local_user(client)
    page = client.get("/").content.decode()
    # Every page loads the script; the home page holds the (hidden, empty) list it fills.
    assert "phylobot/recent.js" in page and f'data-user="{owner.pk}"' in page
    assert 'data-recent-pages data-alive-url="/home/alive/" hidden' in page
    assert "data-page-kind" not in page  # the home page is not somewhere to jump back to
    assert "data-page-kind" not in client.get("/jobs/").content.decode()
    job = Job.objects.create(owner=owner, name="Mine")
    assert 'data-page-kind="Job"' in client.get(f"/jobs/{job.id}/").content.decode()
    # A page the list names must be a real URL name, or a rename would drop it silently.
    names = {getattr(pattern, "name", None) for pattern in urls.urlpatterns}
    assert set(skins.PAGES) <= names
    # Which remembered records still exist: the reader's own, undeleted ones only.
    search = HomologSearch.objects.create(
        owner=owner, service="ebi", database="uniref90", query_length=10
    )
    deleted = Job.objects.create(owner=owner, deleted_at=timezone.now())
    theirs = Job.objects.create(owner=User.objects.create_user("other"))
    asked = [job.id, search.id, deleted.id, theirs.id, "not-a-uuid", "'; drop table"]
    answer = client.get("/home/alive/", {"id": [str(i) for i in asked]})
    assert answer.status_code == 200
    assert answer.json() == {"live": sorted([str(job.id), str(search.id)])}
    assert client.get("/home/alive/").json() == {"live": []}
    assert client.post("/home/alive/").status_code == 405


@pytest.mark.django_db
def test_submit_page_moved_to_jobs_new(client: Client) -> None:
    assert b"<h1>New Job</h1>" in client.get("/jobs/new/").content


@pytest.mark.django_db
def test_footer_on_every_page(client: Client) -> None:
    for url in ("/", "/jobs/new/", "/jobs/", "/homologs/", "/settings/"):
        page = client.get(url).content.decode()
        assert "Copyright © Victor Hanson-Smith" in page, url
        assert "UCSF · HeadwatersDiscovery.org" in page, url


@pytest.mark.django_db
def test_static_files_are_served(client: Client) -> None:
    response = client.get("/static/phylobot/site.css")
    assert response.status_code == 200
    # The High legibility skin's fonts are bundled, not fetched from a font service.
    font = "/static/phylobot/fonts/atkinson-hyperlegible-next-latin-400-normal.woff2"
    assert client.get(font).status_code == 200
    # So are the Broadsheet skin's.
    assert client.get("/static/phylobot/fonts/cinzel-latin-700-normal.woff2").status_code == 200


@pytest.mark.django_db
def test_skin_defaults_to_clean_and_follows_the_cookie(client: Client) -> None:
    page = client.get("/").content.decode()
    assert '<html lang="en" data-skin="clean" data-backdrop="structures"' in page
    assert "phylobot/skins/clean.css" in page and "data-skin-picker" in page
    client.cookies["phylobot_skin"] = "classic"
    assert 'data-skin="classic"' in client.get("/").content.decode()
    client.cookies["phylobot_skin"] = "dark"
    page = client.get("/").content.decode()
    assert 'data-skin="dark"' in page
    assert page.index("skins/clean.css") < page.index("skins/dark.css")  # dark builds on clean
    client.cookies["phylobot_skin"] = "legible"
    page = client.get("/").content.decode()
    assert 'data-skin="legible"' in page and ">High legibility</option>" in page
    client.cookies["phylobot_skin"] = "matrix"
    page = client.get("/").content.decode()
    assert 'data-skin="matrix"' in page and "phylobot/matrix-rain.js" in page
    assert page.index("skins/dark.css") < page.index("skins/matrix.css")  # builds on dark
    client.cookies["phylobot_skin"] = "spectrum"
    assert 'data-skin="spectrum"' in client.get("/").content.decode()
    form = client.get("/jobs/new/").content.decode()  # the steps name their stage
    for stage in ("sequences", "alignments", "trees", "ancestors"):
        assert f'data-stage="{stage}"' in form
    client.cookies["phylobot_skin"] = "lagoon"
    page = client.get("/").content.decode()
    assert 'data-skin="lagoon"' in page and ">Lagoon</option>" in page
    assert page.index("skins/clean.css") < page.index("skins/lagoon.css")  # builds on clean
    client.cookies["phylobot_skin"] = "broadsheet"
    page = client.get("/").content.decode()
    assert 'data-skin="broadsheet"' in page and 'data-backdrop="none"' in page
    assert page.index("skins/clean.css") < page.index("skins/broadsheet.css")
    client.cookies["phylobot_skin"] = "retro95"
    page = client.get("/").content.decode()
    assert 'data-skin="retro95"' in page and ">PhyloBot &#x27;95</option>" in page
    assert page.index("skins/clean.css") < page.index("skins/retro95.css")
    client.cookies["phylobot_skin"] = "sciencefair"
    page = client.get("/").content.decode()
    assert 'data-skin="sciencefair"' in page and 'data-backdrop="none"' in page
    assert page.index("skins/clean.css") < page.index("skins/sciencefair.css")
    for skin in ("dim", "abyss"):  # Dark's remapping of figures and tables, then their own
        client.cookies["phylobot_skin"] = skin
        page = client.get("/").content.decode()
        assert f'data-skin="{skin}"' in page
        assert page.index("skins/dark.css") < page.index(f"skins/{skin}.css")
    client.cookies["phylobot_skin"] = '"><script>'  # unknown skins fall back to the default
    assert 'data-skin="clean"' in client.get("/").content.decode()


@pytest.mark.django_db
def test_only_some_skins_show_a_backdrop(client: Client) -> None:
    for skin in ("clean", "dark", "spectrum", "lagoon", "abyss"):
        client.cookies["phylobot_skin"] = skin
        assert 'data-backdrop="structures"' in client.get("/").content.decode(), skin
    # Quiet whatever the reader chose.
    for skin in ("legible", "matrix", "broadsheet", "dim", "sciencefair", "classic"):
        client.cookies["phylobot_skin"] = skin
        assert 'data-backdrop="none"' in client.get("/").content.decode(), skin
    # The choice is kept while they are on a quiet skin, so coming back restores it, and
    # the Display panel says why nothing is showing.
    page = client.get("/settings/").content.decode()
    assert 'data-backdrop-preference="structures"' in page
    assert '<option value="structures" selected>PDB structures</option>' in page
    assert "data-no-backdrop>" in page
    client.cookies["phylobot_skin"] = "clean"
    assert "data-no-backdrop hidden>" in client.get("/settings/").content.decode()
    # Turning it off is theirs to do, on the skins that do show one.
    client.cookies["phylobot_backdrop"] = "none"
    page = client.get("/").content.decode()
    assert 'data-backdrop="none"' in page and 'data-backdrop-preference="none"' in page
    client.cookies["phylobot_backdrop"] = '"><script>'  # unknown styles fall back
    assert 'data-backdrop-preference="structures"' in client.get("/").content.decode()


@pytest.mark.django_db
def test_backdrop_element_and_its_menu_are_only_where_they_belong(client: Client) -> None:
    menu = 'data-cookie="phylobot_backdrop"'
    for url in ("/", "/about/"):  # draw one, but are not where it is chosen
        page = client.get(url).content.decode()
        assert 'class="backdrop"' in page, url
        assert "/static/phylobot/backdrop/structures.json" in page, url
        assert menu not in page, url
    for url in ("/jobs/new/", "/jobs/", "/homologs/"):  # neither
        page = client.get(url).content.decode()
        assert 'class="backdrop' not in page, url
        assert menu not in page, url
    # Settings draws one as the Display panel's preview, and is the only place to choose.
    page = client.get("/settings/").content.decode()
    assert 'class="backdrop preview"' in page
    assert page.count(menu) == 1 and "<h2>Display</h2>" in page
    assert "<select data-look" not in page.split("<h2>Display</h2>")[0]
    # The structures the "PDB structures" style rains are committed, not fetched live.
    assets = pathlib.Path(phylobot.__file__).parent / "static/phylobot/backdrop"
    entries = json.loads((assets / "structures.json").read_text())
    assert len(entries) >= 8  # enough that the same one rarely falls twice at once
    for entry in entries:
        assert (assets / entry["file"]).exists(), entry["id"]
        assert client.get(f"/static/phylobot/backdrop/{entry['file']}").status_code == 200
    # Editing ENTRIES and re-running scripts/fetch-structures.py is the whole protocol:
    # no image may be left behind that the manifest does not list.
    assert {path.name for path in assets.glob("*.jpeg")} == {e["file"] for e in entries}


@pytest.mark.django_db
def test_backdrop_speed_and_density_are_per_browser_dials(client: Client) -> None:
    page = client.get("/").content.decode()
    assert 'data-backdrop-speed="1" data-backdrop-density="1"' in page  # Normal, Normal
    client.cookies["phylobot_backdrop_speed"] = "2.2"
    client.cookies["phylobot_backdrop_density"] = "0.4"
    page = client.get("/").content.decode()
    assert 'data-backdrop-speed="2.2" data-backdrop-density="0.4"' in page
    settings_page = client.get("/settings/").content.decode()
    assert '<option value="2.2" selected>Fast</option>' in settings_page
    assert '<option value="0.4" selected>Sparse</option>' in settings_page
    # Anything the menus cannot produce falls back to Normal rather than reaching the page.
    client.cookies["phylobot_backdrop_speed"] = "99"
    client.cookies["phylobot_backdrop_density"] = '"><script>'
    assert 'data-backdrop-speed="1" data-backdrop-density="1"' in client.get("/").content.decode()


@pytest.mark.django_db
def test_pages_render_without_leaking_template_comments(client: Client) -> None:
    # A Django {# comment #} is single-line only: a two-line one renders as page text,
    # which no other assertion here would notice.
    for url in ("/", "/about/", "/jobs/new/", "/jobs/", "/homologs/", "/settings/"):
        page = client.get(url).content.decode()
        assert "{#" not in page and "#}" not in page, url


@pytest.mark.django_db
def test_top_bar_marks_the_current_section(client: Client) -> None:
    page = client.get("/jobs/").content.decode()
    assert 'href="/jobs/" aria-current="page"' in page
    assert 'href="/jobs/new/" aria-current' not in page
    assert 'aria-current="page"' not in client.get("/").content.decode()


@pytest.mark.django_db
def test_about_page_cites_the_paper_and_every_method(client: Client) -> None:
    home = client.get("/").content.decode()
    assert 'href="/about/">About</a>' in home  # in the top bar on every page
    page = client.get("/about/").content.decode()
    assert 'href="/about/" aria-current="page">About</a>' in page
    page = page[page.index("<article") :]  # past the contents sidebar
    # The 2016 paper comes first, then what is distinctive, then the methods.
    paper = page.index("doi.org/10.1371/journal.pcbi.1004976")
    assert (
        paper < page.index("What Is Distinctive About PhyloBot 2") < page.index("Included Methods")
    )
    for tool in ("MAFFT", "MUSCLE", "Clustal Omega", "PRANK", "IQ-TREE", "RAxML-NG",
                 "ModelFinder", "GRASP", "Nextflow"):  # fmt: skip
        assert f"<b>{tool}</b>" in page, tool
    assert "doi.org/10.1093/molbev/msag117" in page  # IQ-TREE 3


@pytest.mark.django_db
def test_about_page_draws_the_feature_comparison(client: Client) -> None:
    from phylobot import about

    tools, sections = about.load_comparison()
    assert tools[:2] == ["PhyloBot 2", "PhyloBot 2016"]
    features = [f for _, group in sections for f in group]
    assert all(len(f.marks) == len(tools) for f in features)
    trajectories = next(f for f in features if f.name.startswith("Mutational trajectories"))
    assert trajectories.marks[0].kind == "yes"
    assert trajectories.marks[tools.index("FireProtASR 2.0")].note == "(predicts forward)"
    page = client.get("/about/").content.decode()
    page = page[page.index("<article") :]
    order = ("What Is Distinctive", 'id="comparison"', "Included Methods")
    assert sorted(order, key=page.index) == list(order)
    for tool in tools:
        assert f">{tool}</th>" in page, tool


@pytest.mark.django_db
def test_about_sidebar_links_every_section(client: Client) -> None:
    import re

    page = client.get("/about/").content.decode()
    start = page.index("about-toc")
    sidebar = page[start : page.index("</nav>", start)]
    targets = re.findall(r'href="#([\w-]+)"', sidebar)
    sections = re.findall(r'<section[^>]* id="([\w-]+)"', page)
    assert set(sections) <= set(targets)
    assert all(page.count(f'id="{t}"') == 1 for t in targets)


def test_comparison_marks() -> None:
    from phylobot.about import parse_mark

    assert parse_mark("✓ 4 aligners").note == "4 aligners"
    assert parse_mark("✓?").unsure and parse_mark("✓?").kind == "yes"
    assert parse_mark("?").kind == "unsure" and not parse_mark("?").unsure
    assert parse_mark("◐ MEGA-CC").note == "MEGA-CC" and not parse_mark("◐ MEGA-CC").shown
    assert parse_mark("◐ [cloud not deployed]").shown
    assert parse_mark("◐ [cloud not deployed]").title == "partial: cloud not deployed"
    with pytest.raises(ValueError):
        parse_mark("yes")


def test_every_catalog_model_has_a_citation() -> None:
    from phylobot import about, catalog

    cited = {name for row in about.MATRICES for name in row.names}
    assert cited == {m.name for m in catalog.MATRICES}
    variants = {name for row in about.VARIANTS for name in row.names}
    offered = (*catalog.VARIANTS, *catalog.MIXTURE_VARIANTS)
    assert variants == {v.suffix for v in offered if v.suffix} - {"+F+G4", "+F+R4"}
    assert all(row.refs for row in about.MATRICES)


@pytest.mark.django_db
def test_retro_homepage_furniture(client: Client) -> None:
    # Always in the page (hidden by site.css), so switching skins needs no reload.
    page = client.get("/").content.decode()
    assert 'class="retro-only retro-marquee"' in page and "UNDER CONSTRUCTION" in page
    # The hit counter counts the reader's own jobs.
    assert '<span class="digits">000000</span>' in page
    # The helper, until told to go away.
    assert "data-retro-helper" in page and "phylobot/retro95.js" in page
    client.cookies["phylobot_helper"] = "dismissed"
    page = client.get("/").content.decode()
    assert "data-retro-helper" not in page and "UNDER CONSTRUCTION" in page


@pytest.mark.django_db
def test_skin_menu_has_its_own_order(client: Client) -> None:
    from phylobot import skins

    # The menu lists every skin except Science Fair, which is kept out of the picker.
    menu_keys = {s.key for s in skins.MENU}
    assert "sciencefair" not in menu_keys
    assert menu_keys == {s.key for s in skins.SKINS} - {"sciencefair"}
    page = client.get("/").content.decode()
    menu = page[page.index("data-skin-picker") :]
    assert ">Science Fair</option>" not in menu
    # Clean, Dark, Dark II, Spectrum, High legibility first; then the rest, with Lagoon
    # moved after Broadsheet.
    labels = ["Clean", "Dark", "Dark II", "Spectrum", "High legibility"]
    labels += ["Broadsheet", "Lagoon", "PhyloBot &#x27;95", "Skinless"]
    positions = [menu.index(f">{label}</option>") for label in labels]
    assert positions == sorted(positions)
    # The stylesheets keep the cascade's order: Dark's before Dim's, which build on it.
    assert page.index("skins/dark.css") < page.index("skins/dim.css")
