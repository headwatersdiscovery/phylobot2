# phylobot.com

The public site: one static landing page for people who have not installed PhyloBot 2 yet.

```
make website          # build into website/dist/
make website-serve    # build, then preview at http://localhost:8001
```

`website/dist/` is the whole site: plain files with relative links, about 3 MB. Upload the
folder to any static host (GitHub Pages, Cloudflare Pages, Netlify, an S3 bucket); it
works from the root of a domain, from a folder, or opened straight from disk
(`dist/index.html`). It is not committed: rebuild it before each upload.

## The page, top to bottom

It shows before it tells: a visitor should understand what PhyloBot does from the first two
screens without reading a paragraph.

1. **Title**: the app's own title (logo, taglines) over the falling structures, filling the
   first screen, with Quick Start, the tutorial and GitHub. The skin menu sits in the corner.
2. **From Sequences to Ancestors**: the pipeline, played. An invented family of
   106 sequences goes through the five steps (sequences, alignment, tree, ancestors,
   paths). Every sequence is a row, most a hair tall, with six enlarged as under a lens so
   their residues can be read: a first-time visitor should see that real families are
   large. It plays once it scrolls into view, and each step is a button. This is the
   page's one piece of motion.
3. **Tutorial video** (once it exists; see below).
4. **A Look Inside**: screenshots of the viewers and of an ancestral library, each with
   two or three sentences and a caption.
5. **Why Reconstruct Ancestors**: the three audiences.
6. **Built to Be Checked**: four short claims; the About page's full text and the feature
   comparison table are folded underneath for whoever wants them.
7. **Quick Start** and **Developers**: the two command snippets.
8. **How to Cite**, **Included Methods**, **Substitution Models**, **Related Work**: the
   long lists are folded to a heading and the names in them; a link to anything inside one
   opens it, and "Show every citation" opens them all.

A few links to these sections float at the top once the title has scrolled away. That is
not the app's bar: it only moves around this page.

## The look

Premium, and scientifically rigorous: the page is set like a journal. Serif headings
(Source Serif 4, in the skins that have no typeface of their own), hairline rules instead
of boxes, captioned figures, monospace for sequences and tool names, and a lot
of room. There is **one colour**: everything is ink (the text colour, in tints) except what
the page is pointing at, which takes the skin's accent. In the pipeline figure that is the chosen
ancestor, the path and what changed along it; the residues are tints of ink by chemical
class. So the mascot and the screenshots are the only things on the page in full colour,
and nothing competes with the logo on the first screen. Keep it that way: a second colour
in the page's own drawing is what made the first version of the figure look garish.

No colour is the page's own, so every skin in the menu still dresses it.

## Adding the tutorial video

Paste the video's YouTube ID into `site.toml` (`[tutorial] youtube_id`) and rebuild. The
page gains a "Watch the Tutorial" section after the pipeline, a button in the title, and a
link among the section links; while the ID is empty it has none of them. To see the layout
before the video exists:

```
make website-serve ARGS="--youtube-id XXXXXXXXXXX"
```

Nothing is fetched from YouTube until a reader presses play (the panel is a link until
then, and the player comes from `youtube-nocookie.com`). `poster` in `site.toml` names an
image in `website/static/` to show in the panel; without one it is a plain dark panel with
a play button.

## The screenshots

`site.toml` lists them (`[[shots]]`: title, two sentences, alt text, size); the pictures
are `static/shots/<key>-light.webp` and `<key>-dark.webp` (the dark one shows under the dark
skins). Click one on the page to see it at full size.

**Only ever show public data.** These were taken on 2026-11-11 from a job run on the
bundled Src/Abl example's 20-sequence subset (MAFFT and MUSCLE; LG+G4 and WAG+G4), with
PDB 1QCF attached, in the Clean and Dark skins, with headless Chromium at twice the pixel
density, then scaled to 1600 px wide and saved as WebP (quality 82):

| Key | Page | What is in the picture |
|---|---|---|
| `tree` | `/jobs/<id>/tree/` | From the toolbar down, 1102 x 660, with Node15's panel open |
| `compare` | `/jobs/<id>/ancestors/?node=Node10` | "Node10 in the Other Combinations", 920 x 735 |
| `path` | `/jobs/<id>/path/?path_from=Node3&path_to=Homo_sapiens_P08631` | Path Line and Changes, 1180 x 790 |
| `structure` | `/jobs/<id>/structure/?node=Node15&to=Homo_sapiens_P08631` | The 3D panel with its menus, caption and legend, 870 x 700 |
| `properties` | `/jobs/<id>/ancestors/?node=Node10` | The properties panel, 920 x 640 |
| `library` | `/jobs/<id>/library/<library>/` | Members, 880 x 625: a library of that path's seven ancestors, default design |

Retake them when a viewer changes enough to look different, and update the sizes in
`site.toml` if they change.

## Where the content lives

Little here is a second copy of the app, so the two cannot drift:

| On the page | Comes from |
|---|---|
| The long About text, comparison table, models, related work | `src/phylobot/templates/phylobot/about/` (the app's About page includes the same files) |
| Citations | `phylobot.about` |
| Skins and their menu | `phylobot.skins`, `src/phylobot/static/phylobot/skins/` |
| Styles, fonts, logo, backdrop script | `src/phylobot/static/phylobot/`, copied as they are |
| The backdrop's structures | the app's `backdrop/` folder, written into `static/structures.js` (inline, so they rain from disk too) |
| The page's own words and layout | `templates/index.html` |
| The pipeline demo, section links, zoom, video | `static/website.js` |
| Styles only this page needs | `static/website.css` |
| Address, repository, contact, video, screenshots | `site.toml` |

The page's short copy (the five steps, the three audiences, the four claims) is its own,
condensed from the About page: a claim that changes there should change here too. The tools
named under the pipeline must each be one `phylobot.about` cites, or the build stops.

The skin is remembered in the browser (localStorage; the app uses a cookie because its
server draws the page). `tests/test_website.py` builds the site and checks that every file
it names exists, every in-page link has a target, nothing loads from another site, and the
demo's invented sequences are consistent with the changes it shows.
