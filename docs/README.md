# Maintaining the SpinDynamica Website

This directory (`docs/`) contains the Jekyll source for the SpinDynamica GitHub Pages site at:

`https://christianbengs.github.io/SpinDynamica-pub/`

`README.md` is excluded in `docs/_config.yml` and is not published on the live website.

## Directory Overview

- `docs/_config.yml` — Site metadata, subpath (`baseurl: "/SpinDynamica-pub"`), Kramdown/MathJax settings, and Jekyll collections.
- `docs/_data/navigation.yml` — Top-level navigation bar items.
- `docs/_data/release.yml` — Single source of truth for release version, tag, date, and download links.
- `docs/_layouts/` & `docs/_includes/` — Page templates, header, footer, citation block, download button, and MathJax setup.
- `docs/assets/css/style.css` — Custom academic stylesheet and Rouge syntax-highlighting styles.
- `docs/*.md` — Top-level pages (`index.md`, `download.md`, `installation.md`, `documentation.md`, `examples.md`, `about.md`).
- `docs/_documentation/` — Collection of web documentation pages.
- `docs/_examples/` — Collection of worked web examples.

## 1. Editing Website Content

Edit any Markdown file under `docs/` directly. Keep these conventions in mind:

- **Internal links:** Always use Liquid's `relative_url` filter so links work under the `/SpinDynamica-pub/` subpath (e.g., `[Installation]({{ '/installation/' | relative_url }})`).
- **Mathematica code blocks:** Fence examples with ` ```mathematica ` for syntax highlighting. If a code block contains nested curly braces (`{{` or `}}`), wrap the fenced block in `{% raw %}` and `{% endraw %}` so Jekyll's Liquid parser does not treat `{{` as a template variable.
- **Mathematics (LaTeX / MathJax):** Always use `$$ ... $$` for both inline math (inside a paragraph) and display equations (on its own paragraph surrounded by blank lines). Single `$` delimiters are intentionally disabled so Mathematica symbols such as `$Path` and `$SpinDynamicaVersion` never collide with math rendering.

## 2. Adding New Documentation Pages

Create a new `.md` file inside `docs/_documentation/` (for example, `docs/_documentation/spin-operators.md`) with front matter:

```yaml
---
title: Spin Operators
order: 2
summary: Constructing Cartesian, shift, and spherical tensor operators in SpinDynamica.
---
```

The page is automatically published at `/documentation/spin-operators/` and listed in order on the Documentation page.

## 3. Adding Worked Examples

Create a new `.md` file inside `docs/_examples/` (for example, `docs/_examples/spin-echo.md`) with front matter:

```yaml
---
title: Hahn Spin Echo Simulation
order: 3
summary: Refocusing chemical-shift evolution with a spin-echo pulse sequence.
---
```

The page is automatically published at `/examples/spin-echo/` and listed in order on the Examples page.

## 4. Updating Displayed Release Information

Do not hardcode version numbers across pages. Edit only `docs/_data/release.yml`:

1. Publish the versioned release (and optional archive asset) on GitHub Releases.
2. In `docs/_data/release.yml`:
   - Set `published: true`.
   - Update `version` (e.g. `"3.10.0"`), `tag` (e.g. `"v3.10.0"`), `release_date`, and `asset_name` (if attaching a release archive; leave `""` to link to the release tag page).
   - Optionally add `notes_summary` and move earlier releases into `previous_releases`.

The Home download button, release metadata line, and Download page update automatically from this single file.

## 5. Previewing the Website Locally

### Option A: Using Ruby and Bundler

Install Ruby (3.x) and Bundler, then create a temporary `Gemfile` (excluded by `docs/_config.yml`) inside `docs/` with the GitHub Pages gem:

```ruby
source "https://rubygems.org"
gem "github-pages", group: :jekyll_plugins
```

Run:

```bash
cd docs
bundle install
bundle exec jekyll serve
```

Then open `http://localhost:4000/SpinDynamica-pub/` in your browser.

### Option B: Using Docker (No Local Ruby Required)

From the repository root:

```bash
docker run --rm -it -p 4000:4000 -v "$PWD/docs:/srv/jekyll" jekyll/jekyll:4 jekyll serve --baseurl /SpinDynamica-pub
```

## 6. Deploying Through GitHub Pages

GitHub Pages is configured to build from the `docs/` directory on the `main` branch.

1. Commit your changes inside `docs/`:
   ```bash
   git add docs
   git commit -m "Update website content"
   ```
2. Push to `main`:
   ```bash
   git push origin main
   ```
3. GitHub Pages automatically rebuilds and publishes the site to `https://christianbengs.github.io/SpinDynamica-pub/`.
