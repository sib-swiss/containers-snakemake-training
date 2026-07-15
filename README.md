![GitHub Release Date](https://img.shields.io/github/release-date/sib-swiss/containers-snakemake-training)
[![DOI](https://zenodo.org/badge/691556045.svg)](https://zenodo.org/badge/latestdoi/691556045)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC_BY--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

# Reproducible and Scalable Research with Snakemake and Software Containers

* Website is hosted at <https://sib-swiss.github.io/containers-snakemake-training/latest/>

* Please refer to [issues](https://github.com/sib-swiss/containers-snakemake-training/issues) for improvements/bugs for course material or the website

* Any contribution to this course material is highly appreciated :+1:. Please have a look at the [CONTRIBUTING.md](CONTRIBUTING.md) file to learn more on how to contribute

## Authors

* Antonin Thiébaut [ORCiD](https://orcid.org/0000-0002-7587-5587)
* Rafael Riudavets Puig [ORCiD](https://orcid.org/0000-0002-2855-9952)
* Geert van Geest [ORCiD](https://orcid.org/0000-0002-1561-078X)
* Patricia Palagi [ORCiD](https://orcid.org/0000-0001-9062-6303)

## How to host the website locally?

Once you have cloned the repo, you can host it on your local browser. The website is generated with [Zensical](https://zensical.org/), a modern static site generator built on the [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme, and configured entirely through the existing `mkdocs.yml` (Zensical reads it natively, no config changes required).

* Clone the repository:

    ```bash
    git clone https://github.com/sib-swiss/containers-snakemake-training.git
    ```

* Create and activate a virtual environment (recommended):

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

* Install Zensical:

    ```bash
    pip install zensical
    ```

* Make sure you are in the repository directory and type:

    ```bash
    zensical serve
    ```

The website will be hosted on your local browser at <http://localhost:8000/>.

> **Note:** Zensical is under active development. Plugin support is currently a subset of the MkDocs plugin ecosystem, with full parity planned for upcoming releases. If you rely on custom plugins in `mkdocs.yml`, check [Zensical's compatibility notes](https://zensical.org/docs/community/faqs/?h=compati#compatibility-and-transition) before switching.

## How to generate a github page?

Zensical does not ship a `gh-deploy` shortcut like MkDocs did. Instead, see [Zensical's publishing guide](https://zensical.org/docs/publish-your-site/) for a ready-to-use GitHub Actions workflow that builds the site with `zensical build --clean` and deploys it to GitHub Pages on every push.

This repository publishes multiple versions of the course material (e.g. `latest/`, `2023.10/`), which requires a versioning tool on top of the base publishing setup. Zensical integrates with a [fork of `mike`](https://zensical.org/docs/setup/versioning/) (the versioning tool originally built for MkDocs) as a bridge solution until native versioning support lands. To publish a new version:

```
mike deploy --push --update-aliases <version> latest
mike set-default --push latest
```

Note that `mike`'s Zensical fork isn't on PyPI yet — it must be installed directly from GitHub:

```
pip install git+https://github.com/squidfunk/mike.git
```
