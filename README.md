![GitHub Release Date](https://img.shields.io/github/release-date/sib-swiss/containers-snakemake-training)
[![DOI](https://zenodo.org/badge/331890430.svg)](https://zenodo.org/badge/latestdoi/331890430)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC_BY--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

# Reproducible and Scalable Research with Snakemake and Software Containers

* Website is hosted at [https://sib-swiss.github.io/containers-snakemake-training/latest/](https://sib-swiss.github.io/containers-snakemake-training/latest/)

* Please refer to [issues](https://github.com/sib-swiss/containers-snakemake-training/issues) for improvements/bugs for course material or the website

* Any contribution to this course material is highly appreciated :+1:. Please have a look at the [CONTRIBUTING.md](CONTRIBUTING.md) file to learn more on how to contribute

## Authors

* Antonin Thiébaut [ORCiD](https://orcid.org/0000-0002-7587-5587)
* Rafael Riudavets Puig [ORCiD](https://orcid.org/0000-0002-2855-9952)
* Geert van Geest [ORCiD](https://orcid.org/0000-0002-1561-078X)
* Patricia Palagi [ORCiD](https://orcid.org/0000-0001-9062-6303)

## How to host the website locally?

Once you have cloned the repo, you can host it on your local browser. The website is generated with [MkDocs](https://www.mkdocs.org/), with the [Material](https://squidfunk.github.io/mkdocs-material/) theme.

* Clone the repository:
	```bash
	git clone https://github.com/sib-swiss/containers-snakemake-training.git
	```

* Install MkDocs:
	```bash
	pip install mkdocs
	```

* And Material:
	```bash
	pip install mkdocs-material
	```

* Make sure you are in the repository directory and type:
	```bash
	mkdocs serve
	```

The website will be hosted on your local browser at [http://localhost:8000/](http://localhost:8000/).

## How to generate a github page?

For an automatically generated github page, you can run:

```sh
mkdocs gh-deploy
```
