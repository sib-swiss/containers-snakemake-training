![GitHub Release Date](https://img.shields.io/github/release-date/sib-swiss/containers-snakemake-training)
[![DOI](https://zenodo.org/badge/331890430.svg)](https://zenodo.org/badge/latestdoi/331890430)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC_BY--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

# Introduction to containers and Snakemake course website

* Website is hosted at [https://sib-swiss.github.io/containers-snakemake-training/](https://sib-swiss.github.io/containers-snakemake-training/).

* Please refer to [issues](https://github.com/sib-swiss/containers-snakemake-training/issues) for improvements/bugs for course material or the website.

* Any contribution to this course material is highly appreciated :+1:. Please have a look at the [CONTRIBUTING.md](CONTRIBUTING.md) file to learn more on how to contribute.

## Authors

- Antonin Thiébaut [ORCiD](https://orcid.org/0000-0002-7587-5587)
- Geert van Geest [ORCiD](https://orcid.org/0000-0002-1561-078X)
- Patricia Palagi [ORCiD](https://orcid.org/0000-0001-9062-6303)

## How to reuse this material?

This is a combination of the [introduction to containers course](https://github.com/sib-swiss/containers-introduction-training) with a module on Snakemake. The markdown files of the container course are added as a module.

To clone the entire repository:

* Clone the repository:
	```bash
	git clone https://github.com/sib-swiss/containers-snakemake-training.git
	```

* Initialise the submodules:
	```bash
	cd containers-snakemake-training
	git submodule update --init --recursive
	git submodule update
	```

* Update to the most recent version of the submodules:
	```bash
	git submodule update --remote
	```

### How to host the website locally?

Once you have cloned the repo, you can host it on your local browser. The website is generated with [MkDocs](https://www.mkdocs.org/), with the [Material](https://squidfunk.github.io/mkdocs-material/) theme.

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

### How to generate a github page?

For an automatically generated github page, you can run:

```sh
mkdocs gh-deploy
```
