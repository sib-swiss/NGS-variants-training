![build workflow](https://github.com/sib-swiss/NGS-variants-training/actions/workflows/docker-image.yml/badge.svg)
![GitHub Release Date](https://img.shields.io/github/release-date/sib-swiss/ngs-variants-training)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4668337.svg)](https://doi.org/10.5281/zenodo.4668337)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC_BY--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

Hosted at: [https://sib-swiss.github.io/NGS-variants-training/](https://sib-swiss.github.io/NGS-variants-training/)

# Course website

This website is generated with [Zensical](https://zensical.org), and versioned with [mike (for Zensical)](https://github.com/squidfunk/mike).

## tool installation

To [install Zensical](https://zensical.org/docs/get-started/) you can run:

```
pip install zensical
```

Alternatively you can use `uv` or `pixi`.


To install [mike for Zensical](https://github.com/squidfunk/mike?tab=readme-ov-file#installation), you can use:

```
pip install git+https://github.com/squidfunk/mike.git
```

## host locally

```bash
mike serve
```

Check it out with your browser at [http://localhost:8000/](http://localhost:8000/)

## deploy to gh-page

You can deploy the generated website to gh-pages using:

```
mike deploy -p -u <version tag> latest
```

