# Two-locus selection using [`moments`](https://moments.readthedocs.io/en/latest/)

This repository contains scripts to run the analyses and compile the manuscript
for "Local fitness and epistatic effects lead to distinct patterns of linkage
disequilibrium in protein-coding genes" (formerly "Can we distinguish modes of
selective interactions using linkage disequilibrium?"). The manuscript is built
using [bookdown](https://bookdown.org/), and the repository was first cloned
from kevin thornton's very handy [rmarkdown manuscript
template](https://github.com/thorntonlab/rmarkdown_manuscript_template).

## Changelog

Tagged commits correspond to submissions, revisions, etc, with a short
description of changes between versions.

0.1

- Initial bioRxiv submission, March 25, 2021
  ([doi: 10.1101/2021.03.25.437004](https://www.biorxiv.org/content/10.1101/2021.03.25.437004v1.full))

0.2

- Resubmission with major revisions, April 1, 2022

## Required R packages

Manuscript compilation uses bookdown and rmarkdown. At the very least, you need:

* `bookdown`
* `rmarkdown`
* `knitr`

## Required Python packages

For analyses:

- `moments` >= 1.1.0 (includes `scipy` and `numpy`)

For plotting:

- `matplotlib`
- `bokeh`
