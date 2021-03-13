# Two-locus selection using [`moments`](https://moments.readthedocs.io/en/latest/)

This repository contains scripts to run the analyses and compile the manuscript
for "Can we distinguish modes of selective interactions using linkage
disequilibrium?". The manuscript is built using
[bookdown](https://bookdown.org/), and the repository was first cloned from
Kevin Thornton's very handy [RMarkdown manuscript
template](https://github.com/ThorntonLab/RMarkdown_manuscript_template).

## Changelog

Tagged commits correspond to submissions, revisions, etc, with a short
description of changes between versions.

- Initial BioRxiv submission, March XX, 2021
  ([https://www.biorxiv.org/](https://www.biorxiv.org/))
  - General description
  - other notes?


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
