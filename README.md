# Two-locus selection using `moments`

This repository contains scripts to run the analyses and compile the manuscript
for "Can we distinguish modes of selective interactions using linkage
disequilibrium?". The manuscript is built using
[bookdown](https://bookdown.org/), and the repository was first clone from
Kevin Thornton's very handy [RMarkdown manuscript
template](https://github.com/ThorntonLab/RMarkdown_manuscript_template).

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
