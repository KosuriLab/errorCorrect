A systematic comparison of error correction enzymes by next-generation sequencing
================================================================================

This repository contains everything needed (sans data) to generate our paper. A pre-print of the manuscript is openly availible on [biorXiv](http://biorxiv.org/content/early/2017/01/15/100685). Running `make all` in `analysis` will process the raw data and generate the manuscript. The `Makefile` is also compatible with parallel builds.

Please note that our pipeline will use a significant amount of RAM, and will generate ~8Gb of intermediate files. Much of this comes from running the pipeline in parallel (especially the filtering steps).

## Subdirectories ##
analysis: contains polished scripts, markdown documents, and figures

rawData: location of sequencing reads and scripts to pull them from the relevant repos

explorations: early exploratory data analysis and outdated/dead end scripts

## Necessary Tools ##
- R and (tidyverse, broom, stringr, magrittr, data.table, dtplyr, scale, grid, gridExtra)
- Python 2.7
- [uta-align](https://pypi.python.org/pypi/uta-align) or [equivalently](https://github.com/biocommons/uta-align)
- GNU Make
- latexmk and a reasonable LaTeX installation (ex TexLive)
- BBMap and Bowtie2 (although not critical)




