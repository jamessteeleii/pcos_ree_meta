# Do women with polycystic ovary syndrome (PCOS) have lower resting energy expenditure than women
without PCOS: a systematic review and meta-analysis

This repository contains all the code and data for a project conducting a systematic review and meta-analysis of differences in resting energy expenditure between women with, and without, Polycystic Ovary Syndrome (PCOS).

## Abstract

*To be completed*

## Reproducibility
This repository contains the necessary files and code to reproduce the analyses, figures, and the manuscript. 

## Usage
To reproduce the analyses, you will need to have R (https://cran.r-project.org/) and RStudio (https://www.rstudio.com/products/rstudio/download/#download) installed on your computer.

To help with reproducibility, this project uses the `renv` R package (see https://rstudio.github.io/renv/articles/renv.html). With `renv`, the state of this R project can be easily loaded as `renv` keeps track of the required R packages (including version), and (if known) the external source from which packages were retrieved (e.g., CRAN, Github). With `renv`, packages are installed to a project specific library rather than your user or system library. The `renv` package must be installed on your machine before being able to benefit from its features. The package can be installed using the following command:

``` r
install.packages("renv")
```

Once you have `renv` installed, you can get a copy of this repository on your machine by clicking the green Code button then choose Download zip. Save to your machine and extract. After extraction, double click the `pcos_ree_meta.Rproj` file in the root directory. This will automatically open RStudio. This will ensure all paths work on your system as the working directory will be set to the location of the `.Rproj` file. Upon opening, RStudio will recognize the `renv` files and you will be informed that the project library is out of sync with the lockfile. At shown in the console pane of RStudio, running `renv::restore()` will install the packages recorded in the lockfile. This could take some time depending on your machine and internet connection.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/). Instead of script based pipelines the `targets` package makes use of functions applied to targets specified within the pipeline. The targets can be viewed in the `_targets.R` file, and any user defined functions are available in `R/functions.r`.

You can view the existing targets pipeline by downloading the `targets_pipeline.html` file and opening it in your browser.

Useful console functions:

- `tar_edit()` opens the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

## Software and packages used

The [`grateful`](https://pakillo.github.io/grateful/index.html) package was used to create citations to all software and packages used in the analysis. The `grateful` report can be viewed by downloading the `grateful-report.pdf` file.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg