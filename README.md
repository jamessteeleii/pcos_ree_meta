# Resting energy expenditure of women with and without polycystic ovary syndrome: a systematic review and meta-analysis

This repository contains all the code and data for a project conducting a systematic review and meta-analysis of differences in resting energy expenditure between women with, and without, Polycystic Ovary Syndrome (PCOS).

Preprint, please cite as: Kirwan, R., Peele, L., Nuckols, G., Kohlhoff, G., Cabré, H., Olenick, A., and Steele, J. (2025). Resting energy expenditure of women with and without polycystic ovary syndrome: a systematic review and meta-analysis. medrxiv DOI: [https://doi.org/10.64898/2025.12.03.25341536](https://doi.org/10.64898/2025.12.03.25341536).

## Abstract

Context: Polycystic ovary syndrome (PCOS) is common in reproductive-age women, who often have higher BMI classification. This is assumed to stem from lower resting energy expenditure (REE), influencing lifestyle intervention guidelines. However, evidence for reduced REE in women with PCOS compared with those without is inconsistent. Objective: To systematically search and meta-analyse the existing literature to estimate and describe the difference in REE between women with and without PCOS. Data Sources: A systematic search was conducted using PubMed, Medline and Web of Science databases of published research from January 1990 to January 2025. Study Selection: Studies that measured REE in women living with PCOS, both with and without control arms of women without PCOS, were included. Data Extraction: Bibliometric, demographic, and REE data was extracted by one investigator and checked in triplicate. Data Synthesis: Thirteen studies were included in a Bayesian arm-based multiple condition comparison (i.e., network) type meta-analysis model with informative priors to compare both mean REE, and between person variation in REE, between women with and without PCOS. Mean REE differed between groups by 31 kcal/day [95% quantile interval: -44 to 113 kcal/day] and the contrast ratio for between person standard deviations was 0.98 [95% quantile interval: 0.71 to 1.32]. Conclusions: These findings indicate that REE does not meaningfully differ between women with and without PCOS. Group-level differences in resting energy expenditure are small, insignificant, or not physiologically relevant.

## Reproducibility

The repository contains the source code, canonical analysis dataset, package lockfile, and publication artifacts for the project. It does not contain local package libraries, targets caches, article PDFs, or private journal administration files.

The authoritative dataset is [`data/studies_data.csv`](data/studies_data.csv). Its structure, checksum, limitations, and change protocol are documented in [`data/README.md`](data/README.md). Detailed setup, system requirements, storage policy, and verification commands are in [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Setup

Use R 4.4.3, which matches `renv.lock`. Install `renv`, open `pcos_ree_meta.Rproj`, and restore the recorded packages:

```r
install.packages("renv")
renv::restore()
```

The project also requires a matching Rtools/C++ toolchain, Quarto, XeLaTeX, Latin Modern Roman fonts, and Chrome or Chromium. The Stan models are computationally expensive, so package restoration and a full rebuild may take substantial time and disk space.

Run the fast structural and data checks before fitting models:

```text
Rscript R/check_reproducibility.R
Rscript R/check_reproducibility.R --environment
```

## Targets analysis pipeline

The analysis is defined in [`_targets.R`](_targets.R), with project functions under [`R/functions/`](R/functions/). The committed `_targets.yaml` sets the ignored, project-relative `_targets/` directory as the default store. A machine-specific external store can be selected with an ignored `_targets.yaml.local` file and the `TAR_CONFIG` environment variable.

Inspect the pipeline without fitting models:

```r
targets::tar_manifest()
targets::tar_visnetwork()
```

Prepare every analysis dataset and validate its structure and `brms` formula without compiling or sampling a model:

```r
targets::tar_make(names = analysis_preflight)
```

Validate the exact registered prior specifications against every arm-model formula, also without compiling or sampling:

```r
targets::tar_make(names = model_specification_preflight)
```

Fit models explicitly, then run their corresponding `_diagnostics`, `_trace_plot`, and `_pp_check` targets. Model-dependent estimates and manuscripts are gated: they cannot run until the relevant diagnostics report zero divergences, zero maximum-treedepth hits, R-hat no greater than 1.01, adequate bulk and tail effective sample sizes, and acceptable E-BFMI. After every model has been fitted, verify the combined gate with:

```r
targets::tar_make(names = all_model_diagnostics)
targets::tar_read(all_model_diagnostics)
targets::tar_make(names = all_model_diagnostics_gate)
```

After the model gate passes, rebuild and inspect the post-processing validation target. It checks posterior-summary scales and condition order and verifies that study-level predictions contain exactly one value per posterior draw and study-condition:

```r
targets::tar_make(names = postprocessing_validation)
targets::tar_read(postprocessing_validation)
```

Run individual inexpensive targets with `targets::tar_make(names = ...)`. Run the complete workflow, including figures, tables, and manuscripts, only when the full analysis is intended:

```r
targets::tar_make()
```

Generated files are declared as file targets so deletion or modification invalidates the corresponding target. Publication artifacts may be versioned when deliberately updated and reviewed; development intermediates and caches remain untracked.

## Software and package citations

The [`grateful`](https://pakillo.github.io/grateful/index.html) report is retained in `grateful-report.pdf`. Exact R package versions are recorded in `renv.lock`.

## License

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
