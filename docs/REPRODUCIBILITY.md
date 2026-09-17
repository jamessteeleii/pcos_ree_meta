# Reproducibility and storage

## Workflow architecture

The repository separates durable project state from rebuildable machine state:

- Git stores source code, the canonical CSV, documentation, and deliberately published artifacts.
- `renv.lock` records the R package environment.
- `_targets.R` defines the computational workflow.
- `_targets/` and `renv/library/` are local, rebuildable state and are not committed.
- Private administration, article PDFs, backups, and release archives live outside the repository.

The committed `_targets.yaml` sets the default targets store to the project-relative `_targets/` directory. To place the cache on another local disk, copy that file to the ignored `_targets.yaml.local`, change its `store` value, and set `TAR_CONFIG=_targets.yaml.local` in the local environment. Do not commit the local file or environment setting.

## Required software

- R 4.4.3, matching `renv.lock`.
- A matching Rtools/C++ toolchain for packages that compile C++ and Stan code on Windows.
- Git.
- Quarto CLI and its bundled Pandoc.
- A TeX distribution providing XeLaTeX and the LaTeX packages requested by the manuscripts.
- Latin Modern Roman fonts.
- Chrome or Chromium for `webshot2` when regenerating HTML-derived outputs.

Package restoration and model compilation can take substantial time and disk space. The full pipeline fits multiple Stan models and should not be used as a routine smoke test.

## Restore and validate

From the repository root:

```r
install.packages("renv")
renv::restore()
```

Then run the fast project checks:

```text
Rscript R/check_reproducibility.R
Rscript R/check_reproducibility.R --environment
```

The first command checks repository structure, data integrity, portable paths, locked workflow dependencies, and resolvable manuscript citations. The second also verifies the R version, required packages, and external command-line tools.

Inspect the pipeline without fitting models:

```r
targets::tar_manifest()
targets::tar_visnetwork()
```

Prepare all primary, secondary, and sensitivity-analysis datasets and validate their required fields and `brms` formulas without compilation or sampling:

```r
targets::tar_make(names = analysis_preflight)
```

Validate the exact registered informative priors against each arm-model specification, without compilation or sampling:

```r
targets::tar_make(names = model_specification_preflight)
```

Each fitted model has separate `_diagnostics`, `_trace_plot`, `_pp_check`, and `_diagnostic_gate` targets. Downstream estimates are blocked unless the corresponding gate passes. Once every requested model has been fitted, inspect and enforce the combined diagnostics with:

```r
targets::tar_make(names = all_model_diagnostics)
targets::tar_read(all_model_diagnostics)
targets::tar_make(names = all_model_diagnostics_gate)
```

Run individual inexpensive targets with `targets::tar_make(names = ...)`. Run the full pipeline only when the analysis and system toolchain are ready:

```r
targets::tar_make()
```

## Output policy

Plots, tables, and manuscripts created by the workflow are file targets, so deleting or modifying a generated file invalidates its target. Development intermediates stay ignored. Publication artifacts may be committed only when intentionally updated and reviewed; immutable public versions should be attached to a tagged release or deposited in an archival repository.

The historical Git database still contains previously committed caches, article PDFs, and private forms. This phase intentionally does not rewrite history. Anyone publishing or transferring the repository should account for that historical content separately.
