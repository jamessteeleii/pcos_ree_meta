# Project instructions

## Scope and source of truth

- `data/studies_data.csv` is the authoritative analysis dataset.
- Read `data/README.md` before changing the dataset. Every data change must include a concise rationale, an updated checksum, and a passing reproducibility check.
- Do not recreate or commit the retired extraction workbook, informal extraction notes, article PDFs, journal affirmation forms, or interactive R history.
- Treat `renv.lock` as the authoritative R package environment and `_targets.R` as the authoritative analysis workflow.

## Safe workflow

- Run `Rscript R/check_reproducibility.R` before and after structural changes.
- Do not run the full model pipeline unless the user explicitly requests it. Stan models are computationally expensive.
- For pipeline inspection after restoring the environment, use `targets::tar_manifest()` or `targets::tar_visnetwork()`.
- Keep all committed paths project-relative. A local external targets store may be selected with an ignored `_targets.yaml.local` file and the `TAR_CONFIG` environment variable.
- Do not change model formulas, priors, sampling controls, data transformations, estimands, or reported numerical results outside the requested analytical scope.

## Generated files and storage

- Never commit `_targets/`, `renv/library/`, local caches, logs, lock files, or temporary rendering directories.
- Quarto `.tex` files and `*_files/` directories are build intermediates. Edit the `.qmd` source instead.
- Deliberately published PDFs, plots, and tables may remain versioned when they are part of a documented release. Regenerate them through `targets` and review their diffs before committing.
- Keep private administrative files and third-party article PDFs outside the repository in controlled storage.
- Do not delete the preserved targets cache in the original working copy until a clean rebuild has been verified.

## Verification

- Parse every changed R file.
- Run `Rscript R/check_reproducibility.R` for data, path, dependency, citation, and repository checks.
- After `renv::restore()` on R 4.4.3, run `Rscript R/check_reproducibility.R --environment` to verify the local toolchain.
- When pipeline definitions change, inspect the manifest without fitting models. Run only selected lightweight targets unless a full rebuild is explicitly requested.

## Commit discipline

- Keep cleanup, policy/docs, environment/workflow, and validation changes in focused commits.
- Do not rewrite Git history as part of this project phase.
- Preserve unrelated user changes and report any conflict before proceeding.
