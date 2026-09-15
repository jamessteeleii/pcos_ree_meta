# Analysis data manifest

## Canonical dataset

`studies_data.csv` is the authoritative, version-controlled dataset used by the analysis pipeline. The retired Excel workbook and informal extraction notes are not inputs to the pipeline and are intentionally excluded from the repository.

- Encoding: UTF-8
- Records: 38
- Variables: 89
- Record key: `lab`, `study`, `arm`, and `timepoint`
- SHA-256: `d53f63e600efc6f78b4f6cac94096c2e075042c502c2f85809aa945087478c8c`

The checksum identifies the dataset after restoring the Greek characters in the Kritikou et al. article title. No numeric values were changed during that correction.

## Structure and conventions

Each row represents one study arm at one measurement timepoint. The principal identifiers are followed by bibliographic and condition fields, then repeated groups of summary statistics for age, body mass, fat mass, fat-free mass, height, BMI, and resting energy expenditure.

- Empty CSV fields represent missing values. Do not introduce textual missing-value codes such as `NA`, `N/A`, or `-`.
- Units are recorded in the `units` column and are converted by the scripted data-preparation functions.
- `cond` uses `PCOS` or `Control`.
- `timepoint` uses `baseline`, `mid_intervention`, or `post_intervention`.
- `insulin_resistant_y_n` uses `y`, `n`, or an empty field when unclassified.
- Reported, calculated, and estimated values are retained in the relevant statistic columns. Existing qualifications are recorded in `comments` and the transformation logic in `R/functions/main_data_functions.R`.

The current file does not contain complete cell-level provenance such as source page, extractor, checker, extraction date, and adjudication status. This is a known limitation. Future data collection should add a structured provenance record rather than relying on an informal workbook or notes document.

## Change protocol

1. Edit the CSV as UTF-8 without changing column order or embedded quoted line breaks.
2. Explain the correction and its source in the commit message or accompanying documentation.
3. Update the dimensions or coding conventions above if they change.
4. Recalculate the SHA-256 checksum and update it above.
5. Run `Rscript R/check_reproducibility.R` and resolve every failure.
6. Review the CSV diff before committing. Do not accept a whole-file rewrite caused only by spreadsheet formatting or line-ending changes.
