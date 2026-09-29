# Analysis data manifest

## Canonical dataset

`studies_data.csv` is the authoritative, version-controlled dataset used by the analysis pipeline. The retired Excel workbook and informal extraction notes are not inputs to the pipeline and are intentionally excluded from the repository.

- Encoding: UTF-8
- Records: 38
- Variables: 89
- Record key: `lab`, `study`, `arm`, and `timepoint`
- SHA-256: `d53f63e600efc6f78b4f6cac94096c2e075042c502c2f85809aa945087478c8c`

The checksum identifies the dataset after restoring the Greek characters in the Kritikou et al. article title. No numeric values were changed during that correction.

## Bibliometric impact snapshots

The bibliometric supplement uses four version-controlled data files rather than making live network requests during the reproducible analysis pipeline:

- `bibliometric_footprint_snapshot.csv` contains one row for each of the 17 eligible studies, with identifiers, OpenAlex citation indicators, annualised citation rates, and citing-work type counts retrieved on 21 September 2026.
- `bibliometric_altmetric_snapshot.csv` contains the Altmetric Attention Score and platform-reported source counts retrieved from free Altmetric Details Pages on 22 September 2026. Fourteen studies had scored records, two DOI records explicitly reported no mentions, and one study without a persistent identifier could not be resolved.
- `bibliometric_web_mentions.csv` contains individually verified public webpages and Reddit threads that directly identify an eligible study by link, title/author reference, or distinctive numerical signature. These records are conservative verified minima rather than exhaustive platform totals.
- `bibliometric_web_content_groups.csv` records exact or explicitly documented repost/share relationships. Every page remains auditable, while each content group contributes only once to deduplicated totals.

Missing bibliometric values mean that a metric or index record was unavailable; they must not be recoded as zero. Altmetric source totals are platform-reported counts of users or source accounts and are not deduplicated page counts. Public-web records include a source URL, retrieval date, framing classification, and a concise rationale for inclusion. The pipeline validates study identifiers, citation-rate calculations, citing-type totals, Altmetric statuses and DOI concordance, evidence identifiers, URLs, and duplicate-group consistency before rendering the supplementary bibliometric report.

- `bibliometric_footprint_snapshot.csv` SHA-256: `9b0fcb0bef7f6bac4d5c043c3d321465ebac4464367af23ff072bd1ad21e8369`
- `bibliometric_altmetric_snapshot.csv` SHA-256: `7621335541c934f89baba3967f08ccd16a12fedf90a58318d3099f5f37f61be6`
- `bibliometric_web_mentions.csv` SHA-256: `b369de5d700c9753a2621bf956a4dc439608ae78fb375e33ce7c4725a430f964`
- `bibliometric_web_content_groups.csv` SHA-256: `99815f4d91685324a52834c9eb5b4a4c496530df832f5abf02cb35cc7d9171c7`

These files describe dissemination and attention, not study quality or causal influence. Updates should use a new retrieval date, retain the prior search rules, review every newly added public-web source manually, update all checksums, and rerun `Rscript R/check_reproducibility.R`.

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
