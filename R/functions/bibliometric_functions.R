read_bibliometric_snapshot <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
}

read_bibliometric_web_mentions <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
}

read_bibliometric_web_content_groups <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
}

read_bibliometric_altmetric_snapshot <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, na = c("", "NA"))
}

apply_bibliometric_content_groups <- function(mentions, content_groups) {
  mentions |>
    dplyr::left_join(content_groups, by = "evidence_id") |>
    dplyr::mutate(
      content_group_id = dplyr::coalesce(content_group_id, evidence_id),
      relationship = dplyr::coalesce(relationship, "independent"),
      count_as_unique_item = dplyr::coalesce(count_as_unique_item, TRUE),
      rationale = dplyr::coalesce(rationale, "Independent verified page")
    )
}

validate_bibliometric_snapshot <- function(
  snapshot,
  mentions,
  content_groups,
  altmetrics
) {
  expected_studies <- seq_len(17)
  altmetric_count_columns <- c(
    "news_outlets",
    "blog_mentions",
    "x_users",
    "facebook_pages",
    "reddit_users",
    "youtube_creators",
    "clinical_guideline_sources",
    "policy_document_sources",
    "wikipedia_pages",
    "patents",
    "dimensions_citations",
    "mendeley_readers"
  )
  grouped_mentions <- apply_bibliometric_content_groups(
    mentions,
    content_groups
  )
  duplicate_groups <- grouped_mentions |>
    dplyr::count(content_group_id, name = "records") |>
    dplyr::filter(records > 1L) |>
    dplyr::pull(content_group_id)
  duplicate_group_checks <- grouped_mentions |>
    dplyr::filter(content_group_id %in% duplicate_groups) |>
    dplyr::group_by(content_group_id) |>
    dplyr::summarise(
      canonical_items = sum(count_as_unique_item),
      studies = dplyr::n_distinct(study_id),
      .groups = "drop"
    )

  checks <- tibble::tibble(
    check = c(
      "one row per eligible study",
      "study identifiers are complete",
      "only Kritikou lacks an OpenAlex identifier",
      "citation counts are non-negative",
      "citation rates reproduce from the dated snapshot",
      "citing-type counts do not exceed classified citing works",
      "web evidence uses valid study identifiers",
      "web evidence identifiers are unique",
      "web evidence has HTTP(S) source links",
      "content-group records use known evidence identifiers",
      "duplicate content groups retain one canonical item",
      "duplicate content groups stay within one study",
      "one Altmetric row is present per eligible study",
      "Altmetric statuses use the documented vocabulary",
      "Altmetric records have scores and Details Page URLs",
      "unscored Altmetric rows do not imply zero scores",
      "Altmetric count fields are non-negative",
      "DOI-based Altmetric identifiers match the study snapshot"
    ),
    passed = c(
      nrow(snapshot) == length(expected_studies) && !anyDuplicated(snapshot$study_id),
      length(snapshot$study_id) == length(expected_studies) &&
        setequal(snapshot$study_id, expected_studies),
      length(snapshot$study_id[is.na(snapshot$openalex_id)]) == 1L &&
        snapshot$study_id[is.na(snapshot$openalex_id)] == 11,
      all(snapshot$openalex_citations >= 0, na.rm = TRUE),
      all(
        abs(
          snapshot$citations_per_year -
            snapshot$openalex_citations / snapshot$citation_years
        ) < 0.02,
        na.rm = TRUE
      ),
      all(
        rowSums(
          dplyr::select(
            snapshot,
            citing_articles,
            citing_reviews,
            citing_book_chapters,
            citing_books,
            citing_dissertations,
            citing_preprints,
            citing_paratexts,
            citing_other
          ),
          na.rm = TRUE
        ) <= snapshot$classified_citing_works,
        na.rm = TRUE
      ),
      all(mentions$study_id %in% expected_studies),
      !anyDuplicated(mentions$evidence_id),
      all(grepl("^https?://", mentions$url)),
      all(content_groups$evidence_id %in% mentions$evidence_id) &&
        !anyDuplicated(content_groups$evidence_id),
      nrow(duplicate_group_checks) > 0L &&
        all(duplicate_group_checks$canonical_items == 1L),
      nrow(duplicate_group_checks) > 0L &&
        all(duplicate_group_checks$studies == 1L),
      nrow(altmetrics) == length(expected_studies) &&
        !anyDuplicated(altmetrics$study_id) &&
        setequal(altmetrics$study_id, expected_studies),
      all(
        altmetrics$altmetric_status %in% c(
          "record_found",
          "no_mentions_found",
          "unresolved_no_persistent_identifier"
        )
      ),
      all(
        !is.na(altmetrics$attention_score[altmetrics$altmetric_status == "record_found"])
      ) &&
        all(
          grepl(
            "^https://www\\.altmetric\\.com/details/[0-9]+$",
            altmetrics$altmetric_details_url[
              altmetrics$altmetric_status == "record_found"
            ]
          )
        ),
      all(
        is.na(
          altmetrics$attention_score[
            altmetrics$altmetric_status != "record_found"
          ]
        )
      ),
      all(
        as.matrix(dplyr::select(altmetrics, dplyr::all_of(altmetric_count_columns))) >= 0,
        na.rm = TRUE
      ),
      all(
        tolower(altmetrics$identifier[which(altmetrics$identifier_type == "doi")]) ==
          tolower(
            snapshot$doi[
              match(
                altmetrics$study_id[which(altmetrics$identifier_type == "doi")],
                snapshot$study_id
              )
            ]
          )
      )
    )
  )

  if (!all(checks$passed)) {
    failed <- checks$check[!checks$passed]
    stop(
      "Bibliometric snapshot validation failed: ",
      paste(failed, collapse = "; "),
      call. = FALSE
    )
  }

  checks
}

create_bibliometric_summary <- function(
  snapshot,
  mentions,
  content_groups,
  altmetrics
) {
  grouped_mentions <- apply_bibliometric_content_groups(
    mentions,
    content_groups
  )

  mention_counts <- mentions |>
    dplyr::mutate(
      mention_group = dplyr::case_when(
        source_type == "reddit" ~ "verified_reddit_threads",
        TRUE ~ "verified_other_web_pages"
      )
    ) |>
    dplyr::count(study_id, mention_group, name = "n") |>
    tidyr::pivot_wider(
      names_from = mention_group,
      values_from = n,
      values_fill = 0
    )

  unique_item_counts <- grouped_mentions |>
    dplyr::filter(count_as_unique_item) |>
    dplyr::count(study_id, name = "verified_unique_public_items")

  framing_counts <- grouped_mentions |>
    dplyr::filter(count_as_unique_item) |>
    dplyr::mutate(
      framing_group = dplyr::case_when(
        stringr::str_starts(framing, "supports_") ~
          "verified_supportive_or_amplifying_pages",
        framing %in% c("questions_lower_ree", "mixed_or_corrected") ~
          "verified_critical_or_corrective_pages",
        TRUE ~ "verified_mixed_or_debated_pages"
      )
    ) |>
    dplyr::count(study_id, framing_group, name = "n") |>
    tidyr::pivot_wider(
      names_from = framing_group,
      values_from = n,
      values_fill = 0
    )

  snapshot |>
    dplyr::select(-dplyr::any_of(c("altmetric_score", "altmetric_status"))) |>
    dplyr::left_join(altmetrics, by = "study_id") |>
    dplyr::left_join(mention_counts, by = "study_id") |>
    dplyr::left_join(unique_item_counts, by = "study_id") |>
    dplyr::left_join(framing_counts, by = "study_id") |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("verified_"), ~ tidyr::replace_na(.x, 0L)),
      verified_public_web_total =
        verified_reddit_threads + verified_other_web_pages,
      grey_literature_records = rowSums(
        dplyr::pick(
          citing_dissertations,
          citing_preprints,
          citing_paratexts,
          citing_other
        ),
        na.rm = TRUE
      ),
      extended_scholarly_records = rowSums(
        dplyr::pick(
          citing_reviews,
          citing_book_chapters,
          citing_books,
          citing_dissertations,
          citing_preprints,
          citing_paratexts,
          citing_other
        ),
        na.rm = TRUE
      ),
      author_year = paste0(authors_short, " (", publication_year, ")")
    )
}

render_bibliometric_supplement <- function(
  input,
  snapshot_file,
  mentions_file,
  content_groups_file,
  altmetric_file,
  output
) {
  snapshot_file
  mentions_file
  content_groups_file
  altmetric_file

  quarto::quarto_render(
    input = input,
    output_format = "html",
    output_file = basename(output),
    execute_dir = here::here(),
    quiet = FALSE
  )

  if (!file.exists(output)) {
    stop("Could not publish the bibliometric supplement.", call. = FALSE)
  }

  output
}
