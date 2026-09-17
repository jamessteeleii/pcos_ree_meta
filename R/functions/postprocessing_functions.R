summarise_posterior_draws <- function(draws, by = character(), probability = 0.95) {
  if (!"draw" %in% names(draws)) {
    stop("Posterior draws must contain a draw column.", call. = FALSE)
  }
  missing_groups <- setdiff(by, names(draws))
  if (length(missing_groups)) {
    stop(
      "Posterior draws are missing grouping columns: ",
      paste(missing_groups, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.numeric(draws$draw) || any(!is.finite(draws$draw))) {
    stop("Posterior draws contain non-finite or non-numeric values.", call. = FALSE)
  }
  if (!is.numeric(probability) || length(probability) != 1L ||
      probability <= 0 || probability >= 1) {
    stop("probability must be a single number between zero and one.", call. = FALSE)
  }

  alpha <- (1 - probability) / 2
  grouped_draws <- if (length(by)) {
    dplyr::group_by(draws, dplyr::across(dplyr::all_of(by)))
  } else {
    draws
  }

  summary <- grouped_draws |>
    dplyr::summarise(
      estimate = mean(draw),
      conf.low = unname(stats::quantile(draw, alpha)),
      conf.high = unname(stats::quantile(draw, 1 - alpha)),
      n_draws = dplyr::n(),
      .groups = "drop"
    )

  if ("cond" %in% names(summary)) {
    summary <- summary |>
      dplyr::mutate(
        .condition_order = match(as.character(cond), c("Control", "PCOS"))
      ) |>
      dplyr::arrange(.condition_order) |>
      dplyr::select(-.condition_order)
  }

  summary
}

summarise_fixed_effect <- function(
    model,
    coefficient,
    estimand,
    exponentiate = FALSE) {
  if (!inherits(model, "brmsfit")) {
    stop("The supplied model is not a brmsfit object.", call. = FALSE)
  }

  variable <- paste0("b_", coefficient)
  posterior_draws <- posterior::as_draws_df(model, variable = variable)
  if (!variable %in% names(posterior_draws)) {
    stop("Model does not contain coefficient ", variable, ".", call. = FALSE)
  }

  values <- posterior_draws[[variable]]
  if (isTRUE(exponentiate)) {
    values <- exp(values)
  }

  draws <- tibble::tibble(
    drawid = posterior_draws$.draw,
    draw = values
  )

  summarise_posterior_draws(draws) |>
    dplyr::mutate(
      coefficient = coefficient,
      estimand = estimand,
      scale = if (isTRUE(exponentiate)) "ratio" else "difference",
      .before = estimate
    )
}

summary_table_is_valid <- function(summary, positive = FALSE) {
  required <- c("estimate", "conf.low", "conf.high", "n_draws")
  if (!all(required %in% names(summary)) || !nrow(summary)) {
    return(FALSE)
  }

  values <- unlist(summary[c("estimate", "conf.low", "conf.high")], use.names = FALSE)
  valid <- all(is.finite(values)) &&
    all(summary$conf.low < summary$conf.high) &&
    all(summary$n_draws > 0)
  if (isTRUE(positive)) {
    valid <- valid && all(values > 0)
  }
  valid
}

study_prediction_draws_are_valid <- function(predictions, data, required_data) {
  required_prediction_columns <- c("drawid", "draw", "study", "cond")
  if (!all(required_prediction_columns %in% names(predictions)) || !nrow(predictions)) {
    return(FALSE)
  }
  if (any(!is.finite(predictions$draw))) {
    return(FALSE)
  }
  if (anyDuplicated(predictions[c("drawid", "study", "cond")])) {
    return(FALSE)
  }

  complete <- stats::complete.cases(data[, required_data, drop = FALSE])
  expected_groups <- unique(data[complete, c("study", "cond"), drop = FALSE])
  actual_groups <- unique(predictions[c("study", "cond")])
  expected_keys <- paste(expected_groups$study, expected_groups$cond, sep = "\r")
  actual_keys <- paste(actual_groups$study, actual_groups$cond, sep = "\r")
  if (!setequal(expected_keys, actual_keys)) {
    return(FALSE)
  }

  draw_counts <- table(predictions$drawid)
  length(draw_counts) > 0 && all(draw_counts == nrow(expected_groups))
}

run_postprocessing_validation <- function(
    main_mean_condition_summary,
    main_variance_condition_summary,
    difference_summaries,
    ratio_summaries,
    mean_study_predictions,
    variance_study_predictions,
    main_data) {
  checks <- tibble::tibble(
    check = c(
      "Mean condition summary order and scale",
      "Variance condition summary order and scale",
      "Difference summaries are finite",
      "Ratio summaries are finite and positive",
      "Mean study predictions are unique by draw and study-condition",
      "Variance study predictions are unique by draw and study-condition"
    ),
    passed = c(
      summary_table_is_valid(main_mean_condition_summary) &&
        identical(as.character(main_mean_condition_summary$cond), c("Control", "PCOS")),
      summary_table_is_valid(main_variance_condition_summary, positive = TRUE) &&
        identical(as.character(main_variance_condition_summary$cond), c("Control", "PCOS")),
      length(difference_summaries) > 0 &&
        all(vapply(difference_summaries, summary_table_is_valid, logical(1))),
      length(ratio_summaries) > 0 &&
        all(vapply(
          ratio_summaries,
          summary_table_is_valid,
          logical(1),
          positive = TRUE
        )),
      study_prediction_draws_are_valid(
        mean_study_predictions,
        main_data,
        required_data = "yi_mean"
      ),
      study_prediction_draws_are_valid(
        variance_study_predictions,
        main_data,
        required_data = c("yi_sd", "log_yi_mean", "se_log_yi_mean")
      )
    )
  )

  failed <- checks$check[!checks$passed]
  if (length(failed)) {
    stop(
      "Post-processing validation failed:\n- ",
      paste(failed, collapse = "\n- "),
      call. = FALSE
    )
  }

  checks
}
