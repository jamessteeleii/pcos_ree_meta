analysis_complete_rows <- function(data, columns, label) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      label, " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  complete <- stats::complete.cases(data[, columns, drop = FALSE])
  if (!any(complete)) {
    stop(label, " has no complete rows for the requested model.", call. = FALSE)
  }

  data[complete, , drop = FALSE]
}

analysis_expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

analysis_expect_finite <- function(data, columns, label) {
  for (column in columns) {
    analysis_expect(
      is.numeric(data[[column]]) && all(is.finite(data[[column]])),
      paste0(label, " contains non-finite or non-numeric values in ", column, ".")
    )
  }
}

validate_arm_model_data <- function(
    data,
    outcome = c("mean", "variance"),
    moderator = NULL,
    label = "Arm-based analysis data") {
  outcome <- match.arg(outcome)
  response_columns <- if (identical(outcome, "mean")) {
    c("yi_mean", "vi_mean")
  } else {
    c("yi_sd", "vi_sd", "log_yi_mean", "se_log_yi_mean")
  }
  moderator_columns <- if (is.null(moderator)) {
    character()
  } else {
    c(paste0("m_", moderator), paste0("se_", moderator))
  }
  required_columns <- c(
    "lab", "study", "arm", "effect", "cond", "n",
    response_columns, moderator_columns
  )
  model_data <- analysis_complete_rows(data, required_columns, label)

  analysis_expect(
    setequal(unique(as.character(model_data$cond)), c("Control", "PCOS")),
    paste0(label, " must contain complete Control and PCOS observations.")
  )
  analysis_expect(
    all(model_data$n > 1),
    paste0(label, " contains a modelled observation with n <= 1.")
  )
  analysis_expect_finite(model_data, response_columns, label)
  analysis_expect(
    all(model_data[[if (identical(outcome, "mean")) "vi_mean" else "vi_sd"]] > 0),
    paste0(label, " contains a non-positive sampling variance.")
  )
  if (identical(outcome, "variance")) {
    analysis_expect(
      all(model_data$se_log_yi_mean > 0),
      paste0(label, " contains a non-positive standard error for log mean REE.")
    )
  }
  if (!is.null(moderator)) {
    analysis_expect_finite(model_data, moderator_columns, label)
    analysis_expect(
      all(model_data[[paste0("m_", moderator)]] > 0) &&
        all(model_data[[paste0("se_", moderator)]] > 0),
      paste0(label, " contains a non-positive moderator estimate or standard error.")
    )
  }

  analysis_expect(!anyDuplicated(model_data$effect), paste0(label, " contains duplicate effect identifiers."))
  analysis_expect(
    all(vapply(split(model_data$lab, model_data$study), function(x) length(unique(x)) == 1L, logical(1))),
    paste0(label, " maps at least one study to more than one lab.")
  )
  analysis_expect(
    all(vapply(split(model_data$study, model_data$arm), function(x) length(unique(x)) == 1L, logical(1))),
    paste0(label, " maps at least one arm to more than one study.")
  )
  analysis_expect(
    all(vapply(split(model_data$cond, model_data$arm), function(x) length(unique(x)) == 1L, logical(1))),
    paste0(label, " maps at least one arm to more than one condition.")
  )

  data.frame(
    check = label,
    outcome = outcome,
    moderator = if (is.null(moderator)) "none" else moderator,
    rows = nrow(model_data),
    labs = length(unique(model_data$lab)),
    studies = length(unique(model_data$study)),
    arms = length(unique(model_data$arm)),
    stringsAsFactors = FALSE
  )
}

validate_pairwise_model_data <- function(data, label = "Pairwise analysis data") {
  required_columns <- c(
    "lab", "study", "n_PCOS", "n_Control", "mean_PCOS", "mean_Control",
    "sd_PCOS", "sd_Control", "yi_mean", "vi_mean", "yi_cvr", "vi_cvr"
  )
  model_data <- analysis_complete_rows(data, required_columns, label)
  numeric_columns <- setdiff(required_columns, c("lab", "study"))

  analysis_expect_finite(model_data, numeric_columns, label)
  analysis_expect(
    all(model_data$n_PCOS > 1 & model_data$n_Control > 1),
    paste0(label, " contains a comparison with n <= 1.")
  )
  analysis_expect(
    all(model_data$mean_PCOS > 0 & model_data$mean_Control > 0 &
          model_data$sd_PCOS > 0 & model_data$sd_Control > 0),
    paste0(label, " contains a non-positive mean or standard deviation.")
  )
  analysis_expect(
    all(model_data$vi_mean > 0 & model_data$vi_cvr > 0),
    paste0(label, " contains a non-positive sampling variance.")
  )
  analysis_expect(!anyDuplicated(model_data$study), paste0(label, " contains duplicate study comparisons."))

  data.frame(
    check = label,
    outcome = "mean and variance",
    moderator = "none",
    rows = nrow(model_data),
    labs = length(unique(model_data$lab)),
    studies = length(unique(model_data$study)),
    arms = NA_integer_,
    stringsAsFactors = FALSE
  )
}

validate_brms_formula <- function(formula, data, label) {
  tryCatch(
    brms::get_prior(formula = formula, data = data),
    error = function(error) {
      stop(label, " failed brms formula validation: ", conditionMessage(error), call. = FALSE)
    }
  )
  invisible(TRUE)
}

validate_brms_prior <- function(prior, formula, data, label) {
  tryCatch(
    brms::validate_prior(prior = prior, formula = formula, data = data),
    error = function(error) {
      stop(label, " failed brms prior validation: ", conditionMessage(error), call. = FALSE)
    }
  )
  invisible(TRUE)
}

validate_prior_specification <- function(prior, formula, data, label) {
  validate_brms_prior(prior, formula, data, label)
  data.frame(
    check = label,
    prior = "registered informative priors",
    valid = TRUE,
    stringsAsFactors = FALSE
  )
}

run_model_specification_preflight <- function(
    main_data,
    baseline_data,
    plus_greek_data,
    moderator_data,
    mean_prior,
    variance_prior) {
  checks <- list(
    validate_prior_specification(
      mean_prior,
      arm_mean_formula(),
      main_data,
      "Main arm mean-model priors"
    ),
    validate_prior_specification(
      variance_prior,
      arm_variance_formula(),
      main_data,
      "Main arm variance-model priors"
    ),
    validate_prior_specification(
      mean_prior,
      arm_mean_formula(),
      baseline_data,
      "Baseline mean-model priors"
    ),
    validate_prior_specification(
      variance_prior,
      arm_variance_formula(),
      baseline_data,
      "Baseline variance-model priors"
    ),
    validate_prior_specification(
      mean_prior,
      arm_mean_formula(),
      plus_greek_data,
      "Greek-study mean-model priors"
    ),
    validate_prior_specification(
      variance_prior,
      arm_variance_formula(),
      plus_greek_data,
      "Greek-study variance-model priors"
    ),
    validate_prior_specification(
      mean_prior,
      arm_mean_formula("bmi"),
      moderator_data,
      "BMI-adjusted mean-model priors"
    ),
    validate_prior_specification(
      variance_prior,
      arm_variance_formula("bmi"),
      moderator_data,
      "BMI-adjusted variance-model priors"
    ),
    validate_prior_specification(
      mean_prior,
      arm_mean_formula("fat_free_mass"),
      moderator_data,
      "Fat-free-mass-adjusted mean-model priors"
    ),
    validate_prior_specification(
      variance_prior,
      arm_variance_formula("fat_free_mass"),
      moderator_data,
      "Fat-free-mass-adjusted variance-model priors"
    )
  )

  do.call(rbind, checks)
}

run_analysis_preflight <- function(
    main_data,
    baseline_data,
    plus_greek_data,
    pairwise_data,
    moderator_data) {
  analysis_expect(!any(main_data$lab == 8, na.rm = TRUE), "Main analysis data unexpectedly include lab 8.")
  analysis_expect(any(plus_greek_data$lab == 8, na.rm = TRUE), "Greek-study sensitivity data do not include lab 8.")
  analysis_expect(
    all(baseline_data$timepoint == "baseline"),
    "Baseline sensitivity data contain a non-baseline observation."
  )

  checks <- list(
    validate_arm_model_data(main_data, "mean", label = "Main arm data"),
    validate_arm_model_data(main_data, "variance", label = "Main arm data"),
    validate_arm_model_data(baseline_data, "mean", label = "Baseline sensitivity data"),
    validate_arm_model_data(baseline_data, "variance", label = "Baseline sensitivity data"),
    validate_arm_model_data(plus_greek_data, "mean", label = "Greek-study sensitivity data"),
    validate_arm_model_data(plus_greek_data, "variance", label = "Greek-study sensitivity data"),
    validate_pairwise_model_data(pairwise_data),
    validate_arm_model_data(
      moderator_data, "mean", moderator = "bmi", label = "BMI-adjusted data"
    ),
    validate_arm_model_data(
      moderator_data, "variance", moderator = "bmi", label = "BMI-adjusted data"
    ),
    validate_arm_model_data(
      moderator_data, "mean", moderator = "fat_free_mass",
      label = "Fat-free-mass-adjusted data"
    ),
    validate_arm_model_data(
      moderator_data, "variance", moderator = "fat_free_mass",
      label = "Fat-free-mass-adjusted data"
    )
  )

  validate_brms_formula(arm_mean_formula(), main_data, "Main arm mean model")
  validate_brms_formula(arm_variance_formula(), main_data, "Main arm variance model")
  validate_brms_formula(
    arm_mean_formula(), baseline_data, "Baseline sensitivity arm mean model"
  )
  validate_brms_formula(
    arm_variance_formula(), baseline_data, "Baseline sensitivity arm variance model"
  )
  validate_brms_formula(
    arm_mean_formula(), plus_greek_data, "Greek-study sensitivity arm mean model"
  )
  validate_brms_formula(
    arm_variance_formula(), plus_greek_data, "Greek-study sensitivity arm variance model"
  )
  validate_brms_formula(pairwise_mean_formula(), pairwise_data, "Pairwise mean model")
  validate_brms_formula(pairwise_variance_formula(), pairwise_data, "Pairwise variance model")
  validate_brms_formula(arm_mean_formula("bmi"), moderator_data, "BMI-adjusted mean model")
  validate_brms_formula(arm_variance_formula("bmi"), moderator_data, "BMI-adjusted variance model")
  validate_brms_formula(
    arm_mean_formula("fat_free_mass"), moderator_data, "Fat-free-mass-adjusted mean model"
  )
  validate_brms_formula(
    arm_variance_formula("fat_free_mass"), moderator_data, "Fat-free-mass-adjusted variance model"
  )

  do.call(rbind, checks)
}
