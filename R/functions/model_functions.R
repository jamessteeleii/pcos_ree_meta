# Functions for fitting models

# Setup rstan to run quicker
rstan_setup <- function() {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 4)
}

arm_sampling_control <- function() {
  list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
}

pairwise_mean_sampling_control <- function() {
  list(
    adapt_delta = 0.99,
    max_treedepth = 15
  )
}

pairwise_variance_sampling_control <- function() {
  list(
    adapt_delta = 0.9999,
    max_treedepth = 15
  )
}

# Keep the model formulas in one place so the no-sampling validation target and
# the fitted models always inspect and use exactly the same specifications.
arm_mean_formula <- function(moderator = NULL) {
  moderator_term <- if (is.null(moderator)) {
    ""
  } else {
    paste0(" + me(m_", moderator, ", se_", moderator, ")")
  }

  stats::as.formula(paste0(
    "yi_mean | se(sqrt(vi_mean)) ~ 0 + Intercept + cond",
    moderator_term,
    " + (1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1 | effect)"
  ))
}

arm_variance_formula <- function(moderator = NULL) {
  moderator_term <- if (is.null(moderator)) {
    ""
  } else {
    paste0(" + me(m_", moderator, ", se_", moderator, ")")
  }

  stats::as.formula(paste0(
    "yi_sd | se(sqrt(vi_sd)) ~ 0 + Intercept + cond",
    moderator_term,
    " + me(log_yi_mean, se_log_yi_mean)",
    " + (1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1 | effect)"
  ))
}

pairwise_mean_formula <- function() {
  yi_mean | se(sqrt(vi_mean)) ~ 1 + (1 | lab) + (1 | study)
}

pairwise_variance_formula <- function() {
  yi_cvr | se(sqrt(vi_cvr)) ~ 1 + (1 | lab) + (1 | study)
}

# Arm based models
## Mean effects models
set_prior_arm_mean_effects <- function() {
  # These are the exact priors used in the registered analysis and reported in
  # the manuscript. Define them directly so rebuilding the pipeline does not
  # re-estimate fixed prior constants with diagnostically unstable MCMC fits.
  c(
    brms::set_prior(
      "student_t(3, 1441.81371635515 , 84.5620086616271 )",
      class = "b", coef = "Intercept"
    ),
    brms::set_prior(
      "student_t(3, 149.8864833995 , 82.9076267864978 )",
      class = "sd", coef = "Intercept", group = "study"
    ),
    brms::set_prior(
      "student_t(3, 0, 200)",
      class = "b", coef = "condPCOS"
    )
  )
}

fit_arm_mean_effects_model <- function(data, prior) {
  rstan_setup()
  formula <- arm_mean_formula()
  validate_arm_model_data(data, outcome = "mean", label = "arm mean model")
  validate_brms_prior(prior, formula, data, label = "arm mean model")

  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   control = arm_sampling_control(),
                   iter = 8000)
}

## Variance effects models
set_prior_arm_variance_effects <- function() {
  # Preserve the exact registered and reported prior distributions without
  # rerunning the two-study calibration models that produced these constants.
  c(
    brms::set_prior(
      "student_t(3, 5.54405048313954 , 0.803037850697044 )",
      class = "b", coef = "Intercept"
    ),
    brms::set_prior(
      "student_t(3, 1.08064387735468 , 1.06426381334781 )",
      class = "sd", coef = "Intercept", group = "study"
    ),
    brms::set_prior(
      "student_t(3, 0, 2.5)",
      class = "b", coef = "melog_yi_meanse_log_yi_mean"
    ),
    brms::set_prior(
      "student_t(3, 0, 5.3)",
      class = "b", coef = "condPCOS"
    ),
    brms::set_prior(
      "student_t(3, 7.28248746877373 , 0.624816107345362 )",
      class = "meanme", coef = "melog_yi_mean"
    ),
    brms::set_prior(
      "student_t(3, 0, 5)",
      class = "sdme", coef = "melog_yi_mean"
    )
  )
}

fit_arm_variance_effects_model <- function(data, prior) {
  rstan_setup()
  formula <- arm_variance_formula()
  validate_arm_model_data(data, outcome = "variance", label = "arm variance model")
  validate_brms_prior(prior, formula, data, label = "arm variance model")

  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   control = arm_sampling_control(),
                   iter = 8000)
}

## Get predictions and contrasts

get_mean_preds_condition <- function(model) {
  # pooled preds per condition
  preds <- avg_predictions(
    model,
    by = "cond",
    re_formula = NA
  ) |>
    get_draws()
  
  return(preds)
}

get_mean_preds_study_condition <- function(model, data) {
  # pooled preds per condition
  preds <- predictions(
    model,
    newdata = data |> filter(!is.na(yi_mean)),
    re_formula = NULL,
  ) |>
    get_draws()
  
  return(preds)
}

get_mean_contrast_condition <- function(model) {
  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    re_formula = NA,
    variables = "cond"
  ) |>
    get_draws()
  
  return(preds)
}

get_variance_preds_condition <- function(model, data) {
  reference_values <- get_variance_reference_values(data)

  # pooled preds per condition
  preds <- avg_predictions(
    model,
    newdata = datagrid(
      log_yi_mean = reference_values$log_yi_mean,
      se_log_yi_mean = reference_values$se_log_yi_mean,
      cond = unique(data$cond)
    ),
    by = "cond",
    re_formula = NA
  ) |>
    get_draws() |>
    mutate(draw = exp(draw))
  
  return(preds)
}

get_variance_preds_study_condition <- function(model, data) {
  # pooled preds per condition
  preds <- predictions(
    model,
    newdata = data |>
      filter(!is.na(yi_sd), !is.na(log_yi_mean), !is.na(se_log_yi_mean)),
    re_formula = NULL,
  ) |>
    get_draws() |>
    mutate(draw = exp(draw))
  
  return(preds)
}

get_variance_contrast_condition <- function(model, data) {
  reference_values <- get_variance_reference_values(data)

  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      log_yi_mean = reference_values$log_yi_mean,
      se_log_yi_mean = reference_values$se_log_yi_mean
    ),
    re_formula = NA,
    variables = "cond"
  ) |>
    get_draws() |>
    mutate(draw = exp(draw))
  
  return(preds)
}

get_variance_reference_values <- function(data) {
  complete_data <- data |>
    filter(
      !is.na(yi_mean),
      !is.na(log_yi_mean),
      !is.na(se_log_yi_mean)
    )

  if (nrow(complete_data) == 0) {
    stop("Variance prediction data have no complete mean-effect reference values.", call. = FALSE)
  }

  tibble(
    log_yi_mean = log(median(complete_data$yi_mean)),
    se_log_yi_mean = median(complete_data$se_log_yi_mean)
  )
}


## Create plots

plot_meta_mean_pred <- function(preds) {
  # Meta pred plot
  meta_labels <- preds |>
    group_by(cond) |>
    mean_qi(draw)
  
  meta_pred_plot <- ggplot(preds, aes(x = draw, fill = cond)) +
    stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95) +
    facet_grid(cond~.) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    geom_text(
      data = mutate_if(meta_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{round(draw)} [{round(.lower)}, {round(.upper)}]"),
        x = draw, y = 0.2
      ),
      size = 3
    ) +
    labs(
      x = "Resting Energy Expenditure (kcal)",
      fill = "Condition",
      title = "Global Grand Mean Estimates for Condition"
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank()) +
    theme(title = element_text(size=8))
  
  return(meta_pred_plot)
}

plot_study_mean_pred <- function(preds, data) {
  
  study_pred <- preds |>
    mutate(
      study_label = paste(authors, year)
    ) |>
    group_by(study) |>
    mutate(mean_draw = mean(draw)) |>
    ungroup()
  
  study_labels <- study_pred |>
    group_by(study_label, cond) |>
    mean_qi(draw) 
  
  
  # quick filter for initial incomplete plot
  
  
  study_pred_plot <- ggplot(study_pred, aes(x = draw, 
                                            y = reorder(study_label, mean_draw), 
                                            fill = cond)) +
    stat_halfeye(slab_alpha = .5, point_size = 0.1, linewidth = 0.1, 
                 position = position_dodge(width = 0.5), .width = 0.95) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    geom_text(
      data = mutate_if(study_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{cond}: {round(draw)} [{round(.lower)}, {round(.upper)}]"),
        x = Inf, y =reorder(study_label, draw), group = cond
      ),
      size = 2, position = position_dodge(width = 0.75),
      hjust = 1.1
    ) +
    # Add individual study data
    geom_point(
      data = data |> 
        filter(!is.na(yi_sd)) |>
        mutate(study_label = paste(authors, year)),
      aes(x = yi_mean, y = study_label, color = cond),
      size = 0.25,
      alpha = 0.75,
      position = position_jitterdodge(dodge.width = 1, jitter.height = 0.05) 
    ) +
    scale_x_continuous(limits = c(750,2250)) +
    labs(
      x = "Resting Energy Expenditure (kcal)",
      title = "Conditional Estimates for Condition by Study"
    ) +
    guides(
      fill = "none",
      color = "none"
    ) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    theme(title = element_text(size=8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  return(study_pred_plot)
}

plot_mean_contrast <- function(contrasts) {
  
  contrast_labels <- contrasts |>
    group_by(contrast) |>
    mean_qi(draw)
  
  contrast_plot <- ggplot(contrasts, aes(x = draw)) +
    geom_vline(xintercept = 0, lty = "dashed", size = 0.25, alpha = 0.75) +
    stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95, fill = "black") +
    geom_text(
      data = mutate_if(contrast_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{round(draw)} [{round(.lower)}, {round(.upper)}]"),
        x = draw, y = 0.1
      ),
      size = 3
    ) +
    scale_x_continuous(labels = ~sub("-", "\u2212", .x)) +
    labs(
      x = "Resting Energy Expenditure Contrast (difference in kcal)",
      title = "Contrasts Between Conditions (PCOS - Control)"
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(title = element_text(size=8),
          panel.grid = element_blank())
  
  return(contrast_plot)

}

combine_mean_plots <- function(meta_pred_plot,
                               study_pred_plot,
                               contrast_plot) {
  
  meta_plots <- (study_pred_plot | (meta_pred_plot / contrast_plot)) + 
    plot_annotation(title = "Mean Resting Energy Expenditure",
                    caption = "Point estimates and 95% quantile intervals reported") +
    plot_layout(guides = "collect", axis_titles = "collect",
                # widths = c(2,1,1)
    )  &
    theme(axis.title.x = element_text(size=10),
          legend.position = "bottom")
  
  return(meta_plots)
}


plot_meta_variance_pred <- function(preds) {
  # Meta pred plot
  meta_labels <- preds |>
    group_by(cond) |>
    mean_qi(draw)
  
  meta_pred_plot <- ggplot(preds, aes(x = draw, fill = cond)) +
    stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95) +
    facet_grid(cond~.) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    geom_text(
      data = mutate_if(meta_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{round(draw)} [{round(.lower)}, {round(.upper)}]"),
        x = draw, y = 0.2
      ),
      size = 3
    ) +
    labs(
      x = "Resting Energy Expenditure (kcal)",
      fill = "Condition",
      title = "Global Grand Mean Estimates for Condition"
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank()) +
    theme(title = element_text(size=8))
  
  return(meta_pred_plot)
}

plot_study_variance_pred <- function(preds, data) {
  
  study_pred <- preds |>
    mutate(
      study_label = paste(authors, year)
    ) |>
    group_by(study) |>
    mutate(mean_draw = mean(draw)) |>
    ungroup()
  
  study_labels <- study_pred |>
    group_by(study_label, cond) |>
    mean_qi(draw) 
  
  
  # quick filter for initial incomplete plot
  
  
  study_pred_plot <- ggplot(study_pred, aes(x = draw, 
                                            y = reorder(study_label, mean_draw), 
                                            fill = cond)) +
    stat_halfeye(slab_alpha = .5, point_size = 0.1, linewidth = 0.1, 
                 position = position_dodge(width = 0.5), .width = 0.95) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
    geom_text(
      data = mutate_if(study_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{cond}: {round(draw)} [{round(.lower)}, {round(.upper)}]"),
        x = Inf, y =reorder(study_label, draw), group = cond,
      ),
      size = 2, position = position_dodge(width = 0.75),
      hjust = 1.1
    ) +
    # Add individual study data
    geom_point(
      data = data |> 
        filter(!is.na(yi_sd)) |>
        mutate(study_label = paste(authors, year)),
      aes(x = exp(yi_sd), y = study_label, color = cond),
      size = 0.25,
      alpha = 0.75,
      position = position_jitterdodge(dodge.width = 1, jitter.height = 0.05) 
    ) +
    scale_x_continuous(limits = c(0,750)) +
    labs(
      x = "Resting Energy Expenditure (kcal)",
      title = "Conditional Estimates for Condition by Study"
    ) +
    guides(
      fill = "none",
      color = "none"
    ) +
    theme_bw() +
    theme(axis.title.y = element_blank()) +
    theme(title = element_text(size=8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  return(study_pred_plot)
}

plot_variance_contrast <- function(contrasts) {
  
  contrast_labels <- contrasts |>
    group_by(contrast) |>
    mean_qi(draw)
  
  contrast_plot <- ggplot(contrasts, aes(x = draw)) +
    geom_vline(xintercept = 1, lty = "dashed", size = 0.25, alpha = 0.75) +
    stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95, fill = "black") +
    geom_text(
      data = mutate_if(contrast_labels,
                       is.numeric, round, 2),
      aes(
        label = glue::glue("{round(draw,2)} [{round(.lower,2)}, {round(.upper,2)}]"),
        x = draw, y = 0.1
      ),
      size = 3
    ) +
    scale_x_continuous(labels = ~sub("-", "\u2212", .x)) +
    labs(
      x = "Resting Energy Expenditure Contrast (ratio of standard deviations)",
      title = "Contrasts Between Conditions (PCOS:Control)"
    ) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    theme(title = element_text(size=8),
          panel.grid = element_blank())
  
  return(contrast_plot)
  
}

combine_variance_plots <- function(meta_pred_plot,
                               study_pred_plot,
                               contrast_plot) {
  
  meta_plots <- (study_pred_plot | (meta_pred_plot / contrast_plot)) + 
    plot_annotation(title = "Standard Deviation of Resting Energy Expenditure",
                    caption = "Point estimates and 95% quantile intervals reported") +
    plot_layout(guides = "collect", axis_titles = "collect",
                # widths = c(2,1,1)
    )  &
    theme(axis.title.x = element_text(size=10),
          legend.position = "bottom")
  
  return(meta_plots)
}


# Pairwise models
## Note, we just utilise the weakly regularising default priors from brms for both pairwise models
fit_pairwise_mean_model <- function(data) {
  rstan_setup()
  formula <- pairwise_mean_formula()
  validate_pairwise_model_data(data, label = "pairwise mean model")
  validate_brms_formula(formula, data, label = "pairwise mean model")

  pairwise_model <- brm(formula,
                        data = data,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        control = pairwise_mean_sampling_control(),
                        iter = 8000)
  
  return(pairwise_model)
}

fit_pairwise_variance_model <- function(data) {
  rstan_setup()
  formula <- pairwise_variance_formula()
  validate_pairwise_model_data(data, label = "pairwise variance model")
  validate_brms_formula(formula, data, label = "pairwise variance model")

  pairwise_model <- brm(formula,
                        data = data,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        control = pairwise_variance_sampling_control(),
                        iter = 8000)
  
  return(pairwise_model)
}

# Additional models
get_predictor_medians <- function(data) {
  # Use one observation per control arm to define reference moderator values.
  control_arms <- data |>
    filter(cond == "Control") |>
    group_by(arm) |>
    slice_head(n=1) |>
    ungroup()

  predictor_medians <- control_arms |>
    summarise(
      m_bmi = median(m_bmi, na.rm = TRUE),
      se_bmi = median(se_bmi, na.rm = TRUE),
      m_fat_free_mass = median(m_fat_free_mass, na.rm = TRUE),
      se_fat_free_mass = median(se_fat_free_mass, na.rm = TRUE)
    )
  
  return(predictor_medians)
}
  

fit_arm_mean_effects_model_moderator <- function(data, prior, moderator) {
  rstan_setup()
  formula <- arm_mean_formula(moderator)
  validate_arm_model_data(
    data,
    outcome = "mean",
    moderator = moderator,
    label = paste(moderator, "moderator arm mean model")
  )
  validate_brms_prior(
    prior,
    formula,
    data,
    label = paste(moderator, "moderator arm mean model")
  )
  
  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   control = arm_sampling_control(),
                   iter = 8000)
}

fit_arm_variance_effects_model_moderator <- function(data, prior, moderator) {
  rstan_setup()
  formula <- arm_variance_formula(moderator)
  validate_arm_model_data(
    data,
    outcome = "variance",
    moderator = moderator,
    label = paste(moderator, "moderator arm variance model")
  )
  validate_brms_prior(
    prior,
    formula,
    data,
    label = paste(moderator, "moderator arm variance model")
  )
  
  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   control = arm_sampling_control(),
                   iter = 8000)
}

get_mean_contrast_condition_moderator <- function(model, predictor_medians) {
  
  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      m_bmi = predictor_medians$m_bmi,
      se_bmi = predictor_medians$se_bmi,
      m_fat_free_mass = predictor_medians$m_fat_free_mass,
      se_fat_free_mass = predictor_medians$se_fat_free_mass
    ),
    re_formula = NA,
    variables = "cond"
  ) 
  
  return(preds)
}

get_variance_contrast_condition_moderator <- function(model, data, predictor_medians) {
  reference_values <- get_variance_reference_values(data)

  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      log_yi_mean = reference_values$log_yi_mean,
      se_log_yi_mean = reference_values$se_log_yi_mean,
      m_bmi = predictor_medians$m_bmi,
      se_bmi = predictor_medians$se_bmi,
      m_fat_free_mass = predictor_medians$m_fat_free_mass,
      se_fat_free_mass = predictor_medians$se_fat_free_mass
    ),
    re_formula = NA,
    variables = "cond"
  ) 
  
  return(preds)
}
