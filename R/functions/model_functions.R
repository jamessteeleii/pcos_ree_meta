# Functions for fitting models

# Setup rstan to run quicker
rstan_setup <- function() {
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = 4)
}

# Arm based models
## Mean effects models
set_prior_arm_mean_effects <- function() {
  
  # We take an estimate to set the prior for the CONT arms based on a meta-analysis of two large studies discussed here https://macrofactorapp.com/range-of-bmrs/
  
  data_prior <- tibble(
    study = c(
      "Mifflin, 1990", # DOI: 10.1093/ajcn/51.2.241
      "Pavlidou, 2022" # DOI: 10.3390/metabo13020189
    ),
    n = c(247,549),
    m = c(1349, 1533),
    sd = c(214, 308)
  )
  
  data_prior_mean_effects <- escalc(measure = "MN",
                                   mi = m,
                                   sdi = sd,
                                   ni = n,
                                   data = data_prior)
  
  # We fit a model to estimate with default weakly regularising priors
  estimate_prior_mean_effects <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                                     data = data_prior_mean_effects,
                                     chains = 4,
                                     cores = 4,
                                     seed = 1988,
                                     warmup = 2000,
                                     control = list(adapt_delta = 0.99),
                                     iter = 8000)
  
  
  estimate_prior_mean_effects <- broom.mixed::tidy(estimate_prior_mean_effects)
  
  # The we set the priors taking the estimates from the models (note, df set to 3 to be conservative) 
  prior_arm_mean_effects <-
    c(
      # The prior on the intercept i.e., CONT arms is set from the estimate of the two studies means mentioned above
      set_prior(paste("student_t(3,", estimate_prior_mean_effects$estimate[1],",", estimate_prior_mean_effects$std.error[1],")"),
                class = "b", coef = "Intercept"),
      # The prior on the random effects for intercept i.e., CONT arms is set from the estimate of the two studies variance mentioned above
      set_prior(paste("student_t(3,", estimate_prior_mean_effects$estimate[2],",", estimate_prior_mean_effects$std.error[2],")"),
                class = "sd", coef = "Intercept", group = "study"),
      # The fixed effect coef reflecting the difference between CONT and PCOS is set based on a wide range of possible values
      # This uses the min and max values of ranges reported in the two studies i.e., 2492 - 908 = 1584
      # We then set a student t prior that permits values approximately up to this value with the majority of it's mass centred around zero
      set_prior("student_t(3, 0, 200)", class = "b", coef = "condPCOS")
    )
  
  return(prior_arm_mean_effects)
  
}

fit_arm_mean_effects_model <- function(data, prior) {
  
  arm_model <- brm(yi_mean | se(sqrt(vi_mean)) ~ 0 + Intercept + cond + (1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1|effect),
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   iter = 8000)
}

## Variance effects models
set_prior_arm_variance_effects <- function() {
  
  # We take an estimate to set the prior for the CONT arms based on a meta-analysis of two large studies discussed here https://macrofactorapp.com/range-of-bmrs/
  
  data_prior <- tibble(
    study = c(
      "Mifflin, 1990", # DOI: 10.1093/ajcn/51.2.241
      "Pavlidou, 2022" # DOI: 10.3390/metabo13020189
    ),
    n = c(247,549),
    m = c(1349, 1533),
    sd = c(214, 308)
  )
  
  
  data_prior_variance_effects <- escalc(measure = "SDLN",
                                       mi = m,
                                       sdi = sd,
                                       ni = n,
                                       data = data_prior)
  
  # We fit a model to estimate with default weakly regularising priors
  estimate_prior_variance_effects <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                                         data = data_prior_variance_effects,
                                         chains = 4,
                                         cores = 4,
                                         seed = 1988,
                                         warmup = 2000,
                                         control = list(adapt_delta = 0.99),
                                         iter = 8000)
  
  
  estimate_prior_variance_effects <- broom.mixed::tidy(estimate_prior_variance_effects)
  
  # Re-estimate mean effects so we have a prior for the mean of the means as a covariate in the model (and its se)
  data_prior_mean_effects <- escalc(measure = "MN",
                                   mi = m,
                                   sdi = sd,
                                   ni = n,
                                   data = data_prior)
  
  data_prior_mean_effects <- data_prior_mean_effects |>
    mutate(
      yi_floor = pmax(yi, 1e-6),
      se_log_yi = sqrt(vi) / yi_floor,
      log_yi = log(yi_floor)
    )
  
  # We fit a model to estimate with default weakly regularising priors
  estimate_prior_mean_effects <- brm(log_yi | se(se_log_yi) ~ 1 + (1 | study),
                                     data = data_prior_mean_effects,
                                     chains = 4,
                                     cores = 4,
                                     seed = 1988,
                                     warmup = 2000,
                                     control = list(adapt_delta = 0.99),
                                     iter = 8000)
  
  
  estimate_prior_mean_effects <- broom.mixed::tidy(estimate_prior_mean_effects)
  
  # The we set the priors taking the estimates from the models of 
  prior_arm_variance_effects <-
    c(
      # The prior on the intercept i.e., CONT arms is set from the estimate of the two studies log SDs mentioned above
      set_prior(paste("student_t(3,", estimate_prior_variance_effects$estimate[1],",", estimate_prior_variance_effects$std.error[1],")"),
                class = "b", coef = "Intercept"),
      # The prior on the random effects for intercept i.e., CONT arms is set from the estimate of the two studies variance mentioned above
      set_prior(paste("student_t(3,", estimate_prior_variance_effects$estimate[2],",", estimate_prior_variance_effects$std.error[2],")"),
                class = "sd", coef = "Intercept", group = "study"),
      
      # The fixed effect coef for log(mean) is typically ~1 due to mean-variance relationship being commonplace in other measures 
      # But we set it to be centred there though with a wide scale to indicate uncertainty in this outcome specifically
      set_prior("student_t(3, 0, 2.5)", class = "b", coef = "melog_yi_meanse_log_yi_mean"),
      
      # The fixed effect coef reflecting the difference between CONT and PCOS is set based on a wide range of possible values
      # Given the rough relationship of ~1 for log(mean) on log(sd) in other data we again set it to reflect the range of diffs on the log scale
      # This uses the min and max values of ranges reported in the two studies i.e., 2492 - 908 = 1584
      # We then set a student t prior that permits values approximately up to this value with the majority of it's mass centred around zero
      set_prior("student_t(3, 0, 5.3)", class = "b", coef = "condPCOS"),
      
      # The mean of the log(mean) measurement error has to be positive (as means are positive), as does the sd, so we set these to wide half t distributions
      set_prior(paste("student_t(3,", estimate_prior_mean_effects$estimate[1],",", estimate_prior_mean_effects$std.error[1],")"),
                class = "meanme", coef = "melog_yi_mean"),
      set_prior("student_t(3, 0, 5)", class = "sdme", coef = "melog_yi_mean")
      
    )
  
  return(prior_arm_variance_effects)
  
}

fit_arm_variance_effects_model <- function(data, prior) {
  
  arm_model <- brm(yi_sd | se(sqrt(vi_sd)) ~ 0 + Intercept + cond + me(log_yi_mean, se_log_yi_mean) + (1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1|effect),
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
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
  # pooled preds per condition
  preds <- avg_predictions(
    model,
    newdata = datagrid(
      yi_mean = median(data$yi_mean, na.rm=TRUE),
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
    newdata = data |> filter(!is.na(yi_sd)),
    re_formula = NULL,
  ) |>
    get_draws() |>
    mutate(draw = exp(draw))
  
  return(preds)
}

get_variance_contrast_condition <- function(model, data) {
  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      yi_mean = median(data$yi_mean, na.rm=TRUE)
    ),
    re_formula = NA,
    variables = "cond"
  ) |>
    get_draws() |>
    mutate(draw = exp(draw))
  
  return(preds)
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
  
  pairwise_model <- brm(yi_mean | se(sqrt(vi_mean)) ~ 1 + (1 | lab) + (1 | study),
                        data = data,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        control = list(adapt_delta = 0.99),
                        iter = 8000)
  
  return(pairwise_model)
}

fit_pairwise_variance_model <- function(data) {
  
  pairwise_model <- brm(yi_cvr | se(sqrt(vi_cvr)) ~ 1 + (1 | lab) + (1 | study),
                        data = data,
                        chains = 4,
                        cores = 4,
                        seed = 1988,
                        warmup = 2000,
                        control = list(adapt_delta = 0.99),
                        iter = 8000)
  
  return(pairwise_model)
}

# Additional models
get_predictor_medians <- function(data) {
  # use median BMI in controls
  m_bmi <- data |>
    filter(cond == "Control") |>
    group_by(arm) |>
    slice_head(n=1) |>
    ungroup() |>
    summarise(m_bmi = median(m_bmi, na.rm = TRUE))
  
  # use median BMI in controls
  m_fat_free_mass <- data |>
    filter(cond == "Control") |>
    group_by(arm) |>
    slice_head(n=1) |>
    ungroup() |>
    summarise(m_fat_free_mass = median(m_fat_free_mass, na.rm = TRUE))
  
  predictor_medians <- tibble(
    m_bmi = m_bmi$m_bmi,
    m_fat_free_mass = m_fat_free_mass$m_fat_free_mass
  )
  
  return(predictor_medians)
}
  

fit_arm_mean_effects_model_moderator <- function(data, prior, moderator) {
  
  formula <- as.formula(paste0("yi_mean | se(sqrt(vi_mean)) ~ 0 + Intercept + cond + me(m_",moderator,",se_",moderator,") +", 
                               "(1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1|effect)"))
  
  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   iter = 8000)
}

fit_arm_variance_effects_model_moderator <- function(data, prior, moderator) {
  
  formula <- as.formula(paste0("yi_sd | se(sqrt(vi_sd)) ~ 0 + Intercept + cond + me(m_",moderator,",se_",moderator,") +", 
  "me(log_yi_mean, se_log_yi_mean) + (1 + cond | lab) + (1 + cond | study) + (1 | arm) + (1|effect)"))
  
  arm_model <- brm(formula,
                   data = data,
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   seed = 1988,
                   warmup = 2000,
                   iter = 8000)
}

get_mean_contrast_condition_moderator <- function(model, predictor_medians) {
  
  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      m_bmi = predictor_medians$m_bmi,
      m_fat_free_mass = predictor_medians$m_fat_free_mass
    ),
    re_formula = NA,
    variables = "cond"
  ) 
  
  return(preds)
}

get_variance_contrast_condition_moderator <- function(model, data, predictor_medians) {
  
  # pooled preds per condition
  preds <- avg_comparisons(
    model,
    newdata = datagrid(
      yi_mean = median(data$yi_mean, na.rm=TRUE),
      m_bmi = predictor_medians$m_bmi,
      m_fat_free_mass = predictor_medians$m_fat_free_mass
    ),
    re_formula = NA,
    variables = "cond"
  ) 
  
  return(preds)
}
