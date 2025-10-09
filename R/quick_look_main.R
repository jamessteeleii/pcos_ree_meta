library(tidyverse)
library(patchwork)
library(metafor)
library(brms)
library(marginaleffects)
library(tidybayes)

data <- read_csv("data/studies_data.csv") |>
  
  # convert all to kcal
  mutate(
    across(
      # apply to all numeric columns except 'unit'
      .cols = c(mean:iqr_ffm_adjusted),
      .fns = ~ case_when(
        units == "kj"   ~ .x / 4.184,     # 1 kcal = 4.184 kJ
        units == "MJ"   ~ .x * 1000 / 4.184, # 1 MJ = 1000 kJ
        units == "kcal/kg"  ~ .x * m_body_mass, # note this is just for Pohlmeier et al.
        .default = .x,             # already kcal      
        )
    )
  ) |>
  
  # convert standard errors to standard deviations
  mutate(
    sd = case_when(
     is.na(sd) == TRUE ~ se * sqrt(n),
     .default = sd
    ),
    sd_bm_adjusted  = case_when(
      is.na(sd_bm_adjusted) == TRUE ~ se_bm_adjusted  * sqrt(n),
      .default = sd_bm_adjusted 
    ),
    sd_ffm_adjusted = case_when(
      is.na(sd_ffm_adjusted) == TRUE ~ se_ffm_adjusted * sqrt(n),
      .default = sd_ffm_adjusted
    )
  ) |>
  
  # estimate means and standard deviations from range, iqr, median, and sample size 
  # (see DOI: 10.1186/1471-2288-14-135)
  rowwise() |>
  mutate(
    scenario = case_when(
      # unadjusted REEs
      !is.na(lower_range) & !is.na(upper_range) & is.na(iqr) ~ 1,  # Scenario 1
      is.na(lower_range) & is.na(upper_range) & !is.na(iqr)  ~ 2,  # Scenario 2
      !is.na(lower_range) & !is.na(upper_range) & !is.na(iqr) ~ 3, # Scenario 3
      
      # bm adjusted REEs
      !is.na(lower_range_bm_adjusted) & !is.na(upper_range_bm_adjusted) & is.na(iqr_bm_adjusted) ~ 1,  # Scenario 1
      is.na(lower_range_bm_adjusted) & is.na(upper_range_bm_adjusted) & !is.na(iqr_bm_adjusted)  ~ 2,  # Scenario 2
      !is.na(lower_range_bm_adjusted) & !is.na(upper_range_bm_adjusted) & !is.na(iqr_bm_adjusted) ~ 3, # Scenario 3
      
      # ffm adjusted REEs
      !is.na(lower_range_ffm_adjusted) & !is.na(upper_range_ffm_adjusted) & is.na(iqr_ffm_adjusted) ~ 1,  # Scenario 1
      is.na(lower_range_ffm_adjusted) & is.na(upper_range_ffm_adjusted) & !is.na(iqr_ffm_adjusted)  ~ 2,  # Scenario 2
      !is.na(lower_range_ffm_adjusted) & !is.na(upper_range_ffm_adjusted) & !is.na(iqr_ffm_adjusted) ~ 3, # Scenario 3
      TRUE ~ NA_real_
    ),
    
    # Mean estimates
    mean = case_when(
      is.na(mean) & scenario == 1 ~ (lower_range + 2 * median + upper_range) / 4,
      is.na(mean) & scenario == 2 ~ median,  # Wan: if only IQR given, mean ≈ median
      is.na(mean) & scenario == 3 ~ (lower_range + 2 * median + upper_range) / 4,  # same as scenario 1 for mean
      .default = mean
    ),
    
    mean_bm_adjusted = case_when(
      is.na(mean_bm_adjusted) & scenario == 1 ~ (lower_range_bm_adjusted + 2 * median_bm_adjusted + upper_range_bm_adjusted) / 4,
      is.na(mean_bm_adjusted) & scenario == 2 ~ median_bm_adjusted,  # Wan: if only IQR given, mean ≈ median
      is.na(mean_bm_adjusted) & scenario == 3 ~ (lower_range_bm_adjusted + 2 * median_bm_adjusted + upper_range_bm_adjusted) / 4,  # same as scenario 1 for mean
      .default = mean_bm_adjusted
    ),
    
    mean_ffm_adjusted = case_when(
      is.na(mean_ffm_adjusted) & scenario == 1 ~ (lower_range_ffm_adjusted + 2 * median_ffm_adjusted + upper_range_ffm_adjusted) / 4,
      is.na(mean_ffm_adjusted) & scenario == 2 ~ median_ffm_adjusted,  # Wan: if only IQR given, mean ≈ median
      is.na(mean_ffm_adjusted) & scenario == 3 ~ (lower_range_ffm_adjusted + 2 * median_ffm_adjusted + upper_range_ffm_adjusted) / 4,  # same as scenario 1 for mean
      .default = mean_ffm_adjusted
    ),
    
    # SD estimates
    sd = case_when(
      is.na(sd) & scenario == 1 & n >= 25 ~ (upper_range - lower_range) / 4,
      is.na(sd) & scenario == 1 & n < 25  ~ (upper_range - lower_range) / (2 * qnorm((n - 0.375) / (n + 0.25))),
      is.na(sd) & scenario == 2 ~ iqr / (2 * qnorm(0.75)),  # IQR / 1.349
      is.na(sd) & scenario == 3 & n >= 25 ~ sqrt(((upper_range - lower_range)^2 / 16) + (iqr^2 / (4 * qnorm(0.75)^2))),
      is.na(sd) & scenario == 3 & n < 25  ~ sqrt(((upper_range - lower_range)^2 / (4 * qnorm((n - 0.375) / (n + 0.25)))^2) +
                                       (iqr^2 / (4 * qnorm(0.75)^2))),
      .default = sd
    ),
    
    sd_bm_adjusted = case_when(
      is.na(sd_bm_adjusted) & scenario == 1 & n >= 25 ~ (upper_range_bm_adjusted - lower_range_bm_adjusted) / 4,
      is.na(sd_bm_adjusted) & scenario == 1 & n < 25  ~ (upper_range_bm_adjusted - lower_range_bm_adjusted) / (2 * qnorm((n - 0.375) / (n + 0.25))),
      is.na(sd_bm_adjusted) & scenario == 2 ~ iqr_bm_adjusted / (2 * qnorm(0.75)),  # IQR / 1.349
      is.na(sd_bm_adjusted) & scenario == 3 & n >= 25 ~ sqrt(((upper_range_bm_adjusted - lower_range_bm_adjusted)^2 / 16) + (iqr^2 / (4 * qnorm(0.75)^2))),
      is.na(sd_bm_adjusted) & scenario == 3 & n < 25  ~ sqrt(((upper_range - lower_range_bm_adjusted)^2 / (4 * qnorm((n - 0.375) / (n + 0.25)))^2) +
                                                   (iqr_bm_adjusted^2 / (4 * qnorm(0.75)^2))),
      .default = sd_bm_adjusted
    ),
    
    sd_ffm_adjusted = case_when(
      is.na(sd_ffm_adjusted) & scenario == 1 & n >= 25 ~ (upper_range_ffm_adjusted - lower_range_ffm_adjusted) / 4,
      is.na(sd_ffm_adjusted) & scenario == 1 & n < 25  ~ (upper_range_ffm_adjusted - lower_range_ffm_adjusted) / (2 * qnorm((n - 0.375) / (n + 0.25))),
      is.na(sd_ffm_adjusted) & scenario == 2 ~ iqr_ffm_adjusted / (2 * qnorm(0.75)),  # IQR / 1.349
      is.na(sd_ffm_adjusted) & scenario == 3 & n >= 25 ~ sqrt(((upper_range_ffm_adjusted - lower_range_ffm_adjusted)^2 / 16) + (iqr^2 / (4 * qnorm(0.75)^2))),
      is.na(sd_ffm_adjusted) & scenario == 3 & n < 25  ~ sqrt(((upper_range - lower_range_ffm_adjusted)^2 / (4 * qnorm((n - 0.375) / (n + 0.25)))^2) +
                                                               (iqr_ffm_adjusted^2 / (4 * qnorm(0.75)^2))),
      .default = sd_ffm_adjusted
    )
  ) |>
  ungroup()
  


targets::tar_load(c(main_arm_data_effects, 
                    main_arm_mean_effects_preds_condition,
                    main_arm_mean_effects_preds_study_condition,
                    main_arm_mean_effects_contrast_condition, 
                    main_arm_variance_effects_preds_condition,
                    main_arm_variance_effects_preds_study_condition,
                    main_arm_variance_effects_contrast_condition))


# Meta pred plot
meta_labels <- main_arm_mean_effects_preds_condition |>
  group_by(cond) |>
  mean_qi(draw)

meta_pred_plot <- ggplot(main_arm_mean_effects_preds_condition, aes(x = draw, fill = cond)) +
  stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95) +
  facet_grid(cond~.) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(meta_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{round(draw)} [{round(.lower)}, {round(.upper)}]"),
      x = draw, y = 0.1
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
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))



# Study pred plot

study_pred <- main_arm_mean_effects_preds_study_condition |>
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
      x = 2250, y =reorder(study_label, draw), group = cond
    ),
    size = 2, position = position_dodge(width = 0.75),
    hjust = "inward"
  ) +
  # Add individual study data
  geom_point(
    data = main_arm_data_effects |>
      filter(!is.na(yi_mean)) |>
      mutate(study_label = paste(authors, year)),
    aes(x = yi_mean, y = study_label, color = cond),
    size = 0.25,
    alpha = 0.75,
    position = position_jitterdodge(dodge.width = 1, jitter.height = 0.05) 
  ) +
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
  theme(title = element_text(size=8))





contrast_labels <- meta_slopes |>
  group_by(contrast) |>
  mean_qi(draw)

contrast_plot <- ggplot(meta_slopes, aes(x = draw)) +
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
    title = "Contrasts Between Conditions (Control - PCOS)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))




meta_plots <- (study_pred_plot | (meta_pred_plot / contrast_plot)) + 
  plot_annotation(title = "Mean Resting Energy Expenditure",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect",
              # widths = c(2,1,1)
              )  &
  theme(axis.title.x = element_text(size=10),
        legend.position = "bottom")




# variance model

# pooled preds per condition
meta_pred <- avg_predictions(
  main_arm_variance_effects_model,
  newdata = datagrid(
    yi_mean = median(main_arm_data_effects$yi_mean, na.rm=TRUE),
    cond = unique(main_arm_data_effects$cond)
  ),
  by = "cond",
  re_formula = NA
) |>
  get_draws() |>
  mutate(draw = exp(draw))

# study level preds per condition
study_pred <- predictions(
  main_arm_variance_effects_model,
  re_formula = NULL,
) |>
  get_draws() |>
  mutate(draw = exp(draw))


# Pooled contrast
meta_slopes <- avg_comparisons(
  main_arm_variance_effects_model,
  newdata = datagrid(
    yi_mean = median(main_arm_data_effects$yi_mean, na.rm=TRUE)
  ),
  re_formula = NA,
  variables = "cond"
) |>
  get_draws() |>
  mutate(draw = exp(draw))




# Meta pred plot
meta_labels <- meta_pred |>
  group_by(cond) |>
  mean_qi(draw)

meta_pred_plot <- ggplot(meta_pred, aes(x = draw, fill = cond)) +
  stat_halfeye(slab_alpha = .5, point_size = 0.5, linewidth = 0.5, .width = 0.95) +
  facet_grid(cond~.) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(meta_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{round(draw)} [{round(.lower)}, {round(.upper)}]"),
      x = draw, y = 0.1
    ),
    size = 3
  ) +
  # scale_x_continuous(limits = c(0,500)) +
  labs(
    x = "Resting Energy Expenditure (kcal)",
    fill = "Condition",
    title = "Global Grand Mean Estimates for Condition"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))



# Study pred plot

study_pred <- study_pred |>
  right_join(main_arm_data_effects |>
               filter(!is.na(yi_mean)) |>
               select(study, authors, year),
             by = "study") |>
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
  # # Add individual study data
  # geom_point(
  #   data = main_arm_data_effects |> 
  #     filter(!is.na(yi_mean)) |>
  #     mutate(study_label = paste(authors, year)),
  #   aes(x = yi_mean, y = study_label, color = cond),
  #   position = position_dodge(width = 0.4),
  #   shape = "|"
  # ) +
  stat_halfeye(slab_alpha = .5, point_size = 0.1, linewidth = 0.1, 
               position = position_dodge(width = 0.5), .width = 0.95) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#009E73")) +
  geom_text(
    data = mutate_if(study_labels,
                     is.numeric, round, 2),
    aes(
      label = glue::glue("{cond}: {round(draw)} [{round(.lower)}, {round(.upper)}]"),
      x = 1050, y =reorder(study_label, draw), group = cond
    ),
    size = 2, position = position_dodge(width = 0.75),
    hjust = "inward"
  ) +
  # scale_x_continuous(limits = c(0,500)) +
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
  theme(title = element_text(size=8))





contrast_labels <- meta_slopes |>
  group_by(contrast) |>
  mean_qi(draw)

contrast_plot <- ggplot(meta_slopes, aes(x = draw)) +
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
  # scale_x_continuous(limits = c(0,5)) +
  labs(
    x = "Resting Energy Expenditure Contrast (ratio in kcal)",
    title = "Contrasts Between Conditions (Control:PCOS)"
  ) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(title = element_text(size=8))




meta_plots <- (study_pred_plot | (meta_pred_plot / contrast_plot)) + 
  plot_annotation(title = "Standard Deviation Resting Energy Expenditure",
                  caption = "Point estimates and 95% quantile intervals reported") +
  plot_layout(guides = "collect", axis_titles = "collect",
              # widths = c(2,1,1)
  )  &
  theme(axis.title.x = element_text(size=10),
        legend.position = "bottom")

