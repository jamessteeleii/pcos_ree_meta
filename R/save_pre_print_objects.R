#!/usr/bin/env Rscript

# Export the small, fixed set of objects needed to render the pre-print without
# connecting the manuscript to the targets cache. Run this script only when the
# analysed results or the primary figures have intentionally been regenerated.

object_names <- c(
  "all_model_diagnostics_gate",
  "postprocessing_validation",
  "main_arm_data_imputed_demographics",
  "main_plus_greek_arm_data_imputed_demographics",
  "main_arm_mean_effects_summary_condition",
  "main_arm_mean_effects_summary_contrast",
  "main_arm_variance_effects_summary_condition",
  "main_arm_variance_effects_summary_contrast",
  "pairwise_mean_effects_summary",
  "pairwise_variance_effects_summary",
  "main_plus_greek_arm_mean_effects_summary",
  "main_plus_greek_arm_variance_effects_summary",
  "baseline_arm_mean_effects_summary",
  "baseline_arm_variance_effects_summary",
  "main_arm_mean_effects_contrast_condition_bmi",
  "main_arm_variance_effects_contrast_condition_bmi",
  "main_arm_mean_effects_contrast_condition_fat_free_mass",
  "main_arm_variance_effects_contrast_condition_fat_free_mass"
)

pre_print_objects <- stats::setNames(
  lapply(object_names, targets::tar_read_raw),
  object_names
)

saveRDS(
  pre_print_objects,
  file.path("manuscript", "pre_print_objects.rds"),
  version = 3,
  compress = "xz"
)

ggplot2::ggsave(
  filename = file.path("plots", "combine_mean_plot.pdf"),
  plot = targets::tar_read_raw("combined_mean_plot"),
  device = grDevices::cairo_pdf,
  width = 10,
  height = 5,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path("plots", "combine_variance_plot.pdf"),
  plot = targets::tar_read_raw("combined_variance_plot"),
  device = grDevices::cairo_pdf,
  width = 10,
  height = 5,
  units = "in"
)
