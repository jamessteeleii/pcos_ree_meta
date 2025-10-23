# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "here",
    "metafor",
    "brms",
    "marginaleffects",
    "tidybayes",
    "patchwork",
    "flextable",
    "officer"
  ),
  memory = "transient",
  format = "qs",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
)


tar_source("R/functions/.")
# tar_source("R/other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  #### Miscellaneous ----
  
  tar_target(
    setup_for_rstan,
    rstan_setup()
  ),
  
  #### Setting priors for arm-based models ----
  
  tar_target(
    prior_arm_mean_effects,
    set_prior_arm_mean_effects(),
  ),
  
  tar_target(
    prior_arm_variance_effects,
    set_prior_arm_variance_effects(),
  ),
  
  
  #### Example data and functions for pre-reg pipeline ----
  
  # Create and prepare example data for pairwise contrast and arm based models
  tar_target(
    example_pairwise_data,
    create_pairwise_data_example()
  ),
  
  tar_target(
    example_pairwise_data_effects,
    calculate_pairwise_effects_example(example_pairwise_data)
  ),
  
  tar_target(
    example_arm_data,
    create_arm_data_example(example_pairwise_data)
  ),
  
  tar_target(
    example_arm_data_effects,
    calculate_arm_effects_example(example_arm_data)
  ),
  
  # Fitting example models
  
  tar_target(
    example_arm_mean_effects_model,
    fit_arm_mean_effects_model_example(
      example_arm_data_effects,
      prior_arm_mean_effects
    )
  ),
  
  tar_target(
    example_arm_variance_effects_model,
    fit_arm_variance_effects_model_example(
      example_arm_data_effects,
      prior_arm_variance_effects
    )
  ),
  
  tar_target(
    example_pairwise_mean_effects_model,
    fit_pairwise_mean_model_example(
      example_pairwise_data_effects
    )
  ),
  
  tar_target(
    example_pairwise_variance_effects_model,
    fit_pairwise_variance_model_example(
      example_pairwise_data_effects
    )
  ),
  
  #### Main data and analysis ----
  
  # Read and prepare data for arm-based analysis
  tar_target(
    main_arm_data_file,
    here("data", "studies_data.csv"),
    format = "file"
  ),
  
  tar_target(
    main_arm_data,
    prepare_data(main_arm_data_file)
  ),
  
  tar_target(
    main_arm_data_effects,
    calculate_arm_effects(main_arm_data)
  ),
  
  tar_target(
    main_arm_data_imputed_demographics,
    impute_bmi_estimates(main_arm_data)
  ),
  
  tar_target(
    descriptives_table,
    create_descriptives_table(main_arm_data_imputed_demographics)
  ),
  
  tar_target(
    descriptives_table_docx,
    create_descriptives_table(main_arm_data_imputed_demographics)
  ),
  
  # Fitting main analysis models
  
  tar_target(
    main_arm_mean_effects_model,
    fit_arm_mean_effects_model(
      main_arm_data_effects,
      prior_arm_mean_effects
    )
  ),
  
  tar_target(
    main_arm_variance_effects_model,
    fit_arm_variance_effects_model(
      main_arm_data_effects,
      prior_arm_variance_effects
    )
  ),
  
  # Get predictions and contrasts from main models
  
  tar_target(
    main_arm_mean_effects_preds_condition,
    get_mean_preds_condition(main_arm_mean_effects_model)
  ),
  
  tar_target(
    main_arm_mean_effects_preds_study_condition,
    get_mean_preds_study_condition(main_arm_mean_effects_model, 
                                   main_arm_data_effects)
  ),
  
  tar_target(
    main_arm_mean_effects_contrast_condition,
    get_mean_contrast_condition(main_arm_mean_effects_model)
  ),
  
  tar_target(
    main_arm_variance_effects_preds_condition,
    get_variance_preds_condition(main_arm_variance_effects_model,
                                 main_arm_data_effects)
  ),
  
  tar_target(
    main_arm_variance_effects_preds_study_condition,
    get_variance_preds_study_condition(main_arm_variance_effects_model,
                                       main_arm_data_effects)
  ),
  
  tar_target(
    main_arm_variance_effects_contrast_condition,
    get_variance_contrast_condition(main_arm_variance_effects_model,
                                    main_arm_data_effects)
  ),
  
  # Create plots for main models
  
  tar_target(
    meta_mean_pred_plot,
    plot_meta_mean_pred(main_arm_mean_effects_preds_condition)
  ),
  
  tar_target(
    meta_study_mean_pred_plot,
    plot_study_mean_pred(main_arm_mean_effects_preds_study_condition,
                        main_arm_data_effects)
  ),
  
  tar_target(
    meta_mean_contrast_plot,
    plot_mean_contrast(main_arm_mean_effects_contrast_condition)
  ),
  
  tar_target(
    combined_mean_plot,
    combine_mean_plots(meta_mean_pred_plot,
                       meta_study_mean_pred_plot,
                       meta_mean_contrast_plot)
  ),
  
  tar_target(
    combined_mean_plot_tiff,
    ggsave("plots/combine_mean_plot.tiff",
           plot = combined_mean_plot,
           dpi = 300, device = "tiff",
           w = 10, h = 5)
  ),
  
  tar_target(
    meta_variance_pred_plot,
    plot_meta_variance_pred(main_arm_variance_effects_preds_condition)
  ),
  
  tar_target(
    meta_study_variance_pred_plot,
    plot_study_variance_pred(main_arm_variance_effects_preds_study_condition,
                         main_arm_data_effects)
  ),
  
  tar_target(
    meta_variance_contrast_plot,
    plot_variance_contrast(main_arm_variance_effects_contrast_condition)
  ),
  
  tar_target(
    combined_variance_plot,
    combine_variance_plots(meta_variance_pred_plot,
                           meta_study_variance_pred_plot,
                       meta_variance_contrast_plot)
  ),
  
  tar_target(
    combined_variance_plot_tiff,
    ggsave("plots/combine_variance_plot.tiff",
           plot = combined_variance_plot,
           dpi = 300, device = "tiff",
           w = 10, h = 5)
  ),
  
  #### Pairwise sensitivity analysis ----
  
  tar_target(
    pairwise_data,
    prepare_pairwise_data(main_arm_data)
  ),
  
  tar_target(
    pairwise_data_effects,
    calculate_pairwise_effects(pairwise_data)
  ),
  
  tar_target(
    pairwise_mean_effects_model,
    fit_pairwise_mean_model(
      pairwise_data_effects
    )
  ),
  
  tar_target(
    pairwise_variance_effects_model,
    fit_pairwise_variance_model(
      pairwise_data_effects
    )
  ),
  
  #### Sensitivity analysis with only baseline data ----
  
  # Fitting baseline arm analysis models
  
  tar_target(
    baseline_arm_mean_effects_model,
    fit_arm_mean_effects_model(
      main_arm_data_effects |>
        filter(timepoint == "baseline"), # just take baseline
      prior_arm_mean_effects
    )
  ),
  
  tar_target(
    baseline_arm_variance_effects_model,
    fit_arm_variance_effects_model(
      main_arm_data_effects |>
        filter(timepoint == "baseline"), # just take baseline
      prior_arm_variance_effects
    )
  ),
  
  # Get contrasts from baseline arm models
  
  
  tar_target(
    baseline_arm_mean_effects_contrast_condition,
    get_mean_contrast_condition(baseline_arm_mean_effects_model)
  ),
  
  tar_target(
    baseline_arm_variance_effects_contrast_condition,
    get_variance_contrast_condition(baseline_arm_variance_effects_model,
                                    main_arm_data_effects |>
                                      filter(timepoint == "baseline") # just take baseline
    )
  ),
  
  #### Sensitivity analysis removing Greek lab studies ----
  
  
  # Fitting baseline arm analysis models
  
  tar_target(
    no_greek_arm_mean_effects_model,
    fit_arm_mean_effects_model(
      main_arm_data_effects |>
        filter(lab != 8), # remove Greek studies
      prior_arm_mean_effects
    )
  ),
  
  tar_target(
    no_greek_arm_variance_effects_model,
    fit_arm_variance_effects_model(
      main_arm_data_effects |>
        filter(lab != 8), # remove Greek studies
      prior_arm_variance_effects
    )
  ),
  
  # Get contrasts from no_greek arm models
  
  
  tar_target(
    no_greek_arm_mean_effects_contrast_condition,
    get_mean_contrast_condition(no_greek_arm_mean_effects_model)
  ),
  
  tar_target(
    no_greek_arm_variance_effects_contrast_condition,
    get_variance_contrast_condition(no_greek_arm_variance_effects_model,
                                    main_arm_data_effects |>
                                      filter(lab != 8) # remove Greek studies
    )
  )
  
  
)















