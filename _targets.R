# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)

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
    "officer",
    "webshot2"
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

  #### Setting priors for arm-based models ----
  
  tar_target(
    prior_arm_mean_effects,
    set_prior_arm_mean_effects(),
  ),
  
  tar_target(
    prior_arm_variance_effects,
    set_prior_arm_variance_effects(),
  ),
  #### Main data and analysis ----
  
  # Read and prepare data for arm-based analysis
  tar_target(
    main_arm_data_file,
    here("data", "studies_data.csv"),
    format = "file"
  ),

  tar_target(
    data_manifest_file,
    here("data", "README.md"),
    format = "file"
  ),

  tar_target(
    reproducibility_checks,
    check_project_reproducibility(main_arm_data_file, data_manifest_file)
  ),
  
  tar_target(
    main_arm_data,
    {
      reproducibility_checks
      prepare_data(main_arm_data_file) |>
        filter(lab != 8) # remove Greek studies for main model
    }
  ),
  
  tar_target(
    main_arm_data_effects,
    calculate_arm_effects(main_arm_data)
  ),
  
  tar_target(
    main_arm_data_imputed_demographics,
    impute_bmi_estimates(main_arm_data)
  ),
  
  # Make descriptives table with all studies including Greek ones
  
  tar_target(
    main_plus_greek_arm_data,
    {
      reproducibility_checks
      prepare_data(main_arm_data_file)
    }
  ),
  
  tar_target(
    main_plus_greek_arm_data_effects,
    calculate_arm_effects(main_plus_greek_arm_data)
  ),
  
  tar_target(
    main_plus_greek_arm_data_imputed_demographics,
    impute_bmi_estimates(main_plus_greek_arm_data)
  ),
  
  tar_target(
    descriptives_table,
    create_descriptives_table(main_plus_greek_arm_data_imputed_demographics)
  ),
  
  tar_target(
    descriptives_table_html,
    convert_descriptives_table_to_html(descriptives_table),
    format = "file"
  ),

  tar_target(
    descriptives_table_docx,
    convert_descriptives_table_to_docx(descriptives_table),
    format = "file"
  ),

  # Fitting main analysis models
  
  tar_target(
    main_arm_mean_effects_model,
    {
      model_specification_preflight
      fit_arm_mean_effects_model(
        main_arm_data_effects,
        prior_arm_mean_effects
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_model,
    {
      model_specification_preflight
      fit_arm_variance_effects_model(
        main_arm_data_effects,
        prior_arm_variance_effects
      )
    }
  ),
  
  # Get predictions and contrasts from main models
  
  tar_target(
    main_arm_mean_effects_preds_condition,
    {
      main_arm_mean_effects_diagnostic_gate
      get_mean_preds_condition(main_arm_mean_effects_model)
    }
  ),
  
  tar_target(
    main_arm_mean_effects_preds_study_condition,
    {
      main_arm_mean_effects_diagnostic_gate
      get_mean_preds_study_condition(
        main_arm_mean_effects_model,
        main_arm_data_effects
      )
    }
  ),
  
  tar_target(
    main_arm_mean_effects_contrast_condition,
    {
      main_arm_mean_effects_diagnostic_gate
      get_mean_contrast_condition(main_arm_mean_effects_model)
    }
  ),
  
  tar_target(
    main_arm_variance_effects_preds_condition,
    {
      main_arm_variance_effects_diagnostic_gate
      get_variance_preds_condition(
        main_arm_variance_effects_model,
        main_arm_data_effects
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_preds_study_condition,
    {
      main_arm_variance_effects_diagnostic_gate
      get_variance_preds_study_condition(
        main_arm_variance_effects_model,
        main_arm_data_effects
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_contrast_condition,
    {
      main_arm_variance_effects_diagnostic_gate
      get_variance_contrast_condition(
        main_arm_variance_effects_model,
        main_arm_data_effects
      )
    }
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
    save_plot_tiff(
      plot = combined_mean_plot,
      path = file.path("plots", "combine_mean_plot.tiff"),
      width = 10,
      height = 5
    ),
    format = "file"
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
    save_plot_tiff(
      plot = combined_variance_plot,
      path = file.path("plots", "combine_variance_plot.tiff"),
      width = 10,
      height = 5
    ),
    format = "file"
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
    {
      model_specification_preflight
      fit_pairwise_mean_model(
        pairwise_data_effects
      )
    }
  ),
  
  tar_target(
    pairwise_variance_effects_model,
    {
      model_specification_preflight
      fit_pairwise_variance_model(
        pairwise_data_effects
      )
    }
  ),
  
  #### Sensitivity analysis with only baseline data ----

  tar_target(
    baseline_arm_data_effects,
    main_arm_data_effects |>
      filter(timepoint == "baseline")
  ),

  # Fitting baseline arm analysis models
  
  tar_target(
    baseline_arm_mean_effects_model,
    {
      model_specification_preflight
      fit_arm_mean_effects_model(
        baseline_arm_data_effects,
        prior_arm_mean_effects
      )
    }
  ),
  
  tar_target(
    baseline_arm_variance_effects_model,
    {
      model_specification_preflight
      fit_arm_variance_effects_model(
        baseline_arm_data_effects,
        prior_arm_variance_effects
      )
    }
  ),
  
  #### Sensitivity analysis including Greek lab studies ----
  
  
  # Fitting main arm analysis models including Greek studies
  
  tar_target(
    main_plus_greek_arm_mean_effects_model,
    {
      model_specification_preflight
      fit_arm_mean_effects_model(
        main_plus_greek_arm_data_effects,
        prior_arm_mean_effects
      )
    }
  ),
  
  tar_target(
    main_plus_greek_arm_variance_effects_model,
    {
      model_specification_preflight
      fit_arm_variance_effects_model(
        main_plus_greek_arm_data_effects,
        prior_arm_variance_effects
      )
    }
  ),
  
  #### Additional models including BMI and fat free mass
  
  # Add imputed demographics for moderators
  tar_target(
    main_arm_data_effects_imputed_demographics,
    impute_bmi_estimates(main_arm_data_effects) |>
      mutate(se_bmi = sd_bmi/sqrt(n_bmi),
             se_fat_free_mass = sd_fat_free_mass/sqrt(n_fat_free_mass))
  ),
  
  # Get control medians for predictors
  tar_target(
    predictor_medians,
    get_predictor_medians(main_arm_data_effects_imputed_demographics)
  ),

  # Validate all analysis datasets and brms formulas without compiling or sampling.
  tar_target(
    analysis_preflight,
    run_analysis_preflight(
      main_arm_data_effects,
      baseline_arm_data_effects,
      main_plus_greek_arm_data_effects,
      pairwise_data_effects,
      main_arm_data_effects_imputed_demographics
    )
  ),

  # Validate every registered informative prior against each formula and
  # analysis dataset without compiling or sampling a model.
  tar_target(
    model_specification_preflight,
    {
      analysis_preflight
      run_model_specification_preflight(
        main_arm_data_effects,
        baseline_arm_data_effects,
        main_plus_greek_arm_data_effects,
        main_arm_data_effects_imputed_demographics,
        prior_arm_mean_effects,
        prior_arm_variance_effects
      )
    }
  ),
  
  # BMI
  tar_target(
    main_arm_mean_effects_model_bmi,
    {
      model_specification_preflight
      fit_arm_mean_effects_model_moderator(
        main_arm_data_effects_imputed_demographics,
        prior_arm_mean_effects,
        "bmi"
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_model_bmi,
    {
      model_specification_preflight
      fit_arm_variance_effects_model_moderator(
        main_arm_data_effects_imputed_demographics,
        prior_arm_variance_effects,
        "bmi"
      )
    }
  ),
  
  tar_target(
    main_arm_mean_effects_contrast_condition_bmi,
    {
      main_arm_mean_effects_model_bmi_diagnostic_gate
      get_mean_contrast_condition_moderator(
        main_arm_mean_effects_model_bmi,
        predictor_medians
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_contrast_condition_bmi,
    {
      main_arm_variance_effects_model_bmi_diagnostic_gate
      get_variance_contrast_condition_moderator(
        main_arm_variance_effects_model_bmi,
        main_arm_data_effects_imputed_demographics,
        predictor_medians
      )
    }
  ),
  
  # Fat free mass
  
  tar_target(
    main_arm_mean_effects_model_fat_free_mass,
    {
      model_specification_preflight
      fit_arm_mean_effects_model_moderator(
        main_arm_data_effects_imputed_demographics,
        prior_arm_mean_effects,
        "fat_free_mass"
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_model_fat_free_mass,
    {
      model_specification_preflight
      fit_arm_variance_effects_model_moderator(
        main_arm_data_effects_imputed_demographics,
        prior_arm_variance_effects,
        "fat_free_mass"
      )
    }
  ),
  
  tar_target(
    main_arm_mean_effects_contrast_condition_fat_free_mass,
    {
      main_arm_mean_effects_model_fat_free_mass_diagnostic_gate
      get_mean_contrast_condition_moderator(
        main_arm_mean_effects_model_fat_free_mass,
        predictor_medians
      )
    }
  ),
  
  tar_target(
    main_arm_variance_effects_contrast_condition_fat_free_mass,
    {
      main_arm_variance_effects_model_fat_free_mass_diagnostic_gate
      get_variance_contrast_condition_moderator(
        main_arm_variance_effects_model_fat_free_mass,
        main_arm_data_effects_imputed_demographics,
        predictor_medians
      )
    }
  ),

  #### Model diagnostics ----

  tar_target(
    main_arm_mean_effects_diagnostics,
    model_diagnostic_summary(
      main_arm_mean_effects_model,
      "Primary arm mean model"
    )
  ),
  tar_target(
    main_arm_mean_effects_diagnostic_gate,
    assert_model_diagnostics(main_arm_mean_effects_diagnostics)
  ),
  tar_target(
    main_arm_mean_effects_trace_plot,
    model_trace_plot(main_arm_mean_effects_model, "Primary arm mean model")
  ),
  tar_target(
    main_arm_mean_effects_pp_check,
    model_pp_check(main_arm_mean_effects_model, "Primary arm mean model")
  ),

  tar_target(
    main_arm_variance_effects_diagnostics,
    model_diagnostic_summary(
      main_arm_variance_effects_model,
      "Primary arm variance model"
    )
  ),
  tar_target(
    main_arm_variance_effects_diagnostic_gate,
    assert_model_diagnostics(main_arm_variance_effects_diagnostics)
  ),
  tar_target(
    main_arm_variance_effects_trace_plot,
    model_trace_plot(main_arm_variance_effects_model, "Primary arm variance model")
  ),
  tar_target(
    main_arm_variance_effects_pp_check,
    model_pp_check(main_arm_variance_effects_model, "Primary arm variance model")
  ),

  tar_target(
    baseline_arm_mean_effects_diagnostics,
    model_diagnostic_summary(
      baseline_arm_mean_effects_model,
      "Baseline arm mean model"
    )
  ),
  tar_target(
    baseline_arm_mean_effects_diagnostic_gate,
    assert_model_diagnostics(baseline_arm_mean_effects_diagnostics)
  ),
  tar_target(
    baseline_arm_mean_effects_trace_plot,
    model_trace_plot(baseline_arm_mean_effects_model, "Baseline arm mean model")
  ),
  tar_target(
    baseline_arm_mean_effects_pp_check,
    model_pp_check(baseline_arm_mean_effects_model, "Baseline arm mean model")
  ),

  tar_target(
    baseline_arm_variance_effects_diagnostics,
    model_diagnostic_summary(
      baseline_arm_variance_effects_model,
      "Baseline arm variance model"
    )
  ),
  tar_target(
    baseline_arm_variance_effects_diagnostic_gate,
    assert_model_diagnostics(baseline_arm_variance_effects_diagnostics)
  ),
  tar_target(
    baseline_arm_variance_effects_trace_plot,
    model_trace_plot(baseline_arm_variance_effects_model, "Baseline arm variance model")
  ),
  tar_target(
    baseline_arm_variance_effects_pp_check,
    model_pp_check(baseline_arm_variance_effects_model, "Baseline arm variance model")
  ),

  tar_target(
    main_plus_greek_arm_mean_effects_diagnostics,
    model_diagnostic_summary(
      main_plus_greek_arm_mean_effects_model,
      "Greek-study arm mean model"
    )
  ),
  tar_target(
    main_plus_greek_arm_mean_effects_diagnostic_gate,
    assert_model_diagnostics(main_plus_greek_arm_mean_effects_diagnostics)
  ),
  tar_target(
    main_plus_greek_arm_mean_effects_trace_plot,
    model_trace_plot(
      main_plus_greek_arm_mean_effects_model,
      "Greek-study arm mean model"
    )
  ),
  tar_target(
    main_plus_greek_arm_mean_effects_pp_check,
    model_pp_check(
      main_plus_greek_arm_mean_effects_model,
      "Greek-study arm mean model"
    )
  ),

  tar_target(
    main_plus_greek_arm_variance_effects_diagnostics,
    model_diagnostic_summary(
      main_plus_greek_arm_variance_effects_model,
      "Greek-study arm variance model"
    )
  ),
  tar_target(
    main_plus_greek_arm_variance_effects_diagnostic_gate,
    assert_model_diagnostics(main_plus_greek_arm_variance_effects_diagnostics)
  ),
  tar_target(
    main_plus_greek_arm_variance_effects_trace_plot,
    model_trace_plot(
      main_plus_greek_arm_variance_effects_model,
      "Greek-study arm variance model"
    )
  ),
  tar_target(
    main_plus_greek_arm_variance_effects_pp_check,
    model_pp_check(
      main_plus_greek_arm_variance_effects_model,
      "Greek-study arm variance model"
    )
  ),

  tar_target(
    pairwise_mean_effects_diagnostics,
    model_diagnostic_summary(pairwise_mean_effects_model, "Pairwise mean model")
  ),
  tar_target(
    pairwise_mean_effects_diagnostic_gate,
    assert_model_diagnostics(pairwise_mean_effects_diagnostics)
  ),
  tar_target(
    pairwise_mean_effects_trace_plot,
    model_trace_plot(pairwise_mean_effects_model, "Pairwise mean model")
  ),
  tar_target(
    pairwise_mean_effects_pp_check,
    model_pp_check(pairwise_mean_effects_model, "Pairwise mean model")
  ),

  tar_target(
    pairwise_variance_effects_diagnostics,
    model_diagnostic_summary(
      pairwise_variance_effects_model,
      "Pairwise variance model"
    )
  ),
  tar_target(
    pairwise_variance_effects_diagnostic_gate,
    assert_model_diagnostics(pairwise_variance_effects_diagnostics)
  ),
  tar_target(
    pairwise_variance_effects_trace_plot,
    model_trace_plot(pairwise_variance_effects_model, "Pairwise variance model")
  ),
  tar_target(
    pairwise_variance_effects_pp_check,
    model_pp_check(pairwise_variance_effects_model, "Pairwise variance model")
  ),

  tar_target(
    main_arm_mean_effects_model_bmi_diagnostics,
    model_diagnostic_summary(
      main_arm_mean_effects_model_bmi,
      "BMI-adjusted arm mean model"
    )
  ),
  tar_target(
    main_arm_mean_effects_model_bmi_diagnostic_gate,
    assert_model_diagnostics(main_arm_mean_effects_model_bmi_diagnostics)
  ),
  tar_target(
    main_arm_mean_effects_model_bmi_trace_plot,
    model_trace_plot(main_arm_mean_effects_model_bmi, "BMI-adjusted arm mean model")
  ),
  tar_target(
    main_arm_mean_effects_model_bmi_pp_check,
    model_pp_check(main_arm_mean_effects_model_bmi, "BMI-adjusted arm mean model")
  ),

  tar_target(
    main_arm_variance_effects_model_bmi_diagnostics,
    model_diagnostic_summary(
      main_arm_variance_effects_model_bmi,
      "BMI-adjusted arm variance model"
    )
  ),
  tar_target(
    main_arm_variance_effects_model_bmi_diagnostic_gate,
    assert_model_diagnostics(main_arm_variance_effects_model_bmi_diagnostics)
  ),
  tar_target(
    main_arm_variance_effects_model_bmi_trace_plot,
    model_trace_plot(
      main_arm_variance_effects_model_bmi,
      "BMI-adjusted arm variance model"
    )
  ),
  tar_target(
    main_arm_variance_effects_model_bmi_pp_check,
    model_pp_check(
      main_arm_variance_effects_model_bmi,
      "BMI-adjusted arm variance model"
    )
  ),

  tar_target(
    main_arm_mean_effects_model_fat_free_mass_diagnostics,
    model_diagnostic_summary(
      main_arm_mean_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm mean model"
    )
  ),
  tar_target(
    main_arm_mean_effects_model_fat_free_mass_diagnostic_gate,
    assert_model_diagnostics(
      main_arm_mean_effects_model_fat_free_mass_diagnostics
    )
  ),
  tar_target(
    main_arm_mean_effects_model_fat_free_mass_trace_plot,
    model_trace_plot(
      main_arm_mean_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm mean model"
    )
  ),
  tar_target(
    main_arm_mean_effects_model_fat_free_mass_pp_check,
    model_pp_check(
      main_arm_mean_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm mean model"
    )
  ),

  tar_target(
    main_arm_variance_effects_model_fat_free_mass_diagnostics,
    model_diagnostic_summary(
      main_arm_variance_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm variance model"
    )
  ),
  tar_target(
    main_arm_variance_effects_model_fat_free_mass_diagnostic_gate,
    assert_model_diagnostics(
      main_arm_variance_effects_model_fat_free_mass_diagnostics
    )
  ),
  tar_target(
    main_arm_variance_effects_model_fat_free_mass_trace_plot,
    model_trace_plot(
      main_arm_variance_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm variance model"
    )
  ),
  tar_target(
    main_arm_variance_effects_model_fat_free_mass_pp_check,
    model_pp_check(
      main_arm_variance_effects_model_fat_free_mass,
      "Fat-free-mass-adjusted arm variance model"
    )
  ),

  tar_target(
    all_model_diagnostics,
    dplyr::bind_rows(
      main_arm_mean_effects_diagnostics,
      main_arm_variance_effects_diagnostics,
      baseline_arm_mean_effects_diagnostics,
      baseline_arm_variance_effects_diagnostics,
      main_plus_greek_arm_mean_effects_diagnostics,
      main_plus_greek_arm_variance_effects_diagnostics,
      pairwise_mean_effects_diagnostics,
      pairwise_variance_effects_diagnostics,
      main_arm_mean_effects_model_bmi_diagnostics,
      main_arm_variance_effects_model_bmi_diagnostics,
      main_arm_mean_effects_model_fat_free_mass_diagnostics,
      main_arm_variance_effects_model_fat_free_mass_diagnostics
    )
  ),
  tar_target(
    all_model_diagnostics_gate,
    assert_model_diagnostics(all_model_diagnostics)
  ),

  #### Manuscripts ----

  tarchetypes::tar_quarto(
    pre_print_pdf,
    path = file.path("manuscript", "pre_print.qmd"),
    quiet = FALSE
  ),

  tarchetypes::tar_quarto(
    submission_manuscript_pdf,
    path = file.path("manuscript", "submission_manuscript.qmd"),
    quiet = FALSE
  )
  
)















