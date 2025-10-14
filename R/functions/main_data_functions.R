# Functions for main data preparation
prepare_data <- function(file) {
  
  
  data <- read_csv(file) |>
    
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
    ungroup() |>
    rowid_to_column("effect")
  
  return(data)
}

calculate_arm_effects <- function(data) {
  data <- escalc(measure = "MN",
                 mi = mean,
                 sdi = sd,
                 ni = n,
                 data = data,
                 var.names = c("yi_mean", "vi_mean"))

  data <- escalc(measure = "SDLN",
                 mi = mean,
                 sdi = sd,
                 ni = n,
                 data = data,
                 var.names = c("yi_sd", "vi_sd"))
  
  data <- data |>
    mutate(
      yi_mean_floor = pmax(yi_mean, 1e-6),
      se_log_yi_mean = sqrt(vi_mean) / yi_mean_floor,
      log_yi_mean = log(yi_mean_floor)
    )

  return(data)
}

impute_demographic_estimates <- function(data) {
  data <- data |>
    
    # convert height to cm for reporting
    mutate(m_height = case_when(
      m_height < 100 ~ m_height*100
      ),
      sd_height = case_when(
        m_height < 100 ~ sd_height*100
      )
    ) |>
    
    # impute estimates of missing weight/height/bmi and note where estimate
    mutate(
      m_body_mass_estim = case_when(
        is.na(m_body_mass) & !is.na(m_height) & !is.na(m_bmi) ~ "y",
        .default = "n"
      ),
      m_height_estim = case_when(
        is.na(m_height) & !is.na(m_body_mass) & !is.na(m_bmi) ~ "y",
        .default = "n"
      ),
      m_bmi_estim = case_when(
        is.na(m_bmi) & !is.na(m_body_mass) & !is.na(m_height) ~ "y",
        .default = "n"
      ),
    ) |>
    
    mutate(
      m_body_mass = case_when(
        is.na(m_body_mass) & !is.na(m_height) & !is.na(m_bmi) ~ m_bmi * (m_height/100)^2,
        .default = m_body_mass
      ),
      m_height = case_when(
        is.na(m_height) & !is.na(m_body_mass) & !is.na(m_bmi) ~ (sqrt(m_body_mass/m_bmi))*100,
        .default = m_height
      ),
      m_bmi = case_when(
        is.na(m_bmi) & !is.na(m_body_mass) & !is.na(m_height) ~ m_body_mass/(m_height/100)^2,
        .default = m_bmi
      ),
    )
  
  return(data)
    
}

# Pairwise data preparation for sensitivity analysis

prepare_pairwise_data <- function(data) {
  
  wide_data <- data |> 
    filter(timepoint == "baseline") |> # just take baseline
    select(lab, study, authors, year, cond, n, mean, sd) |> 
    pivot_wider(names_from = "cond", values_from = c("n","mean", "sd")) |>
    
    # just studies with PCOS and Control (no further subgroups)
    filter(map_lgl(n_PCOS, ~ length(.x) == 1) &
             map_lgl(n_Control, ~ length(.x) == 1)) |>
    mutate(across(n_PCOS:sd_Control, as.numeric))
  
  return(wide_data)
}

calculate_pairwise_effects <- function(data) {
  
  data <- escalc(measure = "MD",
                 m1i = mean_PCOS,
                 m2i = mean_Control,
                 sd1i = sd_PCOS,
                 sd2i = sd_Control,
                 n1i = n_PCOS,
                 n2i = n_Control,
                 data = data,
                 var.names = c("yi_mean", "vi_mean"))
  
  data <- escalc(measure = "CVR",
                 m1i = mean_PCOS,
                 m2i = mean_Control,
                 sd1i = sd_PCOS,
                 sd2i = sd_Control,
                 n1i = n_PCOS,
                 n2i = n_Control,
                 data = data,
                 var.names = c("yi_cvr", "vi_cvr"))
  
  return(data)
  
}

