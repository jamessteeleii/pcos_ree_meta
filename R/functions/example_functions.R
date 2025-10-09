# Functions to create initial example data using studies identified in https://macrofactorapp.com/pcos-bmr/ 
# This data will be used for setting up the analysis pipelines and will be replaced by the final extracted data after systematic review
  create_pairwise_data_example <- function() {
    
    # First we create a tibble with the data from the pairwise studies of both PCOS and controls (CONT)
    data <- tibble(
      study = c("Segal, 1990",
                "Robinson, 1992",
                "Cosar, 2008",
                "Georgopoulos, 2009",
                "Graff, 2013",
                "Larsson, 2015",
                "Doh, 2016"),
      n_PCOS = c(10,14,31,91,61,72,14),
      n_CONT = c(11,14,29,48,44,30,10),
      m_PCOS = c(1508,1624,1167,1446,1469,1411,1272),
      m_CONT = c(1558,1633,1046,1841,1453,1325,1240),
      s_PCOS = c(168,137,371,725,227,229,167),
      s_CONT = c(186,215,296,305,249,193,216)
    )
    return(data)
  }
  
  calculate_pairwise_effects_example <- function(data) {
    data <- escalc(measure = "MD",
                  m1i = m_PCOS,
                  m2i = m_CONT,
                  sd1i = s_PCOS,
                  sd2i = s_CONT,
                  n1i = n_PCOS,
                  n2i = n_CONT,
                  data = data,
                  var.names = c("yi_mean", "vi_mean"))
    
    data <- escalc(measure = "CVR",
                  m1i = m_PCOS,
                  m2i = m_CONT,
                  sd1i = s_PCOS,
                  sd2i = s_CONT,
                  n1i = n_PCOS,
                  n2i = n_CONT,
                  data = data,
                  var.names = c("yi_cvr", "vi_cvr"))
    
    return(data)
  }
  

  create_arm_data_example <- function(data) {
    
    # First we pivot the pairwise data to long form (and add the se column)
    data2 <- data |>
      pivot_longer(
        2:7,
        cols_vary = "slowest",
        names_to = c(".value", "cond"),
        names_pattern = "(..)(....)"
      ) |>
      rename(
        n = "n_",
        m = "m_",
        sd = "s_"
      ) |>
      mutate(
        se = NA
      )
    
    # We then bind this with the data for the studies with only PCOS arms
    data_arm <- bind_rows(
      data2, tibble(
        study = c("Bruner, 2006", "Bruner, 2006", "Kritikou, 2006", "Moran, 2006", "Saltamavros, 2007", "Koika, 2009", "Broksey,2017", "Rodrigues, 2017"),
        cond = "PCOS",
        n = c(7,5,63,13,73,156,28,30),
        
        # some of the values are in MJ so converted to kcal multiplying by 239.01
        m = c(1485,1596,1505.6,7.7*239.01,1381.5,1415.7,1689,1677), 
        sd = c(NA,NA,NA,NA,NA,672.9,373,210),
        se = c(177,107,83.2,0.3*239.01,167.9,NA,NA,NA)
      )
    ) |>
      mutate(
        # convert the se to sd
        sd = if_else(is.na(sd), se * sqrt(n), sd),
      )|>
      rowid_to_column("arm")
    
    return(data_arm)
  }
  
  calculate_arm_effects_example <- function(data) {
    data <- escalc(measure = "MN",
                      mi = m,
                      sdi = sd,
                      ni = n,
                      data = data,
                  var.names = c("yi_mean", "vi_mean"))
    
    data <- escalc(measure = "SDLN",
                  mi = m,
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
  

  
# Functions for example arm models (note, main analysis we decided to also include lab as a random effect)
  
  fit_arm_mean_effects_model_example <- function(data, prior) {
    
    arm_model <- brm(yi_mean | se(sqrt(vi_mean)) ~ 0 + Intercept + cond + (1 + cond | study) + (1 | arm),
                     data = data,
                     prior = prior,
                     chains = 4,
                     cores = 4,
                     seed = 1988,
                     warmup = 2000,
                     iter = 8000)
    
    return(arm_model)
  }
  
  fit_arm_variance_effects_model_example <- function(data, prior) {
    
    
    arm_model <- brm(yi_sd | se(sqrt(vi_sd)) ~ 0 + Intercept + cond + me(log_yi_mean, se_log_yi_mean) + (1 + cond | study) + (1 | arm),
                     data = data,
                     prior = prior,
                     chains = 4,
                     cores = 4,
                     seed = 1988,
                     warmup = 2000,
                     iter = 8000)
    
    return(arm_model)
  }