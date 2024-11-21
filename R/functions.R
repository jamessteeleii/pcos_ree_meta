# Functions to create initial example data using studies identified in https://macrofactorapp.com/pcos-bmr/ 
# This data will be used for setting up the analysis pipelines and will be replaced by the final extracted data after systematic review
  create_pairwise_data <- function() {
    
    # First we create a tibble with the data from the pairwise studies of both PCOS and controls (CONT)
    dat <- tibble(
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
    return(dat)
  }
  
  calculate_pairwise_mean_effects <- function(dat) {
    dat <- escalc(measure = "MD",
                  m1i = m_PCOS,
                  m2i = m_CONT,
                  sd1i = s_PCOS,
                  sd2i = s_CONT,
                  n1i = n_PCOS,
                  n2i = n_CONT,
                  data = dat)
    
    return(dat)
  }
  
  calculate_pairwise_variance_effects <- function(dat) {
    dat <- escalc(measure = "CVR",
                  m1i = m_PCOS,
                  m2i = m_CONT,
                  sd1i = s_PCOS,
                  sd2i = s_CONT,
                  n1i = n_PCOS,
                  n2i = n_CONT,
                  data = dat)
    
    return(dat)
  }
  
  create_arm_data <- function(dat) {
    
    # First we pivot the pairwise data to long form (and add the se column)
    dat2 <- dat |>
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
    dat_arm <- bind_rows(
      dat2, tibble(
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
    
    return(dat_arm)
  }
  
  calculate_arm_mean_effects <- function(dat) {
    dat <- escalc(measure = "MN",
                      mi = m,
                      sdi = sd,
                      ni = n,
                      data = dat)
    
    return(dat)
  }
  
  calculate_arm_variance_effects <- function(dat) {
    dat_var <- escalc(measure = "SDLN",
                  mi = m,
                  sdi = sd,
                  ni = n,
                  data = dat)
    
    # to add the mean and it's error for model
    dat_mean <- escalc(measure = "MN",
                       mi = m,
                       sdi = sd,
                       ni = n,
                       data = dat) |>
      rename(m_yi = "yi",
             m_vi = "vi")
    
    dat <- bind_cols(dat_var,
                     dat_mean |>
                       select(m_yi, m_vi))
    
    return(dat)
  }

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
    
    dat_prior <- tibble(
      study = c("Mifflin, 1990", "Pavlidou, 2022"),
      n = c(247,549),
      m = c(1349, 1533),
      sd = c(214, 308)
    )
    
    dat_prior_mean_effects <- escalc(measure = "MN",
                                     mi = m,
                                     sdi = sd,
                                     ni = n,
                                     data = dat_prior)
    
    # We fit a model to estimate with default weakly regularising priors
    estimate_prior_mean_effects <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                                       data = dat_prior_mean_effects,
                                       chains = 4,
                                       cores = 4,
                                       seed = 1988,
                                       warmup = 2000,
                                       control = list(adapt_delta = 0.99),
                                       iter = 8000)
    
    
    estimate_prior_mean_effects <- broom.mixed::tidy(estimate_prior_mean_effects)
    
    # The we set the priors taking the estimates from the models of 
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
  
  fit_arm_mean_effects_model <- function(dat, prior) {
    
    arm_model <- brm(yi | se(sqrt(vi)) ~ 0 + Intercept + cond + (1 + cond | study) + (1 | arm),
                     data = dat,
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
    
    dat_prior <- tibble(
      study = c("Mifflin, 1990", "Pavlidou, 2022"),
      n = c(247,549),
      m = c(1349, 1533),
      sd = c(214, 308)
    )
    
    dat_prior_variance_effects <- escalc(measure = "SDLN",
                                         mi = m,
                                         sdi = sd,
                                         ni = n,
                                         data = dat_prior)
    
    # We fit a model to estimate with default weakly regularising priors
    estimate_prior_variance_effects <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                                           data = dat_prior_variance_effects,
                                           chains = 4,
                                           cores = 4,
                                           seed = 1988,
                                           warmup = 2000,
                                           control = list(adapt_delta = 0.99),
                                           iter = 8000)
    
    
    estimate_prior_variance_effects <- broom.mixed::tidy(estimate_prior_variance_effects)
    
    # Mean effects so we have a prior for the mean of the means as a covariate in the model (and its se)
    dat_prior_mean_effects <- escalc(measure = "MN",
                                     mi = m,
                                     sdi = sd,
                                     ni = n,
                                     data = dat_prior)
    
    # We fit a model to estimate with default weakly regularising priors
    estimate_prior_mean_effects <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                                       data = dat_prior_mean_effects,
                                       chains = 4,
                                       cores = 4,
                                       seed = 1988,
                                       warmup = 2000,
                                       control = list(adapt_delta = 0.99),
                                       iter = 8000)
    
    
    estimate_prior_mean_effects <- broom.mixed::tidy(estimate_prior_mean_effects)
    
    # The we set the priors taking the estimates from the models of 
    prior_arm_mean_effects <-
      c(
        # The prior on the intercept i.e., CONT arms is set from the estimate of the two studies means mentioned above
        set_prior(paste("student_t(3,", estimate_prior_variance_effects$estimate[1],",", estimate_prior_variance_effects$std.error[1],")"),
                  class = "b", coef = "Intercept"),
        # The prior on the random effects for intercept i.e., CONT arms is set from the estimate of the two studies variance mentioned above
        set_prior(paste("student_t(3,", estimate_prior_variance_effects$estimate[2],",", estimate_prior_variance_effects$std.error[2],")"),
                  class = "sd", coef = "Intercept", group = "study"),
        
        # The fixed effect coef for log(mean) is typically ~1 due to mean-variance relationship being commonplace in other measures 
        # But we set it to be centred there though with a wide scale to indicate uncertainty in this outcome specifically
        set_prior("student_t(3, 0, 2.5)", class = "b", coef = "melogm_yisqrtlogm_vi"),
        
        # The fixed effect coef reflecting the difference between CONT and PCOS is set based on a wide range of possible values
        # Given the rough relationship of ~1 for log(mean) on log(sd) in other data we again set it to reflect the range of diffs on the log scale
        # This uses the min and max values of ranges reported in the two studies i.e., 2492 - 908 = 1584
        # We then set a student t prior that permits values approximately up to this value with the majority of it's mass centred around zero
        set_prior("student_t(3, 0, 5.3)", class = "b", coef = "condPCOS"),
        
        # The mean of the log(mean) measurement error has to be positive (as means are positive), as does the sd, so we set these to wide half t distributions
        set_prior(paste("student_t(3,", log(estimate_prior_mean_effects$estimate[1]),",", log(estimate_prior_mean_effects$std.error[1]),")"),
                  class = "meanme", coef = "melogm_yi"),
        set_prior("student_t(3, 0, 5)", class = "sdme", coef = "melogm_yi")
        
      )
    
    return(prior_arm_variance_effects)
    
  }
  
  fit_arm_variance_effects_model <- function(dat, prior) {
    
    arm_model <- brm(yi | se(sqrt(vi)) ~ 0 + Intercept + cond + me(log(m_yi), sqrt(log(m_vi))) + (1 + cond | study) + (1 | arm),
                     data = dat,
                     prior = prior,
                     chains = 4,
                     cores = 4,
                     seed = 1988,
                     warmup = 2000,
                     iter = 8000)
  }
  
  # Pairwise models
  ## Note, we just utilise the weakly regularising default priors from brms for both pairwise models
  fit_pairwise_model <- function(dat) {
    
    pairwise_model <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                          data = dat,
                          chains = 4,
                          cores = 4,
                          seed = 1988,
                          warmup = 2000,
                          control = list(adapt_delta = 0.99),
                          iter = 8000)
    
    return(pairwise_model)
  }
  
  
# Functions for extracting results and generating plots for reporting - TO BE ADDED
  