library(tidyverse)
library(metafor)
library(brms)

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

dat <- escalc(measure = "MD",
              m1i = m_PCOS,
              m2i = m_CONT,
              sd1i = s_PCOS,
              sd2i = s_CONT,
              n1i = n_PCOS,
              n2i = n_CONT,
              data = dat)

# fit brms model

paired_model <- brm(yi | se(sqrt(vi)) ~ 1 + (1 | study),
                    data = dat,
                    chains = 4,
                    cores = 4,
                    seed = 1988,
                    warmup = 2000,
                    iter = 8000)

plot(paired_model)
pp_check(paired_model)


# dat_var <- escalc(measure = "CVR",
#               m1i = m_PCOS,
#               m2i = m_CON,
#               sd1i = sd_PCOS,
#               sd2i = sd_CON,
#               n1i = n_PCOS,
#               n2i = n_CON,
#               data = dat)
# 
# forest(rma(yi, vi, data = dat_var), slab = study)


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
  
dat_arm <- bind_rows(
  dat2, tibble(
    study = c("Bruner, 2006", "Bruner, 2006", "Kritikou, 2006", "Moran, 2006", "Saltamavros, 2007", "Koika, 2009", "Broksey,2017", "Rodrigues, 2017"),
    cond = "PCOS",
    n = c(7,5,63,13,73,156,28,30),
    m = c(1485,1596,1505.6,7.7*239.01,1381.5,1415.7,1689,1677),
    sd = c(NA,NA,NA,NA,NA,672.9,373,210),
    se = c(177,107,83.2,0.3*239.01,167.9,NA,NA,NA)
  )
) |>
  mutate(
    sd = if_else(is.na(sd), se * sqrt(n), sd),
    )|>
  rowid_to_column("arm")

dat_arm <- escalc(measure = "MN",
               mi = m,
               sdi = sd,
               ni = n,
               data = dat_arm)

get_prior(yi | se(sqrt(vi)) ~ 0 + Intercept + cond + me(log(m_yi), sqrt(log(m_vi))) + (1 + cond | study) + (1 | arm),
                 data = dat_arm)

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

arm_model <- brm(yi | se(sqrt(vi)) ~ 0 + Intercept + cond + me(log(m_yi), sqrt(log(m_vi))) + (1 + cond | study) + (1 | arm),
                 data = dat_arm,
                 prior = prior_arm_mean_effects,
                 # sample_prior = "only",
                 chains = 4,
                 cores = 4,
                 seed = 1988,
                 warmup = 2000,
                 save_pars = save_pars(all = TRUE, latent = TRUE),
                 iter = 8000)

plot(arm_model)
pp_check(arm_model, type = "dens_overlay_grouped", group = "cond", ndraws = 100)
pp_check(arm_model, type = "intervals", x = "m_yi", ndraws = 100)

pp_check(arm_model, ndraws = 100)
