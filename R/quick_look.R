library(tidyverse)
library(metafor)

dat <- tibble(
  study = c("Segal, 1990",
            "Robinson, 1992",
            "Cosar, 2008",
            "Georgopoulos, 2009",
            "Graff, 2013",
            "Larsson, 2015",
            "Doh, 2016"),
  n_PCOS = c(10,14,31,91,61,72,14),
  n_CON = c(11,14,29,48,44,30,10),
  m_PCOS = c(1508,1624,1167,1446,1469,1411,1272),
  m_CON = c(1558,1633,1046,1841,1453,1325,1240),
  s_PCOS = c(168,137,371,725,227,229,167),
  s_CON = c(186,215,296,305,249,193,216)
)

dat <- escalc(measure = "MD",
              m1i = m_PCOS,
              m2i = m_CON,
              sd1i = s_PCOS,
              sd2i = s_CON,
              n1i = n_PCOS,
              n2i = n_CON,
              data = dat)

paired_model <- rma(yi, vi, data = dat)

forest(paired_model, slab = study)


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
    names_pattern = "(..)(...)"
  ) |>
  rename(
    n = "n_",
    m = "m_",
    sd = "s_"
    ) |>
  mutate(
    cond = case_when(
      cond == "PCO" ~ "PCOS",
      .default = cond
    ),
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
    effect = 1:22
  )

dat_arm <- escalc(measure = "MN",
               mi = m,
               sdi = sd,
               ni = n,
               data = dat_arm)


arm_model <- rma.mv(yi, vi, 
                    mods = ~ cond, 
                    random = list(~ cond | study, ~ 1 | effect),
                    method="REML", test="t",
                    data = dat_arm)

rob_arm_model <- robust(arm_model, dat_arm$study)

arm_model_preds <- cbind(dat_arm, pred = predict(rob_arm_model)$pred,
                       ci.lb =  predict(rob_arm_model)$ci.lb,
                       ci.ub =  predict(rob_arm_model)$ci.ub) %>%
  mutate(wi = 1/sqrt(vi),
         size = 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))) 

ranef(rob_arm_model)
