finite_diagnostic_values <- function(x) {
  x[is.finite(x)]
}

model_diagnostic_summary <- function(
    model,
    label,
    rhat_threshold = 1.01,
    ess_threshold = 400,
    ebfmi_threshold = 0.30) {
  if (!inherits(model, "brmsfit")) {
    stop(label, " is not a brmsfit object.", call. = FALSE)
  }
  if (!inherits(model$fit, "stanfit")) {
    stop(label, " was not fitted with the rstan backend.", call. = FALSE)
  }

  draws <- posterior::as_draws_array(model)
  draw_diagnostics <- posterior::summarise_draws(
    draws,
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  rhat_values <- finite_diagnostic_values(draw_diagnostics$rhat)
  ess_bulk_values <- finite_diagnostic_values(draw_diagnostics$ess_bulk)
  ess_tail_values <- finite_diagnostic_values(draw_diagnostics$ess_tail)

  sampler_parameters <- rstan::get_sampler_params(
    model$fit,
    inc_warmup = FALSE
  )
  stan_control <- model$fit@stan_args[[1L]]$control
  max_treedepth <- if (is.null(stan_control$max_treedepth)) {
    10
  } else {
    stan_control$max_treedepth
  }
  adapt_delta <- if (is.null(stan_control$adapt_delta)) {
    0.8
  } else {
    stan_control$adapt_delta
  }

  divergences_by_chain <- vapply(
    sampler_parameters,
    function(chain) sum(chain[, "divergent__"]),
    numeric(1)
  )
  treedepth_hits_by_chain <- vapply(
    sampler_parameters,
    function(chain) sum(chain[, "treedepth__"] >= max_treedepth),
    numeric(1)
  )
  ebfmi_by_chain <- vapply(
    sampler_parameters,
    function(chain) {
      energy <- chain[, "energy__"]
      mean(diff(energy)^2) / stats::var(energy)
    },
    numeric(1)
  )
  observation_count <- nrow(model$data)

  diagnostics <- tibble::tibble(
    model = label,
    observations = observation_count,
    chains = length(sampler_parameters),
    post_warmup_draws = sum(vapply(sampler_parameters, nrow, integer(1))),
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    divergences = sum(divergences_by_chain),
    max_treedepth_hits = sum(treedepth_hits_by_chain),
    max_rhat = if (length(rhat_values)) max(rhat_values) else NA_real_,
    min_ess_bulk = if (length(ess_bulk_values)) min(ess_bulk_values) else NA_real_,
    min_ess_tail = if (length(ess_tail_values)) min(ess_tail_values) else NA_real_,
    min_ebfmi = if (all(is.finite(ebfmi_by_chain))) min(ebfmi_by_chain) else NA_real_,
    divergences_by_chain = paste(divergences_by_chain, collapse = ","),
    max_treedepth_hits_by_chain = paste(treedepth_hits_by_chain, collapse = ","),
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold,
    ebfmi_threshold = ebfmi_threshold
  )

  diagnostics |>
    dplyr::mutate(
      diagnostic_ok =
        divergences == 0 &&
        max_treedepth_hits == 0 &&
        is.finite(max_rhat) && max_rhat <= rhat_threshold &&
        is.finite(min_ess_bulk) && min_ess_bulk >= ess_threshold &&
        is.finite(min_ess_tail) && min_ess_tail >= ess_threshold &&
        is.finite(min_ebfmi) && min_ebfmi >= ebfmi_threshold,
      issues = paste(
        c(
          if (divergences > 0) paste(divergences, "divergent transitions") else NULL,
          if (max_treedepth_hits > 0) paste(max_treedepth_hits, "maximum treedepth hits") else NULL,
          if (!is.finite(max_rhat) || max_rhat > rhat_threshold) {
            paste0("maximum R-hat ", signif(max_rhat, 4))
          } else NULL,
          if (!is.finite(min_ess_bulk) || min_ess_bulk < ess_threshold) {
            paste0("minimum bulk ESS ", round(min_ess_bulk))
          } else NULL,
          if (!is.finite(min_ess_tail) || min_ess_tail < ess_threshold) {
            paste0("minimum tail ESS ", round(min_ess_tail))
          } else NULL,
          if (!is.finite(min_ebfmi) || min_ebfmi < ebfmi_threshold) {
            paste0("minimum E-BFMI ", signif(min_ebfmi, 4))
          } else NULL
        ),
        collapse = "; "
      ),
      issues = dplyr::if_else(diagnostic_ok, "none", issues)
    )
}

assert_model_diagnostics <- function(diagnostics) {
  required_columns <- c("model", "diagnostic_ok", "issues")
  missing_columns <- setdiff(required_columns, names(diagnostics))
  if (length(missing_columns)) {
    stop(
      "Model diagnostics are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  failed <- diagnostics[!diagnostics$diagnostic_ok, , drop = FALSE]
  if (nrow(failed)) {
    details <- paste0(failed$model, ": ", failed$issues)
    stop(
      "Model diagnostics failed:\n- ",
      paste(details, collapse = "\n- "),
      call. = FALSE
    )
  }

  TRUE
}

model_trace_plot <- function(model, label) {
  draws <- posterior::as_draws_array(model)
  parameters <- posterior::variables(draws)
  parameters <- grep(
    "^(b_|sd_|cor_|meanme_|sdme_)",
    parameters,
    value = TRUE
  )
  if (!length(parameters)) {
    stop(label, " has no diagnostic parameters available for a trace plot.", call. = FALSE)
  }

  bayesplot::mcmc_trace(
    draws,
    pars = parameters,
    facet_args = list(ncol = 1, strip.position = "left")
  ) +
    ggplot2::labs(title = paste(label, "trace plots"))
}

model_pp_check <- function(model, label, ndraws = 100L) {
  brms::pp_check(
    model,
    type = "dens_overlay",
    ndraws = ndraws
  ) +
    ggplot2::labs(title = paste(label, "posterior predictive check"))
}
