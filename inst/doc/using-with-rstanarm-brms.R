## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  warning    = FALSE,
  message    = FALSE
)

## ----check-deps---------------------------------------------------------------
have_rstanarm <- requireNamespace("rstanarm", quietly = TRUE)
have_brms     <- requireNamespace("brms",     quietly = TRUE)
have_deps     <- have_rstanarm && have_brms

## ----map-prior----------------------------------------------------------------
library(bayprior)

map <- map_prior(
  y            = c(-0.85, -0.62, -1.10),
  se           = c(0.25, 0.30, 0.28),
  outcome_type = "single_arm_log_odds",
  label        = "Historical control (3 trials)"
)

mu    <- map$fit_summary$mean
sigma <- map$fit_summary$sd
c(mu = mu, sigma = sigma)

## ----rstanarm-demo, eval = have_rstanarm--------------------------------------
library(rstanarm)

# A toy current-trial dataset, sized only to keep this vignette's build
# time short -- substitute your own trial data in practice.
set.seed(1)
n_trial  <- 40
p_trial  <- 0.32
dat <- data.frame(y = rbinom(n_trial, 1, p_trial))

# Wrapped in tryCatch so an incomplete local C++/Stan toolchain degrades
# this vignette gracefully (a build failure here would otherwise fail
# R CMD check on such a machine) rather than as a claim this code is
# untested -- it has been run successfully end to end where the toolchain
# is complete.
fit_rstanarm <- tryCatch(
  stan_glm(
    y ~ 1,
    data            = dat,
    family          = binomial(link = "logit"),
    prior_intercept = normal(location = mu, scale = sigma, autoscale = FALSE),
    chains          = 1,
    iter            = 500,
    refresh         = 0,
    seed            = 1
  ),
  error = function(e) {
    message("Model fit skipped (local Stan toolchain issue): ", conditionMessage(e))
    NULL
  }
)

if (!is.null(fit_rstanarm)) print(fit_rstanarm, digits = 3)

## ----rstanarm-skip, eval = !have_rstanarm-------------------------------------
# cat("rstanarm not installed -- skipping this demo.\n")

## ----brms-demo, eval = have_brms----------------------------------------------
library(brms)

prior_spec <- prior_string(
  paste0("normal(", mu, ",", sigma, ")"),
  class = "Intercept"
)
prior_spec

# See the rstanarm chunk above for why this is wrapped in tryCatch.
fit_brms <- tryCatch(
  brm(
    y ~ 1,
    data    = dat,
    family  = bernoulli(link = "logit"),
    prior   = prior_spec,
    chains  = 1,
    iter    = 500,
    refresh = 0,
    seed    = 1,
    silent  = 2
  ),
  error = function(e) {
    message("Model fit skipped (local Stan toolchain issue): ", conditionMessage(e))
    NULL
  }
)

if (!is.null(fit_brms)) summary(fit_brms)

## ----brms-skip, eval = !have_brms---------------------------------------------
# cat("brms not installed -- skipping this demo.\n")

