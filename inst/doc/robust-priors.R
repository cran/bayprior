## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  fig.width  = 6,
  fig.height = 3.8,
  out.width  = "95%",
  dpi        = 120,
  warning    = FALSE,
  message    = FALSE
)
library(bayprior)

## ----robust-------------------------------------------------------------------
informative <- elicit_beta(
  mean      = 0.30,
  sd        = 0.08,
  method    = "moments",
  label     = "Response rate"
)

rob <- robust_prior(
  informative  = informative,
  vague_weight = 0.20,
  label        = "Robust mixture prior"
)
plot(rob)

## ----robust-summary-----------------------------------------------------------
cat("Informative component weight:", 1 - rob$vague_weight, "\n")
cat("Vague component weight:      ", rob$vague_weight, "\n")
cat("Mixture mean:", round(rob$fit_summary$mean, 4), "\n")
cat("Mixture SD:  ", round(rob$fit_summary$sd,   4), "\n")

## ----robust-weights, fig.height = 4-------------------------------------------
oldpar <- par(mfrow = c(1, 1))
weights <- c(0.10, 0.20, 0.30, 0.50)
cols    <- c("#185FA5", "#1D9E75", "#D85A30", "#888780")

x <- seq(0, 0.8, length.out = 300)
plot(x, bayprior:::.eval_density_vec(informative, x),
     type = "l", lwd = 2, col = "#185FA5",
     xlab = "Response rate", ylab = "Density",
     main = "Effect of vague weight on robust prior",
     ylim = c(0, 6))
for (i in seq_along(weights)) {
  r <- robust_prior(informative, vague_weight = weights[i])
  lines(x, bayprior:::.eval_density_vec(r, x),
        col = cols[i], lwd = 1.5, lty = i + 1)
}
legend("topright",
       legend = c("Informative",
                  paste0("w = ", weights)),
       col    = c("#185FA5", cols),
       lwd    = 2,
       lty    = c(1, 2, 3, 4, 5),
       bty    = "n", cex = 0.8)
par(oldpar)       

## ----robust-sd----------------------------------------------------------------
rob_narrow <- robust_prior(informative, vague_weight = 0.20,
                           vague_sd = 2 * informative$fit_summary$sd)
rob_wide   <- robust_prior(informative, vague_weight = 0.20,
                           vague_sd = 20 * informative$fit_summary$sd)

cat("Narrow vague SD: ", round(rob_narrow$components$vague$fit_summary$sd, 3), "\n")
cat("Default vague SD:", round(rob$components$vague$fit_summary$sd, 3), "\n")
cat("Wide vague SD:   ", round(rob_wide$components$vague$fit_summary$sd, 3), "\n")

## ----map-prior----------------------------------------------------------------
# Five historical control-arm trials -- i.e. a single-arm log-odds scale,
# not a two-arm log-odds RATIO. This distinction matters: see below.
y  <- c(0.10, -0.05, 0.22, 0.05, 0.15)
se <- c(0.12, 0.15, 0.10, 0.18, 0.14)

map <- map_prior(
  y, se,
  outcome_type = "single_arm_log_odds",
  label        = "Historical control (log-odds)"
)
print(map)

## ----map-into-robust----------------------------------------------------------
rob_map <- robust_prior(map, vague_weight = 0.20)
plot(rob_map)

## ----historical-effect-sizes, eval=requireNamespace("metafor", quietly = TRUE)----
# The same five trials as above, but as originally reported -- event
# counts (xi out of ni) rather than pre-computed log-odds estimates.
trials <- data.frame(
  xi = c(11, 9, 14, 10, 12),
  ni = c(45, 40, 55, 48, 42)
)
es <- historical_effect_sizes("single_arm_log_odds", trials)
es[, c("xi", "ni", "y", "se")]

map_from_counts <- map_prior(
  y = es$y, se = es$se,
  outcome_type = "single_arm_log_odds",
  label        = "Historical control (from raw counts)"
)
print(map_from_counts)

## ----tau-presets--------------------------------------------------------------
resolve_tau_prior("log_or")               # two-arm log-odds ratio
resolve_tau_prior("single_arm_log_odds")  # e.g. a historical control rate

## ----map-tau-sensitivity------------------------------------------------------
# An explicit tau_prior always overrides the outcome_type preset, and is
# required (map_prior() errors otherwise) for "mean_difference".
map_tight <- map_prior(y, se, tau_prior = list(family = "half_normal", scale = 0.1))
map_wide  <- map_prior(y, se, tau_prior = list(family = "half_normal", scale = 1.0))

cat("Tight tau prior -- posterior SD:", round(map_tight$fit_summary$sd, 4), "\n")
cat("Wide tau prior  -- posterior SD:", round(map_wide$fit_summary$sd, 4), "\n")

## ----map-tau-posterior--------------------------------------------------------
plot(map$tau_posterior$tau, map$tau_posterior$density, type = "l",
     xlab = expression(tau), ylab = "Posterior density",
     main = "Posterior distribution of the heterogeneity SD")

## ----sceptical-normal---------------------------------------------------------
sc_weak <- sceptical_prior(
  null_value = 0, family = "normal", strength = "weak",
  label = "Log OR (weak sceptic)"
)
sc_moderate <- sceptical_prior(
  null_value = 0, family = "normal", strength = "moderate",
  label = "Log OR (moderate sceptic)"
)
sc_strong <- sceptical_prior(
  null_value = 0, family = "normal", strength = "strong",
  label = "Log OR (strong sceptic)"
)

cat("Weak SD:    ", sc_weak$fit_summary$sd, "\n")
cat("Moderate SD:", sc_moderate$fit_summary$sd, "\n")
cat("Strong SD:  ", sc_strong$fit_summary$sd, "\n")

## ----sceptical-normal-plot----------------------------------------------------
plot(sc_moderate)

## ----sceptical-beta-----------------------------------------------------------
# Null response rate of 20%: sceptic believes treatment is no better than 20%
sc_beta <- sceptical_prior(
  null_value = 0.20,
  family     = "beta",
  strength   = "moderate",
  label      = "Response rate (sceptical)"
)
plot(sc_beta)

## ----sceptical-lognormal------------------------------------------------------
sc_hr <- sceptical_prior(
  null_value = 0,         # log(1) = 0, i.e. HR = 1
  family     = "lognormal",
  strength   = "moderate",
  label      = "Hazard ratio (sceptical)"
)
plot(sc_hr)

## ----enthusiastic-sceptical---------------------------------------------------
enthusiastic <- elicit_beta(
  mean = 0.45, sd = 0.08,
  method = "moments", label = "Response rate (enthusiastic)"
)
sceptical <- sceptical_prior(
  null_value = 0.20, family = "beta", strength = "moderate",
  label = "Response rate (sceptical)"
)

data_obs <- list(type = "binary", x = 18, n = 40)

post_enth <- bayprior:::.conjugate_update(enthusiastic, data_obs)
post_scep <- bayprior:::.conjugate_update(sceptical,    data_obs)

cat("Posterior mean (enthusiastic):", round(post_enth$fit_summary$mean, 3), "\n")
cat("Posterior mean (sceptical):   ", round(post_scep$fit_summary$mean, 3), "\n")
cat("Posterior SD (enthusiastic):  ", round(post_enth$fit_summary$sd,   3), "\n")
cat("Posterior SD (sceptical):     ", round(post_scep$fit_summary$sd,   3), "\n")

## ----power-prior, cache = TRUE------------------------------------------------
base <- elicit_beta(
  mean   = 0.50,
  sd     = 0.20,
  method = "moments",
  label  = "Response rate"
)

calib <- calibrate_power_prior(
  historical_data = list(type = "binary", x = 12, n = 40),
  current_data    = list(type = "binary", x = 18, n = 50),
  base_prior      = base,
  target_bf       = 3,
  delta_grid      = seq(0.05, 1.0, by = 0.05),
  method          = "bayes_factor"
)
print(calib)

## ----power-plot---------------------------------------------------------------
plot(calib)

## ----power-compat, cache = TRUE-----------------------------------------------
calib_compat <- calibrate_power_prior(
  historical_data = list(type = "binary", x = 12, n = 40),
  current_data    = list(type = "binary", x = 18, n = 50),
  base_prior      = base,
  method          = "compatibility",
  delta_grid      = seq(0.05, 1.0, by = 0.05)
)
cat("Optimal delta (BF method):           ", calib$delta_opt, "\n")
cat("Optimal delta (compatibility method):", calib_compat$delta_opt, "\n")

## ----power-normal, cache = TRUE-----------------------------------------------
base_norm <- elicit_normal(
  mean = 0.0, sd = 0.5,
  method = "moments", label = "Mean difference"
)

calib_norm <- calibrate_power_prior(
  historical_data = list(type = "continuous", x = 0.35, sd = 0.3, n = 60),
  current_data    = list(type = "continuous", x = 0.42, sd = 0.3, n = 80),
  base_prior      = base_norm,
  target_bf       = 3,
  delta_grid      = seq(0.05, 1.0, by = 0.10),
  method          = "bayes_factor"
)
print(calib_norm)

## ----choice-table-------------------------------------------------------------
library(knitr)
kable(data.frame(
  Situation = c(
    "No conflict, regulatory requirement",
    "Mild conflict detected",
    "Severe conflict detected",
    "Historical data available",
    "Enthusiastic/sceptical pair for sensitivity analysis"
  ),
  `Recommended prior` = c(
    "Robust mixture (w = 0.20)",
    "Robust mixture (w = 0.30-0.40)",
    "Sceptical prior (moderate-strong)",
    "Power prior (calibrated)",
    "Sceptical prior as second arm"
  ),
  check.names = FALSE
), align = "ll")

## ----report-integration, eval = FALSE-----------------------------------------
# # After running the analyses above...
# prior_report(
#   prior           = prior,
#   conflict        = cd,
#   sensitivity     = sa,
#   robust_prior    = rob,        # adds "Robust Mixture" section to report
#   sceptical_prior = scep,       # adds "Sceptical Prior" section to report
#   power_prior     = calib,      # adds "Power Prior" section with calibration table
#   output_format   = "html",
#   output_file     = "prior_justification_report",
#   trial_name      = "TRIAL-001",
#   sponsor         = "BioPharma Ltd",
#   author          = "J. Smith, Biostatistician"
# )

