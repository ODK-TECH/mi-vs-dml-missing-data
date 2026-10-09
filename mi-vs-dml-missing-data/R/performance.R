# Performance measures 
#
# Formulas and Monte Carlo standard errors (MCSE) follow Morris, White and
# Crowther (2019), "Using simulation studies to evaluate statistical methods",
# Statistics in Medicine 38:2074-2102, Table 6.

add_intervals <- function(results, level = 0.95) {
  q <- qt(1 - (1 - level) / 2, df = results$df)
  results$truth <- unname(true_beta[results$term])
  results$lower <- results$est - q * results$se
  results$upper <- results$est + q * results$se
  results
}

performance_one <- function(est, se, lower, upper, truth) {
  ok <- is.finite(est) & is.finite(se)
  est <- est[ok]; se <- se[ok]; lower <- lower[ok]; upper <- upper[ok]
  k <- length(est)
  err <- est - truth
  emp_se <- sd(est)
  mse <- mean(err^2)
  mod_se <- sqrt(mean(se^2))
  cover <- mean(lower <= truth & truth <= upper)
  data.frame(
    n_ok          = k,
    bias          = mean(err),
    bias_mcse     = emp_se / sqrt(k),
    emp_se        = emp_se,
    emp_se_mcse   = emp_se / sqrt(2 * (k - 1)),
    mod_se        = mod_se,
    mod_se_mcse   = sqrt(var(se^2) / (4 * k * mod_se^2)),
    rmse          = sqrt(mse),
    rmse_mcse     = sqrt(var(err^2) / k) / (2 * sqrt(mse)),
    coverage      = cover,
    coverage_mcse = sqrt(cover * (1 - cover) / k)
  )
}

# One row per scenario x method x coefficient
summarise_performance <- function(results) {
  results <- add_intervals(results)
  key_cols <- c("scenario", "mechanism", "rate", "method", "term")
  groups <- split(results, results[, c("scenario", "method", "term")], drop = TRUE)
  rows <- lapply(groups, function(g) {
    cbind(g[1, key_cols], performance_one(g$est, g$se, g$lower, g$upper, g$truth[1]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
