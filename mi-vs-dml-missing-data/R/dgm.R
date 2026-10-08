# Data-generating mechanism ---------------------------------------------------
#
# Full data (no missing values):
#   Z ~ N(0, 1)
#   X = 0.5 Z + sqrt(0.75) E,  E ~ N(0, 1)      so X ~ N(0, 1) and corr(X, Z) = 0.5
#   Y = b0 + b1 X + b2 Z + b3 X Z + U,  U ~ N(0, 1)
#
# Substantive (analysis) model:  lm(y ~ x + z + x:z)
# Estimands: b1 (coefficient of x) and b3 (the x-by-z interaction).
#
# Missing data: X is incomplete; Y and Z are always observed.
#   MCAR: P(X missing) = constant
#   MAR:  P(X missing) = 0.05 + 0.75 * expit(a0 + gamma_y * Y + gamma_z * Z)
#   MNAR: P(X missing) = 0.05 + 0.75 * expit(a0 + gamma_x * X + gamma_z * Z)
# The floor and ceiling keep every unit's chance of being observed at 20% or
# more, so the positivity condition behind weighting estimators holds.
# calibrate_intercept() picks a0 so the marginal proportion missing hits a target.

true_beta <- c("(Intercept)" = 0, "x" = 1, "z" = 1, "x:z" = 0.5)

# Strength of each mechanism on the log-odds scale
mech_coef <- list(
  MCAR = c(y = 0, x = 0,   z = 0),
  MAR  = c(y = 1, x = 0,   z = 0.5),
  MNAR = c(y = 0, x = 1.5, z = 0.5)
)
p_floor <- 0.05
p_ceiling <- 0.80

p_missing <- function(a0, lp) p_floor + (p_ceiling - p_floor) * plogis(a0 + lp)

simulate_full_data <- function(n, beta = true_beta) {
  z <- rnorm(n)
  x <- 0.5 * z + sqrt(0.75) * rnorm(n)
  y <- beta[1] + beta[2] * x + beta[3] * z + beta[4] * x * z + rnorm(n)
  data.frame(y = y, x = x, z = z)
}

# Linear predictor of P(X missing), excluding the intercept
missing_lp <- function(dat, mechanism) {
  g <- mech_coef[[mechanism]]
  g[["y"]] * dat$y + g[["x"]] * dat$x + g[["z"]] * dat$z
}

# Solve for the intercept a0 that gives the target marginal missing rate,
# using one large full-data sample. Called once per scenario, before the
# replicates, with its own seed.
calibrate_intercept <- function(mechanism, target_rate, n_big = 1e6, seed = 1) {
  if (mechanism == "MCAR") return(NA_real_)
  set.seed(seed)
  lp <- missing_lp(simulate_full_data(n_big), mechanism)
  uniroot(function(a0) mean(p_missing(a0, lp)) - target_rate,
          interval = c(-40, 40), tol = 1e-8)$root
}

# Delete values of X according to the mechanism. Returns the incomplete data.
impose_missingness <- function(dat, mechanism, a0, target_rate) {
  p_miss <- if (mechanism == "MCAR") rep(target_rate, nrow(dat)) else p_missing(a0, missing_lp(dat, mechanism))
  miss <- runif(nrow(dat)) < p_miss
  dat$x[miss] <- NA
  dat
}
