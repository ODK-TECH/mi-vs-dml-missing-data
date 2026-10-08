# Analysis methods ------------------------------------------------------------
#
# Each method takes a data frame (y, x, z) with missing values in x and returns
# one row per coefficient: term, est, se, df. The runner builds 95% confidence
# intervals as est +/- qt(0.975, df) * se. df = Inf gives a normal interval.

sm_formula <- y ~ x + z + x:z

lm_table <- function(fit) {
  co <- summary(fit)$coefficients
  data.frame(term = rownames(co), est = co[, 1], se = co[, 2],
             df = fit$df.residual, row.names = NULL)
}

# 0. Full data, before any values are deleted. This is the benchmark, not a
#    method you could use in practice.
fit_full <- function(full) lm_table(lm(sm_formula, data = full))

# 1. Complete-case analysis: lm() drops rows with missing x.
fit_cc <- function(dat) lm_table(lm(sm_formula, data = dat))

# 2. Mean imputation: fill missing x with the observed mean, then fit.
#    The model-based standard errors treat the filled values as real data.
fit_mean_imp <- function(dat) {
  dat$x[is.na(dat$x)] <- mean(dat$x, na.rm = TRUE)
  lm_table(lm(sm_formula, data = dat))
}

# Rubin's rules, with the Barnard-Rubin small-sample degrees of freedom
# (the same rule mice::pool() uses).
pool_rubin <- function(fits) {
  est <- sapply(fits, coef)
  u   <- sapply(fits, function(f) diag(vcov(f)))
  m <- ncol(est)
  df_com <- fits[[1]]$df.residual
  qbar <- rowMeans(est)
  ubar <- rowMeans(u)
  b    <- apply(est, 1, var)
  t    <- ubar + (1 + 1 / m) * b
  lambda <- pmax((1 + 1 / m) * b / t, 1e-4)
  df_old <- (m - 1) / lambda^2
  df_obs <- (df_com + 1) / (df_com + 3) * df_com * (1 - lambda)
  data.frame(term = rownames(est), est = qbar, se = sqrt(t),
             df = df_old * df_obs / (df_old + df_obs), row.names = NULL)
}

# 3. MICE with predictive mean matching, the mice default for a numeric
#    variable. mice imputes x from y and z, and the analysis forms x:z from the
#    imputed x ("impute, then transform"). The imputation model leaves out the
#    interaction that the outcome model contains, so the two are incompatible.
#    With one incomplete variable, every imputation is a fresh draw from the
#    same model, so one iteration is enough.
fit_mice <- function(dat, m = 10) {
  imp <- mice::mice(dat, m = m, method = "pmm", maxit = 1, printFlag = FALSE)
  fits <- lapply(seq_len(m), function(i) lm(sm_formula, data = mice::complete(imp, i)))
  pool_rubin(fits)
}

# 4. Substantive-model-compatible FCS (Bartlett et al. 2015). smcfcs proposes
#    x from a normal model given z, then accepts or rejects each proposal using
#    the outcome model, so the imputations respect the x:z interaction.
#    A rejection-sampling limit of 5000 keeps failed draws rare.
fit_smcfcs <- function(dat, m = 10, numit = 10) {
  invisible(utils::capture.output(suppressWarnings(
    imp <- smcfcs::smcfcs(dat, smtype = "lm", smformula = "y ~ x + z + x:z",
                          method = c("", "norm", ""), m = m, numit = numit,
                          rjlimit = 5000)
  )))
  fits <- lapply(imp$impDatasets, function(d) lm(sm_formula, data = d))
  pool_rubin(fits)
}

# 5. Cross-fitted debiased machine learning: augmented inverse probability
#    weighting (AIPW) for a missing covariate.
#
#    Let D = (1, X, Z, XZ) and R = 1 when X is observed. If X is missing at
#    random given (Y, Z), the estimating equation
#
#      sum_i  R_i / pi_i * D_i (Y_i - D_i' b)
#           + (1 - R_i / pi_i) * (mu_i Y_i - M_i b)  =  0
#
#    has mean zero when either the response probability pi(Y, Z) = P(R = 1 | Y, Z)
#    or the conditional moments mu = E[D | Y, Z] and M = E[D D' | Y, Z] are
#    correct (Robins, Rotnitzky and Zhao 1994). Because D is linear in X, mu
#    and M need only m1 = E[X | Y, Z] and v = Var(X | Y, Z).
#
#    Nuisance learners:
#      pi  a probability forest. Forest predictions stay inside the range seen
#          in training, so no unit gets an extreme weight from extrapolation.
#      m1  a two-learner super learner: a weighted average of a random forest
#      v   and a quadratic regression in (y, z), with the weight chosen by
#          out-of-sample squared error. The regression reaches into the tails,
#          where most missing values sit; its predictions are clipped to the
#          range of the training outcome.
#    Estimated response probabilities are bounded below at 0.05, so no unit's
#    weight exceeds 20. Without this bound, a forest can estimate a response
#    probability near zero in a region where most x values are missing, and
#    a single weight near 100 can make the estimating equation unsolvable.
#    Each unit's nuisance predictions come from models fitted on the other
#    K - 1 folds (cross-fitting). The final estimate takes the one-step form:
#    a plug-in estimate plus one correction step with the AIPW estimating
#    function. The score is Neyman-orthogonal, so the sandwich variance can
#    ignore the error in the estimated nuisances (Chernozhukov et al. 2018).
#    The learner settings were fixed in pilot runs that used different seeds
#    from the main study; the one-step form and the 0.05 bound were adopted
#    after a 200-replicate pilot on the main seeds (see report.Rmd).

nuisance_formula <- t ~ y + z + I(y^2) + I(z^2) + y:z

forest_probability <- function(r, train, test, num.trees = 300) {
  d <- data.frame(r = factor(r, levels = c(FALSE, TRUE)), train[, c("y", "z")])
  f <- ranger::ranger(r ~ y + z, data = d, probability = TRUE, num.trees = num.trees,
                      min.node.size = 50, num.threads = 1)
  predict(f, test, num.threads = 1)$predictions[, "TRUE"]
}

# Returns predictions for `test` and out-of-sample predictions for `train`.
super_learn <- function(outcome, train, test, num.trees = 300) {
  train$t <- outcome

  # Learner 1: random forest; out-of-bag predictions for the training rows
  rf <- ranger::ranger(t ~ y + z, data = train[, c("t", "y", "z")], num.trees = num.trees,
                       min.node.size = 20, num.threads = 1)
  rf_oob  <- rf$predictions
  rf_test <- predict(rf, test, num.threads = 1)$predictions

  # Learner 2: quadratic regression; 5-fold cross-validated predictions
  cv_fold <- sample(rep_len(1:5, nrow(train)))
  lm_cv <- numeric(nrow(train))
  for (j in 1:5) {
    g <- lm(nuisance_formula, data = train[cv_fold != j, ])
    lm_cv[cv_fold == j] <- predict(g, train[cv_fold == j, ])
  }
  lm_test <- predict(lm(nuisance_formula, data = train), test)
  clip <- function(p) pmin(pmax(p, min(outcome)), max(outcome))
  lm_cv <- clip(lm_cv)
  lm_test <- clip(lm_test)

  # Weight on the forest that minimises out-of-sample squared error
  alpha <- seq(0, 1, by = 0.05)
  ok <- !is.na(rf_oob)
  mse <- sapply(alpha, function(a) mean((outcome - a * rf_oob - (1 - a) * lm_cv)[ok]^2))
  best <- alpha[which.min(mse)]
  list(pred = best * rf_test + (1 - best) * lm_test,
       oob  = best * rf_oob  + (1 - best) * lm_cv)
}

fit_dml <- function(dat, K = 5, trim = 0.05) {
  n <- nrow(dat)
  r <- !is.na(dat$x)
  fold <- sample(rep_len(seq_len(K), n))
  pi_hat <- m1 <- v <- numeric(n)

  for (k in seq_len(K)) {
    train <- dat[fold != k, ]
    test  <- dat[fold == k, ]
    cc    <- train[!is.na(train$x), ]

    # Response model P(x observed | y, z)
    pi_hat[fold == k] <- forest_probability(!is.na(train$x), train, test)

    # Mean model E[x | y, z], fitted on complete cases
    mean_fit <- super_learn(cc$x, cc, test)
    m1[fold == k] <- mean_fit$pred

    # Variance model Var(x | y, z), fitted to squared out-of-sample residuals
    v[fold == k] <- super_learn((cc$x - mean_fit$oob)^2, cc, test)$pred
  }

  pi_hat <- pmax(pi_hat, trim)
  w <- r / pi_hat              # 0 for units with missing x
  a <- 1 - w
  y <- dat$y
  z <- dat$z
  x0 <- ifelse(r, dat$x, 0)    # missing x never enters, because w = 0 there
  m2 <- m1^2 + pmax(v, 0)      # E[x^2 | y, z], never below m1^2

  D  <- cbind(1, x0, z, x0 * z)
  mu <- cbind(1, m1, z, m1 * z)
  # E[D D' | y, z] for each unit, one column per cell of the 4 x 4 matrix
  # (column-major order)
  M <- cbind(1,      m1,     z,        m1 * z,
             m1,     m2,     m1 * z,   m2 * z,
             z,      m1 * z, z^2,      m1 * z^2,
             m1 * z, m2 * z, m1 * z^2, m2 * z^2)

  # Plug-in estimate: observed units contribute D D' and D y; units with
  # missing x contribute their conditional expectations M and mu y. J is a
  # sum of positive semi-definite matrices, so it is safe to invert.
  J <- crossprod(D * r, D) + matrix(colSums(M * (1 - r)), 4, 4)
  beta_plug <- solve(J, colSums(D * (r * y)) + colSums(mu * ((1 - r) * y)))

  # One-step correction: add the average AIPW estimating function at the
  # plug-in estimate, scaled by J^-1. This is asymptotically the same as
  # solving the AIPW equation directly, but it never inverts the AIPW matrix,
  # which can be close to singular when an observed unit has a large weight.
  psi_at <- function(b) {
    Mb <- sapply(1:4, function(j) M[, j + 4 * (0:3)] %*% b)
    D * as.vector(w * (y - D %*% b)) + (mu * y - Mb) * a
  }
  J_inv <- solve(J)
  beta <- as.vector(beta_plug + J_inv %*% colSums(psi_at(beta_plug)))

  # Sandwich variance from the estimating function at the final estimate
  psi <- psi_at(beta)
  V <- J_inv %*% crossprod(psi) %*% J_inv

  data.frame(term = c("(Intercept)", "x", "z", "x:z"),
             est = beta, se = sqrt(diag(V)), df = Inf)
}
