# =============================================================================
# functions.R  —  GSCQF core functions only, no simulation code
# Source this instead of paper_replication.R to avoid running the simulation
# =============================================================================

library(quantreg)

# --- GSCQF BFGS version ------------------------------------------------------
gscqf_obj <- function(par, y, X_mat, tau, q) {
  r <- as.vector(y - X_mat %*% par)
  sum(tau * (log1p(exp(-r / tau)) + log1p(exp(r / tau))) + (2*q - 1) * r)
}

gscqf_grad <- function(par, y, X_mat, tau, q) {
  r   <- as.vector(y - X_mat %*% par)
  psi <- 2 / (exp(r / tau) + 1) - 2 * q
  as.vector(-t(X_mat) %*% psi)
}

fit_gscqf <- function(y, X_mat, q, init = NULL) {
  n_obs  <- length(y)
  p      <- ncol(X_mat)
  lm_fit <- lm(y ~ X_mat[, -1])
  tau    <- summary(lm_fit)$sigma / sqrt(n_obs)
  if (is.null(init)) init <- coef(lm_fit)
  fit <- optim(par = init, fn = gscqf_obj, gr = gscqf_grad,
               y = y, X_mat = X_mat, tau = tau, q = q,
               method = "BFGS", control = list(reltol = 1e-10, maxit = 3000))
  beta_hat  <- fit$par
  r_hat     <- y - as.vector(X_mat %*% beta_hat)
  psi_i     <- 2 / (exp(r_hat / tau) + 1) - 2 * q
  dpsi_i    <- (2 / tau) * exp(r_hat / tau) / (1 + exp(r_hat / tau))^2
  A_hat <- B_hat <- matrix(0, p, p)
  for (i in seq_len(n_obs)) {
    xi    <- X_mat[i, , drop = FALSE]
    A_hat <- A_hat + dpsi_i[i] * t(xi) %*% xi
    B_hat <- B_hat + psi_i[i]^2 * t(xi) %*% xi
  }
  A_hat    <- A_hat / n_obs
  B_hat    <- B_hat / n_obs
  cov_sand <- tryCatch(solve(A_hat) %*% B_hat %*% solve(A_hat) / n_obs,
                       error = function(e) matrix(NA_real_, p, p))
  list(beta = beta_hat, se_sand = sqrt(abs(diag(cov_sand))), tau = tau)
}

fit_gscqf_slr <- function(y, x, q) fit_gscqf(y, cbind(1, x), q)

# --- GSCQF NLM version -------------------------------------------------------
gscqf_nlm_fn <- function(par, y, x, tau, q) {
  r      <- y - par[1] - par[2] * x
  er     <- exp(r / tau)
  val    <- sum(tau * (log1p(exp(-r / tau)) + log1p(er)) + (2 * q - 1) * r)
  psi_i  <- 2 / (er + 1) - 2 * q
  dpsi_i <- (2 / tau) * er / (1 + er)^2
  attr(val, "gradient") <- c(sum(psi_i), sum(psi_i * x))
  attr(val, "hessian")  <- matrix(
    c(sum(dpsi_i), sum(dpsi_i * x), sum(dpsi_i * x), sum(dpsi_i * x^2)), 2)
  val
}

fit_gscqf_slr_nlm <- function(y, x, q, tau_fixed = NULL) {
  n_loc   <- length(y)
  lm_fit  <- lm(y ~ x)
  tau_loc <- if (!is.null(tau_fixed)) tau_fixed else summary(lm_fit)$sigma / sqrt(n_loc)
  fit <- nlm(f = gscqf_nlm_fn, p = c(0, 0), y = y, x = x, tau = tau_loc, q = q,
             hessian = TRUE, gradtol = 1e-10, steptol = 1e-10, iterlim = 1000)
  cov_mat <- tryCatch(solve(fit$hessian), error = function(e) matrix(NA_real_, 2, 2))
  list(beta = fit$estimate, se = sqrt(abs(diag(cov_mat))))
}

# --- Bootstrap helpers -------------------------------------------------------
boot_ci_gscqf <- function(y, x, q, B = 50, conf = 0.95) {
  n_loc    <- length(y)
  tau_orig <- summary(lm(y ~ x))$sigma / sqrt(n_loc)
  betas    <- numeric(B)
  for (b in seq_len(B)) {
    idx      <- sample(n_loc, n_loc, replace = TRUE)
    res      <- tryCatch(fit_gscqf_slr_nlm(y[idx], x[idx], q, tau_fixed = tau_orig),
                         error = function(e) NULL)
    betas[b] <- if (!is.null(res)) res$beta[2] else NA_real_
  }
  quantile(betas, probs = c((1-conf)/2, 1-(1-conf)/2), na.rm = TRUE, type = 5)
}

boot_ci_bk <- function(y, x, q, B = 50, conf = 0.95) {
  n_loc <- length(y)
  betas <- numeric(B)
  for (b in seq_len(B)) {
    idx      <- sample(n_loc, n_loc, replace = TRUE)
    betas[b] <- coef(rq(y[idx] ~ x[idx], tau = q))[2]
  }
  quantile(betas, probs = c((1-conf)/2, 1-(1-conf)/2), na.rm = TRUE, type = 5)
}

# --- Permutation helpers -----------------------------------------------------
perm_test_gscqf <- function(y, x, q, n_perm = 50) {
  res0 <- fit_gscqf_slr_nlm(y, x, q)
  t0   <- res0$beta[2] / res0$se[2]
  if (is.na(t0) || is.infinite(t0)) return(NA_real_)
  t_perm <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    res_p     <- tryCatch(fit_gscqf_slr_nlm(sample(y), x, q), error = function(e) NULL)
    t_perm[p] <- if (!is.null(res_p)) res_p$beta[2] / res_p$se[2] else NA_real_
  }
  mean(abs(t_perm) >= abs(t0), na.rm = TRUE)
}

perm_test_bk <- function(y, x, q, n_perm = 50) {
  fit0 <- tryCatch(summary(rq(y ~ x, tau = q), se = "iid"), error = function(e) NULL)
  if (is.null(fit0)) return(NA_real_)
  t0 <- coef(fit0)[2, "Value"] / coef(fit0)[2, "Std. Error"]
  if (is.na(t0) || is.infinite(t0)) return(NA_real_)
  t_perm <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    fit_p     <- tryCatch(summary(rq(sample(y) ~ x, tau = q), se = "iid"),
                          error = function(e) NULL)
    t_perm[p] <- if (!is.null(fit_p)) coef(fit_p)[2, "Value"] / coef(fit_p)[2, "Std. Error"] else NA_real_
  }
  mean(abs(t_perm) >= abs(t0), na.rm = TRUE)
}

perm_test_ls <- function(y, x, n_perm = 50) {
  t0     <- coef(summary(lm(y ~ x)))[2, "t value"]
  t_perm <- vapply(seq_len(n_perm),
                   function(p) coef(summary(lm(sample(y) ~ x)))[2, "t value"],
                   numeric(1))
  mean(abs(t_perm) >= abs(t0), na.rm = TRUE)
}
