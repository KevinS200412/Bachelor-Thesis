# =============================================================================
# 01_utils.R  –  Shared helper functions
# =============================================================================
source("00_config.R")
library(quantreg)
library(numDeriv)

# --- Error generation --------------------------------------------------------
generate_errors <- function(n, dist) {
  switch(dist,
    normal      = rnorm(n, 0, 2),
    exponential = rexp(n, 1) + log(0.5),
    uniform     = runif(n, -sqrt(6), sqrt(6)),
    t3          = rt(n, df = 3),
    stop("Unknown distribution: ", dist)
  )
}

# --- Data generation ---------------------------------------------------------
generate_data <- function(n, dist, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x   <- rnorm(n)
  eps <- generate_errors(n, dist)
  y   <- beta0_true + beta1_true * x + eps
  list(y = y, x = x, X = cbind(1, x))
}

# --- Bandwidth ---------------------------------------------------------------
compute_tau_fixed    <- function(y, X) summary(lm(y ~ X[, -1]))$sigma / sqrt(length(y))
compute_tau_adaptive <- function(y, X, q) compute_tau_fixed(y, X) / (q * (1 - q))

# --- GSCQF objective and gradient (from paper_replication.R) -----------------
gscqf_obj <- function(par, y, X_mat, tau, q) {
  r <- as.vector(y - X_mat %*% par)
  sum(tau * (log1p(exp(-r / tau)) + log1p(exp(r / tau))) + (2 * q - 1) * r)
}

gscqf_grad <- function(par, y, X_mat, tau, q) {
  r   <- as.vector(y - X_mat %*% par)
  psi <- 2 / (exp(r / tau) + 1) - 2 * q
  as.vector(-t(X_mat) %*% psi)
}

# --- Frequentist GSCQF fit ---------------------------------------------------
fit_gscqf <- function(y, X_mat, q) {
  n   <- length(y)
  p   <- ncol(X_mat)
  lm0 <- lm(y ~ X_mat[, -1])
  tau <- summary(lm0)$sigma / sqrt(n)
  fit <- nlminb(
    start     = coef(lm0),
    objective = gscqf_obj,
    gradient  = gscqf_grad,
    y = y, X_mat = X_mat, tau = tau, q = q,
    control = list(rel.tol = 1e-10, iter.max = 3000)
  )
  beta_hat <- fit$par
  r_hat    <- as.vector(y - X_mat %*% beta_hat)
  er       <- exp(r_hat / tau)
  psi_i    <- 2 / (er + 1) - 2 * q
  dpsi_i   <- (2 / tau) * er / (1 + er)^2
  A <- B   <- matrix(0, p, p)
  for (i in seq_len(n)) {
    xi <- X_mat[i, , drop = FALSE]
    A  <- A + dpsi_i[i] * t(xi) %*% xi
    B  <- B + psi_i[i]^2 * t(xi) %*% xi
  }
  A <- A / n;  B <- B / n
  cov_s <- tryCatch(solve(A) %*% B %*% solve(A) / n,
                    error = function(e) matrix(NA_real_, p, p))
  list(beta = beta_hat, se = sqrt(abs(diag(cov_s))), tau = tau)
}

# --- Classical quantile regression (Koenker) ---------------------------------
fit_bk <- function(y, X_mat, q) {
  fit <- tryCatch(
    rq(y ~ X_mat[, -1], tau = q),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(beta = rep(NA, ncol(X_mat)), se = rep(NA, ncol(X_mat))))
  list(beta = coef(fit), se = summary(fit, se = "boot", B = 200)$coefficients[, 2])
}

# --- Bootstrap CI for GSCQF --------------------------------------------------
bootstrap_gscqf <- function(y, X_mat, q, B = B_BOOT, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n       <- length(y)
  tau_fix <- compute_tau_fixed(y, X_mat)
  boots   <- replicate(B, {
    idx <- sample(n, n, replace = TRUE)
    tryCatch(
      fit_gscqf(y[idx], X_mat[idx, ], q)$beta[2],
      error = function(e) NA_real_
    )
  })
  boots <- boots[!is.na(boots)]
  list(
    draws = boots,
    ci    = quantile(boots, c((1 - CI_LEVEL) / 2, 1 - (1 - CI_LEVEL) / 2), na.rm = TRUE)
  )
}

# --- Asymptotic CI for GSCQF -------------------------------------------------
asymp_ci_gscqf <- function(beta, se) {
  z <- qnorm(1 - (1 - CI_LEVEL) / 2)
  c(beta - z * se, beta + z * se)
}

# --- Posterior mode via KDE --------------------------------------------------
posterior_mode_kde <- function(draws) {
  if (length(unique(draws)) < 2) return(draws[1])
  dens <- density(draws)
  dens$x[which.max(dens$y)]
}

# --- File name helpers -------------------------------------------------------
freq_file  <- function(n, dist) file.path(DIR_FREQ,  sprintf("freq_n%d_%s.rds", n, dist))
bayes_file <- function(n, dist, prior_name) {
  pname <- gsub("[^A-Za-z0-9]", "_", prior_name)
  file.path(DIR_BAYES, sprintf("bayes_n%d_%s_%s.rds", n, dist, pname))
}
