###############################################################################
# 0. Packages
###############################################################################

if (!requireNamespace("readxl",   quietly = TRUE))
  install.packages("readxl",   repos = "https://cran.r-project.org/")
if (!requireNamespace("quantreg", quietly = TRUE))
  install.packages("quantreg", repos = "https://cran.r-project.org/")

library(readxl)
library(quantreg)

###############################################################################
# 1. GSCQF CORE FUNCTIONS
###############################################################################

gscqf_obj <- function(par, y, X_mat, tau, q) {
  theta <- X_mat %*% par
  r     <- as.vector(y - theta)
  sum(tau * (log1p(exp(-r / tau)) + log1p(exp(r / tau))) + (2*q - 1) * r)
}

gscqf_grad <- function(par, y, X_mat, tau, q) {
  theta <- X_mat %*% par
  r     <- as.vector(y - theta)
  psi   <- 2 / (exp(r / tau) + 1) - 2 * q
  as.vector(-t(X_mat) %*% psi)
}

fit_gscqf <- function(y, X_mat, q, init = NULL) {
  n_obs <- length(y)
  p     <- ncol(X_mat)

  lm_fit <- lm(y ~ X_mat[, -1])
  tau    <- summary(lm_fit)$sigma / sqrt(n_obs)

  if (is.null(init)) init <- coef(lm_fit)

  fit <- optim(
    par     = init,
    fn      = gscqf_obj,
    gr      = gscqf_grad,
    y       = y,
    X_mat   = X_mat,
    tau     = tau,
    q       = q,
    method  = "BFGS",
    control = list(reltol = 1e-10, maxit = 3000)
  )

  beta_hat  <- fit$par
  theta_hat <- as.vector(X_mat %*% beta_hat)
  r_hat     <- y - theta_hat

  psi_i  <- 2 / (exp(r_hat / tau) + 1) - 2 * q
  dpsi_i <- (2 / tau) * exp(r_hat / tau) / (1 + exp(r_hat / tau))^2

  A_hat <- matrix(0, p, p)
  B_hat <- matrix(0, p, p)
  for (i in seq_len(n_obs)) {
    xi     <- X_mat[i, , drop = FALSE]
    A_hat  <- A_hat + dpsi_i[i] * t(xi) %*% xi
    B_hat  <- B_hat + psi_i[i]^2 * t(xi) %*% xi
  }
  A_hat     <- A_hat / n_obs
  B_hat     <- B_hat / n_obs
  cov_sand  <- tryCatch(solve(A_hat) %*% B_hat %*% solve(A_hat) / n_obs,
                        error = function(e) matrix(NA_real_, p, p))
  se_sand   <- sqrt(abs(diag(cov_sand)))

  list(beta    = beta_hat,
       se_sand = se_sand,
       tau     = tau)
}

fit_gscqf_slr <- function(y, x, q) {
  X_mat <- cbind(1, x)
  fit_gscqf(y, X_mat, q)
}

###############################################################################
# 1b. NLM VERSION OF GSCQF SLR
###############################################################################

gscqf_nlm_fn <- function(par, y, x, tau, q) {
  r      <- y - par[1] - par[2] * x
  er     <- exp(r / tau)

  val    <- sum(tau * (log1p(exp(-r / tau)) + log1p(er)) + (2 * q - 1) * r)

  psi_i  <- 2 / (er + 1) - 2 * q
  dpsi_i <- (2 / tau) * er / (1 + er)^2

  attr(val, "gradient") <- c(sum(psi_i), sum(psi_i * x))

  attr(val, "hessian")  <- matrix(
    c(sum(dpsi_i),        sum(dpsi_i * x),
      sum(dpsi_i * x),    sum(dpsi_i * x^2)),
    nrow = 2
  )

  val
}

fit_gscqf_slr_nlm <- function(y, x, q, tau_fixed = NULL) {
  n_loc   <- length(y)
  lm_fit  <- lm(y ~ x)
  tau_loc <- if (!is.null(tau_fixed)) tau_fixed else summary(lm_fit)$sigma / sqrt(n_loc)
  init    <- c(0, 0)

  fit <- nlm(
    f       = gscqf_nlm_fn,
    p       = init,
    y       = y,
    x       = x,
    tau     = tau_loc,
    q       = q,
    hessian = TRUE,
    gradtol = 1e-10,
    steptol = 1e-10,
    iterlim = 1000
  )

  cov_mat <- tryCatch(solve(fit$hessian),
                      error = function(e) matrix(NA_real_, 2, 2))
  se_vec  <- sqrt(abs(diag(cov_mat)))

  list(beta = fit$estimate, se = se_vec)
}

###############################################################################
# 2. BOOTSTRAP / PERMUTATION HELPERS
###############################################################################

boot_ci_gscqf <- function(y, x, q, B = 50, conf = 0.95) {
  n_loc  <- length(y)
  tau_orig <- summary(lm(y ~ x))$sigma / sqrt(n_loc)
  betas  <- numeric(B)
  for (b in seq_len(B)) {
    idx      <- sample(n_loc, n_loc, replace = TRUE)
    res      <- tryCatch(fit_gscqf_slr_nlm(y[idx], x[idx], q, tau_fixed = tau_orig), error = function(e) NULL)
    betas[b] <- if (!is.null(res)) res$beta[2] else NA_real_
  }
  quantile(betas, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE, type = 5)
}

boot_ci_bk <- function(y, x, q, B = 50, conf = 0.95) {
  n_loc <- length(y)
  betas <- numeric(B)
  for (b in seq_len(B)) {
    idx      <- sample(n_loc, n_loc, replace = TRUE)
    fit_b    <- rq(y[idx] ~ x[idx], tau = q)
    betas[b] <- coef(fit_b)[2]
  }
  quantile(betas, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE, type = 5)
}

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
    fit_p <- tryCatch(summary(rq(sample(y) ~ x, tau = q), se = "iid"),
                      error = function(e) NULL)
    if (is.null(fit_p)) { t_perm[p] <- NA_real_; next }
    t_perm[p] <- coef(fit_p)[2, "Value"] / coef(fit_p)[2, "Std. Error"]
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

gen_errors <- function(n_obs, dist = c("normal", "exponential")) {
  dist <- match.arg(dist)
  if (dist == "normal") rnorm(n_obs, mean = 0, sd = 2)
  else                  rexp(n_obs, rate = 1) + log(0.5)
}

###############################################################################
# 3. SIMULATION PARAMETERS
###############################################################################

set.seed(42)

n_sizes    <- c(10, 30, 50, 250, 1000, 2000)
quantiles  <- c(0.5, 0.60, 0.70, 0.75, 0.9, 0.95, 0.99)
n_sim      <- 1000
n_boot     <- 50
n_perm     <- 50

beta0_true <- 0
beta1_true <- 1
beta1_null <- 0
beta1_alt  <- 0.5

###############################################################################
# 4. SIMULATION STUDY
###############################################################################

if (!exists("RUN_PAPER_SIM") || isTRUE(RUN_PAPER_SIM)) {

cat("\n")
cat(strrep("=", 70), "\n")
cat(" SECTION 5: SIMULATION STUDY\n")
cat(strrep("=", 70), "\n")

##############################################
# TABLES 1 & 2
##############################################

cat("\n--- Tables 1 & 2: Mean and SD of GSCQF slope estimates ---\n")

for (err_dist in c("normal", "exponential")) {
  err_label <- if (err_dist == "normal") "N(0, 2)" else "Exp(1) + log(1/2)"
  cat(sprintf("\nError distribution: %s\n", err_label))
  cat(sprintf("%-10s  %5s  %8s  %8s\n", "Quantile u", "n", "Mean", "Std Dev"))
  cat(strrep("-", 38), "\n")

  for (q in quantiles) {
    for (n_s in n_sizes) {
      slopes  <- numeric(n_sim)
      t_start <- proc.time()[[3]]
      for (sim in seq_len(n_sim)) {
        if (sim %% 500 == 0) {
          cat(sprintf("  [T1&2 | %s | q=%.2f | n=%3d] %d/%d  (%.1fs)\n",
                      err_dist, q, n_s, sim, n_sim,
                      proc.time()[[3]] - t_start))
        }
        x_sim       <- rnorm(n_s, 0, 1)
        y_sim       <- beta0_true + beta1_true * x_sim + gen_errors(n_s, err_dist)
        slopes[sim] <- tryCatch(fit_gscqf_slr_nlm(y_sim, x_sim, q)$beta[2], error = function(e) NA_real_)
      }
      cat(sprintf("%-10.2f  %5d  %8.4f  %8.4f\n",
                  q, n_s, mean(slopes, na.rm = TRUE), sd(slopes, na.rm = TRUE)))
    }
  }
}

##############################################
# TABLES 3 & 4
##############################################

cat("\n\n--- Tables 3 & 4: 95% CI coverage and width for beta1 = 1 ---\n")

for (err_dist in c("normal", "exponential")) {
  err_label <- if (err_dist == "normal") "N(0, 2)" else "Exp(1) + log(1/2)"
  cat(sprintf("\nError distribution: %s\n", err_label))
  cat(sprintf("%-6s %5s  %8s %7s  %11s %10s  %10s %9s\n",
              "u", "n",
              "BK(B)cov", "BK(B)wid",
              "GSCQF(B)cov", "GSCQF(B)wid",
              "GSCQF(A)cov", "GSCQF(A)wid"))
  cat(strrep("-", 80), "\n")

  for (q in quantiles) {
    for (n_s in n_sizes) {
      bk_c   <- bk_w   <- numeric(n_sim)
      gs_bc  <- gs_bw  <- numeric(n_sim)
      gs_ac  <- gs_aw  <- numeric(n_sim)

      t_start <- proc.time()[[3]]
      for (sim in seq_len(n_sim)) {
        if (sim %% 500 == 0) {
          cat(sprintf("  [T3&4 | %s | q=%.2f | n=%3d] %d/%d  (%.1fs)\n",
                      err_dist, q, n_s, sim, n_sim,
                      proc.time()[[3]] - t_start))
        }
        x_sim <- rnorm(n_s, 0, 1)
        y_sim <- beta0_true + beta1_true * x_sim + gen_errors(n_s, err_dist)

        bk_ci      <- tryCatch(boot_ci_bk(y_sim, x_sim, q, B = n_boot),
                               error = function(e) c(NA_real_, NA_real_))
        bk_c[sim]  <- as.numeric(!is.na(bk_ci[1]) &&
                                  bk_ci[1] < beta1_true && beta1_true < bk_ci[2])
        bk_w[sim]  <- bk_ci[2] - bk_ci[1]

        gs_bci      <- tryCatch(boot_ci_gscqf(y_sim, x_sim, q, B = n_boot),
                                error = function(e) c(NA_real_, NA_real_))
        gs_bc[sim]  <- as.numeric(!is.na(gs_bci[1]) &&
                                   gs_bci[1] < beta1_true && beta1_true < gs_bci[2])
        gs_bw[sim]  <- gs_bci[2] - gs_bci[1]

        res_gs      <- tryCatch(fit_gscqf_slr_nlm(y_sim, x_sim, q), error = function(e) NULL)
        gs_alo      <- if (!is.null(res_gs)) res_gs$beta[2] - 1.96 * res_gs$se[2] else NA_real_
        gs_ahi      <- if (!is.null(res_gs)) res_gs$beta[2] + 1.96 * res_gs$se[2] else NA_real_
        gs_ac[sim]  <- as.numeric(!is.na(gs_alo) &&
                                   gs_alo < beta1_true && beta1_true < gs_ahi)
        gs_aw[sim]  <- gs_ahi - gs_alo
      }

      cat(sprintf("%-6.2f %5d  %8.2f %7.2f  %11.2f %10.2f  %10.2f %9.2f\n",
                  q, n_s,
                  mean(bk_c,  na.rm = TRUE), mean(bk_w,  na.rm = TRUE),
                  mean(gs_bc, na.rm = TRUE), mean(gs_bw, na.rm = TRUE),
                  mean(gs_ac, na.rm = TRUE), mean(gs_aw, na.rm = TRUE)))
    }
  }
}

##############################################
# TABLES 5 & 6
##############################################

cat("\n\n--- Tables 5 & 6: Type I error and power (permutation test) ---\n")

n_perm_sizes <- n_sizes[n_sizes <= 50]

for (err_dist in c("normal", "exponential")) {
  err_label <- if (err_dist == "normal") "N(0, 2)" else "Exp(1) + log(1/2)"
  cat(sprintf("\nError distribution: %s\n", err_label))
  cat(sprintf("%-6s %5s  %8s %10s %8s  %8s %10s %8s\n",
              "u", "n",
              "LS(H0)", "GSCQF(H0)", "BK(H0)",
              "LS(H1)", "GSCQF(H1)", "BK(H1)"))
  cat(strrep("-", 78), "\n")

  for (q in quantiles) {
    for (n_s in n_perm_sizes) {
      rej_ls0 <- rej_gs0 <- rej_bk0 <- numeric(n_sim)
      rej_ls1 <- rej_gs1 <- rej_bk1 <- numeric(n_sim)

      t_start <- proc.time()[[3]]
      for (sim in seq_len(n_sim)) {
        if (sim %% 500 == 0) {
          cat(sprintf("  [T5&6 | %s | q=%.2f | n=%3d] %d/%d  (%.1fs)\n",
                      err_dist, q, n_s, sim, n_sim,
                      proc.time()[[3]] - t_start))
        }
        x_sim <- rnorm(n_s, 0, 1)
        eps   <- gen_errors(n_s, err_dist)

        y_h0          <- beta0_true + beta1_null * x_sim + eps
        p_ls0         <- perm_test_ls(y_h0, x_sim, n_perm)
        p_gs0         <- perm_test_gscqf(y_h0, x_sim, q, n_perm)
        p_bk0         <- tryCatch(perm_test_bk(y_h0, x_sim, q, n_perm),
                                  error = function(e) NA_real_)
        rej_ls0[sim]  <- as.numeric(!is.na(p_ls0) && p_ls0 <= 0.05)
        rej_gs0[sim]  <- as.numeric(!is.na(p_gs0) && p_gs0 <= 0.05)
        rej_bk0[sim]  <- as.numeric(!is.na(p_bk0) && p_bk0 <= 0.05)

        y_h1          <- beta0_true + beta1_alt * x_sim + eps
        p_ls1         <- perm_test_ls(y_h1, x_sim, n_perm)
        p_gs1         <- perm_test_gscqf(y_h1, x_sim, q, n_perm)
        p_bk1         <- tryCatch(perm_test_bk(y_h1, x_sim, q, n_perm),
                                  error = function(e) NA_real_)
        rej_ls1[sim]  <- as.numeric(!is.na(p_ls1) && p_ls1 <= 0.05)
        rej_gs1[sim]  <- as.numeric(!is.na(p_gs1) && p_gs1 <= 0.05)
        rej_bk1[sim]  <- as.numeric(!is.na(p_bk1) && p_bk1 <= 0.05)
      }

      ls0_str <- if (q == 0.5) sprintf("%8.2f", mean(rej_ls0, na.rm = TRUE)) else "    --  "
      ls1_str <- if (q == 0.5) sprintf("%8.2f", mean(rej_ls1, na.rm = TRUE)) else "    --  "

      cat(sprintf("%-6.2f %5d  %s %10.2f %8.2f  %s %10.2f %8.2f\n",
                  q, n_s,
                  ls0_str, mean(rej_gs0, na.rm = TRUE), mean(rej_bk0, na.rm = TRUE),
                  ls1_str, mean(rej_gs1, na.rm = TRUE), mean(rej_bk1, na.rm = TRUE)))
    }
  }
}

cat("\nSimulation study complete.\n")

} # end RUN_PAPER_SIM guard

###############################################################################
# 4b. NLM vs BFGS COMPARISON
###############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat(" NLM vs BFGS: Tables 1 & 2")
cat(strrep("=", 70), "\n")

paper_t1 <- list(
  "0.5_10"  = c(0.9830, 0.7696), "0.5_30"  = c(1.0122, 0.3896),
  "0.5_50"  = c(1.0058, 0.3094), "0.5_250" = c(1.0020, 0.1369),
  "0.75_10" = c(1.0145, 0.7526), "0.75_30" = c(0.9889, 0.3959),
  "0.75_50" = c(1.0067, 0.3135), "0.75_250"= c(1.0020, 0.1524),
  "0.9_10"  = c(0.9687, 0.7736), "0.9_30"  = c(1.0198, 0.4327),
  "0.9_50"  = c(0.9939, 0.3530), "0.9_250" = c(1.0110, 0.1738),
  "0.99_10" = c(1.0015, 0.7689), "0.99_30" = c(1.0093, 0.4637),
  "0.99_50" = c(1.0153, 0.4122), "0.99_250"= c(0.9996, 0.2919)
)
paper_t2 <- list(
  "0.5_10"  = c(0.9935, 0.3561), "0.5_30"  = c(1.0060, 0.1853),
  "0.5_50"  = c(1.0017, 0.1358), "0.5_250" = c(0.9971, 0.0518),
  "0.75_10" = c(1.0022, 0.4178), "0.75_30" = c(0.9950, 0.2277),
  "0.75_50" = c(0.9969, 0.1726), "0.75_250"= c(0.9999, 0.0783),
  "0.9_10"  = c(0.9935, 0.4604), "0.9_30"  = c(0.9849, 0.2710),
  "0.9_50"  = c(1.0138, 0.2381), "0.9_250" = c(0.9952, 0.1244),
  "0.99_10" = c(1.0180, 0.4505), "0.99_30" = c(1.0095, 0.3170),
  "0.99_50" = c(0.9894, 0.3470), "0.99_250"= c(0.9815, 0.3663)
)

for (err_dist in c("normal", "exponential")) {
  paper_ref <- if (err_dist == "normal") paper_t1 else paper_t2
  err_label <- if (err_dist == "normal") "N(0, 2)" else "Exp(1) + log(1/2)"
  cat(sprintf("\nError distribution: %s\n", err_label))
  cat(sprintf("%-6s %5s  %9s %8s  %9s %8s  %9s %8s\n",
              "u", "n",
              "BFGS Mean", "BFGS SD",
              "NLM Mean",  "NLM SD",
              "Paper Mean","Paper SD"))
  cat(strrep("-", 78), "\n")

  for (q in quantiles) {
    for (n_s in n_sizes) {
      slopes_bfgs <- numeric(n_sim)
      slopes_nlm  <- numeric(n_sim)
      t_start <- proc.time()[[3]]

      for (sim in seq_len(n_sim)) {
        if (sim %% 500 == 0) {
          cat(sprintf("  [NLM cmp | %s | q=%.2f | n=%3d] %d/%d  (%.1fs)\n",
                      err_dist, q, n_s, sim, n_sim,
                      proc.time()[[3]] - t_start))
        }
        x_sim <- rnorm(n_s, 0, 1)
        y_sim <- beta0_true + beta1_true * x_sim + gen_errors(n_s, err_dist)

        slopes_bfgs[sim] <- fit_gscqf_slr(y_sim, x_sim, q)$beta[2]
        slopes_nlm[sim]  <- tryCatch(
          fit_gscqf_slr_nlm(y_sim, x_sim, q)$beta[2],
          error = function(e) NA_real_
        )
      }

      key <- paste0(q, "_", n_s)
      pv  <- paper_ref[[key]]
      cat(sprintf("%-6.2f %5d  %9.4f %8.4f  %9.4f %8.4f  %9.4f %8.4f\n",
                  q, n_s,
                  mean(slopes_bfgs, na.rm = TRUE), sd(slopes_bfgs, na.rm = TRUE),
                  mean(slopes_nlm,  na.rm = TRUE), sd(slopes_nlm,  na.rm = TRUE),
                  pv[1], pv[2]))
    }
  }
}

cat("\nNLM comparison complete.\n")

###############################################################################
# 5. REAL DATA EXAMPLE
###############################################################################

cat("\n")
cat(strrep("=", 70), "\n")
cat(" SECTION 6: REAL DATA EXAMPLE (Table 7)\n")
cat(strrep("=", 70), "\n")

a <- read_excel("quantregexdata.xlsx")
names(a) <- c("cd19pos_bcells_pct_tot_c", "AGE", "Charlson_CCI")

y_dat  <- a$cd19pos_bcells_pct_tot_c
x1_dat <- a$AGE
x2_dat <- a$Charlson_CCI
n_dat  <- nrow(a)
u_dat  <- 0.5

cat(sprintf("\nDataset: n = %d  |  Outcome: CD19+ B cells (%% lymphocytes)\n",
            n_dat))
cat(sprintf("Covariates: AGE, Charlson Comorbidity Index (CCI)\n"))
cat(sprintf("Quantile: u = %.1f (median)\n\n", u_dat))

fit_bk <- rq(cd19pos_bcells_pct_tot_c ~ AGE + Charlson_CCI,
             tau  = u_dat,
             data = a)
sum_bk  <- summary(fit_bk, se = "iid")
coef_bk <- coef(sum_bk)

perm_p_bk <- function(y, X, j, q, n_perm = 999) {
  fit0   <- summary(rq(y ~ X[, 2] + X[, 3], tau = q), se = "iid")
  t0     <- coef(fit0)[j, "Value"] / coef(fit0)[j, "Std. Error"]
  if (is.na(t0) || is.infinite(t0)) return(NA_real_)
  t_vec  <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    fit_p    <- tryCatch(
      summary(rq(sample(y) ~ X[, 2] + X[, 3], tau = q), se = "iid"),
      error   = function(e) NULL)
    if (is.null(fit_p)) { t_vec[p] <- NA_real_; next }
    t_vec[p] <- coef(fit_p)[j, "Value"] / coef(fit_p)[j, "Std. Error"]
  }
  mean(abs(t_vec) >= abs(t0), na.rm = TRUE)
}

X_dat <- cbind(1, x1_dat, x2_dat)

cat("Computing BK permutation p-values (999 permutations)...\n")
set.seed(123)
p_bk_age <- perm_p_bk(y_dat, X_dat, j = 2, q = u_dat, n_perm = 999)
p_bk_cci <- perm_p_bk(y_dat, X_dat, j = 3, q = u_dat, n_perm = 999)

cat("Fitting GSCQF...\n")
fit_gs <- fit_gscqf(y_dat, X_dat, q = u_dat)

perm_p_gs <- function(y, X, j, q, n_perm = 999) {
  fit0  <- fit_gscqf(y, X, q)
  t0    <- fit0$beta[j] / fit0$se_sand[j]
  if (is.na(t0) || is.infinite(t0)) return(NA_real_)
  t_vec <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    fit_p    <- tryCatch(fit_gscqf(sample(y), X, q), error = function(e) NULL)
    if (is.null(fit_p)) { t_vec[p] <- NA_real_; next }
    t_vec[p] <- fit_p$beta[j] / fit_p$se_sand[j]
  }
  mean(abs(t_vec) >= abs(t0), na.rm = TRUE)
}

cat("Computing GSCQF permutation p-values (999 permutations)...\n")
set.seed(123)
p_gs_age <- perm_p_gs(y_dat, X_dat, j = 2, q = u_dat, n_perm = 999)
p_gs_cci <- perm_p_gs(y_dat, X_dat, j = 3, q = u_dat, n_perm = 999)

zcrit <- qnorm(0.975)

bk_est  <- coef_bk[, "Value"]
bk_se   <- coef_bk[, "Std. Error"]
bk_lo   <- bk_est - zcrit * bk_se
bk_hi   <- bk_est + zcrit * bk_se

gs_est  <- fit_gs$beta
gs_se   <- fit_gs$se_sand
gs_lo   <- gs_est - zcrit * gs_se
gs_hi   <- gs_est + zcrit * gs_se

bk_p  <- c(NA, p_bk_age, p_bk_cci)
gs_p  <- c(NA, p_gs_age, p_gs_cci)
names <- c("Intercept", "AGE", "CCI")

cat("\n\nTable 7. Example results.\n")
cat(strrep("-", 90), "\n")
cat(sprintf("%-12s  %8s  %6s  %7s  %-15s    %8s  %6s  %7s  %-15s\n",
            "Covariate",
            "BK Est", "BK SE", "BK p", "BK 95% CI",
            "GS Est", "GS SE", "GS p", "GS 95% CI"))
cat(strrep("-", 90), "\n")

for (j in seq_along(names)) {
  bk_p_str <- if (is.na(bk_p[j])) "     " else sprintf("%.2f", bk_p[j])
  gs_p_str <- if (is.na(gs_p[j])) "     " else sprintf("%.2f", gs_p[j])
  cat(sprintf("%-12s  %8.2f  %6.2f  %7s  (%.2f, %.2f)    %8.2f  %6.2f  %7s  (%.2f, %.2f)\n",
              names[j],
              bk_est[j], bk_se[j], bk_p_str, bk_lo[j], bk_hi[j],
              gs_est[j], gs_se[j], gs_p_str, gs_lo[j], gs_hi[j]))
}
cat(strrep("-", 90), "\n")
cat("\nNote: p-values from standardised DiCiccio-Romano permutation test (999 permutations).\n")
cat("      Intercept p-value not reported (paper convention).\n")
cat(sprintf("      GSCQF bandwidth tau = %.5f\n\n", fit_gs$tau))

cat("Replication complete.\n")
