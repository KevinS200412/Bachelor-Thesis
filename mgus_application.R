# =============================================================================
# mgus_application.R
# Empirical application: MGUS dataset (Bawek et al. 2024, via Hutson 2026)
# Outcome: CD19+ B cells (% lymphocytes); Covariates: age, CCI
# Bayesian GSCQF + frequentist GSCQF: point estimate + SD per coefficient
# Data loading and objective function sourced from quantregex.R (unmodified).
# quantile_of_interest is set to each q_level in the loop below.
# Output: output/tables/mgus_results.csv
#         output/plots/mgus_est_sd_<coef>.pdf
# =============================================================================

source("00_config.R")

if (!requireNamespace("readxl",   quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr",
                   repos = c("https://stan-dev.r-universe.dev",
                             "https://cloud.r-project.org"))
}
if (!requireNamespace("posterior", quietly = TRUE)) install.packages("posterior")
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")

library(readxl); library(cmdstanr); library(posterior); library(ggplot2)

cat("=== PHASE 8: MGUS Empirical Application ===\n")

# --- Source quantregex.R: loads data, defines objective, tau, init, X --------
# After sourcing, the following are available:
#   y, x1, x2  — outcome and covariates
#   n           — sample size
#   tau         — bandwidth (rmse / sqrt(n))
#   init        — OLS starting values
#   objective   — GSCQF loss function(par, y, x1, x2, tau, q)
#   X           — design matrix cbind(1, x1, x2)
source("quantregex.R")

y_sub      <- y
X_sub      <- X          # cbind(1, x1, x2), defined by quantregex.R
n_sub      <- n
bw         <- tau
K          <- ncol(X_sub)
coef_names <- c("intercept", "age", "cci")
mgus_n_sizes <- c(n_sub)

cat(sprintf("Full dataset: n = %d\n", n_sub))

# --- Compile Stan models ------------------------------------------------------
bsqr_model         <- cmdstan_model("bsqr_sigmoid.stan",         force_recompile = FALSE)
bsqr_model_uniform <- cmdstan_model("bsqr_sigmoid_uniform.stan", force_recompile = FALSE)

# --- Priors to use ------------------------------------------------------------
mgus_prior_names <- c("N(0,10)", "N(0,100)", "t(3,0,2.5)", "Cauchy(0,2.5)",
                      "Uniform(-5,5)", "Uniform(-10,10)")
mgus_prior_specs <- Filter(function(pr) pr$name %in% mgus_prior_names, prior_specs)

results_rows <- list()

cat(sprintf("\n--- n = %d ---\n", n_sub))

# --- Frequentist GSCQF: loop over quantiles, changing quantile_of_interest ---
gs_est_mat <- matrix(NA_real_, nrow = length(q_levels), ncol = K)
gs_sd_mat  <- matrix(NA_real_, nrow = length(q_levels), ncol = K)

for (qi in seq_along(q_levels)) {
  quantile_of_interest <- q_levels[qi]

  fit_nlp_q <- tryCatch(
    optim(par     = init,
          fn      = objective,
          y       = y,
          x1      = x1,
          x2      = x2,
          tau     = tau,
          q       = quantile_of_interest,
          method  = "BFGS",
          control = list(reltol = 1e-10, maxit = 3000)),
    error = function(e) NULL
  )

  if (!is.null(fit_nlp_q)) {
    est_q       <- fit_nlp_q$par
    theta_hat_q <- est_q[1] + est_q[2] * x1 + est_q[3] * x2

    psi_i_q  <- 2 * (1 / (1 + exp((y - theta_hat_q) / tau)) - quantile_of_interest)
    dpsi_i_q <- (2 / tau) * exp((y - theta_hat_q) / tau) /
                (1 + exp((y - theta_hat_q) / tau))^2

    A_q <- matrix(0, K, K)
    B_q <- matrix(0, K, K)
    for (i in seq_len(n_sub)) {
      xi  <- X_sub[i, , drop = FALSE]
      A_q <- A_q + dpsi_i_q[i] * t(xi) %*% xi
      B_q <- B_q + psi_i_q[i]^2 * t(xi) %*% xi
    }
    A_q <- A_q / n_sub
    B_q <- B_q / n_sub

    cov_q <- tryCatch(solve(A_q) %*% B_q %*% solve(A_q) / n_sub,
                      error = function(e) matrix(NA_real_, K, K))

    gs_est_mat[qi, ] <- est_q
    gs_sd_mat[qi, ]  <- sqrt(abs(diag(cov_q)))
  }
}

# --- Bayesian + store results -------------------------------------------------
for (qi in seq_along(q_levels)) {
  q <- q_levels[qi]

  for (pr in mgus_prior_specs) {

    is_uniform <- pr$type == 5L
    model_use  <- if (is_uniform) bsqr_model_uniform else bsqr_model

    stan_data <- if (is_uniform) {
      list(N = n_sub, K = K, X = X_sub, y = as.vector(y_sub),
           q = q, bw = bw, bound = pr$scale)
    } else {
      list(N = n_sub, K = K, X = X_sub, y = as.vector(y_sub),
           q = q, bw = bw, prior_type = pr$type,
           beta_loc   = rep(0.0, K),
           beta_scale = rep(pr$scale, K))
    }

    fit_b <- tryCatch(
      model_use$sample(
        data            = stan_data,
        seed            = MASTER_SEED,
        chains          = MCMC_CHAINS,
        parallel_chains = MCMC_PARALLEL,
        iter_warmup     = MCMC_WARMUP,
        iter_sampling   = MCMC_ITER,
        refresh         = 0,
        show_messages   = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b)) {
      post_mean_vec <- rep(NA_real_, K)
      post_sd_vec   <- rep(NA_real_, K)
    } else {
      draws_all     <- as_draws_matrix(fit_b$draws("beta"))
      post_mean_vec <- colMeans(draws_all)
      post_sd_vec   <- apply(draws_all, 2, sd)
    }

    for (k in seq_len(K)) {
      results_rows[[length(results_rows) + 1]] <- data.frame(
        n_sub     = n_sub,
        q         = q,
        prior     = pr$name,
        coef      = coef_names[k],
        gs_est    = gs_est_mat[qi, k],
        gs_sd     = gs_sd_mat[qi, k],
        post_mean = post_mean_vec[k],
        post_sd   = post_sd_vec[k],
        stringsAsFactors = FALSE
      )
    }
  }
  cat(sprintf("  q=%.2f done\n", q))
}

mgus_df <- do.call(rbind, results_rows)
write.csv(mgus_df, file.path(DIR_TABLES, "mgus_results.csv"), row.names = FALSE)
cat(sprintf("\nSaved -> %s\n", file.path(DIR_TABLES, "mgus_results.csv")))

# --- Print --------------------------------------------------------------------
for (cn in unique(mgus_df$coef)) {
  cat(sprintf("\n=== Coefficient: %s ===\n", cn))
  cat(sprintf("%-6s  %5s  %-18s  %9s  %9s  |  %9s  %9s\n",
              "q", "n", "prior",
              "GS est", "GS sd",
              "Post mean", "Post sd"))
  cat(strrep("-", 72), "\n")
  for (n_sub in mgus_n_sizes) {
    for (q in q_levels) {
      for (pr in mgus_prior_specs) {
        row <- mgus_df[mgus_df$coef == cn & mgus_df$n_sub == n_sub &
                         mgus_df$q == q & mgus_df$prior == pr$name, ]
        if (nrow(row) == 0) next
        cat(sprintf("%-6.2f  %5d  %-18s  %9.3f  %9.3f  |  %9.3f  %9.3f\n",
                    q, n_sub, pr$name,
                    row$gs_est, row$gs_sd,
                    row$post_mean, row$post_sd))
      }
    }
  }
}

# --- Plot: point estimate +/- SD, facet by n, colour by prior ----------------
for (cn in unique(mgus_df$coef)) {
  df_c <- mgus_df[mgus_df$coef == cn, ]
  if (nrow(df_c) == 0) next
  df_c$n_label <- factor(paste0("n=", df_c$n_sub),
                          levels = paste0("n=", mgus_n_sizes))

  df_gs        <- df_c
  df_gs$method <- "GSCQF"
  df_gs$est    <- df_gs$gs_est
  df_gs$sd     <- df_gs$gs_sd

  df_bay        <- df_c
  df_bay$method <- "Bayesian"
  df_bay$est    <- df_bay$post_mean
  df_bay$sd     <- df_bay$post_sd

  df_long <- rbind(
    df_gs[,  c("n_label", "q", "prior", "method", "est", "sd")],
    df_bay[, c("n_label", "q", "prior", "method", "est", "sd")]
  )
  df_long$group <- paste(df_long$prior, df_long$method, sep = " | ")

  p <- ggplot(df_long, aes(x = q, y = est, colour = group, fill = group,
                            linetype = method, group = group)) +
    geom_ribbon(aes(ymin = est - sd, ymax = est + sd),
                alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    facet_wrap(~ n_label, ncol = 2) +
    labs(title = sprintf("MGUS — Point Estimate +/- SD: %s", cn),
         x = "Quantile level", y = cn,
         colour = "Prior | Method", fill = "Prior | Method",
         linetype = "Method") +
    theme_bw()

  ggsave(file.path(DIR_PLOTS, sprintf("mgus_est_sd_%s.pdf", cn)),
         p, width = 11, height = 8)
  cat(sprintf("  Plot saved for coef=%s\n", cn))
}

cat("\nPhase 7 complete.\n")
