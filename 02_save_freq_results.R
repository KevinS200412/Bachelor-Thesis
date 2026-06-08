# =============================================================================
# 02_save_freq_results.R
# Parallelised version using foreach + doParallel (Windows compatible)
# =============================================================================

source("00_config.R")
source("functions.R")

library(quantreg)
library(foreach)
library(doParallel)

N_CORES <- 14
cl <- makeCluster(N_CORES)
registerDoParallel(cl)

cat(sprintf("=== PHASE 1: Saving Frequentist Results (parallel, %d cores) ===\n", N_CORES))
start_all <- proc.time()

gen_errors_ext <- function(n, dist) {
  switch(dist,
    normal      = rnorm(n, 0, 2),
    exponential = rexp(n, 1) + log(0.5),
    uniform     = runif(n, -sqrt(6), sqrt(6)),
    t3          = rt(n, df = 3),
    stop("Unknown dist: ", dist)
  )
}

for (dist in error_dists) {
  for (n in n_sizes) {

    out_file <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n, dist))

    if (file.exists(out_file)) {
      cat(sprintf("SKIP  n=%4d  dist=%-12s (exists)\n", n, dist))
      next
    }

    cat(sprintf("RUN   n=%4d  dist=%-12s\n", n, dist))
    results <- setNames(vector("list", length(q_levels)), sprintf("q%.2f", q_levels))

    for (qi in seq_along(q_levels)) {
      q   <- q_levels[qi]
      key <- sprintf("q%.2f", q)

      # Parallelise over 1000 replications
      reps <- foreach(
        sim = seq_len(N_SIM_FREQ),
        .packages  = c("quantreg"),
        .export    = c("fit_gscqf_slr_nlm", "boot_ci_gscqf", "boot_ci_bk",
                       "gen_errors_ext", "beta0_true", "beta1_true",
                       "MASTER_SEED", "B_BOOT", "CI_LEVEL",
                       "error_dists", "n_sizes", "q_levels"),
        .errorhandling = "pass"
      ) %dopar% {

        seed_sim <- MASTER_SEED * 1e6 +
                    match(dist, error_dists) * 1e5 +
                    match(n, n_sizes)        * 1e4 +
                    qi                       * 1e3 + sim
        set.seed(seed_sim)

        x_sim <- rnorm(n)
        eps   <- gen_errors_ext(n, dist)
        y_sim <- beta0_true + beta1_true * x_sim + eps

        gfit <- tryCatch(fit_gscqf_slr_nlm(y_sim, x_sim, q), error = function(e) NULL)
        if (is.null(gfit)) return(NULL)

        b1_g     <- gfit$beta[2]
        se_g     <- gfit$se[2]
        ci_asymp <- c(b1_g - 1.96 * se_g, b1_g + 1.96 * se_g)

        ci_bg <- tryCatch(
          boot_ci_gscqf(y_sim, x_sim, q, B = B_BOOT),
          error = function(e) c(NA_real_, NA_real_)
        )

        bkfit <- tryCatch(rq(y_sim ~ x_sim, tau = q), error = function(e) NULL)
        b1_bk <- if (!is.null(bkfit)) coef(bkfit)[2] else NA_real_
        ci_bk <- tryCatch(
          boot_ci_bk(y_sim, x_sim, q, B = B_BOOT),
          error = function(e) c(NA_real_, NA_real_)
        )

        # Bootstrap draws for KS test
        set.seed(seed_sim + 1)
        tau_orig   <- summary(lm(y_sim ~ x_sim))$sigma / sqrt(n)
        boot_draws <- replicate(B_BOOT, {
          idx <- sample(n, n, replace = TRUE)
          tryCatch(
            fit_gscqf_slr_nlm(y_sim[idx], x_sim[idx], q,
                               tau_fixed = tau_orig)$beta[2],
            error = function(e) NA_real_
          )
        })

        list(
          beta1_gscqf        = b1_g,
          beta1_bk           = b1_bk,
          boot_draws         = boot_draws,
          ci_boot_gscqf      = ci_bg,
          ci_boot_bk         = ci_bk,
          ci_asymp           = ci_asymp,
          cover_boot_gscqf   = as.integer(!is.na(ci_bg[1]) & ci_bg[1] <= beta1_true & beta1_true <= ci_bg[2]),
          cover_boot_bk      = as.integer(!is.na(ci_bk[1]) & ci_bk[1] <= beta1_true & beta1_true <= ci_bk[2]),
          cover_asymp        = as.integer(!is.na(ci_asymp[1]) & ci_asymp[1] <= beta1_true & beta1_true <= ci_asymp[2]),
          width_boot_gscqf   = ci_bg[2]   - ci_bg[1],
          width_boot_bk      = ci_bk[2]   - ci_bk[1],
          width_asymp        = ci_asymp[2] - ci_asymp[1]
        )
      }

      results[[key]] <- reps
      cat(sprintf("  q=%.2f done\n", q))
    }

    saveRDS(results, out_file)
    cat(sprintf("  Saved -> %s\n", out_file))
  }
}

stopCluster(cl)
cat(sprintf("\nPhase 1 complete in %.1f min.\n",
            (proc.time() - start_all)["elapsed"] / 60))
