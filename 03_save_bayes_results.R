# =============================================================================
# 03_save_bayes_results.R
# Parallelised version using foreach + doParallel (Windows compatible)
# One RDS file per (n, dist, prior, q)
# =============================================================================

sys.source("00_config.R", envir = globalenv())
if (!exists("RUN_Q_LEVELS")) RUN_Q_LEVELS <- q_levels

library(cmdstanr)
library(posterior)
library(foreach)
library(doParallel)

N_CORES <- 7
cl <- makeCluster(N_CORES)
registerDoParallel(cl)
cat(sprintf("Cluster registered with %d workers\n", getDoParWorkers()))

cat("Compiling Stan models...\n")
bsqr_main    <- cmdstan_model("bsqr_sigmoid.stan",         force_recompile = FALSE)
bsqr_uniform <- cmdstan_model("bsqr_sigmoid_uniform.stan", force_recompile = FALSE)
cat("Done.\n\n")

cat(sprintf("=== PHASE 2: Bayesian Simulation (parallel, %d cores) ===\n", N_CORES))
cat(sprintf("Quantile levels: %s\n\n", paste(RUN_Q_LEVELS, collapse = ", ")))
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

posterior_mode_kde <- function(draws) {
  draws <- draws[!is.na(draws)]
  if (length(unique(draws)) < 2) return(draws[1])
  d <- density(draws); d$x[which.max(d$y)]
}

posterior_skewness <- function(draws) {
  draws <- draws[!is.na(draws)]; m <- mean(draws); s <- sd(draws)
  if (s == 0) return(NA_real_)
  mean(((draws - m) / s)^3)
}

posterior_kurtosis <- function(draws) {
  draws <- draws[!is.na(draws)]; m <- mean(draws); s <- sd(draws)
  if (s == 0) return(NA_real_)
  mean(((draws - m) / s)^4) - 3
}

bayes_file_q <- function(n, dist, prior_name, q) {
  pname <- gsub("[^A-Za-z0-9]", "_", prior_name)
  qname <- gsub("\\.", "", sprintf("q%04.2f", q))
  file.path(DIR_BAYES, sprintf("bayes_n%d_%s_%s_%s.rds", n, dist, pname, qname))
}

for (pr in prior_specs) {
  for (dist in error_dists) {
    for (n in N_SIZES_BAYES) {
      for (q in RUN_Q_LEVELS) {

        out_file   <- bayes_file_q(n, dist, pr$name, q)
        is_uniform <- pr$type == 5L

        if (file.exists(out_file)) {
          cat(sprintf("SKIP  prior=%-20s  n=%4d  dist=%-12s  q=%.2f\n",
                      pr$name, n, dist, q))
          next
        }

        cat(sprintf("RUN   prior=%-20s  n=%4d  dist=%-12s  q=%.2f\n",
                    pr$name, n, dist, q))

        qi <- which(q_levels == q)

        reps <- foreach(
          sim            = seq_len(N_SIM_BAYES),
          .packages      = c("cmdstanr", "posterior"),
          .export        = c("gen_errors_ext", "posterior_mode_kde",
                             "posterior_skewness", "posterior_kurtosis",
                             "beta0_true", "beta1_true", "MASTER_SEED",
                             "MCMC_CHAINS", "MCMC_WARMUP", "MCMC_ITER",
                             "error_dists", "n_sizes", "q_levels", "N_SIZES_BAYES",
                             "is_uniform", "pr", "dist", "n", "q", "qi",
                             "DIR_BAYES"),
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
          X_sim <- cbind(1, x_sim)
          bw    <- summary(lm(y_sim ~ x_sim))$sigma / sqrt(n)

          if (is_uniform) {
            model     <- cmdstan_model("bsqr_sigmoid_uniform.stan", force_recompile = FALSE)
            stan_data <- list(
              N = n, K = 2L, X = X_sim,
              y = as.vector(y_sim),
              q = q, bw = bw,
              bound = pr$scale
            )
          } else {
            model     <- cmdstan_model("bsqr_sigmoid.stan", force_recompile = FALSE)
            stan_data <- list(
              N = n, K = 2L, X = X_sim,
              y = as.vector(y_sim),
              q = q, bw = bw,
              prior_type = pr$type,
              beta_loc   = c(0.0, 0.0),
              beta_scale = c(pr$scale, pr$scale)
            )
          }

          fit <- tryCatch(
            model$sample(
              data            = stan_data,
              seed            = seed_sim,
              chains          = MCMC_CHAINS,
              parallel_chains = 1,
              iter_warmup     = MCMC_WARMUP,
              iter_sampling   = MCMC_ITER,
              refresh         = 0,
              show_messages   = FALSE
            ),
            error = function(e) NULL
          )

          if (is.null(fit)) return(NULL)

          draws_all <- as_draws_matrix(fit$draws("beta"))
          b1_draws  <- draws_all[, 2]
          rhat_max  <- max(fit$summary("beta")$rhat, na.rm = TRUE)
          ci        <- quantile(b1_draws, c(0.025, 0.975))

          list(
            draws    = b1_draws,
            mean     = mean(b1_draws),
            median   = median(b1_draws),
            mode     = posterior_mode_kde(b1_draws),
            ci       = ci,
            sd       = sd(b1_draws),
            skewness = posterior_skewness(b1_draws),
            kurtosis = posterior_kurtosis(b1_draws),
            cover    = as.integer(ci[1] <= beta1_true & beta1_true <= ci[2]),
            width    = ci[2] - ci[1],
            rhat_max = rhat_max
          )
        }

        saveRDS(reps, out_file)
        cat(sprintf("  Saved -> %s\n", out_file))
      }
    }
  }
}

stopCluster(cl)
cat(sprintf("\nPhase 2 complete in %.1f min.\n",
            (proc.time() - start_all)["elapsed"] / 60))
