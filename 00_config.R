# =============================================================================
# 00_config.R  –  Central configuration for the BSQR simulation study
# =============================================================================

# --- Quantile levels ---------------------------------------------------------
q_levels <- c(0.50, 0.60, 0.70, 0.75, 0.90, 0.95, 0.99)

# --- QUANTILE SESSION SWITCH -------------------------------------------------
# Set to subset of q_levels to split across sessions.
# Files stored per (n, dist, prior, q) so sessions never overwrite each other.
# Suggested splits:
#   Session 1: c(0.50, 0.60)
#   Session 2: c(0.70, 0.75)
#   Session 3: c(0.90, 0.95)
#   Session 4: c(0.99)
RUN_Q_LEVELS <- q_levels   # run all quantiles at once

# --- Sample sizes ------------------------------------------------------------
n_sizes <- c(10, 30, 50, 250, 1000, 2000)

# --- BAYESIAN SAMPLE SIZE SWITCH ---------------------------------------------
# Run n=10 to n=250 first, evaluate timing, then decide on n=1000 and n=2000
# Session 1: c(10, 30, 50, 250)
# Session 2: c(1000, 2000)
N_SIZES_BAYES <- c(10, 30, 50, 250, 1000)

# --- Error distributions -----------------------------------------------------
error_dists <- c("normal", "exponential", "uniform", "t3")

# --- True DGP ----------------------------------------------------------------
beta0_true <- 0
beta1_true <- 1

# --- Frequentist simulation --------------------------------------------------
N_SIM_FREQ <- 1000
B_BOOT     <- 200

# --- Bayesian simulation -----------------------------------------------------
N_SIM_BAYES   <- 100
MCMC_WARMUP   <- 1000
MCMC_ITER     <- 2000
MCMC_CHAINS   <- 2
MCMC_PARALLEL <- 2

# --- Prior specifications ----------------------------------------------------
prior_specs <- list(
  list(name = "N(0,1)",          type = 1L, scale = 1.0),
  list(name = "N(0,10)",         type = 1L, scale = sqrt(10)),
  list(name = "N(0,100)",        type = 1L, scale = 10.0),
  list(name = "Laplace(0,0.5)",  type = 4L, scale = 0.5),
  list(name = "Laplace(0,1)",    type = 4L, scale = 1.0),
  list(name = "t(3,0,2.5)",      type = 2L, scale = 2.5),
  list(name = "Cauchy(0,2.5)",   type = 3L, scale = 2.5),
  list(name = "Uniform(-1,1)",   type = 5L, scale = 1.0),
  list(name = "Uniform(-5,5)",   type = 5L, scale = 5.0),
  list(name = "Uniform(-10,10)", type = 5L, scale = 10.0)
)

# --- Output directories ------------------------------------------------------
DIR_FREQ   <- "output/phase1_freq"
DIR_BAYES  <- "output/phase2_bayes"
DIR_TABLES <- "output/tables"
DIR_PLOTS  <- "output/plots"

for (d in c(DIR_FREQ, DIR_BAYES, DIR_TABLES, DIR_PLOTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# --- Random seed -------------------------------------------------------------
MASTER_SEED <- 101

# --- CI level ----------------------------------------------------------------
CI_LEVEL <- 0.95
