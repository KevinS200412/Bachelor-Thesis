# =============================================================================
# 00_run_all.R  –  Master controller
# Set each switch to TRUE to run that phase, FALSE to skip
# =============================================================================

# ---- SWITCHES ---------------------------------------------------------------
RUN_FREQ_SIMULATION  <- FALSE   # Phase 1: run frequentist GSCQF simulation & save
RUN_BAYES_SIMULATION <- FALSE   # Phase 2: run Bayesian GSCQF simulation & save
# Analysis, diagnostics and applications are run from standalone scripts:
#   inspect_full.R          — bias / RMSE / coverage / KS / full summary tables
#   posterior_diagnostics.R — skewness / kurtosis / mean-median-mode plots & LaTeX tables
#   mgus_application.R      — empirical MGUS application

# ---- RUN --------------------------------------------------------------------
source("00_config.R")
if (RUN_FREQ_SIMULATION)  source("02_save_freq_results.R")
if (RUN_BAYES_SIMULATION) source("03_save_bayes_results.R")

cat("\nAll selected phases complete.\n")
