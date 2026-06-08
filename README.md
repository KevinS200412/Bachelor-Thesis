# Bayesian Sigmoidal Quantile Regression (BSQR)

Bachelor Thesis - comparing a frequentist and a Bayesian version of the GSCQF estimator through a simulation study and a real data example.

---

## What is this?

This is the R code for my bachelor thesis. The main idea is to take Hutson's (2024, 2026) smooth quantile regression estimator (called GSCQF) and ask: *how does it behave when you put a Bayesian prior on it?*

To answer that I:

1. Coded up the frequentist GSCQF from scratch in R (two optimisers, sandwich SEs, bootstrap CIs, permutation tests).
2. Wrapped the same loss function in Stan to get a Bayesian version, then tried 10 different priors to see how sensitive the results are.
3. Ran a big simulation study: 6 sample sizes, 7 quantile levels, 4 error distributions. In order to compare bias, RMSE, and coverage.
4. Applied everything to a real dataset (MGUS, CD19+ B cells) to see if the conclusions hold in practice.

---

## File overview

```
.
├── 00_config.R               # All settings in one place (sample sizes, priors, paths, seeds)
├── 00_run_freq_bayes.R       # Start here, flip TRUE/FALSE to run each phase
├── 01_utils.R                # Small helper functions used across scripts
├── 02_save_freq_results.R    # Runs the frequentist simulation and saves results
├── 03_save_bayes_results.R   # Runs the Bayesian simulation and saves results
├── functions.R               # Just the GSCQF functions, nothing else
├── paper_replication.R       # Reproduces the original paper tables (Tables 1–7)
├── inspect_full.R            # Summary tables: bias, RMSE, coverage, KS tests, density plots
├── posterior_diagnostics.R   # Posterior shape checks and LaTeX table output
├── mgus_application.R        # Real data analysis (MGUS dataset)
├── bsqr_sigmoid.stan         # Stan model for Normal / t / Cauchy / Laplace priors
├── bsqr_sigmoid_uniform.stan # Stan model for Uniform priors
├── quantregexdata.xlsx       # Contains data for real-life application
├── quantregex.R              # Loads the MGUS data and defines the objective function of GSCQF
└── output/
    ├── phase1_freq/          # Frequentist simulation results (.rds per n × dist)
    ├── phase2_bayes/         # Bayesian simulation results (.rds per n × dist × prior × q)
    ├── tables/               # Tables (CSV / LaTeX)
    └── plots/                # Density comparison plots (PDF)
```

---

## How to run

### 1. Install dependencies

```r
install.packages(c("quantreg", "readxl", "ggplot2", "posterior"))

# cmdstanr (requires a working CmdStan installation)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                        "https://cloud.r-project.org"))
cmdstanr::install_cmdstan()
```

### 2. Run the simulation phases

Open `00_run_freq_bayes.R`, set the switches and source:

```r
RUN_FREQ_SIMULATION  <- TRUE   # Phase 1 — frequentist
RUN_BAYES_SIMULATION <- TRUE   # Phase 2 — Bayesian (very slow, split across sessions)
source("00_run_freq_bayes.R")
```

Results are written to `output/phase1_freq/` and `output/phase2_bayes/` as `.rds` files. Phases are independent: results are stored per `(n, dist, prior, q)` so sessions never overwrite each other.

### 3. Analyse results

After the simulation is complete, source the analysis scripts:

| Script | Purpose |
|---|---|
| `inspect_full.R` | Bias, RMSE, coverage, KS tests, density plots |
| `posterior_diagnostics.R` | Skewness, kurtosis, mean/median/mode plots, LaTeX tables |
| `mgus_application.R` | MGUS empirical application |
| `paper_replication.R` | Replicates paper tables using a self-made simulation |

---

## Configuration (`00_config.R`)

| Parameter | Default | Description |
|---|---|---|
| `q_levels` | `c(0.50, 0.60, 0.70, 0.75, 0.90, 0.95, 0.99)` | Quantile levels |
| `n_sizes` | `c(10, 30, 50, 250, 1000, 2000)` | Sample sizes |
| `error_dists` | `normal, exponential, uniform, t3` | Error distributions |
| `N_SIM_FREQ` | `1000` | Monte Carlo replications (frequentist) |
| `N_SIM_BAYES` | `100` | Monte Carlo replications (Bayesian) |
| `B_BOOT` | `200` | Bootstrap replications |
| `MCMC_ITER` | `2000` | MCMC iterations (post-warmup) |
| `MASTER_SEED` | `101` | Global random seed |

### Priors

Ten priors are compared: Normal (3 scales), Laplace (2 scales), Student-t(3), Cauchy, and Uniform (3 scales).

---

## Stan models

- **`bsqr_sigmoid.stan`** — handles `prior_type` 1–4 (Normal, t, Cauchy, Laplace) via.
- **`bsqr_sigmoid_uniform.stan`** — uniform prior variant with bounded parameter block.

Both implement the GSCQF log-likelihood with bandwidth `bw` and quantile level `q` passed as data.

---

## Output files

- `output/phase1_freq/freq_n{n}_{dist}.rds` — list of 1000 simulation results per `(n, dist)`.
- `output/phase2_bayes/bayes_n{n}_{dist}_{prior}_{q}.rds` — list of 100 posterior summaries per `(n, dist, prior, q)`.
- `output/tables/` — CSV and LaTeX tables produced by the analysis scripts.
- `output/plots/` — PDF density plots comparing bootstrap and posterior distributions.

