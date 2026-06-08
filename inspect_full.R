# =============================================================================
# inspect_full.R  —  Full analysis across all priors, all error distributions
# =============================================================================

sys.source("00_config.R", envir = globalenv())
library(ggplot2)

# --- Helper: find bayes file robustly ----------------------------------------
find_bayes_file <- function(n, dist, prior_name, q) {
  qname  <- gsub("\\.", "", sprintf("q%04.2f", q))
  pname  <- gsub("[^A-Za-z0-9]", "_", prior_name)
  pname  <- gsub("_+$", "", pname)   # strip trailing underscores
  # Try exact match first
  exact  <- file.path(DIR_BAYES, sprintf("bayes_n%d_%s_%s_%s.rds", n, dist, pname, qname))
  if (file.exists(exact)) return(exact)
  # Try with trailing underscore (how files are actually saved)
  exact2 <- file.path(DIR_BAYES, sprintf("bayes_n%d_%s_%s__%s.rds", n, dist, pname, qname))
  if (file.exists(exact2)) return(exact2)
  # Pattern search as fallback
  pattern <- sprintf("bayes_n%d_%s.*%s.*%s\\.rds", n, dist,
                     gsub("[^A-Za-z0-9]", ".", prior_name), qname)
  matches <- list.files(DIR_BAYES, pattern = pattern, full.names = TRUE)
  if (length(matches) > 0) return(matches[1])
  return(NULL)
}

# --- Settings ----------------------------------------------------------------
n_check <- 10 # change depending on sample size for simulation

# =============================================================================
# PART 1: Frequentist summary — all error distributions
# =============================================================================
cat("=== FREQUENTIST ===\n")

freq_all <- data.frame()
for (dist_check in error_dists) {
  fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_check, dist_check))
  if (!file.exists(fname_f)) next
  freq_res <- readRDS(fname_f)

  for (q in q_levels) {
    key  <- sprintf("q%.2f", q)
    reps <- Filter(Negate(is.null), freq_res[[key]])
    if (length(reps) == 0) next

    b1_g  <- sapply(reps, `[[`, "beta1_gscqf")
    b1_bk <- sapply(reps, `[[`, "beta1_bk")
    cov_g <- sapply(reps, `[[`, "cover_boot_gscqf")
    cov_b <- sapply(reps, `[[`, "cover_boot_bk")
    wid_g <- sapply(reps, `[[`, "width_boot_gscqf")
    wid_b <- sapply(reps, `[[`, "width_boot_bk")

    freq_all <- rbind(freq_all, data.frame(
      dist         = dist_check, q = q,
      bias_gscqf   = mean(b1_g,  na.rm = TRUE) - 1,
      rmse_gscqf   = sqrt(mean((b1_g  - 1)^2, na.rm = TRUE)),
      sd_gscqf     = sd(b1_g,  na.rm = TRUE),
      cover_gscqf  = mean(cov_g, na.rm = TRUE),
      width_gscqf  = mean(wid_g, na.rm = TRUE),
      bias_bk      = mean(b1_bk, na.rm = TRUE) - 1,
      rmse_bk      = sqrt(mean((b1_bk - 1)^2, na.rm = TRUE)),
      sd_bk        = sd(b1_bk,  na.rm = TRUE),
      cover_bk     = mean(cov_b, na.rm = TRUE),
      width_bk     = mean(wid_b, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
}
print(freq_all, digits = 3)

# =============================================================================
# PART 2: Bayesian summary — all priors, all error distributions
# =============================================================================
cat("\n=== BAYESIAN (all priors, all distributions) ===\n")

bayes_all <- data.frame()

for (dist_check in error_dists) {
  for (pr in prior_specs) {
    for (q in q_levels) {
      fname <- find_bayes_file(n_check, dist_check, pr$name, q)
      if (is.null(fname)) next
      reps <- Filter(Negate(is.null), readRDS(fname))
      if (length(reps) == 0) next

      means   <- sapply(reps, `[[`, "mean")
      medians <- sapply(reps, `[[`, "median")
      modes   <- sapply(reps, `[[`, "mode")
      covers  <- sapply(reps, `[[`, "cover")
      widths  <- sapply(reps, `[[`, "width")
      skews   <- sapply(reps, `[[`, "skewness")
      kurts   <- sapply(reps, `[[`, "kurtosis")
      rhats   <- sapply(reps, `[[`, "rhat_max")

      bayes_all <- rbind(bayes_all, data.frame(
        dist         = dist_check, prior = pr$name, q = q,
        bias_mean    = mean(means,   na.rm = TRUE) - 1,
        rmse_mean    = sqrt(mean((means   - 1)^2, na.rm = TRUE)),
        sd_mean      = sd(means,   na.rm = TRUE),
        sd_median    = sd(medians, na.rm = TRUE),
        sd_mode      = sd(modes,   na.rm = TRUE),
        bias_median  = mean(medians, na.rm = TRUE) - 1,
        rmse_median  = sqrt(mean((medians - 1)^2, na.rm = TRUE)),
        bias_mode    = mean(modes,   na.rm = TRUE) - 1,
        rmse_mode    = sqrt(mean((modes   - 1)^2, na.rm = TRUE)),
        coverage     = mean(covers,  na.rm = TRUE),
        width        = mean(widths,  na.rm = TRUE),
        skewness     = mean(skews,   na.rm = TRUE),
        kurtosis     = mean(kurts,   na.rm = TRUE),
        pct_rhat_ok  = mean(rhats < 1.05, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}
print(bayes_all, digits = 3)

if (FALSE) { # PART 3: KS test — disabled
# =============================================================================
# PART 3: KS test — all priors, all error distributions
# =============================================================================
cat("\n=== KS TEST: Bootstrap vs Posterior ===\n")

ks_all <- data.frame()

for (dist_check in error_dists) {
  fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_check, dist_check))
  if (!file.exists(fname_f)) next
  freq_res <- readRDS(fname_f)

  for (pr in prior_specs) {
    for (q in q_levels) {
      key   <- sprintf("q%.2f", q)
      fname <- find_bayes_file(n_check, dist_check, pr$name, q)
      if (is.null(fname)) next

      reps_f <- Filter(Negate(is.null), freq_res[[key]])
      reps_b <- Filter(Negate(is.null), readRDS(fname))
      if (length(reps_f) == 0 || length(reps_b) == 0) next

      n_match   <- min(length(reps_f), length(reps_b))
      boot_pool <- unlist(lapply(reps_f[seq_len(n_match)], `[[`, "boot_draws"))
      post_pool <- unlist(lapply(reps_b[seq_len(n_match)], `[[`, "draws"))
      boot_pool <- boot_pool[!is.na(boot_pool)]
      post_pool <- post_pool[!is.na(post_pool)]
      if (length(boot_pool) < 10 || length(post_pool) < 10) next

      ks <- suppressWarnings(ks.test(boot_pool, post_pool))
      ks_all <- rbind(ks_all, data.frame(
        dist      = dist_check, prior = pr$name, q = q,
        ks_stat   = round(ks$statistic, 4),
        ks_pvalue = round(ks$p.value,   4),
        reject    = ks$p.value < 0.05,
        stringsAsFactors = FALSE
      ))
    }
  }
}
print(ks_all, digits = 3)

} # end if (FALSE) PART 3

if (FALSE) { # PART 4: Density plots — disabled
# =============================================================================
# PART 4: Density plots — 3 priors per plot, all dist, all n, all q
# =============================================================================
cat("\n=== DENSITY PLOTS ===\n")

prior_groups <- list(
  list(name = "normal_family",  priors = c("N(0,1)", "N(0,10)", "N(0,100)")),
  list(name = "laplace_cauchy", priors = c("Laplace(0,0.5)", "Laplace(0,1)", "Cauchy(0,2.5)")),
  list(name = "uniform_family", priors = c("Uniform(-1,1)", "Uniform(-5,5)", "Uniform(-10,10)"))
)

for (dist_plot in error_dists) {
  if (dist_plot != "normal") next
  for (n_plot in c(10, 30, 50, 250)) {
    if (n_plot != 50) next
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_plot, dist_plot))
    if (!file.exists(fname_f)) next
    res_f <- readRDS(fname_f)

    for (q in q_levels) {
      if (q != 0.75) next
      key   <- sprintf("q%.2f", q)
      reps_f    <- Filter(Negate(is.null), res_f[[key]])
      boot_pool <- unlist(lapply(reps_f, `[[`, "boot_draws"))
      boot_pool <- boot_pool[!is.na(boot_pool)]
      if (length(boot_pool) < 10) next
      d_boot <- density(boot_pool)

      for (grp in prior_groups) {
        if (grp$name != "uniform_family") next
        dens_list <- list(data.frame(x = d_boot$x, y = d_boot$y,
                                      source = "Bootstrap (GSCQF)"))

        for (pr in prior_specs) {
          if (!pr$name %in% grp$priors) next
          fname_b <- find_bayes_file(n_plot, dist_plot, pr$name, q)
          if (is.null(fname_b)) next
          reps_b    <- Filter(Negate(is.null), readRDS(fname_b))
          post_pool <- unlist(lapply(reps_b, `[[`, "draws"))
          post_pool <- post_pool[!is.na(post_pool)]
          if (length(post_pool) < 10) next
          d_post <- density(post_pool)
          dens_list[[length(dens_list) + 1]] <- data.frame(
            x = d_post$x, y = d_post$y, source = pr$name
          )
        }

        if (length(dens_list) < 2) next
        df <- do.call(rbind, dens_list)
        df <- df[df$x >= -2 & df$x <= 4, ]

        p <- ggplot(df, aes(x = x, y = y, colour = source, linetype = source)) +
          geom_line(linewidth = 0.9) +
          geom_vline(xintercept = 1, linetype = "dashed", colour = "black") +
          coord_cartesian(xlim = c(-2, 4)) +
          labs(x = expression(beta[1]), y = "Density",
               colour = NULL, linetype = NULL) +
          theme_bw() +
          theme(legend.position      = c(0.98, 0.98),
                legend.justification = c("right", "top"),
                legend.background    = element_rect(fill      = alpha("white", 0.85),
                                                    colour    = "grey70",
                                                    linewidth = 0.3),
                legend.key.size      = unit(0.45, "cm"),
                legend.text          = element_text(size = 7),
                legend.title         = element_blank())

        fname_out <- file.path(DIR_PLOTS,
                               sprintf("density_%s_%s_n%d_q%.2f.pdf",
                                       dist_plot, grp$name, n_plot, q))
        ggsave(fname_out, p, width = 9, height = 5)
        cat(sprintf("  Saved: %s dist=%s n=%d q=%.2f\n",
                    grp$name, dist_plot, n_plot, q))
      }
    }
  }
}

} # end if (FALSE) PART 4

if (FALSE) { # PART 5: Posterior shape — disabled
# =============================================================================
# PART 5: Posterior shape
# =============================================================================
cat("\n=== POSTERIOR SHAPE ACROSS QUANTILES ===\n")
print(bayes_all[, c("dist", "prior", "q", "skewness", "kurtosis")], digits = 3)
} # end if (FALSE) PART 5

cat("\nDone.\n")

# =============================================================================
# PART 6: Open all plots
# =============================================================================
cat("\n=== OPENING ALL PLOTS ===\n")
plot_files <- list.files(DIR_PLOTS, pattern = "density_", full.names = TRUE)
cat(sprintf("Found %d plot files\n", length(plot_files)))
for (f in sort(plot_files)) {
  if (.Platform$OS.type == "windows") {
    shell.exec(normalizePath(f))
  } else {
    system(paste("open", shQuote(f)))
  }
  Sys.sleep(0.3)
}
