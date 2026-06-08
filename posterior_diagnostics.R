# =============================================================================
# plot_mean_median_mode_9priors.R
#
# Posterior mean, median and mode across quantile levels
# 3x3 grid of 9 priors (all except Uniform(-1,1))
# One plot per (dist, n) combination
# =============================================================================

sys.source("00_config.R", envir = globalenv())

library(ggplot2)
library(tidyr)

set.seed(MASTER_SEED)

MAKE_PLOTS <- FALSE   # set TRUE to regenerate the 24 PDFs

# ---------------------------------------------------------------------------
# Helper: find Bayesian file
# ---------------------------------------------------------------------------
find_bayes_file <- function(n, dist, prior_name, q) {
  qname <- sprintf("q%03d", as.integer(q * 100))
  pname <- gsub("[^A-Za-z0-9]", "_", prior_name)
  for (sep in c("_", "__")) {
    f <- file.path(DIR_BAYES,
                   sprintf("bayes_n%d_%s_%s%s%s.rds",
                           n, dist, pname, sep, qname))
    if (file.exists(f)) return(f)
  }
  return(NULL)
}

# ---------------------------------------------------------------------------
# Priors: all except Uniform(-1,1)
# Fixed 3x3 order
# ---------------------------------------------------------------------------
priors_plot <- c(
  "N(0,1)",        "N(0,10)",       "N(0,100)",
  "Laplace(0,0.5)","Laplace(0,1)",  "t(3,0,2.5)",
  "Cauchy(0,2.5)", "Uniform(-5,5)", "Uniform(-10,10)"
)

# ---------------------------------------------------------------------------
# Main loop: one plot per (dist, n)
# ---------------------------------------------------------------------------
if (FALSE) {  # mean/median/mode plots already generated — disabled
cat("=== MEAN / MEDIAN / MODE PLOTS ===\n\n")

for (dist in error_dists) {
  for (n in N_SIZES_BAYES) {

    cat(sprintf("Processing dist=%s, n=%d\n", dist, n))
    rows <- list()

    for (prior_name in priors_plot) {
      for (q in q_levels) {

        fname <- find_bayes_file(n, dist, prior_name, q)
        if (is.null(fname)) next

        reps <- tryCatch(
          Filter(Negate(is.null), readRDS(fname)),
          error = function(e) NULL
        )
        if (is.null(reps) || length(reps) == 0) next

        means   <- sapply(reps, `[[`, "mean")
        medians <- sapply(reps, `[[`, "median")
        modes   <- sapply(reps, `[[`, "mode")

        rows[[length(rows) + 1]] <- data.frame(
          prior     = prior_name,
          q         = q,
          estimator = "Mean",
          value     = mean(means,   na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1]] <- data.frame(
          prior     = prior_name,
          q         = q,
          estimator = "Median",
          value     = mean(medians, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1]] <- data.frame(
          prior     = prior_name,
          q         = q,
          estimator = "Mode",
          value     = mean(modes,   na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) == 0) {
      cat(sprintf("  [SKIP] No data for dist=%s n=%d\n", dist, n))
      next
    }

    df <- do.call(rbind, rows)
    df$prior     <- factor(df$prior,     levels = priors_plot)
    df$estimator <- factor(df$estimator, levels = c("Mean", "Median", "Mode"))

    p <- ggplot(df, aes(x        = q,
                        y        = value,
                        colour   = estimator,
                        linetype = estimator)) +
      geom_hline(yintercept = 1, linetype = "dashed",
                 colour = "black", linewidth = 0.5) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~ prior, ncol = 6, scales = "fixed") +
      scale_colour_manual(
        values = c(Mean   = "#2166ac",
                   Median = "#d73027",
                   Mode   = "#4dac26"),
        name   = "Point estimate"
      ) +
      scale_linetype_manual(
        values = c(Mean   = "solid",
                   Median = "dashed",
                   Mode   = "dotted"),
        name   = "Point estimate"
      ) +
      scale_x_continuous(breaks = q_levels,
                         labels = q_levels) +
      labs(
        x = "Quantile level",
        y = expression(hat(beta)[1])
      ) +
      theme_bw(base_size = 10) +
      theme(
        strip.text           = element_text(size = 8),
        legend.position      = c(0.98, 0.84),
        legend.justification = c("right", "top"),
        legend.background    = element_rect(fill      = alpha("white", 0.85),
                                            colour    = "grey70",
                                            linewidth = 0.3),
        legend.key.size      = unit(0.55, "cm"),
        legend.text          = element_text(size = 9),
        legend.title         = element_blank(),
        axis.text.x          = element_text(angle = 45, hjust = 1, size = 7)
      )

    fname_out <- file.path(
      DIR_PLOTS,
      sprintf("mean_median_mode_9priors_%s_n%d.pdf", dist, n)
    )
    ggsave(fname_out, p, width = 10, height = 8)
    cat(sprintf("  Saved: %s\n", basename(fname_out)))
  }
}

cat("\nDone (mean/median/mode).\n")

} # end if (FALSE)

N_SUB <- 200   # subsample size for fair comparison

# =============================================================================
# Grouped skewness & kurtosis – combined facet_grid layout
# Group 1: N(0,10), t(3,0,2.5), Uniform(-10,10)
# Group 2: N(0,100), Cauchy(0,2.5), Uniform(-5,5)
# Layout: 2 rows (Skewness top, Kurtosis bottom) x 3 cols (priors)
# Row/col labels as plain text (no grey strip box)
# 24 PDFs: 2 groups x 3 error dists x 4 n values
# =============================================================================

prior_groups <- list(
  group1 = c("N(0,10)",   "t(3,0,2.5)",    "Uniform(-10,10)"),
  group2 = c("N(0,100)",  "Cauchy(0,2.5)", "Uniform(-5,5)")
)

ref_lines <- data.frame(
  metric = factor(c("Skewness", "Kurtosis"), levels = c("Skewness", "Kurtosis")),
  ref    = c(0, 3)
)

cat("\n=== GROUPED SKEWNESS / KURTOSIS PLOTS (2 groups) ===\n\n")

if (MAKE_PLOTS) {

for (grp_name in names(prior_groups)) {
  grp_priors <- prior_groups[[grp_name]]

  for (dist in error_dists[error_dists != "uniform"]) {
    for (n in N_SIZES_BAYES) {
      cat(sprintf("Processing %s  dist=%s  n=%d\n", grp_name, dist, n))
      rows <- list()

      # --- Posterior ---
      for (prior_name in grp_priors) {
        for (q in q_levels) {
          fname <- find_bayes_file(n, dist, prior_name, q)
          if (is.null(fname)) next
          reps <- tryCatch(
            Filter(Negate(is.null), readRDS(fname)),
            error = function(e) NULL
          )
          if (is.null(reps) || length(reps) == 0) next

          skews <- sapply(reps, function(r) {
            draws <- r$draws
            if (is.null(draws) || length(draws) < 3) return(NA)
            draws <- sample(draws, min(N_SUB, length(draws)), replace = FALSE)
            m <- mean(draws, na.rm = TRUE); s <- sd(draws, na.rm = TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((draws - m) / s)^3, na.rm = TRUE)
          })
          kurts <- sapply(reps, function(r) {
            draws <- r$draws
            if (is.null(draws) || length(draws) < 3) return(NA)
            draws <- sample(draws, min(N_SUB, length(draws)), replace = FALSE)
            m <- mean(draws, na.rm = TRUE); s <- sd(draws, na.rm = TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((draws - m) / s)^4, na.rm = TRUE)
          })
          rows[[length(rows) + 1]] <- data.frame(
            prior    = prior_name, q = q, source = "Posterior",
            skewness = mean(skews, na.rm = TRUE),
            kurtosis = mean(kurts, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }

      # --- GSCQF ---
      fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n, dist))
      if (file.exists(fname_f)) {
        freq_res <- readRDS(fname_f)
        for (q in q_levels) {
          key  <- sprintf("q%.2f", q)
          reps <- Filter(Negate(is.null), freq_res[[key]])
          if (length(reps) == 0) next
          skew_g <- sapply(reps, function(r) {
            bd <- r$boot_draws
            if (is.null(bd) || length(bd) < 3) return(NA)
            m <- mean(bd, na.rm = TRUE); s <- sd(bd, na.rm = TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((bd - m) / s)^3, na.rm = TRUE)
          })
          kurt_g <- sapply(reps, function(r) {
            bd <- r$boot_draws
            if (is.null(bd) || length(bd) < 3) return(NA)
            m <- mean(bd, na.rm = TRUE); s <- sd(bd, na.rm = TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((bd - m) / s)^4, na.rm = TRUE)
          })
          for (prior_name in grp_priors) {
            rows[[length(rows) + 1]] <- data.frame(
              prior    = prior_name, q = q, source = "GSCQF",
              skewness = mean(skew_g, na.rm = TRUE),
              kurtosis = mean(kurt_g, na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          }
        }
      }

      if (length(rows) == 0) {
        cat(sprintf("  [SKIP]\n")); next
      }

      df        <- do.call(rbind, rows)
      df$prior  <- factor(df$prior,  levels = grp_priors)
      df$source <- factor(df$source, levels = c("Posterior", "GSCQF"))

      # Reshape to long format for facet_grid
      df_long <- rbind(
        data.frame(prior = df$prior, q = df$q, source = df$source,
                   metric = "Skewness", value = df$skewness,
                   stringsAsFactors = FALSE),
        data.frame(prior = df$prior, q = df$q, source = df$source,
                   metric = "Kurtosis", value = df$kurtosis,
                   stringsAsFactors = FALSE)
      )
      df_long$metric <- factor(df_long$metric, levels = c("Skewness", "Kurtosis"))
      df_long$prior  <- factor(df_long$prior,  levels = grp_priors)
      df_long$source <- factor(df_long$source, levels = c("Posterior", "GSCQF"))

      p <- ggplot(df_long, aes(x = q, y = value,
                               colour = source, linetype = source)) +
        geom_hline(data          = ref_lines,
                   aes(yintercept = ref),
                   linetype      = "dashed",
                   colour        = "grey50",
                   inherit.aes   = FALSE) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 1.5) +
        facet_grid(metric ~ prior, scales = "free_y") +
        scale_x_continuous(breaks = q_levels, labels = q_levels) +
        scale_colour_manual(values = c(Posterior = "#2166ac", GSCQF = "#d73027")) +
        scale_linetype_manual(values = c(Posterior = "solid", GSCQF = "dashed")) +
        scale_y_continuous(expand = expansion(mult = 0.30)) +
        labs(x = "Quantile level", y = NULL) +
        theme_bw(base_size = 10) +
        theme(
          strip.text           = element_text(size = 9),
          axis.text.x          = element_text(angle = 45, hjust = 1, size = 7),
          legend.position      = c(0.99, 0.99),
          legend.justification = c("right", "top"),
          legend.background    = element_rect(fill      = alpha("white", 0.85),
                                              colour    = "grey70",
                                              linewidth = 0.3),
          legend.key.size      = unit(0.55, "cm"),
          legend.text          = element_text(size = 9),
          legend.title         = element_blank()
        )

      fname_out <- file.path(DIR_PLOTS,
                             sprintf("skew_kurt_%s_%s_n%d.pdf", grp_name, dist, n))
      ggsave(fname_out, p, width = 12, height = 6)
      cat(sprintf("  Saved: %s\n", basename(fname_out)))
    }
  }
}

} # end if (MAKE_PLOTS)

cat("\nAll done.\n")

# ==========================================================================
# Skewness : facet_grid(n ~ prior) — 2 rows × 3 cols, 2 lines each
#            one PDF: skew_group1_normal_n10_n250.pdf
# Kurtosis : facet_wrap(~ prior) — 1 row × 3 cols, 5 lines each
#            Posterior n=10/250 + GSCQF n=10/30/250
#            one PDF: kurt_group1_normal_n10_n250.pdf
# ==========================================================================

grp_priors  <- prior_groups[["group1"]]   # N(0,10), t(3,0,2.5), Uniform(-10,10)
dist_comb   <- "normal"
n_post      <- c(10, 250)          # Posterior loaded for these
n_gscqf     <- c(10, 250)           # GSCQF loaded for these

rows_comb <- list()

# --- Posterior (n = 10 and 250) ---
for (n in n_post) {
  for (prior_name in grp_priors) {
    for (q in q_levels) {
      fname <- find_bayes_file(n, dist_comb, prior_name, q)
      if (is.null(fname)) next
      reps <- tryCatch(Filter(Negate(is.null), readRDS(fname)),
                        error = function(e) NULL)
      if (is.null(reps) || length(reps) == 0) next
      skews <- sapply(reps, function(r) {
        draws <- r$draws
        if (is.null(draws) || length(draws) < 3) return(NA)
        draws <- sample(draws, min(N_SUB, length(draws)), replace = FALSE)
        m <- mean(draws, na.rm = TRUE); s <- sd(draws, na.rm = TRUE)
        if (is.na(s) || s == 0) return(NA)
        mean(((draws - m) / s)^3, na.rm = TRUE)
      })
      kurts <- sapply(reps, function(r) {
        draws <- r$draws
        if (is.null(draws) || length(draws) < 3) return(NA)
        draws <- sample(draws, min(N_SUB, length(draws)), replace = FALSE)
        m <- mean(draws, na.rm = TRUE); s <- sd(draws, na.rm = TRUE)
        if (is.na(s) || s == 0) return(NA)
        mean(((draws - m) / s)^4, na.rm = TRUE)
      })
      rows_comb[[length(rows_comb) + 1]] <- data.frame(
        prior    = prior_name, q = q, n = n, source = "Posterior",
        skewness = mean(skews, na.rm = TRUE),
        kurtosis = mean(kurts, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}

# --- GSCQF (n = 10, 50, 250) ---
for (n in n_gscqf) {
  fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n, dist_comb))
  if (!file.exists(fname_f)) next
  freq_res <- readRDS(fname_f)
  for (q in q_levels) {
    key  <- sprintf("q%.2f", q)
    reps <- Filter(Negate(is.null), freq_res[[key]])
    if (length(reps) == 0) next
    skew_g <- sapply(reps, function(r) {
      bd <- r$boot_draws
      if (is.null(bd) || length(bd) < 3) return(NA)
      m <- mean(bd, na.rm = TRUE); s <- sd(bd, na.rm = TRUE)
      if (is.na(s) || s == 0) return(NA)
      mean(((bd - m) / s)^3, na.rm = TRUE)
    })
    kurt_g <- sapply(reps, function(r) {
      bd <- r$boot_draws
      if (is.null(bd) || length(bd) < 3) return(NA)
      m <- mean(bd, na.rm = TRUE); s <- sd(bd, na.rm = TRUE)
      if (is.na(s) || s == 0) return(NA)
      mean(((bd - m) / s)^4, na.rm = TRUE)
    })
    for (prior_name in grp_priors) {
      rows_comb[[length(rows_comb) + 1]] <- data.frame(
        prior    = prior_name, q = q, n = n, source = "GSCQF",
        skewness = mean(skew_g, na.rm = TRUE),
        kurtosis = mean(kurt_g, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(rows_comb) > 0) {
  df_comb <- do.call(rbind, rows_comb)
  df_comb$prior <- factor(df_comb$prior, levels = grp_priors)

  theme_leg_corner <- theme_bw(base_size = 10) +
    theme(
      strip.text           = element_text(size = 9),
      axis.text.x          = element_text(angle = 45, hjust = 1, size = 7),
      legend.position      = c(0.99, 0.99),
      legend.justification = c("right", "top"),
      legend.background    = element_rect(fill      = alpha("white", 0.85),
                                          colour    = "grey70",
                                          linewidth = 0.3),
      legend.key.size      = unit(0.6, "cm"),
      legend.key.width     = unit(1.0, "cm"),
      legend.text          = element_text(size = 9),
      legend.title         = element_blank()
    )

  # ------------------------------------------------------------------
  # SKEWNESS: 2-row × 3-col facet_grid (rows = n=10 / n=250)
  # Only uses n=10 and n=250; 2 lines per panel (Posterior + GSCQF)
  # ------------------------------------------------------------------
  df_skew <- df_comb[df_comb$n %in% n_post, ]
  df_skew$n_label <- factor(paste0("n = ", df_skew$n),
                              levels = paste0("n = ", n_post))
  df_skew$source  <- factor(df_skew$source, levels = c("Posterior", "GSCQF"))

  p_skew <- ggplot(df_skew,
                    aes(x = q, y = skewness,
                        colour = source, linetype = source, shape = source)) +
    geom_hline(yintercept = 0, linetype = "dashed",
                colour = "grey50", inherit.aes = FALSE) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    facet_grid(n_label ~ prior, scales = "free_y") +
    scale_x_continuous(breaks = q_levels, labels = q_levels) +
    scale_colour_manual(values   = c(Posterior = "#2166ac", GSCQF = "#d73027")) +
    scale_linetype_manual(values = c(Posterior = "solid",   GSCQF = "dashed")) +
    scale_shape_manual(values    = c(Posterior = 16,        GSCQF = 1)) +
    scale_y_continuous(expand = expansion(mult = 0.30)) +
    labs(x = "Quantile level", y = "Skewness") +
    theme_leg_corner

  fname_skew <- file.path(DIR_PLOTS, "skew_group1_normal_n10_n250.pdf")
  ggsave(fname_skew, p_skew, width = 12, height = 6)
  cat(sprintf("  Saved: %s\n", basename(fname_skew)))

  # ------------------------------------------------------------------
  # KURTOSIS: two separate plots — one per n value
  # Each: 1-row × 3-col facet_wrap, 2 lines (Posterior + GSCQF)
  # ------------------------------------------------------------------
  for (n_kurt in c(10, 250)) {
    df_kurt_n <- df_comb[df_comb$n == n_kurt, ]
    df_kurt_n$source <- factor(df_kurt_n$source, levels = c("Posterior", "GSCQF"))

    p_kurt <- ggplot(df_kurt_n,
                      aes(x = q, y = kurtosis,
                          colour   = source,
                          linetype = source,
                          shape    = source)) +
      geom_hline(yintercept = 3, linetype = "dashed",
                  colour = "grey50", inherit.aes = FALSE) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.8) +
      facet_wrap(~ prior, nrow = 1) +
      scale_x_continuous(breaks = q_levels, labels = q_levels) +
      scale_colour_manual(values   = c(Posterior = "#2166ac", GSCQF = "#d73027")) +
      scale_linetype_manual(values = c(Posterior = "solid",   GSCQF = "dashed")) +
      scale_shape_manual(values    = c(Posterior = 16,        GSCQF = 1)) +
      scale_y_continuous(expand = expansion(mult = 0.30)) +
      labs(x = "Quantile level", y = "Kurtosis") +
      theme_leg_corner

    fname_kurt <- file.path(DIR_PLOTS,
                            sprintf("kurt_group1_normal_n%d.pdf", n_kurt))
    ggsave(fname_kurt, p_kurt, width = 12, height = 4)
    cat(sprintf("  Saved: %s\n", basename(fname_kurt)))
  }
}

# =============================================================================
# Appendix summary tables – skewness & kurtosis
# 2 .tex files (table_skewness.tex, table_kurtosis.tex)
# Rows = GSCQF + 6 priors, grouped by error dist
# Cols = 4 sample sizes; each cell = [min; max] across quantile levels
# =============================================================================

if (FALSE) { # skewness/kurtosis tables — disabled
cat("\n=== APPENDIX TABLES ===\n\n")

dist_tex   <- c(normal = "Normal", exponential = "Exponential", t3 = "$t(3)$")
dists_tbl  <- c("normal", "exponential", "t3")
all_priors <- c(prior_groups$group1, prior_groups$group2)

# helper: format a cell as "[min; max]"
fmt_cell <- function(vals) {
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return("--")
  sprintf("[$%.3f$; $%.3f$]", min(vals), max(vals))
}

# --- Load ALL data across n and q into one big frame ---
all_post  <- data.frame()
all_gscqf <- data.frame()

for (n_tbl in N_SIZES_BAYES) {
  cat(sprintf("  Loading n = %d\n", n_tbl))

  for (grp_name in names(prior_groups)) {
    for (prior_name in prior_groups[[grp_name]]) {
      for (dist in dists_tbl) {
        for (q in q_levels) {
          fname <- find_bayes_file(n_tbl, dist, prior_name, q)
          if (is.null(fname)) next
          reps <- tryCatch(Filter(Negate(is.null), readRDS(fname)), error = function(e) NULL)
          if (is.null(reps) || length(reps) == 0) next
          skews <- sapply(reps, function(r) {
            d <- r$draws; if (is.null(d) || length(d) < 3) return(NA)
            d <- sample(d, min(N_SUB, length(d)), replace = FALSE)
            m <- mean(d, na.rm=TRUE); s <- sd(d, na.rm=TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((d - m) / s)^3, na.rm=TRUE)
          })
          kurts <- sapply(reps, function(r) {
            d <- r$draws; if (is.null(d) || length(d) < 3) return(NA)
            d <- sample(d, min(N_SUB, length(d)), replace = FALSE)
            m <- mean(d, na.rm=TRUE); s <- sd(d, na.rm=TRUE)
            if (is.na(s) || s == 0) return(NA)
            mean(((d - m) / s)^4, na.rm=TRUE)
          })
          all_post <- rbind(all_post, data.frame(
            n = n_tbl, prior = prior_name, dist = dist, q = q,
            sk = mean(skews, na.rm=TRUE),
            ku = mean(kurts, na.rm=TRUE),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  for (dist in dists_tbl) {
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_tbl, dist))
    if (!file.exists(fname_f)) next
    freq_res <- readRDS(fname_f)
    for (q in q_levels) {
      key  <- sprintf("q%.2f", q)
      reps <- Filter(Negate(is.null), freq_res[[key]])
      if (length(reps) == 0) next
      sg <- sapply(reps, function(r) {
        b <- r$boot_draws; if (is.null(b) || length(b) < 3) return(NA)
        m <- mean(b, na.rm=TRUE); s <- sd(b, na.rm=TRUE)
        if (is.na(s) || s == 0) return(NA)
        mean(((b - m) / s)^3, na.rm=TRUE)
      })
      kg <- sapply(reps, function(r) {
        b <- r$boot_draws; if (is.null(b) || length(b) < 3) return(NA)
        m <- mean(b, na.rm=TRUE); s <- sd(b, na.rm=TRUE)
        if (is.na(s) || s == 0) return(NA)
        mean(((b - m) / s)^4, na.rm=TRUE)
      })
      all_gscqf <- rbind(all_gscqf, data.frame(
        n = n_tbl, dist = dist, q = q,
        gsk = mean(sg, na.rm=TRUE),
        gku = mean(kg, na.rm=TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# --- Build tables ---
# Columns: Estimator | n=10 | n=30 | n=50 | n=250
n_labels  <- paste0("$n = ", N_SIZES_BAYES, "$")
ncols_tbl <- 1 + length(N_SIZES_BAYES)   # 5
col_spec  <- paste0("l", paste(rep("l", length(N_SIZES_BAYES)), collapse = ""))

header_line <- paste0(
  "Estimator & ",
  paste(n_labels, collapse = " & "),
  " \\\\"
)

for (metric in c("skewness", "kurtosis")) {
  tex <- c(
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    header_line,
    "\\midrule"
  )

  for (i in seq_along(dists_tbl)) {
    dist <- dists_tbl[i]

    tex <- c(tex, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                          ncols_tbl, dist_tex[dist]))

    # GSCQF row (first)
    gcells <- sapply(N_SIZES_BAYES, function(n_tbl) {
      sub <- all_gscqf[all_gscqf$n == n_tbl & all_gscqf$dist == dist, ]
      if (nrow(sub) == 0) return("--")
      fmt_cell(if (metric == "skewness") sub$gsk else sub$gku)
    })
    tex <- c(tex, paste0("GSCQF & ", paste(gcells, collapse = " & "), " \\\\"))

    # One row per prior
    for (pr in all_priors) {
      pr_tex <- gsub("([_#&%${}])", "\\\\\\1", pr)
      pcells <- sapply(N_SIZES_BAYES, function(n_tbl) {
        sub <- all_post[all_post$n == n_tbl & all_post$prior == pr &
                          all_post$dist == dist, ]
        if (nrow(sub) == 0) return("--")
        fmt_cell(if (metric == "skewness") sub$sk else sub$ku)
      })
      tex <- c(tex, paste0(pr_tex, " & ", paste(pcells, collapse = " & "), " \\\\"))
    }

    if (i < length(dists_tbl)) tex <- c(tex, "\\addlinespace[4pt]")
  }

  tex <- c(tex, "\\bottomrule", "\\end{tabular}")

  fname_tex <- file.path(DIR_TABLES, sprintf("table_%s.tex", metric))
  writeLines(tex, fname_tex)
  cat(sprintf("  Saved: %s\n", basename(fname_tex)))
}

cat("\nTables done.\n")

} # end if (FALSE) skewness/kurtosis tables

if (FALSE) { # mean estimator table — disabled
# =============================================================================
# Mean estimator table — all 9 priors × 4 sample sizes × 4 error dists
# Rows = 9 priors (priors_plot order) grouped by dist
# Cols = n=10, n=30, n=50, n=250; cell = [min; max] across quantile levels
# =============================================================================

cat("\n=== MEAN ESTIMATOR TABLE ===\n\n")

dists_mean <- c("normal", "exponential", "uniform", "t3")
dist_tex_mean <- c(normal = "Normal", exponential = "Exponential",
                   uniform = "Uniform", t3 = "$t(3)$")

# Load mean values for all combinations
mean_tbl <- data.frame()
for (n_tbl in N_SIZES_BAYES) {
  cat(sprintf("  Loading n = %d\n", n_tbl))
  for (prior_name in priors_plot) {
    for (dist in dists_mean) {
      for (q in q_levels) {
        fname <- find_bayes_file(n_tbl, dist, prior_name, q)
        if (is.null(fname)) next
        reps <- tryCatch(Filter(Negate(is.null), readRDS(fname)), error = function(e) NULL)
        if (is.null(reps) || length(reps) == 0) next
        vals <- sapply(reps, `[[`, "mean")
        mean_tbl <- rbind(mean_tbl, data.frame(
          n     = n_tbl,
          prior = prior_name,
          dist  = dist,
          q     = q,
          val   = mean(vals, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Build table
ncols_m     <- 1 + length(N_SIZES_BAYES)
col_spec_m  <- paste0("l", paste(rep("l", length(N_SIZES_BAYES)), collapse = ""))
n_labels_m  <- paste0("$n = ", N_SIZES_BAYES, "$")
header_m    <- paste0("Prior & ", paste(n_labels_m, collapse = " & "), " \\\\")

tex <- c(
  paste0("\\begin{tabular}{", col_spec_m, "}"),
  "\\toprule",
  header_m,
  "\\midrule"
)

for (i in seq_along(dists_mean)) {
  dist <- dists_mean[i]
  tex <- c(tex, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                        ncols_m, dist_tex_mean[dist]))
  for (pr in priors_plot) {
    pr_tex <- gsub("([_#&%${}])", "\\\\\\1", pr)
    cells <- sapply(N_SIZES_BAYES, function(n_tbl) {
      sub <- mean_tbl[mean_tbl$n == n_tbl & mean_tbl$prior == pr &
                        mean_tbl$dist == dist, ]
      if (nrow(sub) == 0) return("--")
      fmt_cell(sub$val)
    })
    tex <- c(tex, paste0(pr_tex, " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  if (i < length(dists_mean)) tex <- c(tex, "\\addlinespace[4pt]")
}

tex <- c(tex, "\\bottomrule", "\\end{tabular}")

fname_tex <- file.path(DIR_TABLES, "table_mean.tex")
writeLines(tex, fname_tex)
cat(sprintf("  Saved: %s\n", basename(fname_tex)))

cat("\nMean table done.\n")

} # end if (FALSE) mean estimator table

if (FALSE) { # comparison tables — disabled
# =============================================================================
# Comparison tables — Bias / SD / RMSE for n = 10, 30, 50, 250
# Load freq + bayes data, then write table_compare_{bias,sd,rmse}_n{n}.tex
# (BK row and Uniform distribution excluded)
# =============================================================================

# --- Table layout (rows = estimators, cols = q levels) ---
prior_names_ord <- setdiff(
  sapply(prior_specs, `[[`, "name"),
  c("N(0,1)", "Laplace(0,0.5)", "Laplace(0,1)", "Uniform(-1,1)")
)
error_dists_cmp <- setdiff(error_dists, "uniform")
dist_tex_cmp    <- c(normal = "Normal", exponential = "Exponential",
                     t3 = "$t(3)$")
prior_short_cmp <- c(
  "N(0,1)"          = "$\\mathcal{N}(0,1)$",
  "N(0,10)"         = "$\\mathcal{N}(0,10)$",
  "N(0,100)"        = "$\\mathcal{N}(0,100)$",
  "Laplace(0,0.5)"  = "Lap$(\\frac{1}{2})$",
  "Laplace(0,1)"    = "Lap$(1)$",
  "t(3,0,2.5)"      = "$t_3(2.5)$",
  "Cauchy(0,2.5)"   = "Ca$(2.5)$",
  "Uniform(-5,5)"   = "U$(-5,5)$",
  "Uniform(-10,10)" = "U$(-10,10)$"
)
col_count_x <- 1 + length(q_levels)
col_spec_x  <- paste0("l", paste(rep("r", length(q_levels)), collapse = ""))
header_x    <- paste0(
  "Estimator & ",
  paste(sprintf("$%.2f$", q_levels), collapse = " & "),
  " \\\\"
)

cat("\n=== WRITING COMPARISON TABLES (n = 10, 30, 50, 250) ===\n")
for (n_cmp in c(10, 30, 50, 250)) {
  cat(sprintf("\n--- n = %d ---\n", n_cmp))

  # --- Frequentist ---
  freq_cmp <- data.frame()
  for (dist_c in error_dists_cmp) {
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_cmp, dist_c))
    if (!file.exists(fname_f)) next
    freq_res <- readRDS(fname_f)
    for (q in q_levels) {
      key  <- sprintf("q%.2f", q)
      reps <- Filter(Negate(is.null), freq_res[[key]])
      if (length(reps) == 0) next
      b1_g  <- sapply(reps, `[[`, "beta1_gscqf")
      freq_cmp <- rbind(freq_cmp, data.frame(
        dist       = dist_c, q = q,
        bias_gscqf = mean(b1_g,  na.rm = TRUE) - 1,
        rmse_gscqf = sqrt(mean((b1_g  - 1)^2, na.rm = TRUE)),
        sd_gscqf   = sd(b1_g,  na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  # --- Bayesian (posterior mean only) ---
  bayes_cmp <- data.frame()
  for (dist_c in error_dists_cmp) {
    for (pr in prior_specs) {
      if (pr$name %in% c("N(0,1)", "Laplace(0,0.5)", "Laplace(0,1)", "Uniform(-1,1)")) next
      for (q in q_levels) {
        fname <- find_bayes_file(n_cmp, dist_c, pr$name, q)
        if (is.null(fname)) next
        reps <- Filter(Negate(is.null), readRDS(fname))
        if (length(reps) == 0) next
        means <- sapply(reps, `[[`, "mean")
        bayes_cmp <- rbind(bayes_cmp, data.frame(
          dist      = dist_c, prior = pr$name, q = q,
          bias_mean = mean(means, na.rm = TRUE) - 1,
          rmse_mean = sqrt(mean((means - 1)^2, na.rm = TRUE)),
          sd_mean   = sd(means, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  for (metric in c("bias", "sd", "rmse")) {
    fc_g <- paste0(metric, "_gscqf")
    bc   <- paste0(metric, "_mean")

    tex <- c(
      paste0("\\begin{tabular}{", col_spec_x, "}"),
      "\\toprule",
      header_x,
      "\\midrule"
    )
    for (i in seq_along(error_dists_cmp)) {
      dist <- error_dists_cmp[i]
      tex <- c(tex, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\",
                            col_count_x, dist_tex_cmp[dist]))

      # GSCQF row
      g_vals <- sapply(q_levels, function(q) {
        frow <- freq_cmp[freq_cmp$dist == dist & freq_cmp$q == q, ]
        if (nrow(frow) > 0) sprintf("$%.4f$", frow[[fc_g]][1]) else "--"
      })
      tex <- c(tex, paste0("GSCQF & ", paste(g_vals, collapse = " & "), " \\\\"))

      # One row per prior
      for (pr in prior_names_ord) {
        pr_vals <- sapply(q_levels, function(q) {
          brow <- bayes_cmp[bayes_cmp$dist == dist & bayes_cmp$prior == pr &
                              bayes_cmp$q == q, ]
          if (nrow(brow) == 0) return("--")
          sprintf("$%.4f$", brow[[bc]][1])
        })
        tex <- c(tex, paste0(
          prior_short_cmp[pr], " & ", paste(pr_vals, collapse = " & "), " \\\\"
        ))
      }
      if (i < length(error_dists_cmp)) tex <- c(tex, "\\addlinespace[4pt]")
    }
    tex <- c(tex, "\\bottomrule", "\\end{tabular}")
    fname_tex <- file.path(DIR_TABLES,
                           sprintf("table_compare_%s_n%d.tex", metric, n_cmp))
    writeLines(tex, fname_tex)
    cat(sprintf("  Saved: %s\n", basename(fname_tex)))
  }
}
cat("\nComparison tables done.\n")

} # end if (FALSE) comparison tables

if (FALSE) { # per-distribution tables — disabled
# =============================================================================
# Per-distribution tables: one .tex per error dist
# Rows: metric (RMSE/Bias/SD) > n (10/30/50/250) > estimators
# Cols: 7 q levels
# Output: table_dist_normal.tex, table_dist_exponential.tex, table_dist_t3.tex
# =============================================================================
cat("\n=== PER-DISTRIBUTION TABLES ===\n")

metric_labels <- c(rmse = "RMSE", bias = "Bias", sd = "SD")
n_sizes_dist  <- c(10, 30, 50, 250)

for (dist in error_dists_cmp) {

  # Load all data for this distribution
  freq_d <- data.frame()
  for (n_d in n_sizes_dist) {
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_d, dist))
    if (!file.exists(fname_f)) next
    freq_res <- readRDS(fname_f)
    for (q in q_levels) {
      key  <- sprintf("q%.2f", q)
      reps <- Filter(Negate(is.null), freq_res[[key]])
      if (length(reps) == 0) next
      b1_g <- sapply(reps, `[[`, "beta1_gscqf")
      freq_d <- rbind(freq_d, data.frame(
        n = n_d, q = q,
        bias_gscqf = mean(b1_g, na.rm = TRUE) - 1,
        rmse_gscqf = sqrt(mean((b1_g - 1)^2, na.rm = TRUE)),
        sd_gscqf   = sd(b1_g, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  bayes_d <- data.frame()
  for (n_d in n_sizes_dist) {
    for (pr in prior_specs) {
      if (pr$name %in% c("N(0,1)", "Laplace(0,0.5)", "Laplace(0,1)", "Uniform(-1,1)")) next
      for (q in q_levels) {
        fname <- find_bayes_file(n_d, dist, pr$name, q)
        if (is.null(fname)) next
        reps <- Filter(Negate(is.null), readRDS(fname))
        if (length(reps) == 0) next
        means <- sapply(reps, `[[`, "mean")
        bayes_d <- rbind(bayes_d, data.frame(
          n = n_d, prior = pr$name, q = q,
          bias_mean = mean(means, na.rm = TRUE) - 1,
          rmse_mean = sqrt(mean((means - 1)^2, na.rm = TRUE)),
          sd_mean   = sd(means, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # One table per n value
  for (n_d in n_sizes_dist) {
    tex <- c(
      paste0("\\begin{tabular}{", col_spec_x, "}"),
      "\\toprule",
      header_x,
      "\\midrule"
    )

    first_metric <- TRUE
    for (metric in c("rmse", "bias", "sd")) {
      fc_g <- paste0(metric, "_gscqf")
      bc   <- paste0(metric, "_mean")

      if (!first_metric) tex <- c(tex, "\\addlinespace[6pt]")
      first_metric <- FALSE

      tex <- c(tex, sprintf("\\multicolumn{%d}{l}{\\textbf{%s}} \\\\",
                            col_count_x, metric_labels[metric]))

      # GSCQF row
      g_vals <- sapply(q_levels, function(q) {
        frow <- freq_d[freq_d$n == n_d & freq_d$q == q, ]
        if (nrow(frow) > 0) sprintf("$%.4f$", frow[[fc_g]][1]) else "--"
      })
      tex <- c(tex, paste0("GSCQF & ", paste(g_vals, collapse = " & "), " \\\\"))

      # Prior rows
      for (pr in prior_names_ord) {
        pr_vals <- sapply(q_levels, function(q) {
          brow <- bayes_d[bayes_d$n == n_d & bayes_d$prior == pr & bayes_d$q == q, ]
          if (nrow(brow) == 0) return("--")
          sprintf("$%.4f$", brow[[bc]][1])
        })
        tex <- c(tex, paste0(prior_short_cmp[pr], " & ",
                             paste(pr_vals, collapse = " & "), " \\\\"))
      }
    }

    tex <- c(tex, "\\bottomrule", "\\end{tabular}")
    fname_tex <- file.path(DIR_TABLES, sprintf("table_dist_%s_n%d.tex", dist, n_d))
    writeLines(tex, fname_tex)
    cat(sprintf("  Saved: %s\n", basename(fname_tex)))
  }
}
cat("\nPer-distribution tables done.\n")
} # end if (FALSE) per-distribution tables

# =============================================================================
# Coverage tables — one per error distribution
# Rows: n (section) > GSCQF + all 9 priors; Cols: 7 q levels
# Output: table_coverage_{normal,exponential,t3}.tex
# =============================================================================
cat("\n=== COVERAGE TABLES ===\n")

cov_dists    <- setdiff(error_dists, "uniform")   # normal, exponential, t3
n_sizes_cov  <- c(10, 30, 50, 250)
col_count_cv <- 1 + length(q_levels)
col_spec_cv  <- paste0("l", paste(rep("r", length(q_levels)), collapse = ""))
header_cv    <- paste0(
  "Estimator & ",
  paste(sprintf("$%.2f$", q_levels), collapse = " & "),
  " \\\\"
)
all_priors_cov <- c(priors_plot, "Uniform(-1,1)")
prior_lbl_cv <- c(
  "N(0,1)"          = "$\\mathcal{N}(0,1)$",
  "N(0,10)"         = "$\\mathcal{N}(0,10)$",
  "N(0,100)"        = "$\\mathcal{N}(0,100)$",
  "Laplace(0,0.5)"  = "Lap$(0,\\frac{1}{2})$",
  "Laplace(0,1)"    = "Lap$(0,1)$",
  "t(3,0,2.5)"      = "$t_3(0,2.5)$",
  "Cauchy(0,2.5)"   = "Cauchy$(0,2.5)$",
  "Uniform(-5,5)"   = "U$(-5,5)$",
  "Uniform(-10,10)" = "U$(-10,10)$",
  "Uniform(-1,1)"   = "U$(-1,1)$"
)

for (dist in cov_dists) {
  cat(sprintf("  Loading %s...\n", dist))

  freq_cov <- data.frame()
  for (n_d in n_sizes_cov) {
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n_d, dist))
    if (!file.exists(fname_f)) next
    freq_res <- readRDS(fname_f)
    for (q in q_levels) {
      key  <- sprintf("q%.2f", q)
      reps <- Filter(Negate(is.null), freq_res[[key]])
      if (length(reps) == 0) next
      cov_g <- sapply(reps, `[[`, "cover_boot_gscqf")
      freq_cov <- rbind(freq_cov, data.frame(
        n = n_d, q = q,
        coverage = mean(cov_g, na.rm = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  bayes_cov <- data.frame()
  for (n_d in n_sizes_cov) {
    for (pr in all_priors_cov) {
      for (q in q_levels) {
        fname <- find_bayes_file(n_d, dist, pr, q)
        if (is.null(fname)) next
        reps <- Filter(Negate(is.null), readRDS(fname))
        if (length(reps) == 0) next
        covers <- sapply(reps, `[[`, "cover")
        bayes_cov <- rbind(bayes_cov, data.frame(
          n = n_d, prior = pr, q = q,
          coverage = mean(covers, na.rm = TRUE),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  tex <- c(
    paste0("\\begin{tabular}{", col_spec_cv, "}"),
    "\\toprule",
    header_cv,
    "\\midrule"
  )

  for (i_n in seq_along(n_sizes_cov)) {
    n_d <- n_sizes_cov[i_n]
    if (i_n > 1) tex <- c(tex, "\\addlinespace[6pt]")
    tex <- c(tex, sprintf("\\multicolumn{%d}{l}{\\textbf{$n = %d$}} \\\\",
                          col_count_cv, n_d))

    g_vals <- sapply(q_levels, function(q) {
      frow <- freq_cov[freq_cov$n == n_d & freq_cov$q == q, ]
      if (nrow(frow) > 0) sprintf("$%.4f$", frow$coverage[1]) else "--"
    })
    tex <- c(tex, paste0("GSCQF & ", paste(g_vals, collapse = " & "), " \\\\"))

    for (pr in all_priors_cov) {
      pr_vals <- sapply(q_levels, function(q) {
        brow <- bayes_cov[bayes_cov$n == n_d & bayes_cov$prior == pr & bayes_cov$q == q, ]
        if (nrow(brow) == 0) return("--")
        sprintf("$%.4f$", brow$coverage[1])
      })
      tex <- c(tex, paste0(prior_lbl_cv[pr], " & ",
                           paste(pr_vals, collapse = " & "), " \\\\"))
    }
  }

  tex <- c(tex, "\\bottomrule", "\\end{tabular}")
  fname_tex <- file.path(DIR_TABLES, sprintf("table_coverage_%s.tex", dist))
  writeLines(tex, fname_tex)
  cat(sprintf("  Saved: %s\n", basename(fname_tex)))
}
cat("\nCoverage tables done.\n")

# =============================================================================
# Coverage plots: one PDF per (dist, n)
# Layout: 2 rows x 3 cols — 6 priors matching the two prior_groups
#   Row 1: N(0,10), t(3,0,2.5), Uniform(-10,10)
#   Row 2: N(0,100), Cauchy(0,2.5), Uniform(-5,5)
# X-axis: q level, Y-axis: coverage
# Lines: Posterior (solid blue) + GSCQF (dashed red)
# Grey dashed diagonal = perfect calibration (y = x)
# =============================================================================

cat("\n=== COVERAGE PLOTS ===\n\n")

if (FALSE) {

cov_priors <- c(
  "N(0,10)",        "t(3,0,2.5)",    "Uniform(-10,10)",
  "N(0,100)",       "Cauchy(0,2.5)", "Uniform(-5,5)"
)

cov_prior_labels <- c(
  "N(0,10)"          = "N(0, 10)",
  "N(0,100)"         = "N(0, 100)",
  "t(3,0,2.5)"       = "t(3, 0, 2.5)",
  "Cauchy(0,2.5)"    = "Cauchy(0, 2.5)",
  "Uniform(-5,5)"    = "Uniform(-5, 5)",
  "Uniform(-10,10)"  = "Uniform(-10, 10)"
)

for (dist in error_dists_cmp) {
  for (n in N_SIZES_BAYES) {
    cat(sprintf("  Coverage plot: dist=%s  n=%d\n", dist, n))

    rows <- list()

    # --- Bayesian coverage per prior ---
    for (prior_name in cov_priors) {
      for (q in q_levels) {
        fname <- find_bayes_file(n, dist, prior_name, q)
        if (is.null(fname)) next
        reps <- tryCatch(
          Filter(Negate(is.null), readRDS(fname)),
          error = function(e) NULL
        )
        if (is.null(reps) || length(reps) == 0) next
        covers <- sapply(reps, `[[`, "cover")
        rows[[length(rows) + 1]] <- data.frame(
          prior    = prior_name,
          q        = q,
          source   = "Posterior",
          coverage = mean(covers, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }

    # --- GSCQF coverage ---
    fname_f <- file.path(DIR_FREQ, sprintf("freq_n%d_%s.rds", n, dist))
    if (file.exists(fname_f)) {
      freq_res <- readRDS(fname_f)
      for (q in q_levels) {
        key  <- sprintf("q%.2f", q)
        reps <- Filter(Negate(is.null), freq_res[[key]])
        if (length(reps) == 0) next
        cov_g <- mean(sapply(reps, `[[`, "cover_boot_gscqf"), na.rm = TRUE)
        for (prior_name in cov_priors) {
          rows[[length(rows) + 1]] <- data.frame(
            prior    = prior_name,
            q        = q,
            source   = "GSCQF",
            coverage = cov_g,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(rows) == 0) { cat("  [SKIP]\n"); next }

    df        <- do.call(rbind, rows)
    df$prior  <- factor(df$prior,  levels = cov_priors)
    df$source <- factor(df$source, levels = c("Posterior", "GSCQF"))

    p <- ggplot(df, aes(x = q, y = coverage, colour = source, linetype = source)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      facet_wrap(~ prior, nrow = 2, ncol = 3,
                 labeller = as_labeller(cov_prior_labels)) +
      scale_x_continuous(breaks = q_levels, labels = q_levels) +
      scale_y_continuous(limits = c(NA, 1.02), expand = expansion(mult = 0.05)) +
      scale_colour_manual(values = c(Posterior = "#2166ac", GSCQF = "#d73027")) +
      scale_linetype_manual(values = c(Posterior = "solid", GSCQF = "dashed")) +
      labs(x = "Quantile level", y = "Coverage") +
      theme_bw(base_size = 10) +
      theme(
        strip.text           = element_text(size = 9),
        axis.text.x          = element_text(angle = 45, hjust = 1, size = 7),
        legend.position      = c(0.99, 0.01),
        legend.justification = c("right", "bottom"),
        legend.background    = element_rect(fill      = alpha("white", 0.85),
                                            colour    = "grey70",
                                            linewidth = 0.3),
        legend.key.size      = unit(0.55, "cm"),
        legend.text          = element_text(size = 9),
        legend.title         = element_blank()
      )

    fname_out <- file.path(DIR_PLOTS, sprintf("coverage_%s_n%d.pdf", dist, n))
    ggsave(fname_out, p, width = 12, height = 6)
    cat(sprintf("  Saved: %s\n", basename(fname_out)))
  }
}

} # end if (MAKE_PLOTS) coverage plots

# =============================================================================
# Combined n=10 vs n=250 — mean/median/mode — normal distribution
# 3x3 grid of 9 priors; 6 lines per panel, each with a unique colour,
# linetype AND point shape:
#   Mean   n=10  | Mean   n=250
#   Median n=10  | Median n=250
#   Mode   n=10  | Mode   n=250
# =============================================================================

cat("\n=== COMBINED n=10 vs n=250 MEAN/MEDIAN/MODE (normal) ===\n\n")

dist_combined <- "normal"
ns_combined   <- c(10, 250)

rows_combined <- list()

for (n_c in ns_combined) {
  for (prior_name in priors_plot) {
    for (q in q_levels) {
      fname <- find_bayes_file(n_c, dist_combined, prior_name, q)
      if (is.null(fname)) next
      reps <- tryCatch(
        Filter(Negate(is.null), readRDS(fname)),
        error = function(e) NULL
      )
      if (is.null(reps) || length(reps) == 0) next

      means   <- sapply(reps, `[[`, "mean")
      medians <- sapply(reps, `[[`, "median")
      modes   <- sapply(reps, `[[`, "mode")

      for (est_name in c("Mean", "Median", "Mode")) {
        vals <- switch(est_name,
                       Mean   = means,
                       Median = medians,
                       Mode   = modes)
        rows_combined[[length(rows_combined) + 1]] <- data.frame(
          prior  = prior_name,
          q      = q,
          group  = paste0(est_name, " (n=", n_c, ")"),
          value  = mean(vals, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(rows_combined) > 0) {
  df_c <- do.call(rbind, rows_combined)
  df_c$prior <- factor(df_c$prior, levels = priors_plot)

  group_levels <- c(
    "Mean (n=10)",   "Median (n=10)",   "Mode (n=10)",
    "Mean (n=250)",  "Median (n=250)",  "Mode (n=250)"
  )
  df_c$group <- factor(df_c$group, levels = group_levels)

  # 6 distinct colours: paired shades per estimator
  grp_colours <- c(
    "Mean (n=10)"   = "#08519c",   # dark blue
    "Mean (n=250)"  = "#6baed6",   # light blue
    "Median (n=10)" = "#a50f15",   # dark red
    "Median (n=250)"= "#fc8d59",   # orange-red
    "Mode (n=10)"   = "#006d2c",   # dark green
    "Mode (n=250)"  = "#74c476"    # light green
  )

  # 6 distinct linetypes
  grp_ltypes <- c(
    "Mean (n=10)"   = "solid",
    "Mean (n=250)"  = "dashed",
    "Median (n=10)" = "dotted",
    "Median (n=250)"= "dotdash",
    "Mode (n=10)"   = "longdash",
    "Mode (n=250)"  = "twodash"
  )

  # 6 distinct point shapes
  grp_shapes <- c(
    "Mean (n=10)"   = 16,   # filled circle
    "Mean (n=250)"  = 1,    # open circle
    "Median (n=10)" = 17,   # filled triangle
    "Median (n=250)"= 2,    # open triangle
    "Mode (n=10)"   = 15,   # filled square
    "Mode (n=250)"  = 0     # open square
  )

  p_combined <- ggplot(df_c,
                       aes(x        = q,
                           y        = value,
                           colour   = group,
                           linetype = group,
                           shape    = group)) +
    geom_hline(yintercept = 1, linetype = "longdash",
               colour = "black", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ prior, ncol = 3, scales = "fixed") +
    scale_colour_manual(values = grp_colours,   name = NULL) +
    scale_linetype_manual(values = grp_ltypes,  name = NULL) +
    scale_shape_manual(values = grp_shapes,     name = NULL) +
    guides(colour   = guide_legend(ncol = 2),
           linetype = guide_legend(ncol = 2),
           shape    = guide_legend(ncol = 2)) +
    scale_x_continuous(breaks = q_levels, labels = q_levels) +
    labs(
      x     = "Quantile level",
      y     = expression(hat(beta)[1])
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text           = element_text(size = 8),
      legend.background    = element_rect(fill      = alpha("white", 0.85),
                                          colour    = "grey70",
                                          linewidth = 0.3),
      legend.position      = c(0.99, 0.70),
      legend.justification = c("right", "bottom"),
      legend.key.size      = unit(0.7, "cm"),
      legend.key.width     = unit(1.2, "cm"),
      legend.text          = element_text(size = 9),
      axis.text.x          = element_text(angle = 45, hjust = 1, size = 7)
    )

  fname_out <- file.path(DIR_PLOTS,
                         "mean_median_mode_9priors_normal_n10_n250_combined.pdf")
  ggsave(fname_out, p_combined, width = 12, height = 9)
  cat(sprintf("  Saved: %s\n", basename(fname_out)))
} else {
  cat("  [SKIP] No data found.\n")
}

cat("\nAll done.\n")
