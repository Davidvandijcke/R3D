#!/usr/bin/env Rscript
#
# 03_mc_compare.R
# Compare R vs Stata Monte Carlo results.
# Produces summary tables and LaTeX output.
#
# Usage: Rscript 03_mc_compare.R

cat("=== R3D Monte Carlo Comparison: R vs Stata ===\n\n")

outdir <- file.path(getwd(), "output")
stata_dir <- file.path(outdir, "results_Stata")
r_dir <- file.path(outdir, "results_R")
comp_dir <- file.path(outdir, "comparison")
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

# Load R results
R_all <- readRDS(file.path(r_dir, "mc_R_all.rds"))

NQ <- 9
Q_GRID <- (1:9) / 10
N_STATA <- 100

# ============================================================================
# CONFIGURATION: cells to compare
# ============================================================================
dgps <- c("dgp1", "dgp2")
sample_sizes <- c(200, 500)
deltas <- c(0, 2)
methods <- c("simple", "frechet")

# Pre-compute DGP2 true tau
true_tau_dgp2 <- list()
for (delta in c(0, 1, 2)) {
  set.seed(99999)
  large_n <- 1e6
  y_above <- rnorm(large_n, mean = delta, sd = 1) + 2 * rexp(large_n, rate = 1)
  y_below <- rnorm(large_n, mean = 0, sd = 1) + 2 * rexp(large_n, rate = 1)
  true_tau_dgp2[[as.character(delta)]] <-
    quantile(y_above, probs = Q_GRID, names = FALSE) -
    quantile(y_below, probs = Q_GRID, names = FALSE)
}

# ============================================================================
# A. DIRECT AGREEMENT (shared bandwidths, same xi)
# ============================================================================
cat("--- A. Direct Agreement (shared bw) ---\n")

agreement_rows <- list()

for (dgp in dgps) {
  for (n_obs in sample_sizes) {
    for (delta in deltas) {
      for (method in methods) {
        cell <- sprintf("%s_sharp_%d_d%s_%s", dgp, n_obs, delta, method)
        r_key <- cell

        # Check R results exist
        if (!(r_key %in% names(R_all))) {
          cat("  Skipping", cell, "- no R results\n")
          next
        }

        tau_diffs <- c()
        cb_diffs <- c()
        pval_null_diffs <- c()
        pval_homo_diffs <- c()
        n_matched <- 0

        for (s in 1:N_STATA) {
          # Load Stata shared-bw result
          stata_file <- file.path(stata_dir,
                                  sprintf("mc_res_%s_sharedbw_%d.csv", cell, s))
          if (!file.exists(stata_file)) next

          stata_res <- read.csv(stata_file)
          if (nrow(stata_res) != NQ) next

          # Load R per-rep result
          r_file <- file.path(r_dir, sprintf("mc_res_%s_%d.csv", r_key, s))
          if (!file.exists(r_file)) next

          r_res <- read.csv(r_file)
          if (nrow(r_res) != NQ) next

          n_matched <- n_matched + 1

          # Tau difference
          tau_diff <- max(abs(r_res$tau - stata_res$tau), na.rm = TRUE)
          tau_diffs <- c(tau_diffs, tau_diff)

          # CB difference
          cb_diff_lo <- max(abs(r_res$cb_lower - stata_res$cb_lower), na.rm = TRUE)
          cb_diff_hi <- max(abs(r_res$cb_upper - stata_res$cb_upper), na.rm = TRUE)
          cb_diffs <- c(cb_diffs, max(cb_diff_lo, cb_diff_hi))

          # P-value differences
          stata_pv_file <- file.path(stata_dir,
                                     sprintf("mc_pval_%s_sharedbw_%d.csv", cell, s))
          r_pv <- R_all[[r_key]]
          if (file.exists(stata_pv_file) && !is.null(r_pv$pval_nullity)) {
            stata_pv <- read.csv(stata_pv_file)
            if ("pval_nullity" %in% names(stata_pv) && !is.na(r_pv$pval_nullity[s])) {
              pval_null_diffs <- c(pval_null_diffs,
                                   abs(r_pv$pval_nullity[s] - stata_pv$pval_nullity[1]))
            }
            if ("pval_homogeneity" %in% names(stata_pv) &&
                !is.null(r_pv$pval_homogeneity) && !is.na(r_pv$pval_homogeneity[s])) {
              pval_homo_diffs <- c(pval_homo_diffs,
                                   abs(r_pv$pval_homogeneity[s] - stata_pv$pval_homogeneity[1]))
            }
          }
        }

        if (n_matched > 0) {
          agreement_rows[[length(agreement_rows) + 1]] <- data.frame(
            cell = cell, dgp = dgp, method = method, n = n_obs, delta = delta,
            n_matched = n_matched,
            tau_diff_mean = mean(tau_diffs),
            tau_diff_median = median(tau_diffs),
            tau_diff_p95 = quantile(tau_diffs, 0.95),
            cb_diff_mean = mean(cb_diffs),
            cb_diff_median = median(cb_diffs),
            cb_diff_p95 = quantile(cb_diffs, 0.95),
            pval_null_diff_mean = if (length(pval_null_diffs) > 0) mean(pval_null_diffs) else NA,
            pval_homo_diff_mean = if (length(pval_homo_diffs) > 0) mean(pval_homo_diffs) else NA,
            stringsAsFactors = FALSE, row.names = NULL
          )
          cat(sprintf("  %s: %d matched, tau_diff median=%.6f, p95=%.6f\n",
                      cell, n_matched, median(tau_diffs), quantile(tau_diffs, 0.95)))
        }
      }
    }
  }
}

agreement_df <- if (length(agreement_rows) > 0) do.call(rbind, agreement_rows) else data.frame()
if (nrow(agreement_df) > 0) {
  write.csv(agreement_df, file.path(comp_dir, "mc_agreement.csv"), row.names = FALSE)
}

# ============================================================================
# B-E. COVERAGE, BIAS, SIZE, POWER
# ============================================================================
cat("\n--- B-E. Coverage, Bias, Size, Power ---\n")

summary_rows <- list()

for (dgp in dgps) {
  for (n_obs in sample_sizes) {
    for (delta in c(0, 1, 2)) {
      for (method in methods) {
        cell <- sprintf("%s_sharp_%d_d%s_%s", dgp, n_obs, delta, method)

        # True tau
        if (dgp == "dgp1") {
          true_tau <- rep(delta, NQ)
        } else {
          true_tau <- true_tau_dgp2[[as.character(delta)]]
        }

        # --- R results (full N_SIM) ---
        r_key <- cell
        r_res <- R_all[[r_key]]
        if (!is.null(r_res)) {
          valid_r <- !is.na(r_res$tau_mat[, 1])
          n_r <- sum(valid_r)

          if (n_r > 0) {
            bias_r <- colMeans(r_res$tau_mat[valid_r, , drop = FALSE]) - true_tau
            rmse_r <- sqrt(colMeans((r_res$tau_mat[valid_r, , drop = FALSE] -
                                      matrix(true_tau, nrow = n_r, ncol = NQ, byrow = TRUE))^2))
            cov_r <- mean(apply(
              r_res$cb_lower_mat[valid_r, , drop = FALSE] <=
                matrix(true_tau, nrow = n_r, ncol = NQ, byrow = TRUE) &
              r_res$cb_upper_mat[valid_r, , drop = FALSE] >=
                matrix(true_tau, nrow = n_r, ncol = NQ, byrow = TRUE),
              1, all), na.rm = TRUE)
            rej_null_r <- mean(r_res$pval_nullity[valid_r] < 0.05, na.rm = TRUE)
            rej_homo_r <- if (!is.null(r_res$pval_homogeneity))
              mean(r_res$pval_homogeneity[valid_r] < 0.05, na.rm = TRUE) else NA
          }
        } else {
          n_r <- 0
        }

        # --- Stata results (shared-bw, up to N_STATA sims) ---
        stata_taus <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_cb_lo <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_cb_hi <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_pn <- rep(NA, N_STATA)
        stata_ph <- rep(NA, N_STATA)

        for (s in 1:N_STATA) {
          f <- file.path(stata_dir, sprintf("mc_res_%s_sharedbw_%d.csv", cell, s))
          if (!file.exists(f)) next
          d <- read.csv(f)
          if (nrow(d) != NQ) next
          stata_taus[s, ] <- d$tau
          stata_cb_lo[s, ] <- d$cb_lower
          stata_cb_hi[s, ] <- d$cb_upper

          pf <- file.path(stata_dir, sprintf("mc_pval_%s_sharedbw_%d.csv", cell, s))
          if (file.exists(pf)) {
            pd <- read.csv(pf)
            if ("pval_nullity" %in% names(pd)) stata_pn[s] <- pd$pval_nullity[1]
            if ("pval_homogeneity" %in% names(pd)) stata_ph[s] <- pd$pval_homogeneity[1]
          }
        }

        valid_s <- !is.na(stata_taus[, 1])
        n_s <- sum(valid_s)

        # --- Stata own-bw results ---
        stata_own_taus <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_own_cb_lo <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_own_cb_hi <- matrix(NA, nrow = N_STATA, ncol = NQ)
        stata_own_pn <- rep(NA, N_STATA)

        for (s in 1:N_STATA) {
          f <- file.path(stata_dir, sprintf("mc_res_%s_ownbw_%d.csv", cell, s))
          if (!file.exists(f)) next
          d <- read.csv(f)
          if (nrow(d) != NQ) next
          stata_own_taus[s, ] <- d$tau
          stata_own_cb_lo[s, ] <- d$cb_lower
          stata_own_cb_hi[s, ] <- d$cb_upper

          pf <- file.path(stata_dir, sprintf("mc_pval_%s_ownbw_%d.csv", cell, s))
          if (file.exists(pf)) {
            pd <- read.csv(pf)
            if ("pval_nullity" %in% names(pd)) stata_own_pn[s] <- pd$pval_nullity[1]
          }
        }

        valid_so <- !is.na(stata_own_taus[, 1])
        n_so <- sum(valid_so)

        # Compute metrics for Stata shared-bw
        if (n_s > 0) {
          bias_s <- colMeans(stata_taus[valid_s, , drop = FALSE]) - true_tau
          rmse_s <- sqrt(colMeans((stata_taus[valid_s, , drop = FALSE] -
                                    matrix(true_tau, nrow = n_s, ncol = NQ, byrow = TRUE))^2))
          cov_s <- mean(apply(
            stata_cb_lo[valid_s, , drop = FALSE] <=
              matrix(true_tau, nrow = n_s, ncol = NQ, byrow = TRUE) &
            stata_cb_hi[valid_s, , drop = FALSE] >=
              matrix(true_tau, nrow = n_s, ncol = NQ, byrow = TRUE),
            1, all), na.rm = TRUE)
          rej_null_s <- mean(stata_pn[valid_s] < 0.05, na.rm = TRUE)
          rej_homo_s <- mean(stata_ph[valid_s] < 0.05, na.rm = TRUE)
        }

        # Compute metrics for Stata own-bw
        if (n_so > 0) {
          bias_so <- colMeans(stata_own_taus[valid_so, , drop = FALSE]) - true_tau
          cov_so <- mean(apply(
            stata_own_cb_lo[valid_so, , drop = FALSE] <=
              matrix(true_tau, nrow = n_so, ncol = NQ, byrow = TRUE) &
            stata_own_cb_hi[valid_so, , drop = FALSE] >=
              matrix(true_tau, nrow = n_so, ncol = NQ, byrow = TRUE),
            1, all), na.rm = TRUE)
          rej_null_so <- mean(stata_own_pn[valid_so] < 0.05, na.rm = TRUE)
        }

        row <- data.frame(
          cell = cell, dgp = dgp, method = method, n = n_obs, delta = delta,
          # R metrics
          n_R = n_r,
          bias_R = if (n_r > 0) mean(abs(bias_r)) else NA,
          rmse_R = if (n_r > 0) mean(rmse_r) else NA,
          coverage_R = if (n_r > 0) cov_r else NA,
          rej_null_R = if (n_r > 0) rej_null_r else NA,
          rej_homo_R = if (n_r > 0) rej_homo_r else NA,
          # Stata shared-bw metrics
          n_Stata_shared = n_s,
          bias_Stata_shared = if (n_s > 0) mean(abs(bias_s)) else NA,
          rmse_Stata_shared = if (n_s > 0) mean(rmse_s) else NA,
          coverage_Stata_shared = if (n_s > 0) cov_s else NA,
          rej_null_Stata_shared = if (n_s > 0) rej_null_s else NA,
          rej_homo_Stata_shared = if (n_s > 0) rej_homo_s else NA,
          # Stata own-bw metrics
          n_Stata_own = n_so,
          bias_Stata_own = if (n_so > 0) mean(abs(bias_so)) else NA,
          coverage_Stata_own = if (n_so > 0) cov_so else NA,
          rej_null_Stata_own = if (n_so > 0) rej_null_so else NA,
          stringsAsFactors = FALSE, row.names = NULL
        )

        summary_rows[[length(summary_rows) + 1]] <- row

        if (n_r > 0 || n_s > 0) {
          cat(sprintf("  %s: R(%d) ucov=%.3f rej=%.3f | Stata-shared(%d) ucov=%.3f rej=%.3f | Stata-own(%d) ucov=%.3f rej=%.3f\n",
                      cell, n_r,
                      if (n_r > 0) cov_r else NA,
                      if (n_r > 0) rej_null_r else NA,
                      n_s,
                      if (n_s > 0) cov_s else NA,
                      if (n_s > 0) rej_null_s else NA,
                      n_so,
                      if (n_so > 0) cov_so else NA,
                      if (n_so > 0) rej_null_so else NA))
        }
      }
    }
  }
}

summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
if (nrow(summary_df) > 0) {
  write.csv(summary_df, file.path(comp_dir, "mc_summary.csv"), row.names = FALSE)
}

# ============================================================================
# LaTeX TABLE
# ============================================================================
cat("\n--- Generating LaTeX tables ---\n")

if (nrow(summary_df) > 0) {
  tex_file <- file.path(comp_dir, "mc_tables.tex")
  sink(tex_file)

  cat("\\begin{table}[htbp]\n")
  cat("\\centering\n")
  cat("\\caption{Monte Carlo Comparison: R vs Stata (shared bandwidths)}\n")
  cat("\\label{tab:mc_comparison}\n")
  cat("\\small\n")
  cat("\\begin{tabular}{llcccccccc}\n")
  cat("\\hline\\hline\n")
  cat("& & \\multicolumn{4}{c}{R} & \\multicolumn{4}{c}{Stata (shared bw)} \\\\\n")
  cat("\\cmidrule(lr){3-6} \\cmidrule(lr){7-10}\n")
  cat("Cell & $n$ & $|$Bias$|$ & RMSE & Cov. & Rej. & $|$Bias$|$ & RMSE & Cov. & Rej. \\\\\n")
  cat("\\hline\n")

  for (i in seq_len(nrow(summary_df))) {
    r <- summary_df[i, ]
    # Choose the relevant rejection rate based on delta
    rej_label <- if (r$delta == 0) "Size" else "Power"
    cat(sprintf("%s & %d & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n",
                gsub("_", "\\\\_", r$cell), r$n,
                ifelse(is.na(r$bias_R), 0, r$bias_R),
                ifelse(is.na(r$rmse_R), 0, r$rmse_R),
                ifelse(is.na(r$coverage_R), 0, r$coverage_R),
                ifelse(is.na(r$rej_null_R), 0, r$rej_null_R),
                ifelse(is.na(r$bias_Stata_shared), 0, r$bias_Stata_shared),
                ifelse(is.na(r$rmse_Stata_shared), 0, r$rmse_Stata_shared),
                ifelse(is.na(r$coverage_Stata_shared), 0, r$coverage_Stata_shared),
                ifelse(is.na(r$rej_null_Stata_shared), 0, r$rej_null_Stata_shared)))
  }

  cat("\\hline\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\begin{tablenotes}\\small\n")
  cat("\\item Notes: Cov. = uniform coverage at $\\alpha=0.05$. Rej. = rejection rate of nullity test at $\\alpha=0.05$ (size when $\\Delta=0$, power otherwise). Shared bw = Stata uses R bandwidths.\n")
  cat("\\end{tablenotes}\n")
  cat("\\end{table}\n")

  sink()
  cat("  Saved:", tex_file, "\n")
}

# ============================================================================
# PASS/FAIL ASSESSMENT
# ============================================================================
cat("\n=== PASS/FAIL Assessment ===\n")

if (nrow(agreement_df) > 0) {
  # Direct agreement thresholds
  frechet_ok <- all(agreement_df$tau_diff_p95[agreement_df$method == "frechet"] < 0.01,
                     na.rm = TRUE)
  simple_ok <- all(agreement_df$tau_diff_p95[agreement_df$method == "simple"] < 0.05,
                    na.rm = TRUE)

  cat(sprintf("  Frechet tau agreement (p95 < 0.01): %s\n",
              if (frechet_ok) "PASS" else "FAIL"))
  cat(sprintf("  Simple tau agreement (p95 < 0.05): %s\n",
              if (simple_ok) "PASS" else "FAIL"))
} else {
  cat("  No agreement data (Stata results not yet available)\n")
  frechet_ok <- NA; simple_ok <- NA
}

if (nrow(summary_df) > 0) {
  # Coverage check (for delta > 0 cells)
  cov_cells <- summary_df[summary_df$delta > 0, ]

  if (nrow(cov_cells) > 0 && any(!is.na(cov_cells$coverage_R))) {
    r_cov_ok <- all(cov_cells$coverage_R[!is.na(cov_cells$coverage_R)] > 0.80, na.rm = TRUE)
    cat(sprintf("  R coverage > 0.80 (delta>0): %s\n",
                if (r_cov_ok) "PASS" else "FAIL"))
  }

  if (nrow(cov_cells) > 0 && any(!is.na(cov_cells$coverage_Stata_shared))) {
    s_cov_ok <- all(cov_cells$coverage_Stata_shared[!is.na(cov_cells$coverage_Stata_shared)] > 0.80,
                     na.rm = TRUE)
    cat(sprintf("  Stata coverage > 0.80 (delta>0): %s\n",
                if (s_cov_ok) "PASS" else "FAIL"))
  }

  # Size check (delta == 0)
  size_cells <- summary_df[summary_df$delta == 0, ]
  if (nrow(size_cells) > 0 && any(!is.na(size_cells$rej_null_R))) {
    r_size_ok <- all(size_cells$rej_null_R[!is.na(size_cells$rej_null_R)] < 0.15, na.rm = TRUE)
    cat(sprintf("  R size < 0.15 (delta=0): %s\n",
                if (r_size_ok) "PASS" else "FAIL"))
  }

  # R vs Stata within 5pp for power
  power_cells <- summary_df[summary_df$delta > 0 &
                              !is.na(summary_df$rej_null_R) &
                              !is.na(summary_df$rej_null_Stata_shared), ]
  if (nrow(power_cells) > 0) {
    power_diff <- abs(power_cells$rej_null_R - power_cells$rej_null_Stata_shared)
    power_ok <- all(power_diff < 0.10, na.rm = TRUE)
    cat(sprintf("  R vs Stata power within 10pp: %s (max diff=%.3f)\n",
                if (power_ok) "PASS" else "FAIL", max(power_diff)))
  }
} else {
  cat("  No summary data available\n")
}

cat("\n=== Comparison complete ===\n")
cat("Output saved to:", comp_dir, "\n")
