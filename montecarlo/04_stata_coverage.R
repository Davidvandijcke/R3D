#!/usr/bin/env Rscript
#
# 04_stata_coverage.R
# Compute uniform coverage from Stata Monte Carlo results.
#
# Usage: Rscript 04_stata_coverage.R

cat("=== R3D Monte Carlo: Stata Uniform Coverage ===\n\n")

outdir <- file.path(getwd(), "output")
stata_dir <- file.path(outdir, "results_Stata")

NQ <- 9
Q_GRID <- (1:9) / 10
N_SIMS <- 100

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

# Cell definitions
cells <- expand.grid(
  dgp = c("dgp1", "dgp2"),
  n = c(200, 500),
  delta = c(0, 1, 2),
  method = c("simple", "frechet"),
  stringsAsFactors = FALSE
)

results <- list()

for (i in seq_len(nrow(cells))) {
  dgp <- cells$dgp[i]
  n_obs <- cells$n[i]
  delta <- cells$delta[i]
  method <- cells$method[i]

  cell <- sprintf("%s_sharp_%d_d%d_%s", dgp, n_obs, delta, method)

  # True tau
  if (dgp == "dgp1") {
    true_tau <- rep(delta, NQ)
  } else {
    true_tau <- true_tau_dgp2[[as.character(delta)]]
  }

  for (tag in c("sharedbw", "ownbw")) {
    # Collect results across sims
    tau_mat <- matrix(NA, nrow = N_SIMS, ncol = NQ)
    cb_lo_mat <- matrix(NA, nrow = N_SIMS, ncol = NQ)
    cb_hi_mat <- matrix(NA, nrow = N_SIMS, ncol = NQ)
    pval_null <- rep(NA, N_SIMS)
    pval_homo <- rep(NA, N_SIMS)

    for (s in 1:N_SIMS) {
      f <- file.path(stata_dir, sprintf("mc_res_%s_%s_%d.csv", cell, tag, s))
      if (!file.exists(f)) next
      d <- read.csv(f)
      if (nrow(d) != NQ) next
      tau_mat[s, ] <- d$tau
      cb_lo_mat[s, ] <- d$cb_lower
      cb_hi_mat[s, ] <- d$cb_upper

      pf <- file.path(stata_dir, sprintf("mc_pval_%s_%s_%d.csv", cell, tag, s))
      if (file.exists(pf)) {
        pd <- read.csv(pf)
        if ("pval_nullity" %in% names(pd)) pval_null[s] <- pd$pval_nullity[1]
        if ("pval_homogeneity" %in% names(pd)) pval_homo[s] <- pd$pval_homogeneity[1]
      }
    }

    valid <- !is.na(tau_mat[, 1])
    n_valid <- sum(valid)

    if (n_valid == 0) next

    # Uniform coverage: true tau within [cb_lo, cb_hi] at ALL quantiles simultaneously
    covered_mat <- cb_lo_mat[valid, , drop = FALSE] <=
      matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE) &
      cb_hi_mat[valid, , drop = FALSE] >=
      matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE)
    uniform_covered <- apply(covered_mat, 1, all)
    uniform_coverage <- mean(uniform_covered, na.rm = TRUE)

    # Also compute pointwise for comparison
    pointwise_coverage <- colMeans(covered_mat, na.rm = TRUE)

    # Bias
    bias <- colMeans(tau_mat[valid, , drop = FALSE]) - true_tau
    mean_abs_bias <- mean(abs(bias))

    # RMSE
    rmse <- sqrt(colMeans((tau_mat[valid, , drop = FALSE] -
                             matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE))^2))
    mean_rmse <- mean(rmse)

    # Rejection rates
    rej_null <- mean(pval_null[valid] < 0.05, na.rm = TRUE)
    rej_homo <- mean(pval_homo[valid] < 0.05, na.rm = TRUE)

    results[[length(results) + 1]] <- data.frame(
      cell = cell, dgp = dgp, method = method, n = n_obs, delta = delta,
      bw_type = tag, n_valid = n_valid,
      uniform_coverage = uniform_coverage,
      min_pointwise_cov = min(pointwise_coverage),
      mean_pointwise_cov = mean(pointwise_coverage),
      mean_abs_bias = mean_abs_bias,
      mean_rmse = mean_rmse,
      rej_nullity = rej_null,
      rej_homogeneity = rej_homo,
      stringsAsFactors = FALSE, row.names = NULL
    )
  }
}

df <- do.call(rbind, results)

# Print summary table
cat(sprintf("%-45s  %4s  %6s  %6s  %6s  %6s  %6s\n",
            "Cell", "N", "UCov", "MnPCov", "Bias", "RMSE", "RejH0"))
cat(paste(rep("-", 95), collapse = ""), "\n")

for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  cat(sprintf("%-45s  %4d  %6.3f  %6.3f  %6.4f  %6.4f  %6.3f\n",
              paste0(r$cell, " (", r$bw_type, ")"),
              r$n_valid, r$uniform_coverage,
              r$mean_pointwise_cov, r$mean_abs_bias,
              r$mean_rmse, r$rej_nullity))
}

# Save CSV
comp_dir <- file.path(outdir, "comparison")
dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(df, file.path(comp_dir, "stata_coverage.csv"), row.names = FALSE)

cat("\n\n=== Coverage Summary ===\n")
cat("Nominal level: 95% (alpha = 0.05)\n\n")

# Group by delta for interpretation
for (d in unique(df$delta)) {
  sub <- df[df$delta == d, ]
  cat(sprintf("Delta = %d:\n", d))
  if (d == 0) {
    cat(sprintf("  Uniform coverage range: [%.3f, %.3f]  (target: >= 0.95)\n",
                min(sub$uniform_coverage), max(sub$uniform_coverage)))
    cat(sprintf("  Nullity rejection rate (size): [%.3f, %.3f]  (target: <= 0.05)\n",
                min(sub$rej_nullity), max(sub$rej_nullity)))
  } else {
    cat(sprintf("  Uniform coverage range: [%.3f, %.3f]  (target: >= 0.95)\n",
                min(sub$uniform_coverage), max(sub$uniform_coverage)))
    cat(sprintf("  Nullity rejection rate (power): [%.3f, %.3f]  (target: high)\n",
                min(sub$rej_nullity), max(sub$rej_nullity)))
  }
  cat("\n")
}

cat("Results saved to:", file.path(comp_dir, "stata_coverage.csv"), "\n")
