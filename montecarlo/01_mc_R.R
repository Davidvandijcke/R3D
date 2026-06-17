#!/usr/bin/env Rscript
#
# 01_mc_R.R
# Monte Carlo simulations for R3D: R-side
# Generates DGP data, runs R3D, exports data/xi/bw for Stata comparison.
#
# Usage: Rscript 01_mc_R.R [--subset-only]
#   --subset-only: Only generate the Stata subset (100 sims), skip full 500-sim R run

library(R3D)

args <- commandArgs(trailingOnly = TRUE)
subset_only <- "--subset-only" %in% args

# ============================================================================
# CONFIGURATION
# ============================================================================
N_SIM         <- 500     # Total R simulations
N_STATA       <- 100     # Subset exported for Stata
SAMPLE_SIZES  <- c(200, 500)
N_I           <- 50      # Draws per unit
Q_GRID        <- (1:9) / 10  # 9 deciles
BOOT_REPS     <- 200
ALPHA         <- 0.05
DELTAS        <- c(0, 1, 2)
METHODS       <- c("simple", "frechet")
KERNEL        <- "epanechnikov"
P_ORDER       <- 2
NCORES        <- parallel::detectCores() - 1
if (NCORES < 1) NCORES <- 1
if (subset_only) N_SIM <- N_STATA

NQ <- length(Q_GRID)

outdir <- file.path(getwd(), "output")
dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "xi"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "results_R"), recursive = TRUE, showWarnings = FALSE)

cat("=== R3D Monte Carlo Simulation (R side) ===\n")
cat("N_SIM:", N_SIM, " N_STATA:", N_STATA, " BOOT_REPS:", BOOT_REPS, "\n")
cat("Sample sizes:", SAMPLE_SIZES, "\n")
cat("Deltas:", DELTAS, "\n")
cat("Methods:", METHODS, "\n")
cat("Cores:", NCORES, "\n\n")

# ============================================================================
# DGP FUNCTIONS
# ============================================================================

#' DGP 1: Normal-Normal (constant treatment effect)
#' True tau(q) = Delta for all q
generate_dgp1 <- function(n, delta, n_i, seed) {
  set.seed(seed)
  X <- runif(n, -1, 1)
  above <- as.numeric(X >= 0)
  Y_list <- lapply(seq_len(n), function(i) {
    mu_i <- 5 + 5 * X[i] + delta * above[i]
    sigma_i <- abs(1 + X[i])
    if (sigma_i < 0.1) sigma_i <- 0.1
    rnorm(n_i, mean = mu_i, sd = sigma_i)
  })
  true_tau <- rep(delta, length(Q_GRID))
  list(X = X, Y_list = Y_list, true_tau = true_tau)
}

#' DGP 2: Normal-Exponential mixture (heterogeneous effect)
#' True tau(q) computed numerically
generate_dgp2 <- function(n, delta, n_i, seed) {
  set.seed(seed)
  X <- runif(n, -1, 1)
  above <- as.numeric(X >= 0)
  Y_list <- lapply(seq_len(n), function(i) {
    mu_i <- 2 * X[i] + delta * above[i]
    lambda_i <- max(0.5, min(1.5, 1 + 0.5 * X[i]))
    rnorm(n_i, mean = mu_i, sd = 1) + 2 * rexp(n_i, rate = lambda_i)
  })
  list(X = X, Y_list = Y_list)
}

#' Fuzzy variant of DGP 1
generate_dgp1_fuzzy <- function(n, delta, n_i, seed) {
  set.seed(seed)
  X <- runif(n, -1, 1)
  above <- as.numeric(X >= 0)
  p_treat <- 0.1 + 0.8 * above
  T_vec <- rbinom(n, 1, p_treat)
  Y_list <- lapply(seq_len(n), function(i) {
    mu_i <- 5 + 5 * X[i] + delta * T_vec[i]
    sigma_i <- abs(1 + X[i])
    if (sigma_i < 0.1) sigma_i <- 0.1
    rnorm(n_i, mean = mu_i, sd = sigma_i)
  })
  # True LATE tau(q) = delta (constant) — same as sharp since effect is homogeneous
  true_tau <- rep(delta, length(Q_GRID))
  list(X = X, Y_list = Y_list, T_vec = T_vec, true_tau = true_tau)
}

# ============================================================================
# HELPER: Save data CSV for Stata
# ============================================================================
save_mc_data <- function(X, Y_list, q_grid, filepath, T_vec = NULL) {
  Qmat <- t(sapply(Y_list, quantile, probs = q_grid))
  df <- data.frame(X = X)
  if (!is.null(T_vec)) df$T_treat <- T_vec
  for (j in seq_along(q_grid)) df[[paste0("Q", j)]] <- Qmat[, j]
  write.csv(df, filepath, row.names = FALSE)
}

# ============================================================================
# MAIN SIMULATION LOOP
# ============================================================================

# Pre-compute true tau for DGP2 at each delta (these don't depend on n or seed)
true_tau_dgp2 <- list()
for (delta in DELTAS) {
  key <- as.character(delta)
  set.seed(99999)
  large_n <- 1e6
  y_above <- rnorm(large_n, mean = delta, sd = 1) + 2 * rexp(large_n, rate = 1)
  y_below <- rnorm(large_n, mean = 0, sd = 1) + 2 * rexp(large_n, rate = 1)
  true_tau_dgp2[[key]] <- quantile(y_above, probs = Q_GRID, names = FALSE) -
                           quantile(y_below, probs = Q_GRID, names = FALSE)
}

# Storage for R results
all_R_results <- list()

dgp_names <- c("dgp1", "dgp2")

for (n_obs in SAMPLE_SIZES) {
  cat("\n===== n =", n_obs, "=====\n")

  # Generate xi matrices for all sims (shared across DGPs/deltas for same n)
  cat("Generating xi matrices...\n")
  for (s in 1:N_STATA) {
    set.seed(2000000 + s)
    xi <- matrix(rnorm(n_obs * BOOT_REPS), nrow = n_obs, ncol = BOOT_REPS)
    write.csv(xi, file.path(outdir, "xi", sprintf("mc_xi_%d_%d.csv", n_obs, s)),
              row.names = FALSE)
  }

  for (dgp_name in dgp_names) {
    for (delta in DELTAS) {
      cell_key <- sprintf("%s_%s_%d_d%s", dgp_name, "sharp", n_obs, delta)
      cat("\n--- Cell:", cell_key, "---\n")

      # True tau
      if (dgp_name == "dgp1") {
        true_tau <- rep(delta, NQ)
      } else {
        true_tau <- true_tau_dgp2[[as.character(delta)]]
      }

      for (method in METHODS) {
        result_key <- sprintf("%s_%s_%d_d%s_%s", dgp_name, "sharp", n_obs, delta, method)
        cat("  Method:", method, "\n")

        # Storage for this cell
        tau_mat <- matrix(NA, nrow = N_SIM, ncol = NQ)
        cb_lower_mat <- matrix(NA, nrow = N_SIM, ncol = NQ)
        cb_upper_mat <- matrix(NA, nrow = N_SIM, ncol = NQ)
        bw_mat <- matrix(NA, nrow = N_SIM, ncol = NQ)
        pval_nullity <- rep(NA, N_SIM)
        pval_homogeneity <- rep(NA, N_SIM)

        sim_run <- function(s) {
          set.seed(3000000 + s)
          data_seed <- 1000 + s
          if (dgp_name == "dgp1") {
            dat <- generate_dgp1(n_obs, delta, N_I, data_seed)
          } else {
            dat <- generate_dgp2(n_obs, delta, N_I, data_seed)
          }

          # Regenerate xi from seed rather than reading CSV to avoid file I/O in parallel workers;
          # seed matches CSV generation above so values are identical.
          xi <- NULL
          if (s <= N_STATA) {
            set.seed(2000000 + s)
            xi <- matrix(rnorm(n_obs * BOOT_REPS), nrow = n_obs, ncol = BOOT_REPS)
          }

          # Run r3d
          fit <- tryCatch({
            r3d(dat$X, dat$Y_list, cutoff = 0, method = method, p = P_ORDER,
                q_grid = Q_GRID, kernel_fun = KERNEL,
                boot = TRUE, boot_reps = BOOT_REPS, alpha = ALPHA,
                test = c("nullity", "homogeneity"),
                xi_mat = xi)
          }, error = function(e) {
            warning(sprintf("Sim %d failed: %s", s, e$message))
            NULL
          })

          if (is.null(fit)) return(NULL)

          # Extract bandwidths
          bw <- fit$bandwidths$h_star_num
          if (is.null(bw) || length(bw) == 0) { warning('NULL bandwidth returned'); return(NULL) }
          if (length(bw) == 1) bw <- rep(bw, NQ)

          # Extract p-values
          pn <- NA; ph <- NA
          if (!is.null(fit$boot_out$test_results)) {
            if (!is.null(fit$boot_out$test_results$nullity)) {
              first_null <- fit$boot_out$test_results$nullity[[1]]
              pn <- first_null$p_value
            }
            if (!is.null(fit$boot_out$test_results$homogeneity)) {
              first_homo <- fit$boot_out$test_results$homogeneity[[1]]
              ph <- first_homo$p_value
            }
          }

          c(
            list(
              tau = as.numeric(fit$tau),
              cb_lower = as.numeric(fit$boot_out$cb_lower),
              cb_upper = as.numeric(fit$boot_out$cb_upper),
              bw = bw,
              pval_nullity = pn,
              pval_homogeneity = ph,
              s = s
            ),
            if (s <= N_STATA) list(X = dat$X, Y_list = dat$Y_list)
            else list(X = NULL, Y_list = NULL)
          )
        }

        # Run simulations (parallel for speed)
        cat("    Running", N_SIM, "simulations...\n")
        t0 <- proc.time()

        if (NCORES > 1) {
          results_list <- parallel::mclapply(1:N_SIM, sim_run, mc.cores = NCORES)
        } else {
          results_list <- lapply(1:N_SIM, sim_run)
        }

        elapsed <- (proc.time() - t0)[3]
        cat("    Done in", round(elapsed, 1), "seconds\n")

        # Collect results
        for (res in results_list) {
          if (is.null(res)) next
          s <- res$s
          tau_mat[s, ] <- res$tau
          cb_lower_mat[s, ] <- res$cb_lower
          cb_upper_mat[s, ] <- res$cb_upper
          bw_mat[s, ] <- res$bw
          pval_nullity[s] <- res$pval_nullity
          pval_homogeneity[s] <- res$pval_homogeneity

          # Save data CSVs for Stata subset
          if (s <= N_STATA) {
            save_mc_data(res$X, res$Y_list, Q_GRID,
                         file.path(outdir, "data",
                                   sprintf("mc_%s_%d_%d_%d.csv", dgp_name, n_obs, delta, s)))

            # Save R bandwidths for shared-bw Stata runs
            write.csv(data.frame(q = Q_GRID, bw = res$bw),
                      file.path(outdir, "results_R",
                                sprintf("mc_bw_%s_%s_%d_%d_%d.csv",
                                        dgp_name, method, n_obs, delta, s)),
                      row.names = FALSE)
          }
        }

        # Store in master list
        all_R_results[[result_key]] <- list(
          tau_mat = tau_mat,
          cb_lower_mat = cb_lower_mat,
          cb_upper_mat = cb_upper_mat,
          bw_mat = bw_mat,
          pval_nullity = pval_nullity,
          pval_homogeneity = pval_homogeneity,
          true_tau = true_tau,
          dgp = dgp_name,
          method = method,
          n = n_obs,
          delta = delta,
          design = "sharp"
        )

        # Print quick summary
        valid <- !is.na(tau_mat[, 1])
        n_valid <- sum(valid)
        if (n_valid > 0) {
          bias <- colMeans(tau_mat[valid, , drop = FALSE]) - true_tau
          rmse <- sqrt(colMeans((tau_mat[valid, , drop = FALSE] -
                                  matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE))^2))
          uniform_covered <- apply(
            cb_lower_mat[valid, , drop = FALSE] <=
              matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE) &
            cb_upper_mat[valid, , drop = FALSE] >=
              matrix(true_tau, nrow = n_valid, ncol = NQ, byrow = TRUE),
            1, all)
          coverage_uniform <- mean(uniform_covered, na.rm = TRUE)
          cat("    Valid sims:", n_valid, "/", N_SIM, "\n")
          cat("    Mean |bias|:", round(mean(abs(bias)), 4), "\n")
          cat("    Mean RMSE:", round(mean(rmse), 4), "\n")
          cat("    Uniform coverage:", round(coverage_uniform, 4), "\n")
          if (delta == 0) {
            cat("    Rejection rate (nullity, alpha=0.05):",
                round(mean(pval_nullity[valid] < 0.05, na.rm = TRUE), 4), "\n")
          } else {
            cat("    Power (nullity, alpha=0.05):",
                round(mean(pval_nullity[valid] < 0.05, na.rm = TRUE), 4), "\n")
          }
        }
      }
    }
  }

  # Fuzzy: 1 cell — DGP1, n=200, delta=2, simple method only
  if (n_obs == 200) {
    delta_f <- 2
    method_f <- "simple"
    cell_key_f <- sprintf("dgp1_fuzzy_%d_d%d_%s", n_obs, delta_f, method_f)
    cat("\n--- Fuzzy cell:", cell_key_f, "---\n")

    true_tau_f <- rep(delta_f, NQ)
    tau_mat_f <- matrix(NA, nrow = N_STATA, ncol = NQ)
    cb_lower_mat_f <- matrix(NA, nrow = N_STATA, ncol = NQ)
    cb_upper_mat_f <- matrix(NA, nrow = N_STATA, ncol = NQ)
    bw_mat_f <- matrix(NA, nrow = N_STATA, ncol = NQ)
    pval_nullity_f <- rep(NA, N_STATA)

    sim_run_fuzzy <- function(s) {
      data_seed <- 1000 + s
      dat <- generate_dgp1_fuzzy(n_obs, delta_f, N_I, data_seed)

      set.seed(2000000 + s)
      xi <- matrix(rnorm(n_obs * BOOT_REPS), nrow = n_obs, ncol = BOOT_REPS)

      fit <- tryCatch({
        r3d(dat$X, dat$Y_list, T = dat$T_vec, cutoff = 0, method = method_f,
            p = P_ORDER, q_grid = Q_GRID, kernel_fun = KERNEL, fuzzy = TRUE,
            boot = TRUE, boot_reps = BOOT_REPS, alpha = ALPHA,
            test = c("nullity"), xi_mat = xi)
      }, error = function(e) {
        warning(sprintf("Fuzzy sim %d failed: %s", s, e$message))
        NULL
      })

      if (is.null(fit)) return(NULL)

      bw <- fit$bandwidths$h_star_num
      if (is.null(bw) || length(bw) == 0) { warning('NULL bandwidth returned'); return(NULL) }
      if (length(bw) == 1) bw <- rep(bw, NQ)

      pn <- NA
      if (!is.null(fit$boot_out$test_results$nullity)) {
        pn <- fit$boot_out$test_results$nullity[[1]]$p_value
      }

      # Save data for Stata
      save_mc_data(dat$X, dat$Y_list, Q_GRID,
                   file.path(outdir, "data",
                             sprintf("mc_dgp1fuzzy_%d_%d_%d.csv", n_obs, delta_f, s)),
                   T_vec = dat$T_vec)
      write.csv(data.frame(q = Q_GRID, bw = bw),
                file.path(outdir, "results_R",
                          sprintf("mc_bw_dgp1fuzzy_%s_%d_%d_%d.csv",
                                  method_f, n_obs, delta_f, s)),
                row.names = FALSE)
      # Save denominator bandwidth
      bw_den <- fit$bandwidths$h_star_den
      write.csv(data.frame(bw_den = bw_den),
                file.path(outdir, "results_R",
                          sprintf("mc_bwden_dgp1fuzzy_%s_%d_%d_%d.csv",
                                  method_f, n_obs, delta_f, s)),
                row.names = FALSE)

      list(
        tau = as.numeric(fit$tau),
        cb_lower = as.numeric(fit$boot_out$cb_lower),
        cb_upper = as.numeric(fit$boot_out$cb_upper),
        bw = bw,
        pval_nullity = pn,
        s = s
      )
    }

    cat("    Running", N_STATA, "fuzzy simulations...\n")
    t0 <- proc.time()
    if (NCORES > 1) {
      fuzz_list <- parallel::mclapply(1:N_STATA, sim_run_fuzzy, mc.cores = NCORES)
    } else {
      fuzz_list <- lapply(1:N_STATA, sim_run_fuzzy)
    }
    elapsed <- (proc.time() - t0)[3]
    cat("    Done in", round(elapsed, 1), "seconds\n")

    for (res in fuzz_list) {
      if (is.null(res)) next
      s <- res$s
      tau_mat_f[s, ] <- res$tau
      cb_lower_mat_f[s, ] <- res$cb_lower
      cb_upper_mat_f[s, ] <- res$cb_upper
      bw_mat_f[s, ] <- res$bw
      pval_nullity_f[s] <- res$pval_nullity
    }

    all_R_results[[cell_key_f]] <- list(
      tau_mat = tau_mat_f,
      cb_lower_mat = cb_lower_mat_f,
      cb_upper_mat = cb_upper_mat_f,
      bw_mat = bw_mat_f,
      pval_nullity = pval_nullity_f,
      true_tau = true_tau_f,
      dgp = "dgp1", method = method_f, n = n_obs,
      delta = delta_f, design = "fuzzy"
    )

    valid_f <- !is.na(tau_mat_f[, 1])
    if (sum(valid_f) > 0) {
      cat("    Valid:", sum(valid_f), "/", N_STATA, "\n")
      cat("    Power:", round(mean(pval_nullity_f[valid_f] < 0.05, na.rm = TRUE), 4), "\n")
    }
  }
}

# ============================================================================
# SAVE ALL R RESULTS
# ============================================================================
cat("\n=== Saving R results ===\n")
saveRDS(all_R_results, file.path(outdir, "results_R", "mc_R_all.rds"))

# Also save per-replication tau/cb/pvals as CSV for the Stata-overlap subset
for (key in names(all_R_results)) {
  res <- all_R_results[[key]]
  n_save <- min(N_STATA, nrow(res$tau_mat))

  # Per-replication results (first N_STATA sims)
  for (s in 1:n_save) {
    if (is.na(res$tau_mat[s, 1])) next
    df <- data.frame(
      q = Q_GRID,
      tau = res$tau_mat[s, ],
      cb_lower = res$cb_lower_mat[s, ],
      cb_upper = res$cb_upper_mat[s, ],
      bw = res$bw_mat[s, ]
    )
    write.csv(df,
              file.path(outdir, "results_R",
                        sprintf("mc_res_%s_%d.csv", key, s)),
              row.names = FALSE)
  }

  # Per-replication p-values
  pv_df <- data.frame(sim = 1:n_save)
  if (!is.null(res$pval_nullity)) pv_df$pval_nullity <- res$pval_nullity[1:n_save]
  if (!is.null(res$pval_homogeneity)) pv_df$pval_homogeneity <- res$pval_homogeneity[1:n_save]
  write.csv(pv_df,
            file.path(outdir, "results_R", sprintf("mc_pvals_%s.csv", key)),
            row.names = FALSE)
}

# Save summary table
summary_rows <- list()
for (key in names(all_R_results)) {
  res <- all_R_results[[key]]
  valid <- !is.na(res$tau_mat[, 1])
  n_valid <- sum(valid)
  if (n_valid == 0) next

  tt <- res$true_tau
  bias <- colMeans(res$tau_mat[valid, , drop = FALSE]) - tt
  rmse <- sqrt(colMeans((res$tau_mat[valid, , drop = FALSE] -
                          matrix(tt, nrow = n_valid, ncol = NQ, byrow = TRUE))^2))
  uniform_covered <- apply(
    res$cb_lower_mat[valid, , drop = FALSE] <=
      matrix(tt, nrow = n_valid, ncol = NQ, byrow = TRUE) &
    res$cb_upper_mat[valid, , drop = FALSE] >=
      matrix(tt, nrow = n_valid, ncol = NQ, byrow = TRUE),
    1, all)
  coverage_uniform <- mean(uniform_covered, na.rm = TRUE)

  rej_null <- if (!is.null(res$pval_nullity))
    mean(res$pval_nullity[valid] < 0.05, na.rm = TRUE) else NA
  rej_homo <- if (!is.null(res$pval_homogeneity))
    mean(res$pval_homogeneity[valid] < 0.05, na.rm = TRUE) else NA

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    cell = key,
    dgp = res$dgp, method = res$method, n = res$n,
    delta = res$delta, design = res$design,
    n_valid = n_valid,
    mean_abs_bias = mean(abs(bias)),
    max_abs_bias = max(abs(bias)),
    mean_rmse = mean(rmse),
    uniform_coverage = coverage_uniform,
    rej_nullity = rej_null,
    rej_homogeneity = rej_homo,
    stringsAsFactors = FALSE
  )
}
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(outdir, "results_R", "mc_R_summary.csv"), row.names = FALSE)

cat("\n=== R Monte Carlo complete ===\n")
cat("Summary:\n")
print(summary_df[, c("cell", "n_valid", "mean_abs_bias", "uniform_coverage",
                      "rej_nullity", "rej_homogeneity")])
cat("\nAll results saved to:", file.path(outdir, "results_R"), "\n")
