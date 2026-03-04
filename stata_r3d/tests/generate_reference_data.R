#!/usr/bin/env Rscript
#
# generate_reference_data.R
# Generates comprehensive reference datasets and R3D outputs for Stata equivalence testing.
# Run from tests/ directory: Rscript generate_reference_data.R
#

library(R3D)

cat("=== Generating R3D reference data for Stata equivalence tests ===\n")

# ============================================================================
# Helpers
# ============================================================================
save_data_for_stata <- function(X, Y_list, q_grid, filename, T_vec = NULL) {
  n <- length(X)
  nq <- length(q_grid)
  Qmat <- R3D:::.compute_empirical_qmat(Y_list, q_grid)
  df <- data.frame(X = X)
  if (!is.null(T_vec)) df$T_treat <- T_vec
  for (j in seq_len(nq)) df[[paste0("Q", j)]] <- Qmat[, j]
  write.csv(df, filename, row.names = FALSE)
  cat("  Saved:", filename, "(", n, "obs x", ncol(df), "vars)\n")
  return(Qmat)
}

save_results <- function(fit, filename) {
  nq <- length(fit$q_grid)
  df <- data.frame(
    q = fit$q_grid,
    tau = as.numeric(fit$tau),
    bw_num = as.numeric(fit$bandwidths$h_star_num),
    int_plus = as.numeric(fit$int_plus),
    int_minus = as.numeric(fit$int_minus)
  )
  if (!is.null(fit$bandwidths$h_star_den) && length(fit$bandwidths$h_star_den) > 0) {
    df$bw_den <- rep(fit$bandwidths$h_star_den, nq)
  }
  if (!is.null(fit$boot_out)) {
    df$cb_lower <- as.numeric(fit$boot_out$cb_lower)
    df$cb_upper <- as.numeric(fit$boot_out$cb_upper)
  }
  write.csv(df, filename, row.names = FALSE)
  cat("  Saved:", filename, "\n")
}

save_alpha_matrix <- function(mat, filename) {
  write.csv(mat, filename, row.names = FALSE)
  cat("  Saved:", filename, "\n")
}

# ============================================================================
# GENERATE BASE DATA (shared across multiple tests)
# ============================================================================
cat("\n--- Generating base datasets ---\n")

# Sharp RDD data (n=200)
set.seed(42)
n_sharp <- 200
X_sharp <- runif(n_sharp, -1, 1)
Y_list_sharp <- lapply(seq_len(n_sharp), function(i) {
  mu <- 5 + 0.5 * X_sharp[i] + 3 * (X_sharp[i] >= 0)
  rnorm(50, mean = mu, sd = 0.5)
})
q_grid_20 <- (1:20) / 21
q_grid_99 <- (1:99) / 100

# Fuzzy RDD data (n=200)
set.seed(42)
n_fuzzy <- 200
X_fuzzy <- runif(n_fuzzy, -1, 1)
p_treat <- 0.1 + 0.8 * (X_fuzzy >= 0)
T_fuzzy <- rbinom(n_fuzzy, 1, p_treat)
Y_list_fuzzy <- lapply(seq_len(n_fuzzy), function(i) {
  mu <- 5 + 0.5 * X_fuzzy[i] + 3 * T_fuzzy[i]
  rnorm(50, mean = mu, sd = 0.5)
})

# Save data files
Qmat_sharp_20 <- save_data_for_stata(X_sharp, Y_list_sharp, q_grid_20,
                                      "ref_data_sharp_simple_20.csv")
Qmat_sharp_99 <- save_data_for_stata(X_sharp, Y_list_sharp, q_grid_99,
                                      "ref_data_sharp_simple_99.csv")
Qmat_fuzzy_20 <- save_data_for_stata(X_fuzzy, Y_list_fuzzy, q_grid_20,
                                      "ref_data_fuzzy_simple_20.csv", T_vec = T_fuzzy)


# ############################################################################
# LAYER 1: INTERMEDIATE COMPUTATIONS
# ############################################################################
cat("\n\n=== LAYER 1: Intermediate Computations ===\n")

# --------------------------------------------------------------------------
# Test 1.1: Quantile computation (R type-7)
# --------------------------------------------------------------------------
cat("\n--- Test 1.1: Quantile computation ---\n")

# Use first few observations from sharp data with different grid sizes
# Case A: 5 obs, single quantile
set.seed(42)
Y_q_test <- list(
  rnorm(5),
  rnorm(50),
  rnorm(500)
)
q_grids_test <- list(
  c(0.5),
  c(0.1, 0.5, 0.9),
  (1:99)/100
)

q_results <- data.frame()
for (ci in 1:3) {
  y_vec <- Y_q_test[[ci]]
  for (gi in 1:3) {
    qg <- q_grids_test[[gi]]
    qvals <- quantile(y_vec, probs = qg, type = 7, names = FALSE)
    for (qi in seq_along(qg)) {
      q_results <- rbind(q_results, data.frame(
        case = ci, grid = gi, q = qg[qi], value = qvals[qi]
      ))
    }
  }
}
write.csv(q_results, "ref_quantiles.csv", row.names = FALSE)
# Also save the raw Y vectors for Stata to load
for (ci in 1:3) {
  write.csv(data.frame(y = Y_q_test[[ci]]),
            paste0("ref_quantile_y_", ci, ".csv"), row.names = FALSE)
}
cat("  Saved: ref_quantiles.csv + ref_quantile_y_{1,2,3}.csv\n")

# --------------------------------------------------------------------------
# Test 1.2: Density estimation
# --------------------------------------------------------------------------
cat("\n--- Test 1.2: Density estimation ---\n")

Xc <- X_sharp - 0
sigma_X <- sd(Xc)
h_bw <- 1.06 * sigma_X * n_sharp^(-1/5)

kernel_tri <- function(u) pmax(0, 1 - abs(u))
kernel_epa <- function(u) 0.75 * pmax(0, 1 - u^2)
kernel_uni <- function(u) 0.5 * (abs(u) <= 1)

f_X_tri <- mean(kernel_tri(Xc / h_bw)) / h_bw
f_X_epa <- mean(kernel_epa(Xc / h_bw)) / h_bw
f_X_uni <- mean(kernel_uni(Xc / h_bw)) / h_bw

density_df <- data.frame(
  kernel = c("triangular", "epanechnikov", "uniform"),
  f_X_hat = c(f_X_tri, f_X_epa, f_X_uni),
  h_bw = h_bw,
  sigma_X = sigma_X
)
write.csv(density_df, "ref_density_all.csv", row.names = FALSE)
cat("  Saved: ref_density_all.csv\n")
cat("  f_X_hat: tri=", f_X_tri, " epa=", f_X_epa, " uni=", f_X_uni, "\n")

# --------------------------------------------------------------------------
# Test 1.3: Isotonic regression (PAVA)
# --------------------------------------------------------------------------
cat("\n--- Test 1.3: Isotonic regression ---\n")

iso_cases <- list(
  monotone = c(1, 2, 3, 4, 5),              # already monotone
  single_viol = c(1, 3, 2, 4, 5),           # single violation
  multi_viol = c(5, 3, 4, 1, 2, 6),         # multiple violations
  all_equal = c(3, 3, 3, 3, 3)              # all equal
)

iso_results <- data.frame()
for (nm in names(iso_cases)) {
  y_in <- iso_cases[[nm]]
  iso_out <- as.numeric(isoreg(seq_along(y_in), y_in)$yf)
  for (i in seq_along(y_in)) {
    iso_results <- rbind(iso_results, data.frame(
      case = nm, idx = i, input = y_in[i], output = iso_out[i]
    ))
  }
}
write.csv(iso_results, "ref_isotonic.csv", row.names = FALSE)
cat("  Saved: ref_isotonic.csv\n")

# --------------------------------------------------------------------------
# Test 1.4: Gini from quantile function
# --------------------------------------------------------------------------
cat("\n--- Test 1.4: Gini coefficient ---\n")

# Uniform distribution: Q(p) = p, Gini = 1/3
q_uniform <- (1:99)/100
qf_uniform <- q_uniform

# Exponential(1): Q(p) = -ln(1-p), Gini ≈ 0.5
qf_exponential <- -log(1 - q_uniform)

# Degenerate (constant): Q(p) = c, Gini = 0
qf_degenerate <- rep(5, 99)

gini_df <- data.frame(
  case = c("uniform", "exponential", "degenerate"),
  gini = c(
    R3D:::calculate_gini_from_quantile(q_uniform, qf_uniform),
    R3D:::calculate_gini_from_quantile(q_uniform, qf_exponential),
    R3D:::calculate_gini_from_quantile(q_uniform, qf_degenerate)
  )
)
write.csv(gini_df, "ref_gini.csv", row.names = FALSE)
# Also save the quantile functions for Stata
for (cs in c("uniform", "exponential", "degenerate")) {
  qf <- switch(cs, uniform = qf_uniform, exponential = qf_exponential,
                degenerate = qf_degenerate)
  write.csv(data.frame(q = q_uniform, qf = qf),
            paste0("ref_gini_qf_", cs, ".csv"), row.names = FALSE)
}
cat("  Saved: ref_gini.csv + ref_gini_qf_{uniform,exponential,degenerate}.csv\n")
cat("  Gini values:", gini_df$gini, "\n")


# ############################################################################
# LAYER 2: LOCAL POLYNOMIAL REGRESSION
# ############################################################################
cat("\n\n=== LAYER 2: Local Polynomial Regression ===\n")

# --------------------------------------------------------------------------
# Test 2.1: locpoly coefficients + intercept weights
# --------------------------------------------------------------------------
for (kname in c("triangular", "epanechnikov", "uniform")) {
  cat("\n--- Test 2.1: locpoly with", kname, "kernel ---\n")

  kernel_type_int <- switch(kname, triangular = 1L, epanechnikov = 2L, uniform = 3L)
  p_order <- 2L
  nq <- 20L

  # Use automatic bandwidth for this kernel
  kernel_fn <- switch(kname,
    triangular = function(u) pmax(0, 1 - abs(u)),
    epanechnikov = function(u) 0.75 * pmax(0, 1 - u^2),
    uniform = function(u) 0.5 * (abs(u) <= 1))

  fit_temp <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
                  q_grid = q_grid_20, kernel_fun = kname)

  h_vec <- fit_temp$bandwidths$h_star_num

  Xc <- X_sharp - 0
  N <- as.integer(n_sharp)
  NQ <- as.integer(nq)

  outPlus <- .Fortran("locweights",
    X = as.double(Xc), YMAT = as.double(Qmat_sharp_20),
    N = N, P = p_order,
    H = as.double(if (length(h_vec) == 1) rep(h_vec, nq) else h_vec),
    SIDE = as.integer(1), KERNEL_TYPE = kernel_type_int,
    ALPHA = double((p_order + 1) * nq), WINT = double(n_sharp * nq),
    INFO = integer(1), NQ = NQ, PACKAGE = "R3D")

  outMinus <- .Fortran("locweights",
    X = as.double(Xc), YMAT = as.double(Qmat_sharp_20),
    N = N, P = p_order,
    H = as.double(if (length(h_vec) == 1) rep(h_vec, nq) else h_vec),
    SIDE = as.integer(0), KERNEL_TYPE = kernel_type_int,
    ALPHA = double((p_order + 1) * nq), WINT = double(n_sharp * nq),
    INFO = integer(1), NQ = NQ, PACKAGE = "R3D")

  alpha_plus <- matrix(outPlus$ALPHA, nrow = p_order + 1, ncol = nq)
  alpha_minus <- matrix(outMinus$ALPHA, nrow = p_order + 1, ncol = nq)
  w_plus <- matrix(outPlus$WINT, nrow = n_sharp, ncol = nq)
  w_minus <- matrix(outMinus$WINT, nrow = n_sharp, ncol = nq)

  save_alpha_matrix(alpha_plus, paste0("ref_locpoly_alpha_plus_", kname, ".csv"))
  save_alpha_matrix(alpha_minus, paste0("ref_locpoly_alpha_minus_", kname, ".csv"))
  save_alpha_matrix(w_plus, paste0("ref_locpoly_w_plus_", kname, ".csv"))
  save_alpha_matrix(w_minus, paste0("ref_locpoly_w_minus_", kname, ".csv"))

  # Also save the bandwidths used
  write.csv(data.frame(q = q_grid_20, h = h_vec),
            paste0("ref_locpoly_bw_", kname, ".csv"), row.names = FALSE)
}


# ############################################################################
# LAYER 3: FULL ESTIMATION PIPELINE
# ############################################################################
cat("\n\n=== LAYER 3: Full Estimation Pipeline ===\n")

# --------------------------------------------------------------------------
# Test 3.1: Sharp-Simple — all kernels
# --------------------------------------------------------------------------
for (kname in c("epanechnikov", "triangular", "uniform")) {
  cat("\n--- Test 3.1: Sharp-Simple,", kname, "---\n")
  set.seed(42)
  fit <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
             q_grid = q_grid_20, kernel_fun = kname)
  save_results(fit, paste0("ref_results_sharp_simple_", kname, ".csv"))

  # Save alpha matrices for deeper comparison
  save_alpha_matrix(fit$alpha_plus, paste0("ref_est_alpha_plus_sharp_simple_", kname, ".csv"))
  save_alpha_matrix(fit$alpha_minus, paste0("ref_est_alpha_minus_sharp_simple_", kname, ".csv"))
}

# Sharp-Simple, 99 quantiles, triangular
cat("\n--- Test 3.1b: Sharp-Simple, nq=99, triangular ---\n")
set.seed(42)
fit_99 <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
              q_grid = q_grid_99, kernel_fun = "triangular")
save_results(fit_99, "ref_results_sharp_simple_99.csv")

# --------------------------------------------------------------------------
# Test 3.2: Sharp-Frechet
# --------------------------------------------------------------------------
cat("\n--- Test 3.2: Sharp-Frechet ---\n")
set.seed(42)
fit_frechet <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "frechet", p = 2,
                   q_grid = q_grid_20, kernel_fun = "epanechnikov")
save_results(fit_frechet, "ref_results_sharp_frechet_20.csv")

# --------------------------------------------------------------------------
# Test 3.3: Fuzzy-Simple
# --------------------------------------------------------------------------
cat("\n--- Test 3.3: Fuzzy-Simple ---\n")
set.seed(42)
fit_fuzzy <- r3d(X_fuzzy, Y_list_fuzzy, T = T_fuzzy, cutoff = 0,
                 method = "simple", p = 2, q_grid = q_grid_20,
                 kernel_fun = "epanechnikov", fuzzy = TRUE)
save_results(fit_fuzzy, "ref_results_fuzzy_simple_20.csv")

# --------------------------------------------------------------------------
# Test 3.4: Fuzzy-Frechet
# --------------------------------------------------------------------------
cat("\n--- Test 3.4: Fuzzy-Frechet ---\n")
set.seed(42)
fit_fuzzy_fr <- r3d(X_fuzzy, Y_list_fuzzy, T = T_fuzzy, cutoff = 0,
                    method = "frechet", p = 2, q_grid = q_grid_20,
                    kernel_fun = "epanechnikov", fuzzy = TRUE)
save_results(fit_fuzzy_fr, "ref_results_fuzzy_frechet_20.csv")

# --------------------------------------------------------------------------
# Test 3.5: User-supplied bandwidth
# --------------------------------------------------------------------------
cat("\n--- Test 3.5: User-supplied bandwidth ---\n")
set.seed(42)
fit_fixbw <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
                 q_grid = q_grid_20, kernel_fun = "epanechnikov",
                 bandwidths = 0.5)
save_results(fit_fixbw, "ref_results_sharp_simple_fixedbw.csv")


# ############################################################################
# LAYER 4: BANDWIDTH SELECTION
# ############################################################################
cat("\n\n=== LAYER 4: Bandwidth Selection ===\n")

# --------------------------------------------------------------------------
# Test 4.1: Sharp bandwidth selection
# --------------------------------------------------------------------------
cat("\n--- Test 4.1: Sharp bandwidth selection ---\n")
kernel_tri_fn <- function(u) pmax(0, 1 - abs(u))

bw_sharp <- r3d_bwselect(X = X_sharp, Y_list = Y_list_sharp,
                          q_grid = q_grid_20, method = "simple",
                          s = 1, p = 2, kernel = kernel_tri_fn, cutoff = 0)

bw_sharp_df <- data.frame(
  q = q_grid_20,
  h_star_num = bw_sharp$h_star_num,
  pilot_h_num = bw_sharp$pilot_h_num,
  B_plus = bw_sharp$B_plus,
  B_minus = bw_sharp$B_minus,
  V_plus = bw_sharp$V_plus,
  V_minus = bw_sharp$V_minus
)
bw_sharp_df$f_X_hat <- bw_sharp$f_X_hat
write.csv(bw_sharp_df, "ref_bwselect_sharp.csv", row.names = FALSE)
cat("  Saved: ref_bwselect_sharp.csv\n")

# Frechet bandwidth selection
bw_frechet <- r3d_bwselect(X = X_sharp, Y_list = Y_list_sharp,
                            q_grid = q_grid_20, method = "frechet",
                            s = 1, p = 2, kernel = kernel_tri_fn, cutoff = 0)
bw_frechet_df <- data.frame(
  q = q_grid_20,
  h_star_num = bw_frechet$h_star_num,
  pilot_h_num = bw_frechet$pilot_h_num,
  B_plus = bw_frechet$B_plus,
  B_minus = bw_frechet$B_minus,
  V_plus = bw_frechet$V_plus,
  V_minus = bw_frechet$V_minus
)
bw_frechet_df$f_X_hat <- bw_frechet$f_X_hat
write.csv(bw_frechet_df, "ref_bwselect_frechet.csv", row.names = FALSE)
cat("  Saved: ref_bwselect_frechet.csv\n")

# --------------------------------------------------------------------------
# Test 4.2: Fuzzy bandwidth selection
# --------------------------------------------------------------------------
cat("\n--- Test 4.2: Fuzzy bandwidth selection ---\n")
bw_fuzzy <- r3d_bwselect(X = X_fuzzy, Y_list = Y_list_fuzzy, T = T_fuzzy,
                          q_grid = q_grid_20, method = "simple",
                          s = 1, p = 2, kernel = kernel_tri_fn, cutoff = 0,
                          fuzzy = TRUE)
bw_fuzzy_df <- data.frame(
  q = q_grid_20,
  h_star_num = bw_fuzzy$h_star_num,
  pilot_h_num = bw_fuzzy$pilot_h_num,
  B_plus = bw_fuzzy$B_plus,
  B_minus = bw_fuzzy$B_minus,
  V_plus = bw_fuzzy$V_plus,
  V_minus = bw_fuzzy$V_minus
)
bw_fuzzy_df$f_X_hat <- bw_fuzzy$f_X_hat
bw_fuzzy_df$h_star_den <- bw_fuzzy$h_star_den
bw_fuzzy_df$pilot_h_den <- bw_fuzzy$pilot_h_den
write.csv(bw_fuzzy_df, "ref_bwselect_fuzzy.csv", row.names = FALSE)
cat("  Saved: ref_bwselect_fuzzy.csv\n")


# ############################################################################
# LAYER 5: BOOTSTRAP & HYPOTHESIS TESTS
# ############################################################################
cat("\n\n=== LAYER 5: Bootstrap & Hypothesis Tests ===\n")

# --------------------------------------------------------------------------
# Test 5.1: Deterministic bootstrap comparison
# --------------------------------------------------------------------------
cat("\n--- Test 5.1: Deterministic bootstrap ---\n")

# First run r3d to get the fitted object
set.seed(42)
fit_boot <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
                q_grid = q_grid_20, kernel_fun = "epanechnikov")

# Generate xi matrix (n x B)
B_det <- 200
set.seed(123)
xi_mat <- matrix(rnorm(n_sharp * B_det), nrow = n_sharp, ncol = B_det)

# Save xi matrix for Stata
write.csv(xi_mat, "ref_xi_matrix.csv", row.names = FALSE)
cat("  Saved: ref_xi_matrix.csv (", n_sharp, "x", B_det, ")\n")

# Run deterministic bootstrap
boot_det <- r3d_bootstrap(
  object = fit_boot, X = X_sharp, Y_list = Y_list_sharp,
  B = B_det, alpha = 0.05,
  test = c("nullity", "homogeneity"),
  xi_mat = xi_mat
)

# Save deterministic bootstrap results
det_df <- data.frame(
  q = q_grid_20,
  cb_lower = boot_det$cb_lower,
  cb_upper = boot_det$cb_upper,
  crit_val = boot_det$crit_val
)
write.csv(det_df, "ref_bootstrap_deterministic.csv", row.names = FALSE)
cat("  Saved: ref_bootstrap_deterministic.csv\n")

# Save bootstrap tau matrix (B x nq)
write.csv(boot_det$boot_taus, "ref_bootstrap_taus.csv", row.names = FALSE)
cat("  Saved: ref_bootstrap_taus.csv\n")

# --------------------------------------------------------------------------
# Test 5.2: P-value comparison (deterministic)
# --------------------------------------------------------------------------
cat("\n--- Test 5.2: P-values (deterministic) ---\n")

pval_df <- data.frame(
  test = character(), p_value = numeric(), stringsAsFactors = FALSE
)

if (!is.null(boot_det$test_results)) {
  for (tname in names(boot_det$test_results)) {
    for (rname in names(boot_det$test_results[[tname]])) {
      tr <- boot_det$test_results[[tname]][[rname]]
      pval_df <- rbind(pval_df, data.frame(
        test = paste0(tname, ":", rname),
        p_value = tr$p_value,
        test_stat = tr$test_stat,
        stringsAsFactors = FALSE
      ))
    }
  }
}
write.csv(pval_df, "ref_pvalues_deterministic.csv", row.names = FALSE)
cat("  Saved: ref_pvalues_deterministic.csv\n")

# --------------------------------------------------------------------------
# Test 5.3: Stochastic bootstrap (for CI width comparison)
# --------------------------------------------------------------------------
cat("\n--- Test 5.3: Stochastic bootstrap ---\n")
set.seed(42)
fit_boot_stoch <- r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
                      q_grid = q_grid_20, kernel_fun = "epanechnikov",
                      boot = TRUE, boot_reps = 500, alpha = 0.05,
                      test = c("nullity", "homogeneity"))
save_results(fit_boot_stoch, "ref_results_bootstrap_sharp.csv")

# Save p-values
if (!is.null(fit_boot_stoch$boot_out$test_results)) {
  test_names <- c()
  p_vals <- c()
  for (tname in names(fit_boot_stoch$boot_out$test_results)) {
    for (rname in names(fit_boot_stoch$boot_out$test_results[[tname]])) {
      tr <- fit_boot_stoch$boot_out$test_results[[tname]][[rname]]
      test_names <- c(test_names, paste0(tname, ":", rname))
      p_vals <- c(p_vals, tr$p_value)
    }
  }
  pv <- data.frame(test = test_names, p_value = p_vals)
  write.csv(pv, "ref_pvalues_bootstrap_sharp.csv", row.names = FALSE)
  cat("  Saved: ref_pvalues_bootstrap_sharp.csv\n")
}


# ############################################################################
# LAYER 6: EDGE CASES & ROBUSTNESS
# ############################################################################
cat("\n\n=== LAYER 6: Edge Cases & Robustness ===\n")

# --------------------------------------------------------------------------
# Test 6.1: Small sample (n=30, nq=5)
# --------------------------------------------------------------------------
cat("\n--- Test 6.1: Small sample ---\n")
set.seed(42)
n_small <- 30
X_small <- runif(n_small, -1, 1)
Y_list_small <- lapply(seq_len(n_small), function(i) {
  mu <- 5 + 0.5 * X_small[i] + 3 * (X_small[i] >= 0)
  rnorm(20, mean = mu, sd = 0.5)
})
q_grid_5 <- (1:5) / 6

Qmat_small <- save_data_for_stata(X_small, Y_list_small, q_grid_5,
                                   "ref_data_small_sample.csv")

fit_small <- tryCatch({
  r3d(X_small, Y_list_small, cutoff = 0, method = "simple", p = 2,
      q_grid = q_grid_5, kernel_fun = "epanechnikov")
}, error = function(e) {
  cat("  Small sample error:", e$message, "\n")
  NULL
})
if (!is.null(fit_small)) {
  save_results(fit_small, "ref_results_small_sample.csv")
}

# --------------------------------------------------------------------------
# Test 6.2: Large nq (nq=200)
# --------------------------------------------------------------------------
cat("\n--- Test 6.2: Large nq ---\n")
set.seed(42)
q_grid_200 <- (1:200) / 201

Qmat_200 <- save_data_for_stata(X_sharp, Y_list_sharp, q_grid_200,
                                 "ref_data_large_nq.csv")

fit_200 <- tryCatch({
  r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
      q_grid = q_grid_200, kernel_fun = "epanechnikov")
}, error = function(e) {
  cat("  Large nq error:", e$message, "\n")
  NULL
})
if (!is.null(fit_200)) {
  save_results(fit_200, "ref_results_large_nq.csv")
}

# --------------------------------------------------------------------------
# Test 6.3: Polynomial order variation
# --------------------------------------------------------------------------
cat("\n--- Test 6.3: Polynomial order variation ---\n")

for (p_order in c(1, 3)) {
  cat("  p =", p_order, "\n")
  set.seed(42)
  fit_p <- tryCatch({
    r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = p_order,
        q_grid = q_grid_20, kernel_fun = "epanechnikov")
  }, error = function(e) {
    cat("    Error:", e$message, "\n")
    NULL
  })
  if (!is.null(fit_p)) {
    save_results(fit_p, paste0("ref_results_sharp_simple_p", p_order, ".csv"))
  }
}

# --------------------------------------------------------------------------
# Test 6.6: Coverage correction
# --------------------------------------------------------------------------
cat("\n--- Test 6.6: Coverage correction ---\n")
set.seed(42)
fit_cov <- tryCatch({
  r3d(X_sharp, Y_list_sharp, cutoff = 0, method = "simple", p = 2,
      q_grid = q_grid_20, kernel_fun = "epanechnikov", coverage = TRUE)
}, error = function(e) {
  cat("  Coverage error:", e$message, "\n")
  NULL
})
if (!is.null(fit_cov)) {
  save_results(fit_cov, "ref_results_coverage.csv")
}

# --------------------------------------------------------------------------
# Test 6.7: Paper DGP 1 (Normal-Normal)
# --------------------------------------------------------------------------
cat("\n--- Test 6.7: Paper DGP 1 (Normal-Normal) ---\n")
set.seed(42)
n_dgp <- 500
X_dgp1 <- runif(n_dgp, -1, 1)
Delta <- 2.27
Y_list_dgp1 <- lapply(seq_len(n_dgp), function(i) {
  above <- as.numeric(X_dgp1[i] >= 0)
  mu_i <- 5 + 5 * X_dgp1[i] + Delta * above
  sigma_i <- abs(1 + X_dgp1[i])
  if (sigma_i < 0.1) sigma_i <- 0.1
  rnorm(50, mean = mu_i, sd = sigma_i)
})
q_grid_dgp <- c(0.10, 0.25, 0.50, 0.75, 0.90)

Qmat_dgp1 <- save_data_for_stata(X_dgp1, Y_list_dgp1, q_grid_dgp,
                                   "ref_data_dgp1.csv")

fit_dgp1 <- r3d(X_dgp1, Y_list_dgp1, cutoff = 0, method = "simple", p = 2,
                 q_grid = q_grid_dgp, kernel_fun = "epanechnikov")
save_results(fit_dgp1, "ref_results_dgp1.csv")

# --------------------------------------------------------------------------
# Test 6.8: Paper DGP 2 (Normal-Exponential mixture)
# --------------------------------------------------------------------------
cat("\n--- Test 6.8: Paper DGP 2 (Normal-Exponential mixture) ---\n")
set.seed(42)
X_dgp2 <- runif(n_dgp, -1, 1)
Delta2 <- 1.86
Y_list_dgp2 <- lapply(seq_len(n_dgp), function(i) {
  above <- as.numeric(X_dgp2[i] >= 0)
  mu_i <- 2 * X_dgp2[i] + Delta2 * above
  lambda_i <- max(0.5, min(1.5, 1 + 0.5 * X_dgp2[i]))
  rnorm(25, mean = mu_i, sd = 1) + 2 * rexp(25, rate = lambda_i)
})

Qmat_dgp2 <- save_data_for_stata(X_dgp2, Y_list_dgp2, q_grid_dgp,
                                   "ref_data_dgp2.csv")

fit_dgp2 <- r3d(X_dgp2, Y_list_dgp2, cutoff = 0, method = "simple", p = 2,
                 q_grid = q_grid_dgp, kernel_fun = "epanechnikov")
save_results(fit_dgp2, "ref_results_dgp2.csv")


cat("\n\n=== Reference data generation complete ===\n")
cat("Total CSV files generated in tests/ directory.\n")
