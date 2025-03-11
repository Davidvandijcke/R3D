# File: tests/testthat/test-r3d-features.R
context("R3D Package Feature Tests")

library(testthat)
library(R3D)
library(dplyr, warn.conflicts = FALSE)

###############################################################################
# 1) Functions to generate fake data for testing
###############################################################################

# (A) Simple test data with known constant effect = 2
#    Below:   mean =  5 + 0.5*x
#    Above:   mean = (5 + 0.5*x) + 2
create_test_data <- function(n = 500, seed = 123) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  
  # Known, uniform effect of +2 for x>=0
  effect_size <- 2
  y_list <- vector("list", n)
  for(i in seq_len(n)) {
    if(x[i] < 0) {
      y_list[[i]] <- rnorm(200, mean = 5 + 0.5*x[i], sd = 0.5)
    } else {
      y_list[[i]] <- rnorm(200, mean = 5 + 0.5*x[i] + effect_size, sd = 0.5)
    }
  }
  
  # Optional fuzzy RD for testing
  p_t <- 0.1 + 0.8*(x>0)  # e.g. 10% below, 90% above
  p_t <- pmax(0, pmin(1, p_t))
  t   <- rbinom(n, 1, p_t)
  
  list(
    x = x,
    y_list = y_list,
    t = t,
    true_effect = effect_size  # exactly 2 at all quantiles
  )
}


# (B) Heterogeneous effect data
#    We do a 3-part mixture for x>=0, and a simple N(...) for x<0.
#    Then we approximate the distribution-level effect by a large “internal” sample.

create_hetero_effect_data <- function(n = 1000, seed = 456, big_reps = 50000) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  
  # Actually generate the sample data
  y_list <- vector("list", n)
  for(i in seq_len(n)) {
    if(x[i] < 0) {
      # Control group
      y_list[[i]] <- rnorm(100, mean = 2 + x[i], sd = 0.5)
    } else {
      # Treatment group => mixture of 3 subgroups
      y_list[[i]] <- c(
        rnorm(30, mean = 2 + x[i] + 4, sd = 0.8),  # big shift on lower quantiles
        rnorm(40, mean = 2 + x[i] + 2, sd = 0.5),
        rnorm(30, mean = 2 + x[i] + 1, sd = 1)
      )
    }
  }
  
  # -----------------------------------------------------------------------
  # Approximate the *true distributional effect* by a large replicate.
  # We'll simulate from the design *as if* x=0 for below vs. x=0 for above,
  # i.e. ignoring the slope to highlight "pure" effect differences.
  # If you prefer, you can incorporate x's distribution for a partial effect,
  # but here's a simpler approach so that "q-> effect" is clearer.
  # -----------------------------------------------------------------------
  
  # Big replicate for "below" design:
  # => rnorm(...) with mean=2 + xBelow + ... or a single, simpler assumption
  # Actually let's do the same logic: (x<0) => normal(2,1). 
  # Then (x>=0) => mixture of means: 2+4, 2+2, 2+1, etc.
  
  set.seed(9999)
  # "Below" distribution:
  big_below <- rnorm(big_reps, mean = 2, sd = 1)
  # "Above" distribution => 3-part mixture
  #   30% => mean=6, sd=0.8
  #   40% => mean=4, sd=1
  #   30% => mean=3, sd=1.2
  # We'll sample them in proportion and combine.
  
  idxs <- sample.int(3, size=big_reps, replace=TRUE, prob=c(0.3,0.4,0.3))
  big_above <- numeric(big_reps)
  big_above[idxs==1] <- rnorm(sum(idxs==1), mean=6, sd=0.8)
  big_above[idxs==2] <- rnorm(sum(idxs==2), mean=4, sd=1)
  big_above[idxs==3] <- rnorm(sum(idxs==3), mean=3, sd=1.2)
  
  # Approximate difference of average quantiles from these large samples
  # We'll create a grid of quantiles and store the difference.
  q_grid <- seq(0.01, 0.99, 0.01)
  
  # Empirical quantiles
  q_below  <- quantile(big_below,  probs=q_grid)
  q_above  <- quantile(big_above,  probs=q_grid)
  diff_q    <- q_above - q_below   # "true" difference in distribution
  
  # We'll store this in a named vector for easy reference
  names(diff_q) <- as.character(q_grid)
  
  list(
    x          = x,
    y_list     = y_list,
    true_effects = diff_q   # A named numeric vector for many quantiles
  )
}

################################################################################
# 1) Test r3d_bwselect() basic usage
################################################################################

test_that("r3d_bwselect() works for method='simple' and 'frechet'", {
  set.seed(101)
  # Example data
  n <- 50
  x <- runif(n, -1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, mean=2+0.5*x[i]))
  
  # Simple
  bwres_simp <- r3d_bwselect(x, y_list, method="simple", p=1)
  expect_equal(bwres_simp$method, "simple")
  expect_true(is.numeric(bwres_simp$h_star))
  expect_true(length(bwres_simp$h_star) > 1)  # per-quantile => vector
  
  # Frechet
  bwres_frech <- r3d_bwselect(x, y_list, method="frechet", p=1)
  expect_equal(bwres_frech$method, "frechet")
  expect_true(is.numeric(bwres_frech$h_star))
  expect_length(bwres_frech$h_star, 1)     # single IMSE bandwidth
})

test_that("r3d_bwselect() properly handles edge cases", {
  set.seed(1010)
  n <- 30
  x <- runif(n, -1, 1)
  
  # Case 1: All observations on one side of cutoff
  x_one_side <- abs(x) # all positive
  y_list <- lapply(seq_len(n), function(i) rnorm(20, mean=3+x_one_side[i]))
  
  # Should give warning but still return reasonable bandwidth
  expect_error(
    r3d_bwselect(x_one_side, y_list, method="simple")
  )
  
  # Case 2: Very large variance in distributions
  y_list_var <- lapply(seq_len(n), function(i) {
    if(x[i] < 0) {
      rnorm(20, mean=2, sd=1)
    } else {
      rnorm(20, mean=2, sd=10) # very large variance
    }
  })
  
  # Should handle the heteroskedasticity properly
  bw_var <- r3d_bwselect(x, y_list_var, method="frechet")
  expect_true(is.numeric(bw_var$h_star))
  expect_false(is.na(bw_var$h_star))
})

################################################################################
# 2) Test r3d() usage for method='simple' and 'frechet' (sharp vs fuzzy)
################################################################################

test_that("r3d() runs for simple (sharp) with no errors", {
  set.seed(102)
  n <- 50
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, 2+0.5*x[i]))
  out_simp <- r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1, boot=FALSE)
  
  # Basic class and structure checks
  expect_s3_class(out_simp, "r3d")
  expect_equal(out_simp$method, "simple")
  expect_false("boot_out" %in% names(out_simp$bootstrap))  # no bootstrap
  
  # Check dimensions
  expect_equal(length(out_simp$tau), length(out_simp$q_grid))
  
  # Check that weights arrays were correctly computed
  expect_true("w_plus" %in% names(out_simp))
  expect_true("w_minus" %in% names(out_simp))
  expect_equal(dim(out_simp$w_plus), c(n, length(out_simp$q_grid)))
  
  # Check nested structure
  expect_true("results" %in% names(out_simp))
  expect_true("coefficients" %in% names(out_simp))
  expect_true("bootstrap" %in% names(out_simp))
  expect_true("inputs" %in% names(out_simp))
  expect_true("diagnostics" %in% names(out_simp))
  
  # Check key results are accessible in both ways (top level and nested)
  expect_equal(out_simp$tau, out_simp$results$tau)
  expect_equal(out_simp$q_grid, out_simp$results$q_grid)
})

test_that("r3d() runs for simple (fuzzy) with no errors", {
  set.seed(103)
  n <- 60
  x <- runif(n, -1, 1)
  Tvar <- as.integer(x>=0)  # simple fuzzy
  y_list <- lapply(seq_len(n), function(i) rnorm(40, 2 + 0.7*x[i]))
  out_simp_fuzzy <- r3d(X=x, Y_list=y_list, T=Tvar, fuzzy=TRUE, method="simple", p=2)
  
  # Basic checks
  expect_s3_class(out_simp_fuzzy, "r3d")
  expect_true(out_simp_fuzzy$fuzzy)
  
  # Check that alphaT_plus and alphaT_minus are computed
  expect_true("alphaT_plus" %in% names(out_simp_fuzzy))
  expect_true("alphaT_minus" %in% names(out_simp_fuzzy))
  expect_false(is.null(out_simp_fuzzy$alphaT_plus))
  expect_false(is.null(out_simp_fuzzy$alphaT_minus))
  
  # Check fuzzy-specific elements
  expect_true(!is.null(out_simp_fuzzy$diagnostics$denominator))
  expect_true(!is.null(out_simp_fuzzy$bootstrap$e2_mat))
})

test_that("r3d() runs for frechet (sharp) with no errors", {
  set.seed(104)
  n <- 70
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, 3 + 0.2*x[i]))
  out_frech <- r3d(X=x, Y_list=y_list, method="frechet", p=1, boot=FALSE)
  
  # Basic checks
  expect_s3_class(out_frech, "r3d")
  expect_equal(out_frech$method, "frechet")
  
  # Check IMSE bandwidth (should be a single value)
  expect_length(out_frech$bandwidths, 1)
  
  # For Frechet method, check that tau_vec is monotonic (due to isotonic regression)
  if(length(out_frech$int_minus) > 1) {
    is_monotonic <- all(diff(out_frech$int_minus) >= -1e-10) # allow for numerical imprecision
    expect_true(is_monotonic)
  }
})

test_that("r3d() runs for frechet (fuzzy) with no errors", {
  set.seed(105)
  n <- 60
  x <- runif(n, -1, 1)
  Tvar <- (x>0)*1
  y_list <- lapply(seq_len(n), function(i) rexp(50, rate=0.4 + 0.2*(x[i]>0)))
  out_frech_fuzzy <- r3d(X=x, Y_list=y_list, T=Tvar, fuzzy=TRUE, method="frechet", p=1)
  
  # Basic checks
  expect_s3_class(out_frech_fuzzy, "r3d")
  expect_true(out_frech_fuzzy$fuzzy)
  
  # For Frechet method, check that tau_vec is monotonic (due to isotonic regression)
  if(length(out_frech_fuzzy$int_plus) > 1) {
    is_monotonic <- all(diff(out_frech_fuzzy$int_plus) >= -1e-10) # allow for numerical imprecision
    expect_true(is_monotonic)
  }
})

################################################################################
# 3) Test bootstrap (boot=TRUE) with method='simple', check we have boot_out
################################################################################

test_that("r3d() boot=TRUE yields boot_out with confidence bands & tests", {
  set.seed(106)
  n <- 40
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(20, 1+ x[i]))
  
  # Small example => few reps => alpha=0.2, test="nullity"
  out_boot <- r3d(X=x, Y_list=y_list, method="simple", p=1,
                  boot=TRUE, boot_reps=10, test="nullity", alpha=0.2)
  
  # Basic checks
  expect_s3_class(out_boot, "r3d")
  expect_true("boot_out" %in% names(out_boot))
  expect_true("boot_out" %in% names(out_boot$bootstrap))
  expect_identical(out_boot$boot_out, out_boot$bootstrap$boot_out)
  
  # Check bootstrap results structure
  bo <- out_boot$boot_out
  expect_true(all(c("cb_lower","cb_upper","test_stat","test_crit_val","p_value") %in% names(bo)))
  expect_equal(length(bo$cb_lower), length(out_boot$q_grid))
  expect_equal(length(bo$cb_upper), length(out_boot$q_grid))
  
  # Test bootstrap matrix dimensions
  expect_equal(dim(bo$boot_taus), c(length(out_boot$q_grid), 10)) # 10 bootstrap reps
  
  # Test confidence band properties
  expect_true(all(bo$cb_upper >= out_boot$tau))
  expect_true(all(bo$cb_lower <= out_boot$tau))
})

test_that("r3d() with boot=TRUE and test='homogeneity' works properly", {
  set.seed(1066)
  
  # Create data with heterogeneous treatment effects
  hetero_data <- create_hetero_effect_data(n=60)
  q_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  # Run r3d with homogeneity test
  out_test <- r3d(X=hetero_data$x, Y_list=hetero_data$y_list, cutoff=0, 
                  method="simple", p=1, q_grid=q_grid,
                  boot=TRUE, boot_reps=20, test="homogeneity")
  
  # Check test results
  expect_true(!is.na(out_test$boot_out$test_stat))
  expect_true(!is.na(out_test$boot_out$p_value))
  
  # With heterogeneous effects, we expect a low p-value
  # (though with only 20 bootstrap reps, it might not be significant)
  # This is more of a directional check than a strict test
  expect_true(out_test$boot_out$test_stat > 0)
})

# Test direct bootstrap function
test_that("r3d_bootstrap() can be called directly with r3d object", {
  set.seed(107)
  n <- 40
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(20, 1+ x[i]))
  
  # Create r3d object without bootstrap
  r3d_obj <- r3d(X=x, Y_list=y_list, method="simple", p=1, boot=FALSE)
  
  # Call bootstrap directly
  bo <- r3d_bootstrap(object=r3d_obj, X=x, Y_list=y_list, 
                      B=5, alpha=0.1, test="homogeneity", cores=1)
  
  # Check structure
  expect_true(all(c("cb_lower", "cb_upper", "p_value") %in% names(bo)))
  expect_equal(length(bo$cb_lower), length(r3d_obj$q_grid))
  expect_true(!is.na(bo$p_value))
})

################################################################################
# 4) Scenario-based test with random distributions
################################################################################

delta_mu <- 2
delta_sigma <- 0.5
base_mu_below <- 5
slope_mu_below<- 0.5
base_sigma_below<-1
slope_sigma_below<-0.2
c_cutoff <- 0

simulate_data_scenario1 <- function(N, n_obs) {
  x <- runif(N, -1, 1)
  y <- vector("list", N)
  mu_at_c <- base_mu_below + slope_mu_below*c_cutoff
  sigma_at_c<- base_sigma_below + slope_sigma_below*c_cutoff
  for(i in seq_len(N)){
    if(x[i]< c_cutoff){
      mu    <- rnorm(1, mean=base_mu_below + slope_mu_below*x[i], sd=1)
      sigma <- abs(rnorm(1, mean= base_sigma_below + slope_sigma_below*x[i], sd=0.5))
    } else {
      mu    <- rnorm(1, mean= mu_at_c + delta_mu, sd=1)
      sigma <- abs(rnorm(1, mean= sigma_at_c + delta_sigma, sd=0.5))
    }
    y[[i]]<- rnorm(n_obs, mean=mu, sd=sigma)
  }
  list(x=x, y=y)
}

true_treatment_effect_scenario1 <- function(q_levels) {
  mu_at_c <- base_mu_below + slope_mu_below*c_cutoff
  sigma_at_c<- base_sigma_below + slope_sigma_below*c_cutoff
  mu_below <- mu_at_c
  sigma_below<-sigma_at_c
  mu_above <- mu_below + delta_mu
  sigma_above<- sigma_below + delta_sigma
  q_below <- qnorm(q_levels, mean=mu_below, sd=sigma_below)
  q_above <- qnorm(q_levels, mean=mu_above, sd=sigma_above)
  q_above - q_below
}

test_that("Scenario 1 distribution test => simple & frechet have finite bias", {
  set.seed(108)
  n_rep <- 5  # smaller for test speed
  N     <- 80
  n_obs <- 30
  q_levels<- seq(0.1, 0.9, by=0.1)
  true_te <- true_treatment_effect_scenario1(q_levels)
  
  results_simp <- list()
  results_frech <- list()
  for(r in seq_len(n_rep)){
    dat <- simulate_data_scenario1(N, n_obs)
    # run r3d simple
    out_simp <- r3d(X=dat$x, Y_list=dat$y, cutoff=c_cutoff, 
                    method="simple", p=1, q_grid=q_levels, boot=FALSE)
    # run r3d frechet
    out_frech <- r3d(X=dat$x, Y_list=dat$y, cutoff=c_cutoff,
                     method="frechet", p=1, q_grid=q_levels, boot=FALSE)
    results_simp[[r]] <- out_simp$tau
    results_frech[[r]] <- out_frech$tau
  }
  # convert to matrix => row=rep, col=quantile
  simp_mat  <- do.call(rbind, results_simp)
  frech_mat <- do.call(rbind, results_frech)
  
  # compute average error
  avg_bias_simp <- mean(rowMeans(simp_mat - matrix(true_te, nrow=n_rep, ncol=length(q_levels), byrow=TRUE)), na.rm=TRUE)
  avg_bias_frech <- mean(rowMeans(frech_mat - matrix(true_te, nrow=n_rep, ncol=length(q_levels), byrow=TRUE)), na.rm=TRUE)
  
  # check that the absolute bias is not huge
  expect_lt(abs(avg_bias_simp), 2)
  expect_lt(abs(avg_bias_frech), 2)
})

test_that("Known treatment effect is recovered by both methods", {
  set.seed(1088)
  
  # Create data with known constant treatment effect
  test_data <- create_test_data(n=100)
  q_grid <- seq(0.1, 0.9, by=0.1)
  
  # Run both methods
  out_simp <- r3d(X=test_data$x, Y_list=test_data$y_list, cutoff=0, 
                  method="simple", p=1, q_grid=q_grid, boot=FALSE)
  
  out_frech <- r3d(X=test_data$x, Y_list=test_data$y_list, cutoff=0, 
                   method="frechet", p=1, q_grid=q_grid, boot=FALSE)
  
  # Check recovery of treatment effect (allowing modest error)
  expect_true(mean(abs(out_simp$tau - test_data$true_effect)) < 0.5)
  expect_true(mean(abs(out_frech$tau - test_data$true_effect)) < 0.5)
})

################################################################################
# 5) Test edge cases and error handling
################################################################################

test_that("r3d handles very small samples with warning", {
  # Create tiny dataset
  x <- c(-0.5, 0.5, -0.2, 0.2)
  y_list <- list(
    rnorm(5, 1),
    rnorm(5, 3),
    rnorm(5, 1.5),
    rnorm(5, 2.5)
  )
  
  # Should give warning but not error
  expect_warning(
    r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1),
    "small sample"
  )
})

test_that("r3d errors correctly when fuzzy=TRUE but T=NULL", {
  set.seed(110)
  n <- 20
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(10, 2 + 0.5*x[i]))
  
  # Should error
  expect_error(
    r3d(X=x, Y_list=y_list, cutoff=0, method="simple", fuzzy=TRUE, T=NULL),
    "fuzzy.*requires T"
  )
})

test_that("r3d handles empty distributions in Y_list", {
  set.seed(111)
  n <- 20
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(10, 2 + 0.5*x[i]))
  
  # Add an empty distribution
  y_list[[5]] <- numeric(0)
  
  # Should give warning but not error
  expect_warning(
    result <- r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1),
    "Empty"
  )
  
  # Result should still be valid
  expect_s3_class(result, "r3d")
})

################################################################################
# 6) Test summary.r3d() and plot.r3d() methods
################################################################################

test_that("summary.r3d() and plot.r3d() work properly", {
  skip_on_cran()  # Skip on CRAN since these are visual tests
  
  set.seed(111)
  n <- 40
  x <- runif(n, -1, 1)
  y_list <- lapply(seq_len(n), function(i) rnorm(20, 2 + 0.3*x[i]))
  
  # Do a small bootstrap too
  out <- r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1,
             q_grid=seq(0.25, 0.75, by=0.25),
             boot=TRUE, boot_reps=5, test="none")
  
  # Check summary method
  expect_output(summary(out), "Method")
  expect_output(summary(out), "Polynomial order")
  expect_output(summary(out), "Uniform Confidence Bands")
  
  # Check print method
  expect_output(print(out), "R3D: Regression Discontinuity with Distributional Outcomes")
  expect_output(print(out), "Treatment Effect Summary")
  
  # Check plot method visually (this just verifies it runs without error)
  pdf(NULL)  # Redirect plot output to null device 
  expect_silent(plot(out, main="Test Plot R3D"))
  dev.off()
})


###############################################################################
# 2) Visual validation checks for numeric testing
###############################################################################

test_that("Visual validation with known treatment effects", {
  skip_on_cran()  # Skip on CRAN - purely for local/manual checking
  
  set.seed(777)
  
  #---------------------------------------------------------------------------
  # 1) CONSTANT-EFFECT DATA
  #---------------------------------------------------------------------------
  
  test_data <- create_test_data(n = 500)
  
  # We want a 100-point grid
  q_grid <- seq(0.1, 0.9, by=0.01)
  
  # r3d with Frechet approach
  out_constant <- r3d(
    X        = test_data$x,
    Y_list   = test_data$y_list,
    cutoff   = 0,
    method   = "simple",
    p        = 2,
    s        = 1,
    q_grid   = q_grid,
    boot     = TRUE,
    boot_reps= 1000,
    test     = "homogeneity"
  )
  
  # We know the true effect is exactly 2 at all quantiles:
  expected_constant <- rep(test_data$true_effect, length(q_grid))
  
  #---------------------------------------------------------------------------
  # 2) HETEROGENEOUS-EFFECT DATA
  #---------------------------------------------------------------------------
  
  hetero_data <- create_hetero_effect_data(n = 500)
  
  # We have the same q_grid = seq(0.01, 0.99, 0.01).
  # We'll compare our estimate to hetero_data$true_effects, which is named for 0.01..0.99
  out_hetero <- r3d(
    X         = hetero_data$x,
    Y_list    = hetero_data$y_list,
    cutoff    = 0,
    method    = "frechet",
    p         = 2,
    s         = 1,
    q_grid    = q_grid,
    boot      = TRUE,
    boot_reps = 1000,
    test      = "homogeneity"
  )
  
  # Now for each q in q_grid, we'll see the difference we expect from the big replicate
  # stored in hetero_data$true_effects (which is a named vector).
  # If any exact matching fails because of floating float in names, we can do:
  q_strs <- sprintf("%.2f", q_grid) # character
  expected_hetero <- as.numeric(hetero_data$true_effects[q_strs])
  
  #---------------------------------------------------------------------------
  # Print results for manual checking
  #---------------------------------------------------------------------------
  cat("\n\n==== VISUAL VALIDATION ====\n\n")
  
  # (A) Constant Effect
  cat("CONSTANT EFFECT MODEL:\n")
  cat("The effect is exactly:", test_data$true_effect, "for all q.\n\n")
  
  res_const <- data.frame(
    Quantile  = q_grid,
    Estimated = round(out_constant$tau, 2),
    Expected  = round(expected_constant, 2),
    Bias      = round(out_constant$tau - expected_constant, 2),
    LowerCB   = round(out_constant$boot_out$cb_lower, 2),
    UpperCB   = round(out_constant$boot_out$cb_upper, 2)
  )
  print(head(res_const, 10))   # Print just the first few lines
  cat("...\n")
  cat("Homogeneity test p-value:", out_constant$boot_out$p_value, "\n")
  cat("(Should be > 0.05 if effect is truly constant.)\n\n")
  
  # (B) Heterogeneous Effect
  cat("HETEROGENEOUS EFFECT MODEL:\n")
  cat("(We used a big replicate to approximate these differences)\n\n")
  
  res_hetero <- data.frame(
    Quantile  = q_grid,
    Estimated = round(out_hetero$tau, 2),
    Expected  = round(expected_hetero, 2),
    Bias      = round(out_hetero$tau - expected_hetero, 2),
    LowerCB   = round(out_hetero$boot_out$cb_lower, 2),
    UpperCB   = round(out_hetero$boot_out$cb_upper, 2)
  )
  print(head(res_hetero, 10))
  cat("...\n")
  cat("Homogeneity test p-value:", out_hetero$boot_out$p_value, "\n")
  cat("(Should be < 0.05 since effect varies by quantile.)\n\n")
  
  #---------------------------------------------------------------------------
  # Quick Plots
  #---------------------------------------------------------------------------
  par(mfrow=c(1,2))
  
  # Plot 1: Constant effect
  plot(out_constant, main="Constant Treatment Effect (Frechet)", 
       col="blue", lwd=2)
  abline(h=test_data$true_effect, col="red", lty=2, lwd=2)
  legend("topright", c("Estimated", "True (2)"), 
         col=c("blue", "red"), lty=c(1,2), lwd=2)
  
  # Plot 2: Heterogeneous effect
  plot(out_hetero, main="Heterogeneous Treatment Effect (Frechet)", 
       col="blue", lwd=2)
  lines(q_grid, expected_hetero, col="red", lty=2, lwd=2)
  legend("topright", c("Estimated", "Approx. True"), 
         col=c("blue", "red"), lty=c(1,2), lwd=2)
  
  par(mfrow=c(1,1))
  
  # No formal expectations, just no errors
  expect_true(TRUE)
})