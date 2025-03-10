# File: tests/testthat/test-r3d-features.R
context("R3D Package Feature Tests")

library(testthat)
library(R3D)  # your R3D package
library(dplyr)
library(data.table)

################################################################################
# 1) Test bwselect() basic usage
################################################################################

test_that("bwselect() works for method='simple' and 'frechet'", {
  set.seed(101)
  # Example data
  n <- 50
  x <- runif(n, -1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, mean=2+0.5*x[i]))
  
  # Simple
  bwres_simp <- bwselect(x, y_list, method="simple", p=1)
  expect_equal(bwres_simp$method, "simple")
  expect_true(is.numeric(bwres_simp$h_star))
  expect_true(length(bwres_simp$h_star)>1)  # per-quantile => vector
  
  # Frechet
  bwres_frech <- bwselect(x, y_list, method="frechet", p=1)
  expect_equal(bwres_frech$method, "frechet")
  expect_true(is.numeric(bwres_frech$h_star))
  expect_length(bwres_frech$h_star, 1)     # single IMSE bandwidth
})


################################################################################
# 2) Test r3d() usage for method='simple' and 'frechet' (sharp vs fuzzy)
################################################################################

test_that("r3d() runs for simple (sharp) with no errors", {
  set.seed(102)
  n <- 50
  x <- runif(n,-1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, 2+0.5*x[i]))
  out_simp <- r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1, boot=FALSE)
  expect_s3_class(out_simp, "r3d")
  expect_equal(out_simp$method, "simple")
  expect_false("boot_out" %in% names(out_simp))  # no bootstrap
  # check dims
  expect_equal(length(out_simp$tau), length(out_simp$q_grid))
})


test_that("r3d() runs for simple (fuzzy) with no errors", {
  set.seed(103)
  n <- 60
  x <- runif(n, -1,1)
  Tvar <- as.integer(x>=0)  # simple fuzzy
  y_list <- lapply(seq_len(n), function(i) rnorm(40, 2 + 0.7*x[i]))
  out_simp_fuzzy <- r3d(X=x, Y_list=y_list, T=Tvar, fuzzy=TRUE, method="simple", p=2)
  expect_s3_class(out_simp_fuzzy, "r3d")
  expect_true(out_simp_fuzzy$fuzzy)
})


test_that("r3d() runs for frechet (sharp) with no errors", {
  set.seed(104)
  n <- 70
  x <- runif(n,-1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(30, 3 + 0.2*x[i]))
  out_frech <- r3d(X=x, Y_list=y_list, method="frechet", p=1, boot=FALSE)
  expect_s3_class(out_frech, "r3d")
  expect_equal(out_frech$method, "frechet")
})


test_that("r3d() runs for frechet (fuzzy) with no errors", {
  set.seed(105)
  n <- 60
  x <- runif(n, -1,1)
  Tvar <- (x>0)*1
  y_list <- lapply(seq_len(n), function(i) rexp(50, rate=0.4 + 0.2*(x[i]>0)))
  out_frech_fuzzy <- r3d(X=x, Y_list=y_list, T=Tvar, fuzzy=TRUE, method="frechet", p=2)
  expect_s3_class(out_frech_fuzzy, "r3d")
  expect_true(out_frech_fuzzy$fuzzy)
})


################################################################################
# 3) Test bootstrap (boot=TRUE) with method='simple', check we have boot_out
################################################################################

test_that("r3d() boot=TRUE yields boot_out with confidence bands & tests", {
  set.seed(106)
  n <- 40
  x <- runif(n,-1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(20, 1+ x[i]))
  # small example => few reps => alpha=0.2, test="nullity"
  out_boot <- r3d(X=x, Y_list=y_list, method="simple", p=1,
                  boot=TRUE, boot_reps=10, test="nullity", alpha=0.2)
  expect_s3_class(out_boot, "r3d")
  expect_true("boot_out" %in% names(out_boot))
  bo <- out_boot$boot_out
  expect_true(all(c("cb_lower","cb_upper","test_stat","test_crit_val","p_value") %in% names(bo)))
  expect_equal(length(bo$cb_lower), length(out_boot$q_grid))
  expect_equal(length(bo$cb_upper), length(out_boot$q_grid))
})

################################################################################
# 4) Scenario-based test with random distributions (like your prior approach)
#    We'll do just scenario 1 for brevity. You can replicate for scenario 2, 3.
################################################################################

delta_mu <- 2
delta_sigma <- 0.5
base_mu_below <- 5
slope_mu_below<- 0.5
base_sigma_below<-1
slope_sigma_below<-0.2
c_cutoff <- 0

simulate_data_scenario1 <- function(N, n_obs) {
  x <- runif(N, -1,1)
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
  set.seed(107)
  n_rep <- 20  # smaller for speed
  N     <- 80
  n_obs <- 30
  q_levels<- seq(0.1,0.9,by=0.1)
  true_te <- true_treatment_effect_scenario1(q_levels)
  
  results_simp<- list()
  results_frech<- list()
  for(r in seq_len(n_rep)){
    dat <- simulate_data_scenario1(N, n_obs)
    # run r3d simple
    out_simp<- r3d(X=dat$x, Y_list=dat$y, cutoff=c_cutoff, 
                   method="simple", p=1, q_grid=q_levels, boot=FALSE)
    # run r3d frechet
    out_frech<- r3d(X=dat$x, Y_list=dat$y, cutoff=c_cutoff,
                    method="frechet", p=1, q_grid=q_levels, boot=FALSE)
    results_simp[[r]] <- out_simp$tau
    results_frech[[r]]<- out_frech$tau
  }
  # convert to matrix => row=rep, col=quantile
  simp_mat  <- do.call(rbind, results_simp)
  frech_mat <- do.call(rbind, results_frech)
  # compute average error
  avg_bias_simp<- mean(rowMeans(simp_mat - matrix(true_te, nrow=n_rep, ncol=length(q_levels), byrow=TRUE)))
  avg_bias_frech<-mean(rowMeans(frech_mat- matrix(true_te, nrow=n_rep, ncol=length(q_levels), byrow=TRUE)))
  # check that the absolute bias is not huge
  expect_lt(abs(avg_bias_simp), 2)
  expect_lt(abs(avg_bias_frech),2)
})

################################################################################
# 5) Test summary.r3d() and plot.r3d() run without error and produce output
################################################################################

test_that("summary.r3d() and plot.r3d() do not crash", {
  set.seed(108)
  n <- 40
  x <- runif(n,-1,1)
  y_list <- lapply(seq_len(n), function(i) rnorm(20, 2 + 0.3*x[i]))
  # do a small bootstrap too
  out <- r3d(X=x, Y_list=y_list, cutoff=0, method="simple", p=1,
             q_grid=seq(0.25,0.75,by=0.25),
             boot=TRUE, boot_reps=5, test="none")
  # check summary
  expect_silent({
    sm <- summary(out, samples=c(0.25,0.5))
    print(sm)
  })
  # check plot
  expect_silent({
    plot(out, main="Test Plot R3D")
  })
})

