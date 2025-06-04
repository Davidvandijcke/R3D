# R3D: Regression Discontinuity Design with Distribution-Valued Outcomes <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/dvdijcke/R3D/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/R3D/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/R3D)](https://CRAN.R-project.org/package=R3D)
[![R-CMD-check](https://github.com/Davidvandijcke/r3d/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Davidvandijcke/r3d/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Adding a New Dimension to RDD! 🚀

R3D brings your regression discontinuity designs into the third dimension by handling *distributional outcomes* instead of just scalar responses. When your treatment effect isn't just a single number but an entire distribution shift, R3D has you covered!

### Why "R3D"? 🤔

- **3 D's**: **D**iscontinuity **D**esign with **D**istributions
- **3D**: Because we're adding a new dimension (quantiles) to traditional RDD
- **R-3D**: It's like watching your results in 3D, but with R! 🥽

## Installation 📦

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("Davidvandijcke/R3D")

# Or once on CRAN:
# install.packages("R3D")
```

## Features ✨

- 📏 Estimates treatment effects across the entire distribution
- 🎯 Handles both sharp and fuzzy RD designs
- 🧮 Choose between simple quantile-by-quantile or global Fréchet approaches
- 📊 Beautiful visualizations of distributional treatment effects
- 🔍 Uniform inference with multiplier bootstrap
- 🧪 Tests of effect nullity and homogeneity

## Quick Example 📋

```r
library(R3D)

# Simulate data
set.seed(123)
n <- 100
X <- runif(n, -1, 1)  # Running variable
Y_list <- lapply(seq_len(n), function(i) {
  # Distribution is Normal with mean depending on X
  rnorm(sample(30:50, 1), mean = 2 + 2 * (X[i] >= 0))
})

# Sharp RDD with distributional outcome
fit <- r3d(X, Y_list, cutoff = 0, 
           method = "frechet", p = 2,
           boot = TRUE, boot_reps = 200)

# Visualize the results
plot(fit)

# Examine detailed results
summary(fit)
```

## How It Works 🔧

R3D uses advanced local polynomial and Fréchet regression techniques to estimate how an entire distribution changes at a discontinuity threshold when treatment is at a higher level of aggregation than the outcome variable. Instead of just estimating E[Y|X] at the cutoff, we estimate the entire *average* conditional quantile function E[Q_Y(τ|X)] at the threshold. Note the word average! Unlike traditional quantile RDD, we are now sampling *random* distributions, which
requires a different approach to estimation and inference. For more details on this fascinating subject, check out the supporting paper by Van Dijcke (that's me!) (2025 (that's now!)).

## Citation 📄

If you use R3D in your research, please cite:

```
Van Dijcke, D. (2025). Regression Discontinuity Design with Distributional Outcomes.
Working paper.
```

## Learn More 📚

Check out the [full documentation](https://Davidvandijcke.github.io/R3D/) for tutorials, examples, and detailed function references.

## Contributing 🤝

Contributions welcome! Feel free to submit issues or pull requests.
