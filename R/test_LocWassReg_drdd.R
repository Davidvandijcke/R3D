# Load necessary libraries
library(plotly)
library(osqp)
library(Matrix)

# Simulation parameters
set.seed(123)
threshold <- 0
x_below <- seq(-1, 0, 0.01)
x_above <- seq(0, 1, 0.01)

# Generate data below and above threshold
y_below <- lapply(x_below, function(x) rnorm(1000, mean = 5 * x^3, sd = 1))
y_above <- lapply(x_above, function(x) rnorm(1000, mean = 5 * x^3 + 7, sd = 2))

# Combine data
x <- c(x_below, x_above)
y <- c(y_below, y_above)

# Define function to compute quantile functions
compute_quantile <- function(y_values, probs = seq(0, 1, length.out = 100)) {
  quantile(y_values, probs = probs)
}

# Convert `y` to quantiles
q_below <- lapply(y_below, compute_quantile)
q_above <- lapply(y_above, compute_quantile)
q_all <- c(q_below, q_above)
qin <- do.call(rbind, q_all)

# Prepare input for LocWassRegDRDD
xin <- matrix(x, ncol = 1) # Predictor matrix with one column (x values)
xout <- xin                 # Use same x values for output predictions

# Run local Wasserstein regression with D-RDD estimator
# Set options, including bandwidth and kernel type
optns <- list(bw = 0.05, ker = "gauss", c = threshold, qSup = seq(0, 1, length.out = 100))

# Run the modified LocWassRegDRDD function, estimating from both sides
result <- LocWassRegDRDD(xin = xin, qin = qin, xout = xout, optns = optns, compute_density = FALSE)

# Extract results
qout_left <- result$left$qout   # Quantile outputs for "from the left" estimation
qout_right <- result$right$qout # Quantile outputs for "from the right" estimation

# Visualization (3D plot) of quantile regression results
fig <- plot_ly()

# Plot quantile curves from the left
for (i in 1:length(x_below)) {
  fig <- fig %>%
    add_trace(x = rep(x_below[i], ncol(qout_left)), 
              y = seq(0, 1, length.out = ncol(qout_left)), 
              z = qout_left[i, ], 
              type = 'scatter3d', mode = 'lines', name = 'Left')
}

# Plot quantile curves from the right
for (i in 1:length(x_above)) {
  fig <- fig %>%
    add_trace(x = rep(x_above[i], ncol(qout_right)), 
              y = seq(0, 1, length.out = ncol(qout_right)), 
              z = qout_right[i, ], 
              type = 'scatter3d', mode = 'lines', name = 'Right')
}

# Customize 3D plot
fig <- fig %>%
  layout(scene = list(
    xaxis = list(title = 'Running Variable (X)'),
    yaxis = list(title = 'Quantile Level (t)'),
    zaxis = list(title = "Quantile Value (Q)")
  ),
  showlegend = FALSE)

# Display plot
fig

# Optionally, save the plot
# fig_path <- "quantile_regression_plot.png"
# save_image(fig, fig_path, scale = 3)
