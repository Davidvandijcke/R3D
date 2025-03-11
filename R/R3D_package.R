#' @title R3D: Regression Discontinuity with Distributional Outcomes
#' @name R3D-package
#' @docType package
#'
#' @description
#' The \pkg{R3D} package provides methods to estimate, infer, and visualize 
#' a new class of regression discontinuity designs where the outcome is 
#' a distribution rather than a single scalar. It includes:
#' \itemize{
#'   \item \code{\link{r3d}}: main estimation function
#'   \item \code{\link{r3d_bwselect}}: bandwidth selection
#'   \item \code{\link{r3d_bootstrap}}: multiplier bootstrap for uniform inference
#'   \item S3 methods: \code{\link{summary.r3d}}, \code{\link{plot.r3d}}, \code{\link{print.r3d}}
#' }
#'
#' @author 
#' Your Name (David Van Dijcke) \email{dvdijcke@umich.edu}
#'
#' @references
#' Van Dijcke, D. (2025). \emph{Regression Discontinuity Design with Distributional Outcomes.}
#' Working paper. 
#'
#' @keywords package
NULL