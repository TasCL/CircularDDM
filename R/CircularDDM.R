#' Circular Drift-diffusion Model
#'
#' Circular drift-diffusion model for continuous report. This package
#' implements the density function and the random number generator for circular
#' drift-diffusion model described in Smith (2016). For circular
#' drift-diffusion model, the random number refers to response time and
#' response angles. Both are continuous. Most functions are implemented in
#' Armadillo C++ via RcppArmadillo.
#'
#' @keywords CircularDDM
#' @name CircularDDM
#' @docType package
#' @author  Yi-Shin Lin <yishin.lin@utas.edu.au> \cr
#' Andrew Heathcote <andrew.heathcote@utas.edu.au> \cr
#' Peter Kvam <kvam.peter@gmail.com> \cr
#' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
#' Continuous Report, Psychological Review, 123(4), 425-451. http://dx.doi.org/10.1037/rev0000023
#' @importFrom Rcpp evalCpp
#' @useDynLib CircularDDM
NULL
