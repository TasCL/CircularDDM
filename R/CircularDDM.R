#' Circular Drift-diffusion Model
#'
#' \code{CircularDDM} provides density functions for circular
#' decision-diffusion/decision-diffusion model for continuous report. It
#' implements equations and convenient functions
#' for calculating likelihood and generating random deviates. For circular
#' DDM, random deviates are continuous response times and response angles.
#'
#' @keywords CircularDDM
#' @name CircularDDM
#' @docType package
#' @author  Yi-Shin Lin <yishin.lin@utas.edu.au> \cr
#' Andrew Heathcote <andrew.heathcote@utas.edu.au> \cr
#' Peter Kvam <kvam.peter@gmail.com> \cr
#' @references
#' \enumerate{
#'       \item Smith, P. L. (2016). Diffusion Theory of Decision Making in
#'       Continuous Report, \emph{Psychological Review}, 123 (4), 425--451.
#'       \url{http://dx.doi.org/10.1037/rev0000023}.
#'       \item Hamana, Y., & Matsumoto, H. (2013). The probability distributions
#'       of the first hitting times of Bessel processes, \emph{Transactions}
#'       \emph{of the American Mathematical Society}, 365, 5237--5257.
#'       \item Agostinelli, C. and Lund, U (2017). R package 'circular':
#'       Circular Statistics (version 0.4-93). \url{https://r-forge.r-project.org/projects/circular/}.
#'       \item Von Winckel, G. besselzero.m
#'       \url{http://au.mathworks.com/matlabcentral/fileexchange/6794-bessel-function-zeros/content/besselzero.m}
#' }
#' @importFrom Rcpp evalCpp
#' @useDynLib CircularDDM
NULL

#' The Circular Decision-diffusion Model
#'
#' \code{dcddm0} is the zero-drift first passage of the Bessel process through
#' the boundary \code{a}. It is equation (22) on page 432.
#' \code{dcddm} is equation (23) on page 433 (Smith, 2016), the nonzero-drift
#' model. \code{rcddm} gives random deviates for the circular DDM.
#'
#' @param x a vector storing RTs or a matrix storing RTs and responses.
#' Each row represents a trial.
#' @param n number of observations. This must be a scalar.
#' @param a decision boundary. This must be a scalar.
#' @param mu1 x axis drift rate. This must be a scalar.
#' @param mu2 y axis drift rate. This must be a scalar.
#' @param t0 non-decision time. This must be a scalar.
#' @param sigma square root of the diffusion coefficient. This is the
#' within-trial, moment-to-moment, standard deviation. This must be a scalar.
#' @param k maximum number of eigenvalues. This is the parameter
#' controling the precision of the bessel function. The larger the k is, the
#' higher the precision (the slower the computation). Note we use a default
#' value, 141. Smith (2016) used 50.
#' @param ntheta number of response angles. \code{rcddm} draws a range of
#' response angles between pi and -pi. \code{ntheta} the numbers of response
#' angles between pi and -pi to draw.
#' @param upper,lower upper and lower bounds of RT. Similar to \code{ntheta}
#'
#' @return \code{dcddm0} and \code{dcccm} give the density. \code{rcddm}
#' generates random deviates (RT and responses).
#' @examples
#' ## Example 1
#' a <- 1
#' sigma <- 1
#' kmax <- 50
#' h <- 1e-4 # second
#' tmax <- 2 # second
#' t0 <- 0
#'
#' x <- cbind(RT = seq(0, tmax, by = h), R  = seq(-pi, pi, length.out = tmax/h + 1))
#' den1 <- dcddm0(x[,1], a, t0, sigma, kmax)
#' plot(x[,1], den1, type = "l", lwd = 2, xlab = "RT (s)", ylab ="Density")
#'
#' ## Example 2
#' x <- cbind(
#' RT = c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
#' R  = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
#' )
#'
#' den1 <- dcddm(x, a = 2.45, mu1 = 1.5, mu2 = 1.25, t0 = .1, s = 1)
#'
#' ## Example 3
#' a <- 1
#' t0 <- .15
#' sigma <- sqrt(2)
#' kmax <- 121
#' driftMag <- .5
#' driftDir <- -pi/3
#' mu1 <- driftMag * cos(driftDir)
#' mu2 <- driftMag * sin(driftDir)
#' n <- 15000
#' dat1 <- rcddm(n, a, mu1, mu2, t0, sigma, kmax)
#' @rdname dcddm
#' @export
dcddm <- function(x, a = 1, mu1 = 0, mu2 = 0, t0 = 0, sigma = 1, k = 141) {
  if(!is.matrix(x) | ncol(x) != 2) stop("x must be a n x 2 matrix")
  return(dcircle(x, a, mu1, mu2, t0, sigma, k)[,1])
}

#' @rdname dcddm
#' @export
dcddm0 <- function(x, a = 1, t0 = 0, sigma = 1, k = 141) {
  if(!is.vector(x)) stop("x must be a vector")
  return(dhamana(x, a, t0, sigma, k)[,1])
}

#' @rdname dcddm
#' @export
rcddm <- function(n, a = 1, mu1 = 0, mu2 = 0, t0 = 0, sigma = 1, k = 141,
  lower = 0.01, upper = 2.5, ntheta = 41) {
  if (lower < 0 | upper < 0) stop("lower and upper must be positive.")
  if (ntheta < 0) stop("ntheta must be a positive integer")
  out <- data.frame(rcircle(n, a, mu1, mu2, t0, sigma, k, lower, upper ,ntheta))
  names(out) <- c("RT", "R")
  return(out)
}

#' Find First k Zeros of a Bessel Function
#'
#' The function modifies from Von Winckel's \code{besselzero.m}, using Halley's
#' method to find the first k zeros of a Bessel function.
#' The Bessel function is either the first kind, \code{R::bessel_j} or the
#' second kind, \code{R::bessel_y}.
#'
#' @param nu The order of the corresponding Bessel function.
#' @param k an integer for the first k positive zeros.
#' @param kind 0, 1, or 2. A switch selects \link{besselI}, \link{besselJ} or
#' \link{besselY}
#'
#' @return a vector
#' @export
#' @examples
#' nu <- seq(0, 5, length.out=10)
#' output <- matrix(numeric(5*length(nu)), nrow=5)
#' for (i in 1:length(nu)) { output[,i] <- bessel_zero(nu[i], 5, 1) }
#' print(output)
#'
#' output <- matrix(numeric(5*length(nu)), nrow=5)
#' for(i in 1:length(nu)) { output[,i] <- bessel_zero(nu[i], 5, 2) }
#' print(output)
bessel_zero <- function(nu, k, kind) {
  out <- besselzero(nu, k, kind)
  return(out[,1])
}

#' Convert Dirft Rate between Magnitude-Direction and X-Y coordinate
#'
#' The function convert dirft rate from magnitude-direction to X-Y coordinate or
#' vice versa.
#'
#' @param mean_v a \code{c(magnitude, direction)} or \code{c(mu_x mu_y)} vector
#' @param mag2xy a boolean switch indicating converting from \code{c(magnitude,
#' direction)} (TRUE) to \code{c(mu_x, mu_y)} or the other way (FALSE).
#'
#' @return a vector
#' @export
#' @examples
#' conv_drift(c(.5, -pi/3), TRUE)
#' conv_drift(c(.25, -0.4330127), FALSE)
#'
conv_drift <- function(mean_v, mag2xy) {
  out <- convDrift(mean_v, mag2xy)
  return(out[,1])
}
