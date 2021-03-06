% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CircularDDM.R
\name{dcddm}
\alias{dcddm}
\alias{dcddm0}
\alias{rcddm}
\title{The Density Function of the Circular Decision-diffusion Model}
\usage{
dcddm(x, a = 1, mu1 = 0, mu2 = 0, t0 = 0, sigma = 1, k = 141)

dcddm0(x, a = 1, t0 = 0, sigma = 1, k = 141)

rcddm(n, a = 1, mu1 = 0, mu2 = 0, t0 = 0, sigma = 1, k = 141,
  lower = 0.01, upper = 2.5, ntheta = 41)
}
\arguments{
\item{x}{a vector storing RTs or a \emph{matrix} storing RTs and responses.
Each row is a response.}

\item{a}{decision boundary. This must be a scalar.}

\item{mu1}{x axis drift rate. This must be a scalar.}

\item{mu2}{y axis drift rate. This must must be a scalar.}

\item{t0}{non-decision time. This must be a scalar.}

\item{sigma}{square root of the diffusion coefficient. This is the
within-trial, moment-to-moment, standard deviation. This must be a scalar.}

\item{k}{maximum number of eigenvalues. The tuning parameter controls
the precision of the bessel function. The larger the k is, the higher
the precision. Default is 141. Note Smith (2016) used 50.}

\item{n}{number of observations. This must be a scalar.}

\item{upper, lower}{upper and lower bounds of RT for; These three options are
for \code{rcddm}.}

\item{ntheta}{number of response angles.}
}
\value{
\code{dcddm0} and \code{dcccm} give the density. \code{rcddm}
generates random deviates (RT and responses).
}
\description{
\code{dcddm0} is the zero-drift first passage of the Bessel process through
the boundary \code{a}. It is equation (22) on page 432.
\code{dcddm} is equation (23) on page 433 (Smith, 2016).
}
\examples{
## Example 1
a <- 1
sigma <- 1
kmax <- 50
h <- 1e-4 # second
tmax <- 2 # second
t0 <- 0

x <- cbind(RT = seq(0, tmax, by = h), R  = seq(-pi, pi, length.out = tmax/h + 1))
den1 <- dcddm0(x[,1], a, t0, sigma, kmax)
plot(x[,1], den1, type = "l", lwd = 2, xlab = "RT (s)", ylab ="Density")

## Example 2
x <- cbind(
RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
)

den1 <- dcddm(x, a = 2.45, mu1 = 1.5, mu2 = 1.25, t0 = .1, s = 1)

## Example 3
a <- 1
t0 <- .15
sigma <- sqrt(2)
kmax <- 121
driftMag <- .5
driftDir <- -pi/3
mu1 <- driftMag * cos(driftDir)
mu2 <- driftMag * sin(driftDir)
n <- 15000
dat1 <- rcddm(n, a, mu1, mu2, t0, sigma, kmax)
}
