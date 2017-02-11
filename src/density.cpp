//   Copyright (C) <2017>  <Yi-Shin Lin>
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; version 2
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License along
//   with this program; if not, write to the Free Software Foundation, Inc.,
//   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <CircularDDM.hpp>

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

arma::vec besselJ(arma::vec x, double nu=1) {
  arma::vec out(x.n_elem);
  for(arma::vec::iterator i=x.begin(); i!=x.end(); ++i)
  {
    int idx  = std::distance(x.begin(), i);
    out[idx] = R::bessel_j(*i, nu);
  }
  return out;
}

//' Von Mises Distribution: Random Number Generator
//'
//' Generate random numbers for von Mises circular distribution.
//'
//' @param n number of observations.
//' @param mu mean direction of the distribution.
//' @param k non-negative numeric value for the concentration parameter of the
//' distribution
//'
//' @return a vector
//' @examples
//' n  <- 100
//' mu <- 0
//' k  <- 10
//' vm3_de <- rvm(n, mu, k)       ## in degree unit
//' vm3_pi <- vm3_de %% (2 * pi)  ## in radian unit
//' @export
// [[Rcpp::export]]
arma::vec rvm(int n, double mu, double k) {
  double U, a, b, r, z, f, c;

  a = 1 + std::sqrt(1+4 * k * k);
  b = (a - std::sqrt(2*a))/(2* k);
  r = (1 + b*b)/(2*b);

  arma::vec out(n);
  arma::vec::iterator i = out.begin() ;
  do {
    z  = std::cos(M_PI * R::runif(0,1));
    f  = (1. + r * z)/(r + z);
    c  = k * (r - f);

    U = R::runif(0, 1);
    if(c * (2 - c) > U) {
      *i = (R::runif(0,1) > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
      if(k == 0) { *i = R::runif(0, 2*M_PI); }
      i++;
    } else {
      if(std::log(c/U) + 1 >= c) {
        *i = (R::runif(0,1) > .50) ? std::acos(f) + mu : -std::acos(f) + mu;
        if(k == 0) { *i = R::runif(0, 2*M_PI); }
        i++;
      }
    }
  } while(i < out.end());

  return out;
}

inline double findzero(double n, double x0, int kind, double tol=1e-12,
  int MAXIT=100, double err=1) {
  // Tolerance; Maximum number of times to iterate; Initial error

  double a, b, x, n1=n+1;
  int iter=0;

  do {
    switch (kind) {
    case 0 :
      a = R::bessel_i(x0, n,  0);
      b = R::bessel_i(x0, n1, 0);
      break;
    case 1:
      a = R::bessel_j(x0, n);
      b = R::bessel_j(x0, n1);
      break;
    case 2:
      a = R::bessel_y(x0, n);
      b = R::bessel_y(x0, n1);
      break;
    default:
      a = 0;
      b = 0;
    }
    err  = (2*a*x0*(n*a - b*x0) ) /
      ( 2*b*b*x0*x0 - a*b*x0*(4*n1) + (n*n1+x0*x0)*a*a );
    x    = x0 - err;
    x0   = x;
    iter = iter + 1;
  } while ( (std::abs(err) > tol) & (iter < MAXIT) );

  if (iter > (MAXIT - 1)) {
    std::cout << "Failed to converge to within tolerance.\n" <<
      "Try a different initial guess";
    x=INFINITY ;
  }

  return x;

}

//' Find First k Positive Zeros for the Bessel Functions
//'
//' Find first k positive zeros of the Bessel function J(n,x) or Y(n,x)
//' using Halley's method.
//'
//' @param nu The order of the corresponding Bessel function.
//' @param k an integer for first k positive zeros.
//' @param kind 0, 1, or 2. A switch selects \link{besselI}, \link{besselJ} or
//' \link{besselY}
//'
//' @return a vector
//' @references
//' \href{http://au.mathworks.com/matlabcentral/fileexchange/6794-bessel-function-zeros/content/besselzero.m}{besselzero.m}
//' @examples
//' nu <- seq(0, 5, length.out=10)
//' output <- matrix(numeric(5*length(nu)), nrow=5)
//'   for(i in 1:length(nu)) {
//'     output[,i] <- besselzero(nu[i], 5, 1)
//'   }
//' output
//'
//' output <- matrix(numeric(5*length(nu)), nrow=5)
//' for(i in 1:length(nu)) {
//'     output[,i] <- besselzero(nu[i], 5, 2)
//' }
//' output
//' @export
// [[Rcpp::export]]
arma::vec besselzero(double nu, int k, int kind) {
  double x0;

  arma::vec x = arma::zeros<arma::vec>(3*k);
  for (int j=1; j<=3*k; j++) {
    // Initial guess of zeros
    x0     = 1 + std::sqrt(2) + (j-1) * M_PI + nu + std::pow(nu, 0.4);
    x(j-1) = findzero(nu, x0, kind);     // Halley's method
    if (x(j-1) == INFINITY) {Rcpp::stop("Bad guess.");}
  }
  if(!x.is_sorted()) { x = sort(x); };
  arma::vec onevec = arma::ones<arma::vec>(1);
  arma::vec dx     = arma::join_vert(onevec, arma::abs(arma::diff(x)));
  arma::vec out    = x(arma::find(dx > 1e-8));
  if( out.has_nan() ) {Rcpp::stop("NA found.");}
  return out.rows(0, k-1);
}

//' Log-Likelihood for Continuous Reports
//'
//' Calculate log-likelihood of the continuous reports, using
//' part part in equation (23) on p 433.
//'
//' @param x a matrix storing a first column as RT and a second column of
//' continuous responses/reports/outcomes. Each row is a trial.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, s],
//' or [thresh, mu1, mu2, ndt, sigmasq], using alternative names.
//' @param log a switch to return log-likelihood or just likelihood. Default
//' is to return log-likelihood
//'
//' @return a vector
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' x <- cbind(
//' RT=c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R =c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' den  <- logLik_resp(x, pVec=pVec); den
//' @export
// [[Rcpp::export]]
arma::vec logLik_resp(arma::mat x, arma::vec pVec) {
  // This is the first part of equation (23) with exp in Smith (2016)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  arma::vec rts = x.col(0);
  arma::vec choices = x.col(1);
  int n = choices.n_elem;

  // Each row of pMat is a replicates of eg pVec[0] of the number of choices
  arma::mat pMat = arma::repmat(pVec, 1, n);
  arma::vec term0, term1, term2, term3, term4, t0_vec, term0_vec, term3_vec;
  term0 = pVec[0] / pVec[4];
  term1 = trans(pMat.row(1)) % arma::cos(choices);
  term2 = trans(pMat.row(2)) % arma::sin(choices);
  term3 = (0.5 * pVec[4]) * ( std::pow(pVec[1], 2) + std::pow(pVec[2], 2) );

  term0_vec = arma::repmat(term0, n, 1);
  term3_vec = arma::repmat(term3, n, 1);
  t0_vec    = arma::repmat(pVec.row(3), n, 1);
  term4     = rts - t0_vec;
  arma::vec out = term0_vec % (term1 + term2) - (term3_vec % term4) ;
  return out;
}

//' Log-Likelihood for Circular First Passage Time
//'
//' Calculate circular log-likelihood of the first passage time, using
//' equation (22) on p 432.
//'
//' @param x a matrix storing a first column as RT and a second column of
//' continuous responses/reports/outcomes. Each row is a trial.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, s],
//' or [thresh, mu1, mu2, ndt, sigmasq].
//' @param k a precision for bessel function. The larger the k is, the larger
//' the memory space is required. Default is 141.
//' @param log a switch to return log-likelihood or just likelihood. Default
//' is to return log-likelihood
//'
//' @return a vector
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' x <- cbind(
//' RT=c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R =c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' den  <- logLik_dt(x, pVec=pVec); den
//' @export
// [[Rcpp::export]]
arma::vec logLik_dt(arma::mat x, arma::vec pVec, int k=141) {
  // This is the first part of equation (23) with exp in Smith (2016)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  int n, idx; double dt, tmp;
  arma::vec rts, j0k, j0k2, J1, scalar, scalar_vec;
  rts     = x.col(0);
  n       = rts.n_elem;
  j0k     = besselzero(0, k, 1);
  j0k2    = arma::pow(j0k, 2); // squared j0k
  J1      = besselJ(j0k, 1);
  J1(k-1) = J1(k-1) / 2; // replace the last element

  arma::vec out(n);
  for(arma::vec::iterator i=rts.begin(); i!=rts.end(); ++i)
  {
    idx        = std::distance(rts.begin(), i);
    dt         = *i - pVec[3];
    tmp        = -0.5 * pVec[4] * dt / std::pow(pVec[0] , 2);
    scalar     = tmp;
    scalar_vec = arma::repmat(scalar, k, 1);
    out[idx]   = pVec[4] / std::pow(pVec[0], 2) *
      arma::accu(j0k / J1 % arma::exp(scalar_vec % j0k2));
    // When RT and t0 is almost identical, we deem it unlikely.
    if ((*i - pVec[3]) < 0.01) {out[idx] = 1e-10;} //
  }

  return arma::log(out);

}

//' The Circular Drift-diffusion Distribution
//'
//' Density function and random generation for the circular drift-diffusion
//' model with theta vector equal to \code{pVec}.  \code{dddm} is the
//' equation (23) in p 433 in Smith (2016).
//'
//' @param x a matrix storing a first column as RT and a second column of
//' continuous responses/reports/outcomes. Each row is a trial.
//' @param n number of observations.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, s],
//' or [thresh, mu1, mu2, ndt, sigmasq]. The order matters.
//' @param k a precision for calculate the infinite series in \code{dddm}. The
//' larger the k is, the larger the memory space is required. Default is 141.
//' @param p a precision for random walk step in \code{rddm}. Default is 0.15
//' second
//' @return \code{dddm} gives a log-likelihood vector. \code{rddm} generates
//' random deviates, returning a n x 3 matrix with the columns: RTs, choices
//' and then angles.
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' ## dddm example
//' x <- cbind(
//' RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' dddm(x, pVec)
//'
//' ## rddm example
//' pVec <- c(a=2, vx=1.5, vy=1.25, t0=.25, s=1)
//' den  <- rddm(1e3, pVec);
//' hist(den[,1], breaks = "fd", xlab="Response Time",  main="Density")
//' hist(den[,3], breaks = "fd", xlab="Response Angle", main="Density")
//' @export
// [[Rcpp::export]]
arma::vec dddm(arma::mat x, arma::vec pVec, int k=141) {
  arma::vec LL_dt   = logLik_dt(x, pVec, k);
  arma::vec LL_resp = logLik_resp(x, pVec);
  return LL_dt + LL_resp;
}

//' @rdname dddm
//' @export
// [[Rcpp::export]]
arma::mat rddm(int n, arma::vec pVec, double p=.15) {
  int step;   // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  double rPos, xPos, yPos, thPos, theta; // thPos stands for theta position
  arma::vec RT(n), R(n), A(n); // R for responses, A for angle
  for (int i = 0; i < n; i++) {
    step = 0; rPos = 0; xPos = 0; yPos = 0;
    do {
      theta = arma::as_scalar(rvm(1, 1, 1));
      xPos  = xPos + std::cos(theta);
      yPos  = yPos + std::sin(theta);
      rPos  = std::sqrt(std::pow(xPos, 2) + std::pow(yPos, 2));
      thPos = std::atan2(yPos, xPos);
      step++;
    } while (std::abs(rPos) < pVec[0]);

    // dt = 0; // rexp take scale==mean==mu
    // for(int j=0; j<step; j++) { dt = dt + R::rexp(p); }
    // rts[i] = pVec[3] + dt; // gamma a=shape b=scale=1/rate
    RT[i] = pVec[3] + R::rgamma(step, p); // gamma a=shape b=scale=1/rate
    R[i]  = thPos/2;
    A[i]  = ((0.5 - R[i]) > 0.5*M_PI) ? M_PI - (0.5 - R[i]) : 0.5 - R[i];
  }

  arma::mat tmp = arma::join_horiz(RT, R);
  arma::mat out = arma::join_horiz(tmp, A);
  return out;
}


