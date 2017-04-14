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
// #include <cstdlib>      // std::rand, std::srand
// #include <ctime>        // std::time
// #include <algorithm>    // std::random_shuffle

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

arma::vec besselJ(arma::vec x, double nu=1.0) {
  arma::vec out(x.n_elem);
  for(arma::vec::iterator i=x.begin(); i!=x.end(); ++i)
  {
    int idx  = std::distance(x.begin(), i);
    out[idx] = R::bessel_j(*i, nu);
  }
  return out;
}

//' Generate random deviates from a von Mises distribution
//'
//' This function generates random deviates from a von Mises distribution with
//' parameters, mu and k (ie kappa).
//'
//' A random deviate for a circular normal distribution has the form:\cr
//' \deqn{f(theta; mu, kappa) = 1 / (2*pi*I0(kappa)) * exp(kappa*cos(theta-mu))}
//' theta is within 0 and 2*pi.
//'
//' \code{I0(kappa)} in the normalizing constant is the modified Bessel
//' function of the first kind and order zero.
//'
//' @param n number of observations.
//' @param mu mean direction of the distribution.
//' @param k kappa parameter. A positive numeric value for the concentration
//' parameter of the distribution.
//'
//' @return a column vector
//' @examples
//' n  <- 100
//' mu <- 0
//' k  <- 10
//' vm_de <- rvm(n, mu, k)      ## in degree unit
//' vm_pi <- vm_de %% (2 * pi)  ## in radian unit
//' @export
// [[Rcpp::export]]
arma::vec rvm(int n, double mu, double k) {
  double U, a, b, r, z, f, c;

  a = 1.0 + std::sqrt(1.0 + 4.0 * k * k);
  b = (a - std::sqrt(2.0 * a)) / (2.0 * k);
  r = (1.0 + b * b) / (2.0 * b);

  arma::vec out(n);
  arma::vec::iterator i = out.begin() ;
  do {
    z  = std::cos(M_PI * R::runif(0.0, 1.0));
    f  = (1.0 + r * z) / (r + z);
    c  = k * (r - f);

    U = R::runif(0.0, 1.0);
    if(c * (2.0 - c) > U) {
      *i = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
      if(k == 0) { *i = R::runif(0.0, 2.0*M_PI); }
      i++;
    } else {
      if(std::log(c/U) + 1.0 >= c) {
        *i = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
        if(k == 0) { *i = R::runif(0.0, 2.0*M_PI); }
        i++;
      }
    }
  } while(i < out.end());

  return out;
}

inline double findzero(double n, double x0, int kind, double tol=1e-12,
  int MAXIT=100, double err=1.0) {
  // Tolerance; Maximum number of times to iterate; Initial error
  double a, b, x, n1 = n + 1.0;
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
    err  = (2.0*a*x0*(n*a - b*x0) ) /
      ( 2.0*b*b*x0*x0 - a*b*x0*(4.0 * n1) + (n*n1+x0*x0)*a*a );
    x    = x0 - err;
    x0   = x;
    iter = iter + 1;
  } while ( (std::abs(err) > tol) & (iter < MAXIT) );

  if (iter > (MAXIT - 1)) {
    Rcpp::Rcout << "Failed to converge within tolerance.\n" <<
      "Try a different initial guess";
    x=INFINITY;
  }
  return x;
}

//' Find First k Positive Zeros for the Bessel Functions
//'
//' This function finds the first k positive zeros of the Bessel function
//' J(n, x) or Y(n, x) using Halley's method.
//'
//' @param nu The order of the corresponding Bessel function.
//' @param k an integer for the first k positive zeros.
//' @param kind 0, 1, or 2. A switch selects \link{besselI}, \link{besselJ} or
//' \link{besselY}
//'
//' @return a column vector
//' @references
//' \href{http://au.mathworks.com/matlabcentral/fileexchange/6794-bessel-function-zeros/content/besselzero.m}{besselzero.m}
//' @examples
//' nu <- seq(0, 5, length.out=10)
//' output <- matrix(numeric(5*length(nu)), nrow=5)
//' for (i in 1:length(nu)) { output[,i] <- besselzero(nu[i], 5, 1) }
//' print(output)
//'
//' output <- matrix(numeric(5*length(nu)), nrow=5)
//' for(i in 1:length(nu)) { output[,i] <- besselzero(nu[i], 5, 2) }
//' print(output)
//' @export
// [[Rcpp::export]]
arma::vec besselzero(double nu, int k, int kind) {
  double x0;

  arma::vec x = arma::zeros<arma::vec>(3.0 * (double)k);
  for (int j=1; j<=3*k; j++) {
    // Initial guess of zeros
    x0     = 1.0 + std::sqrt(2.0) + ((double)j - 1.0) * M_PI + nu + std::pow(nu, 0.4);
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
//' equation (23) on p 433 in Smith (2016).
//'
//' This function does not do 'exp' in Smith (2016, p433), so it calculates
//' log-likelihood, instead of just likelihood.
//'
//' @param x a matrix storing a first column of RT and a second column of
//' continuous responses/reports. Each row is a trial.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, sv_sq], which
//' are decision threshold, drift rate along x axis,  drift rate along y axis,
//' non-decision time, and within trial variance.
//'
//' @return a column vector
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' ## An example data, showing the format of input data to logLik_resp
//' x <- cbind (
//' RT=c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R =c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' den  <- logLik_resp(x, pVec=pVec); den
//' @export
// [[Rcpp::export]]
arma::vec logLik_resp(arma::mat x, arma::vec pVec) {
  /* This is the first part of equation (23) without exp in Smith (2016, p433),
   * so it is log-likelihood, instead of likelihood.
   * pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]             */
  arma::vec RT = x.col(0);
  arma::vec R  = x.col(1);
  int n        = R.n_elem;

  // Each row in pMat is a replicate of eg pVec[0] of the number of choices
  arma::mat pMat = arma::repmat(pVec, 1, n);
  arma::vec term0, term1, term2, term3, term4, t0_vec, term0_vec, term3_vec;
  term0 = pVec[0] / pVec[4];
  term1 = trans(pMat.row(1)) % arma::cos(R);
  term2 = trans(pMat.row(2)) % arma::sin(R);
  term3 = (0.5 * pVec[4]) * (pVec[1] * pVec[1] + pVec[2] * pVec[2]);

  term0_vec = arma::repmat(term0, n, 1);
  term3_vec = arma::repmat(term3, n, 1);
  t0_vec    = arma::repmat(pVec.row(3), n, 1);
  term4     = RT - t0_vec; // substracting out non-decision time from RT
  return (term0_vec % (term1 + term2) - (term3_vec % term4));
}

//' Log-Likelihood for Circular First Passage Time
//'
//' Calculate circular log-likelihood of the first passage time, using
//' equation (22) on page 432. Note this function return log density so as to
//' be compatiable with \code{logLik_resp}. Equation (22) in Smith (2016)
//' describes the density function.
//'
//' @param x a matrix storing a first column of RT and a second column of
//' continuous responses/reports. Each row is a trial.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, s].
//' a stands for response threshold, vx is the drift rate along x axis,
//' vy is the drift rate along y axis, t0 is the non-decision time, and s
//' is the within-trial variance.
//' @param k a precision for bessel function. The larger the k is, the larger
//' the memory space is required. Default is 141.
//'
//' @return a column vector
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' x <- cbind(
//' RT=c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R =c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' den  <- logLik_dt(x, pVec=pVec);
//' den
//' @export
// [[Rcpp::export]]
arma::vec logLik_dt(arma::mat x, arma::vec pVec, int k=141) {
  // This is the first part of equation (22) in Smith (2016, page 432)
  // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  int n, idx;
  double dt, a2;         // tmp is sigma^2 / (2*a^2) * t
  arma::vec RT, j0k, j0k2, J1, tmp(k);
  RT      = x.col(0);
  n       = RT.n_elem;
  j0k     = besselzero(0.0, k, 1);
  j0k2    = j0k % j0k;   // element-wise j0k^2
  J1      = besselJ(j0k, 1);
  J1(k-1) = J1(k-1) / 2; // replace the last element
  a2      = pVec[0] *  pVec[0];

  arma::vec out(n);
  for(arma::vec::iterator i=RT.begin(); i!=RT.end(); ++i)
  {
    idx        = std::distance(RT.begin(), i);
    dt         = *i - pVec[3];  // *i extracts individual RT from RT vector.
    tmp.fill(-0.5 * pVec[4] * dt / a2);
    out[idx]   = pVec[4] / a2 * arma::accu(j0k / J1 % arma::exp(tmp % j0k2));
    // When RT and t0 is almost identical, we deem it unlikely.
    if ((*i - pVec[3]) < 0.01) {out[idx] = 1e-10;}
  }
  return arma::log(out);
}

//' The Circular Drift-diffusion Distribution
//'
//' Density function for the circular drift-diffusion model with theta
//' vector equal to \code{pVec}.  \code{dcddm} is the
//' equation (23) on page 433 in Smith (2016).
//'
//' @param x a matrix storing a first column of RT and a second column of
//' continuous responses/reports. Each row is a trial.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, s].
//' a stands for response threshold, vx is the drift rate along x axis,
//' vy is the drift rate along y axis, t0 is the non-decision time, and s
//' is the within-trial variance. The order matters.
//' @param k a precision for calculating the infinite series in \code{dcddm}.
//' The larger the k is, the larger the memory space is required. Default is
//' 141.
//' @return \code{dcddm} gives a likelihood column vector.
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' ## dcddm example
//' x <- cbind(
//' RT= c(1.2595272, 0.8693937, 0.8009044, 1.0018933, 2.3640007, 1.0521304),
//' R = c(1.9217430, 1.7844653, 0.2662521, 2.1569724, 1.7277440, 0.8607271)
//' )
//' pVec <- c(a=2.45, vx=1.5, vy=1.25, t0=.1, s=1)
//' dcddm(x, pVec)
//' @export
// [[Rcpp::export]]
arma::vec dcddm(arma::mat x, arma::vec pVec, int k=141) {
  arma::vec LL_dt   = logLik_dt(x, pVec, k);
  arma::vec LL_resp = logLik_resp(x, pVec);
  return arma::exp(LL_dt + LL_resp);
}

//' The Circular Drift-diffusion Distribution
//'
//' Generate random deviates for the circular drift-diffusion
//' model with a theta vector, namely \code{pVec}. \code{rcddm1} is a canonical
//' form of the random number generator for circular drift-diffusion model.
//' \code{rcddm2} is an extension, allowing one to supply a threshold vector,
//' an angle vector, a starting point matrix [xPos, yPos] and a non-decision
//' time.
//'
//' @param n number of observations.
//' @param pVec a parameter vector with the order [a, vx, vy, t0, sv_sq].
//' a stands for response threshold, vx is the drift rate along x axis,
//' vy is the drift rate along y axis, t0 is the non-decision time, and sv_sq
//' is the within-trial variance. The order matters.
//' @param p a precision for random walk step in \code{rcddm}. Default is 0.01
//' second.
//' @param angle an angle vector, allowing one to supply a vector of random
//' draws from an arbitary, instead of von Mise, distribution.
//' @param threshold a threshold vector.
//' @param sp a starting point matrix.  First column must be the x coordinate
//' of the starting point and the second column must be the y coordinate of
//' the starting point. Each row is a random deviate.
//' @param t0 non-decision time. Must be a scalar.
//'
//' @return a n x 2 matrix with two columns: RTs and angles.
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' ## rcddm1 example
//' pVec <- c(a=2, vx=1.5, vy=1.25, t0=.25, s=1)
//' den  <- rcddm1(1e3, pVec);
//' hist(den[,1], breaks = "fd", xlab="Response Time",  main="Density")
//' hist(den[,2], breaks = "fd", xlab="Response Angle", main="Density")
//'
//' ## rcddm2 example
//' t0 <- .25
//' threshold <- runif(10, 0, 2)
//' angle <- rvm(1e3, 0, 10)
//' startpoint <- cbind(runif(10, 0, 1), runif(10, 0,1))
//' den  <- rcddm2_internal(1e3, threshold, angle, startpoint, t0)
//' head(den)
//'
//' ## If starting points only has one row, this is equilvalent to fixing
//' ## the starting points at constant values. This is because internally the
//' ## function will always use the first row when 'sp' is a row vector.
//' ## Similarly, one can apply the same trick on threshold and angle vectors.
//' startpoint <- cbind(runif(1, 0, 1),  runif(1, 0,1));
//' den  <- rcddm2(1e3, threshold, angle, startpoint, t0)
//'
//' ## If the user enters a R's vector, which is not a matrix,
//' ## rcddm2 wrapper helps the user to convert it to a matrix (ie row vector)
//' ## Note this flexibility requires checking steps.
//' den  <- rcddm2(1e3, threshold, angle, startpoint, t0)
//'
//' @export
// [[Rcpp::export]]
arma::mat rcddm1(int n, arma::vec pVec, double p=0.01) {
  int step;   // pVec [a, vx, vy, t0, s] == [thresh, mu1, mu2, ndt, sigmasq]
  double rPos, xPos, yPos, thPos, theta; // thPos stands for theta position
  arma::vec RT(n), A(n); // R for responses, A for angle

  // page 435 in Smith (2016) equation (29)
  double mu = std::atan2(pVec[2], pVec[1]); // drift direction
  // drift magnitude == std::sqrt(pVec[2]*pVec[2] + pVec[1]*pVec[1])
  // update step size from 1.0 to p
  double k  = p * std::sqrt(pVec[2]*pVec[2] + pVec[1]*pVec[1]) / pVec[4];

  for (int i = 0; i < n; i++) {
    step = 0; // If the starting point is (0, 0), rPos=0. So I did not calculate
    xPos = 0; // initial rPos here.
    yPos = 0;
    do { // should be drift direction
      theta = arma::as_scalar(rvm(1, mu, k)); // get 1 von miss random number
      xPos  += std::cos(theta);               // with parameter mu and kappa
      yPos  += std::sin(theta);
      rPos  = std::sqrt(xPos*xPos + yPos*yPos); // should be drift magnitude
      thPos = std::atan2(yPos, xPos);
      step++;
      if(step > (10.0/p)) {
        // If a trial has taken more than 10 seconds, but yet reached the
        // threshold, we stop this trial.
        Rcpp::Rcout << "Trial " << i << " has taken more than 10 seconds, " <<
          "but yet reached the threshold. Please check the threshold\n";
        break;
      }
    } while (std::abs(rPos) < pVec[0]);

    // dt = 0; // rexp take scale==mean==mu
    // for(int j=0; j<step; j++) { dt = dt + R::rexp(p); }
    // rts[i] = pVec[3] + dt; // gamma a=shape b=scale=1/rate
    RT[i] = pVec[3] + R::rgamma(step, p); // gamma a=shape b=scale=1/rate
    A[i]  = fmod(thPos + 2.0*M_PI, 2*M_PI);
  }
  return arma::join_horiz(RT, A);
}

//' @rdname rcddm1
//' @export
// [[Rcpp::export]]
arma::mat rcddm2_internal(int n, arma::vec threshold, arma::vec angle,
  arma::mat sp, double t0, double p=0.01) {
    int step, idx_a, idx_ang, idx_sp, nang, na, nsp;
    double rPos, xPos, yPos, thPos, theta;
    arma::vec RT(n), A(n);
    nang = angle.n_elem;
    na   = threshold.n_elem;
    nsp  = sp.n_rows;

    // No mu and kappa for von mise, because the user should supply an angle vector
    for (int i = 0; i < n; i++) { // each i is a trial
      step   = 0;                 // '-1' corrects for R-C indexing
      idx_sp = (nsp==1) ? 0 : (std::ceil((double)nsp * R::runif(0.0, 1.0)) - 1);
      idx_a  = (na ==1) ? 0 : (std::ceil((double)na  * R::runif(0.0, 1.0)) - 1);
      arma::rowvec spPos = sp.row(idx_sp);
      xPos = spPos[0]; // sp_x;
      yPos = spPos[1]; // sp_y;
      rPos = std::sqrt(xPos*xPos + yPos*yPos);

      do { // within trial accumulation process
        idx_ang = (nang==1) ? 0 : (std::ceil((double)nang * R::runif(0.0, 1.0)) - 1);
        theta   = angle[idx_ang];
        xPos += std::cos(theta);
        yPos += std::sin(theta);
        rPos  = std::sqrt(xPos*xPos + yPos*yPos);
        thPos = std::atan2(yPos, xPos);
        step++;
        if(step > (10.0/p)) {
          // If a trial has taken more than 10 seconds, but yet reached the
          // threshold, we stop this trial.
          Rcpp::Rcout << "Trial " << i << " has taken more than 10 seconds, " <<
            "but yet reached the threshold. Please check the threshold\n";
          break;
        }
      } while (std::abs(rPos) < threshold[idx_a]);

      RT[i] = t0 + R::rgamma(step, p); // gamma a=shape b=scale=1/rate
      A[i]  = fmod(thPos + 2.0*M_PI, 2.0*M_PI);
    }
    return arma::join_horiz(RT, A);
}

struct generateIntVec {
    int x ;
    generateIntVec() {x=0;}
    int operator()() {return x++;}
} ; // integer generator:


/*
std::vector<int> genInt (int n) {
    // generate an integer vector based on C index from 0 to n-1
    std::vector<int> out(n); // a empty vector with n integer mem space.
    generateIntVec generator; //
    std::generate(out.begin(), out.end(), generator);
    return out;
}

std::vector<int> genInt2 (int n) {
    // generate an integer vector based on C index from 0 to n-1
    std::vector<int> out(n); // a empty vector with n integer mem space.
    generateIntVec generator; //
    std::generate(out.begin(), out.end(), generator);

    std::srand (unsigned (std::time(0)));
    std::random_shuffle(out.begin(), out.end());
    return out;
}

std::vector<int> genInt3 (int n) {
    std::vector<int> out(n); // a empty vector with n integer mem space.
    for(std::vector<int>::iterator i=out.begin(); i!=out.end(); ++i)
    {
        *i = std::ceil((double)n * R::runif(0.0, 1.0)); // idx == 0
    }
    return out;
}

*/

// #include <RcppArmadilloExtensions/sample.h>
// arma::vec test_sample(int a) {
//   arma::vec x   = arma::randn<arma::vec>(5);
//   int size = 1;
//   bool replace = 1;
//   arma::vec prob = arma::randu<arma::vec>(5);
//   arma::vec out = Rcpp::RcppArmadillo::sample(x, size, replace, prob);
//   return out;
// }

