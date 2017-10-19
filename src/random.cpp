#include <CircularDDM.hpp>

arma::vec fmod(arma::vec dividend, double divisor) {
  // a float point modulus operator taking armadillo vector and a double
  return dividend - arma::floor(dividend / divisor)*divisor;
}

//' Generate random deviates from a von Mises distribution
//'
//' This function generates random numbers in radian unit from a von Mises
//' distribution using the location (ie mean) parameter, mu and the
//' concentration (ie precision) parameter kappa.
//'
//' A random number for a circular normal distribution has the form:\cr
//' \deqn{f(theta; mu, kappa) = 1 / (2*pi*I0(kappa)) * exp(kappa*cos(theta-mu))}
//' theta is between 0 and 2*pi.
//'
//' \code{I0(kappa)} in the normalizing constant is the modified Bessel
//' function of the first kind and order zero.
//'
//' @param n number of observations. Must be a scalar.
//' @param mu mean direction of the distribution. Must be a scalar.
//' @param kappa concentration parameter. That is, the kappa. A positive value
//' for the concentration parameter of the distribution. Must be a scalar.
//'
//' @return a column vector
//' @examples
//' n  <- 1e2
//' mu <- 0
//' k  <- 10
//'
//' \dontrun{
//' vm1 <- circular:::RvonmisesRad(n, mu, k)
//' vm2 <- rvm(n, mu, k)
//' vm3 <- circular:::conversion.circular(circular:::circular(vm1))
//' vm4 <- circular:::conversion.circular(circular:::circular(vm2))
//' plot(vm3)
//' plot(vm4)
//' }
//' @export
// [[Rcpp::export]]
arma::vec rvm(int n, double mu, double kappa) {
  /* This function is modified from Ulric Lund, Claudio Agostinelli, et al's
   * rvm C code and RvonmisesRad R code in 'circular' package, Version 0.4-91.
   * Ported to RcppArmadillo by Yi-Shin Lin (2017)
   */

  double U, a, b, r, z, f, c, tmp, twopi = 2*M_PI;
  arma::vec out(n);
  arma::vec::iterator i = out.begin() ;

  // If kappa is small, sample angles from an uniform distribution [0 2*pi]
  if (kappa < 1e-10 ) {
    do { *i =  R::runif(0.0, twopi); i++; } while (i < out.end());
  } else {
    a = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    b = (a - std::sqrt(2.0 * a)) / (2.0 * kappa);
    r = (1.0 + b * b) / (2.0 * b);

    do {
      z  = std::cos(M_PI * R::runif(0.0, 1.0));
      f  = (1.0 + r*z) / (r + z);
      c  = kappa*(r - f);
      U = R::runif(0.0, 1.0);

      if (c * (2.0 - c) > U) {
        tmp = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
        *i = tmp - std::floor(tmp/twopi)*twopi;
        i++;
      } else {
        if (std::log(c/U) + 1.0 >= c) {
          tmp = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
          *i = tmp - std::floor(tmp/twopi)*twopi;
          i++;
        }
      }
    } while(i < out.end());
  }

  return out;
}


//' The Circular Decision-diffusion Distribution
//'
//' Generate random deviates from a circular decision-diffusion
//' model, using parameters, a, vx, vy, t0, and s. \code{rcddm1} is
//' canonical form of the random number generator for the circular
//' decision-diffusion model. \code{rcddm2} is an extension, allowing one to
//' supply a threshold vector, an angle vector, a starting point matrix
//' [xPos, yPos] and a non-decision time. The angle vector carries the
//' information of the drift magnitude and the drift direction.
//'
//' @param n numbers of observation.
//' @param a decision threshold. Must be a scalar.
//' @param mu1 x axis drift rate. This must be a scalar.
//' @param mu2 y axis drift rate. This must must be a scalar.
//' @param t0 non-decision time. This must be a scalar.
//' @param sigma square root of the diffusion coefficient. This is the
//' within-trial, moment-to-moment, standard deviation. This must be a scalar.
//' @param stepSize the size of a random walk step.
//' @param stepTime the unit time for a random walk step.
//' @param angle an angle vector, allowing one to supply a vector of randomly
//' drawn numbers from any distributions.
//' @param threshold a threshold vector.
//' @param sp a starting point matrix.  First column must be the x coordinate
//' of the starting point and the second column must be the y coordinate of
//' the starting point. Each row is a random number.
//'
//' @return a n x 2 matrix with two columns: RTs and angles.
//' @references Smith, P. L. (2016). Diffusion Theory of Decision Making in
//' Continuous Report, Psychological Review, 123 (4), 425--451.
//' @examples
//' ## rcddm1 example
//' den  <- rcddm1(1024, a = 2, mu1 = 1.5, mu2 = 1.25, t0 = 0.25, sigma = 1);
//' par(mfrow = c(1, 2))
//' hist(den[,1], breaks = "fd", xlab="RT (s)",  main="")
//' hist(den[,2], breaks = "fd", xlab="Angle (Rad)", main="")
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
//' den <- rcddm2_internal(1e3, threshold, angle, startpoint, t0)
//'
//' ## If the user enters a R's vector, which is not a matrix,
//' ## rcddm2 wrapper helps the user to convert it to a matrix (ie row vector)
//' ## Note this flexibility requires checking steps.
//' den <- rcddm2_internal(1e3, threshold, angle, startpoint, t0)
//'
//' @export
// [[Rcpp::export]]
arma::mat rcddm1(unsigned int n, double a, double mu1, double mu2, double t0,
  double sigma, double stepSize = 1, double stepTime = 1e-4) {
  // thPos = theta position; thStep = theta step; R = response angles,
  unsigned int nstep;
  double rPos, xPos, yPos, thStep;
  arma::vec RT(n), R(n);
  // the displacement of a Brownian motion is a function of the square root of
  // the time
  // double stepTime = stepSize * stepSize;

  // page 435 in Smith (2016) equation (29)
  // drift magnitude == boundary * std::sqrt(pVec[2]*pVec[2] + pVec[1]*pVec[1])
  // step size * the magnitude of Euclidean norm / sigmasq
  double mu     = std::atan2(mu2, mu1); // drift direction
  double kappa  = stepSize * std::sqrt(mu1*mu1 + mu2*mu2) / sigma * sigma;

  for (size_t i = 0; i < n; i++) {
    nstep = 0; // Assume starting from (0, 0),
    xPos  = 0;
    yPos  = 0;
    do {
      thStep = arma::as_scalar(rvm(1, mu, kappa)); // get 1 von miss random number in radian
      xPos  += stepSize*std::cos(thStep);   // with parameter mu and kappa
      yPos  += stepSize*std::sin(thStep);
      rPos   = std::sqrt(xPos*xPos + yPos*yPos); // drift magnitude
      nstep++;
    } while (rPos < a);

    RT[i] = t0 + R::rgamma(nstep, stepTime); // gamma a=shape b=scale=1/rate
    R[i]  = std::atan2(yPos, xPos);
    // R[i]  = thPos - std::floor(thPos/ 2.0*M_PI)*2.0*M_PI;
    // R[i]  = std::fmod(thPos + 2.0*M_PI, 2*M_PI);
    // R[i]  = thPos;
  }
  return arma::join_horiz(RT, R);
}

//' @rdname rcddm1
//' @export
// [[Rcpp::export]]
arma::mat rcddm2_internal(unsigned int n, arma::vec threshold, arma::vec angle,
  arma::mat sp, double t0, double stepSize = .1, double stepTime = .005) {

  int step, idx_a, idx_ang, idx_sp, nang, na, nsp;
  double rPos, xPos, yPos, thPos, theta;
  arma::vec RT(n), A(n);
  nang = angle.n_elem;
  na   = threshold.n_elem;
  nsp  = sp.n_rows;

  // No mu and kappa for von mise, because the user should supply an angle vector
  for (size_t i = 0; i < n; i++) { // each i is a trial
    step   = 0;                    // '-1' corrects for R-C indexing
    idx_sp = (nsp==1) ? 0 : (std::ceil((double)nsp * R::runif(0.0, 1.0)) - 1);
    idx_a  = (na ==1) ? 0 : (std::ceil((double)na  * R::runif(0.0, 1.0)) - 1);
    arma::rowvec spPos = sp.row(idx_sp);
    xPos = spPos[0]; // sp_x;
    yPos = spPos[1]; // sp_y;
    // Check initial position does not suppass threshold
    rPos = std::sqrt(xPos*xPos + yPos*yPos);

    do { // within trial accumulation process
      idx_ang = (nang==1) ? 0 : (std::ceil((double)nang * R::runif(0.0, 1.0)) - 1);
      theta   = angle[idx_ang];
      xPos += stepSize*std::cos(theta);
      yPos += stepSize*std::sin(theta);
      rPos = std::sqrt(xPos*xPos + yPos*yPos);
      step++;
    } while (rPos < threshold[idx_a]);

    thPos = std::atan2(yPos, xPos);
    RT[i] = t0 + R::rgamma(step, stepTime); // gamma a=shape b=scale=1/rate
    A[i]  = std::fmod(thPos + 2.0*M_PI, 2.0*M_PI);

  }
  return arma::join_horiz(RT, A);
}


