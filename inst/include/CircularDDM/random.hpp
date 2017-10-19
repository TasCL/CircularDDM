#include <RcppArmadillo.h>
arma::vec mod(arma::vec dividend, double divisor);
arma::vec rvm(int n, double mu, double kappa);
arma::mat rcddm1(int n, double a, double vx, double vy, double t0, double s,
  double stepSize, int maxStep);
arma::mat rcddm2_internal(int n, arma::vec threshold, arma::vec angle,
  arma::mat sp, double t0, double p);


