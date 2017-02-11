#include <RcppArmadillo.h>
arma::vec rvm(int n, double mu, double k);
arma::vec besselzero(double nu, int k, int kind);
arma::vec logLik_resp(arma::mat x, arma::vec pVec);
arma::vec logLik_dt(arma::mat x, arma::vec pVec, int k);
arma::vec dddm(arma::mat x, arma::vec pVec, int k);
arma::mat rddm(int n, arma::vec pVec, double p);
