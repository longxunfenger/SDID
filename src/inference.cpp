#include <RcppArmadillo.h>
#include "unit_weights.h"
#include "time_weights.h"
#include "att.h"
// [[Rcpp::depends("RcppArmadillo")]]

//' Inference module. Obtains the empirical distribution of ATT estimate
//' @param X - N*T x k data matrix
//' @param Y - N*T vector of outcomes
//' @param Y00 - T0 x N0 matrix of control unit outcomes before the treatment
//' @param Y10 - T0 x N - N0 matrix of control unit outcomes after the treatment
//' @param Y01 - T-T0 x N - N0 matrix of control unit outcomes after the treatment
//' @param dzeta - Ridge penalty parameter
//' @param rho - ADMM parameter
//' @param tol - ADMM convergence tolerance
//' @return List containing vector atts of estimates
//' @export
// [[Rcpp::export]]
Rcpp::List inference(arma::mat &X, arma::colvec &Y, arma::mat &Y00, arma::mat &Y01, arma::mat &Y10, double dzeta = 10,
                     double rho = 1, double tol = 1e-12)
{
  //first we construct our matrices 
  int N0 = Y00.n_cols;
  int N1 = Y10.n_cols;
  int T0 = Y00.n_rows;
  int T1 = Y01.n_rows;
  int k = 1;//X.n_cols;
  
  
  arma::mat w_unit(N0 + N1, 1);
  arma::mat w_time(T0 + T1, 1);
  
  
  //get the weights
  w_unit = unit_weights(Y00, Y10, N1, dzeta, rho, tol, "SDID");
  w_time = time_weights(Y00, Y01, dzeta, rho, tol, "SDID");
  
  
  arma::mat att_ests(1, N0);
  att_ests.zeros();
  

  //now we leave one unit out and estimate att
  for (int i = 0; i < N0; i++)
  {
    att_ests.col(i) = att(X, Y, i, w_unit, w_time);
  }


  return Rcpp::List::create(Rcpp::Named("atts") = att_ests);
}