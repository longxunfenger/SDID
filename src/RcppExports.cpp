// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// admm_solver
arma::mat admm_solver(arma::mat& Dmat, arma::colvec& dvec, double rho, int max_iter, double tol);
RcppExport SEXP _sdid_admm_solver(SEXP DmatSEXP, SEXP dvecSEXP, SEXP rhoSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Dmat(DmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type dvec(dvecSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(admm_solver(Dmat, dvec, rho, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// admm_time_solver
arma::mat admm_time_solver(arma::mat& Dmat, arma::colvec& dvec, double rho, int max_iter, double tol);
RcppExport SEXP _sdid_admm_time_solver(SEXP DmatSEXP, SEXP dvecSEXP, SEXP rhoSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Dmat(DmatSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type dvec(dvecSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(admm_time_solver(Dmat, dvec, rho, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// att
arma::mat att(arma::mat& X, arma::colvec& y, int unit_exclude, arma::mat& unit_weights, arma::mat& time_weights);
RcppExport SEXP _sdid_att(SEXP XSEXP, SEXP ySEXP, SEXP unit_excludeSEXP, SEXP unit_weightsSEXP, SEXP time_weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type unit_exclude(unit_excludeSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type unit_weights(unit_weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type time_weights(time_weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(att(X, y, unit_exclude, unit_weights, time_weights));
    return rcpp_result_gen;
END_RCPP
}
// inference
Rcpp::List inference(arma::mat& X, arma::colvec& Y, arma::mat& Y00, arma::mat& Y01, arma::mat& Y10, double dzeta, double rho, double tol, std::string weight_type, std::string resampling_procedure);
RcppExport SEXP _sdid_inference(SEXP XSEXP, SEXP YSEXP, SEXP Y00SEXP, SEXP Y01SEXP, SEXP Y10SEXP, SEXP dzetaSEXP, SEXP rhoSEXP, SEXP tolSEXP, SEXP weight_typeSEXP, SEXP resampling_procedureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y00(Y00SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y01(Y01SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y10(Y10SEXP);
    Rcpp::traits::input_parameter< double >::type dzeta(dzetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type weight_type(weight_typeSEXP);
    Rcpp::traits::input_parameter< std::string >::type resampling_procedure(resampling_procedureSEXP);
    rcpp_result_gen = Rcpp::wrap(inference(X, Y, Y00, Y01, Y10, dzeta, rho, tol, weight_type, resampling_procedure));
    return rcpp_result_gen;
END_RCPP
}
// time_weights
arma::mat time_weights(arma::mat& Y00, arma::mat& Y01, double dzeta, double rho, double tol, std::string weight_type);
RcppExport SEXP _sdid_time_weights(SEXP Y00SEXP, SEXP Y01SEXP, SEXP dzetaSEXP, SEXP rhoSEXP, SEXP tolSEXP, SEXP weight_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Y00(Y00SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y01(Y01SEXP);
    Rcpp::traits::input_parameter< double >::type dzeta(dzetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type weight_type(weight_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(time_weights(Y00, Y01, dzeta, rho, tol, weight_type));
    return rcpp_result_gen;
END_RCPP
}
// unit_weights
arma::mat unit_weights(arma::mat& Y00, arma::mat& Y10, double dzeta, double rho, double tol, std::string weight_type);
RcppExport SEXP _sdid_unit_weights(SEXP Y00SEXP, SEXP Y10SEXP, SEXP dzetaSEXP, SEXP rhoSEXP, SEXP tolSEXP, SEXP weight_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Y00(Y00SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Y10(Y10SEXP);
    Rcpp::traits::input_parameter< double >::type dzeta(dzetaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type weight_type(weight_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(unit_weights(Y00, Y10, dzeta, rho, tol, weight_type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sdid_admm_solver", (DL_FUNC) &_sdid_admm_solver, 5},
    {"_sdid_admm_time_solver", (DL_FUNC) &_sdid_admm_time_solver, 5},
    {"_sdid_att", (DL_FUNC) &_sdid_att, 5},
    {"_sdid_inference", (DL_FUNC) &_sdid_inference, 10},
    {"_sdid_time_weights", (DL_FUNC) &_sdid_time_weights, 6},
    {"_sdid_unit_weights", (DL_FUNC) &_sdid_unit_weights, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_sdid(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
