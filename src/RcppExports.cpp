// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cumulantDerivatives
NumericVector cumulantDerivatives(double t, IntegerVector x, NumericVector p, NumericVector s);
RcppExport SEXP ncdDetect2_cumulantDerivatives(SEXP tSEXP, SEXP xSEXP, SEXP pSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(cumulantDerivatives(t, x, p, s));
    return rcpp_result_gen;
END_RCPP
}
