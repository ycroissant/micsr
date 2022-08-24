// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// unorm
double unorm(double x);
RcppExport SEXP _micsr_unorm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(unorm(x));
    return rcpp_result_gen;
END_RCPP
}
// bnorm
double bnorm(double h1, double h2, double r);
RcppExport SEXP _micsr_bnorm(SEXP h1SEXP, SEXP h2SEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< double >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(bnorm(h1, h2, r));
    return rcpp_result_gen;
END_RCPP
}
// punorm
NumericVector punorm(NumericVector x);
RcppExport SEXP _micsr_punorm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(punorm(x));
    return rcpp_result_gen;
END_RCPP
}
// pbnorm
NumericVector pbnorm(NumericVector z1, NumericVector z2, NumericVector rho);
RcppExport SEXP _micsr_pbnorm(SEXP z1SEXP, SEXP z2SEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z2(z2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(pbnorm(z1, z2, rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_micsr_unorm", (DL_FUNC) &_micsr_unorm, 1},
    {"_micsr_bnorm", (DL_FUNC) &_micsr_bnorm, 3},
    {"_micsr_punorm", (DL_FUNC) &_micsr_punorm, 1},
    {"_micsr_pbnorm", (DL_FUNC) &_micsr_pbnorm, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_micsr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}