// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// uirtestimate
Rcpp::List uirtestimate(Rcpp::NumericMatrix data, int model_);
RcppExport SEXP IRTpp_uirtestimate(SEXP dataSEXP, SEXP model_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type model_(model_SEXP);
    __result = Rcpp::wrap(uirtestimate(data, model_));
    return __result;
END_RCPP
}
// abilityinterface
Rcpp::List abilityinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix data, int model_, int method, bool matrix_flag, Rcpp::NumericVector prob_matrix);
RcppExport SEXP IRTpp_abilityinterface(SEXP zita_parSEXP, SEXP dataSEXP, SEXP model_SEXP, SEXP methodSEXP, SEXP matrix_flagSEXP, SEXP prob_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type zita_par(zita_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type model_(model_SEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< bool >::type matrix_flag(matrix_flagSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_matrix(prob_matrixSEXP);
    __result = Rcpp::wrap(abilityinterface(zita_par, data, model_, method, matrix_flag, prob_matrix));
    return __result;
END_RCPP
}
// mapinterface
Rcpp::List mapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, bool matrix_flag, Rcpp::NumericVector prob_matrix);
RcppExport SEXP IRTpp_mapinterface(SEXP zita_parSEXP, SEXP datSEXP, SEXP e_modelSEXP, SEXP matrix_flagSEXP, SEXP prob_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type zita_par(zita_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type matrix_flag(matrix_flagSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_matrix(prob_matrixSEXP);
    __result = Rcpp::wrap(mapinterface(zita_par, dat, e_model, matrix_flag, prob_matrix));
    return __result;
END_RCPP
}
// eapinterface
Rcpp::List eapinterface(Rcpp::NumericMatrix zita_par, Rcpp::NumericMatrix dat, int e_model, bool matrix_flag, Rcpp::NumericVector prob_matrix);
RcppExport SEXP IRTpp_eapinterface(SEXP zita_parSEXP, SEXP datSEXP, SEXP e_modelSEXP, SEXP matrix_flagSEXP, SEXP prob_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type zita_par(zita_parSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type e_model(e_modelSEXP);
    Rcpp::traits::input_parameter< bool >::type matrix_flag(matrix_flagSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type prob_matrix(prob_matrixSEXP);
    __result = Rcpp::wrap(eapinterface(zita_par, dat, e_model, matrix_flag, prob_matrix));
    return __result;
END_RCPP
}
