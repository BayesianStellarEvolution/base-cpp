// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// listFilters
std::vector<std::string> listFilters();
RcppExport SEXP _BASE9_listFilters() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(listFilters());
    return rcpp_result_gen;
END_RCPP
}
// getAGBt_zmass
double getAGBt_zmass();
RcppExport SEXP _BASE9_getAGBt_zmass() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getAGBt_zmass());
    return rcpp_result_gen;
END_RCPP
}
// initBase
void initBase(std::string modelDir, std::string msModel, int wdModel, int ifmr, bool distModIsParallax);
RcppExport SEXP _BASE9_initBase(SEXP modelDirSEXP, SEXP msModelSEXP, SEXP wdModelSEXP, SEXP ifmrSEXP, SEXP distModIsParallaxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type modelDir(modelDirSEXP);
    Rcpp::traits::input_parameter< std::string >::type msModel(msModelSEXP);
    Rcpp::traits::input_parameter< int >::type wdModel(wdModelSEXP);
    Rcpp::traits::input_parameter< int >::type ifmr(ifmrSEXP);
    Rcpp::traits::input_parameter< bool >::type distModIsParallax(distModIsParallaxSEXP);
    initBase(modelDir, msModel, wdModel, ifmr, distModIsParallax);
    return R_NilValue;
END_RCPP
}
// setClusterParameters
void setClusterParameters(double age, double feh, double distMod, double av, double y, double carbonicity);
RcppExport SEXP _BASE9_setClusterParameters(SEXP ageSEXP, SEXP fehSEXP, SEXP distModSEXP, SEXP avSEXP, SEXP ySEXP, SEXP carbonicitySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type age(ageSEXP);
    Rcpp::traits::input_parameter< double >::type feh(fehSEXP);
    Rcpp::traits::input_parameter< double >::type distMod(distModSEXP);
    Rcpp::traits::input_parameter< double >::type av(avSEXP);
    Rcpp::traits::input_parameter< double >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type carbonicity(carbonicitySEXP);
    setClusterParameters(age, feh, distMod, av, y, carbonicity);
    return R_NilValue;
END_RCPP
}
// changeModels
void changeModels(std::string msModel, int wdModel, int ifmr, bool distModIsParallax);
RcppExport SEXP _BASE9_changeModels(SEXP msModelSEXP, SEXP wdModelSEXP, SEXP ifmrSEXP, SEXP distModIsParallaxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type msModel(msModelSEXP);
    Rcpp::traits::input_parameter< int >::type wdModel(wdModelSEXP);
    Rcpp::traits::input_parameter< int >::type ifmr(ifmrSEXP);
    Rcpp::traits::input_parameter< bool >::type distModIsParallax(distModIsParallaxSEXP);
    changeModels(msModel, wdModel, ifmr, distModIsParallax);
    return R_NilValue;
END_RCPP
}
// setIFMRParameters
void setIFMRParameters(double intercept, double slope, double quadCoef);
RcppExport SEXP _BASE9_setIFMRParameters(SEXP interceptSEXP, SEXP slopeSEXP, SEXP quadCoefSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< double >::type slope(slopeSEXP);
    Rcpp::traits::input_parameter< double >::type quadCoef(quadCoefSEXP);
    setIFMRParameters(intercept, slope, quadCoef);
    return R_NilValue;
END_RCPP
}
// evolve
std::vector<double> evolve(double mass1, double mass2);
RcppExport SEXP _BASE9_evolve(SEXP mass1SEXP, SEXP mass2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mass1(mass1SEXP);
    Rcpp::traits::input_parameter< double >::type mass2(mass2SEXP);
    rcpp_result_gen = Rcpp::wrap(evolve(mass1, mass2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BASE9_listFilters", (DL_FUNC) &_BASE9_listFilters, 0},
    {"_BASE9_getAGBt_zmass", (DL_FUNC) &_BASE9_getAGBt_zmass, 0},
    {"_BASE9_initBase", (DL_FUNC) &_BASE9_initBase, 5},
    {"_BASE9_setClusterParameters", (DL_FUNC) &_BASE9_setClusterParameters, 6},
    {"_BASE9_changeModels", (DL_FUNC) &_BASE9_changeModels, 4},
    {"_BASE9_setIFMRParameters", (DL_FUNC) &_BASE9_setIFMRParameters, 3},
    {"_BASE9_evolve", (DL_FUNC) &_BASE9_evolve, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_BASE9(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}