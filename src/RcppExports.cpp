// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// computeBOAEigen
void computeBOAEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, Eigen::Map<Eigen::VectorXd> wc, Eigen::Map<Eigen::VectorXd> w0c, Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> Regc, Eigen::Map<Eigen::VectorXd> Bc, Eigen::Map<Eigen::VectorXd> Vc, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeBOAEigen(SEXP awakeSEXP, SEXP etaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP wcSEXP, SEXP w0cSEXP, SEXP RcSEXP, SEXP RegcSEXP, SEXP BcSEXP, SEXP VcSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type wc(wcSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type w0c(w0cSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Rc(RcSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Regc(RegcSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Bc(BcSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Vc(VcSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    computeBOAEigen(awake, eta, experts, weights, y, predictions, wc, w0c, Rc, Regc, Bc, Vc, loss_name, loss_tau, loss_gradient);
    return R_NilValue;
END_RCPP
}
// computeEWAEigen
double computeEWAEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, Eigen::Map<Eigen::VectorXd> w0c, double eta, double cumulativeLoss, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeEWAEigen(SEXP awakeSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP w0cSEXP, SEXP etaSEXP, SEXP cumulativeLossSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type w0c(w0cSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type cumulativeLoss(cumulativeLossSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEWAEigen(awake, experts, weights, y, predictions, w0c, eta, cumulativeLoss, loss_name, loss_tau, loss_gradient));
    return rcpp_result_gen;
END_RCPP
}
// computeMLPolCPP
double computeMLPolCPP(NumericMatrix awake, NumericMatrix eta, NumericMatrix experts, NumericMatrix weights, NumericVector y, NumericVector predictions, NumericVector R, NumericVector w, double B, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeMLPolCPP(SEXP awakeSEXP, SEXP etaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP RSEXP, SEXP wSEXP, SEXP BSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMLPolCPP(awake, eta, experts, weights, y, predictions, R, w, B, loss_name, loss_tau, loss_gradient));
    return rcpp_result_gen;
END_RCPP
}
// computeMLPolEigen
double computeMLPolEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, Eigen::Map<Eigen::VectorXd> R, Eigen::Map<Eigen::VectorXd> w, double B, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeMLPolEigen(SEXP awakeSEXP, SEXP etaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP RSEXP, SEXP wSEXP, SEXP BSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type R(RSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMLPolEigen(awake, eta, experts, weights, y, predictions, R, w, B, loss_name, loss_tau, loss_gradient));
    return rcpp_result_gen;
END_RCPP
}
// computeMLPolEigenSimpleLoss
double computeMLPolEigenSimpleLoss(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, Eigen::Map<Eigen::VectorXd> Rc, Eigen::Map<Eigen::VectorXd> wc, double B, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeMLPolEigenSimpleLoss(SEXP awakeSEXP, SEXP etaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP RcSEXP, SEXP wcSEXP, SEXP BSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Rc(RcSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type wc(wcSEXP);
    Rcpp::traits::input_parameter< double >::type B(BSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMLPolEigenSimpleLoss(awake, eta, experts, weights, y, predictions, Rc, wc, B, loss_name, loss_tau, loss_gradient));
    return rcpp_result_gen;
END_RCPP
}
// computeMLProdEigen
void computeMLProdEigen(Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> eta, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> predictions, Eigen::Map<Eigen::VectorXd> R, Eigen::Map<Eigen::VectorXd> L, double maxloss, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeMLProdEigen(SEXP awakeSEXP, SEXP etaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP ySEXP, SEXP predictionsSEXP, SEXP RSEXP, SEXP LSEXP, SEXP maxlossSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type predictions(predictionsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type R(RSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type maxloss(maxlossSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    computeMLProdEigen(awake, eta, experts, weights, y, predictions, R, L, maxloss, loss_name, loss_tau, loss_gradient);
    return R_NilValue;
END_RCPP
}
// RidgeCalibStep1
size_t RidgeCalibStep1(size_t tp1, double dbestlambda, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::MatrixXd> wlambda, Eigen::Map<Eigen::MatrixXd> w0, Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, Eigen::Map<Eigen::MatrixXd> gridlambda, Eigen::Map<Eigen::MatrixXd> predlambda, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> lambda, Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction);
RcppExport SEXP _opera_RidgeCalibStep1(SEXP tp1SEXP, SEXP dbestlambdaSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP wlambdaSEXP, SEXP w0SEXP, SEXP AtSEXP, SEXP btSEXP, SEXP gridlambdaSEXP, SEXP predlambdaSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP cumulativelossSEXP, SEXP predictionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type tp1(tp1SEXP);
    Rcpp::traits::input_parameter< double >::type dbestlambda(dbestlambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type wlambda(wlambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type At(AtSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type bt(btSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type gridlambda(gridlambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type predlambda(predlambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type cumulativeloss(cumulativelossSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type prediction(predictionSEXP);
    rcpp_result_gen = Rcpp::wrap(RidgeCalibStep1(tp1, dbestlambda, experts, weights, wlambda, w0, At, bt, gridlambda, predlambda, y, lambda, cumulativeloss, prediction));
    return rcpp_result_gen;
END_RCPP
}
// RidgeCalibStep2
void RidgeCalibStep2(Eigen::Map<Eigen::MatrixXd> wlambda, Eigen::Map<Eigen::VectorXd> w0, Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, Eigen::Map<Eigen::VectorXd> gridlambda);
RcppExport SEXP _opera_RidgeCalibStep2(SEXP wlambdaSEXP, SEXP w0SEXP, SEXP AtSEXP, SEXP btSEXP, SEXP gridlambdaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type wlambda(wlambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type At(AtSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type bt(btSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type gridlambda(gridlambdaSEXP);
    RidgeCalibStep2(wlambda, w0, At, bt, gridlambda);
    return R_NilValue;
END_RCPP
}
// computeRidgeCPP
size_t computeRidgeCPP(Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> w, Eigen::Map<Eigen::MatrixXd> At, Eigen::Map<Eigen::MatrixXd> bt, Eigen::Map<Eigen::VectorXd> y);
RcppExport SEXP _opera_computeRidgeCPP(SEXP expertsSEXP, SEXP wSEXP, SEXP AtSEXP, SEXP btSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type At(AtSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type bt(btSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(computeRidgeCPP(experts, w, At, bt, y));
    return rcpp_result_gen;
END_RCPP
}
// computeEWACalib
size_t computeEWACalib(size_t tp1, double dbesteta, Eigen::Map<Eigen::MatrixXd> awake, Eigen::Map<Eigen::MatrixXd> experts, Eigen::Map<Eigen::MatrixXd> weights, Eigen::Map<Eigen::MatrixXd> weta, Eigen::Map<Eigen::MatrixXd> w0, Eigen::Map<Eigen::MatrixXd> grideta, Eigen::Map<Eigen::VectorXd> y, Eigen::Map<Eigen::VectorXd> eta, Eigen::Map<Eigen::VectorXd> cumulativeloss, Eigen::Map<Eigen::VectorXd> prediction, String loss_name, double loss_tau, bool loss_gradient);
RcppExport SEXP _opera_computeEWACalib(SEXP tp1SEXP, SEXP dbestetaSEXP, SEXP awakeSEXP, SEXP expertsSEXP, SEXP weightsSEXP, SEXP wetaSEXP, SEXP w0SEXP, SEXP gridetaSEXP, SEXP ySEXP, SEXP etaSEXP, SEXP cumulativelossSEXP, SEXP predictionSEXP, SEXP loss_nameSEXP, SEXP loss_tauSEXP, SEXP loss_gradientSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< size_t >::type tp1(tp1SEXP);
    Rcpp::traits::input_parameter< double >::type dbesteta(dbestetaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type awake(awakeSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type experts(expertsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weta(wetaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type grideta(gridetaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type cumulativeloss(cumulativelossSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< String >::type loss_name(loss_nameSEXP);
    Rcpp::traits::input_parameter< double >::type loss_tau(loss_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type loss_gradient(loss_gradientSEXP);
    rcpp_result_gen = Rcpp::wrap(computeEWACalib(tp1, dbesteta, awake, experts, weights, weta, w0, grideta, y, eta, cumulativeloss, prediction, loss_name, loss_tau, loss_gradient));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_opera_computeBOAEigen", (DL_FUNC) &_opera_computeBOAEigen, 15},
    {"_opera_computeEWAEigen", (DL_FUNC) &_opera_computeEWAEigen, 11},
    {"_opera_computeMLPolCPP", (DL_FUNC) &_opera_computeMLPolCPP, 12},
    {"_opera_computeMLPolEigen", (DL_FUNC) &_opera_computeMLPolEigen, 12},
    {"_opera_computeMLPolEigenSimpleLoss", (DL_FUNC) &_opera_computeMLPolEigenSimpleLoss, 12},
    {"_opera_computeMLProdEigen", (DL_FUNC) &_opera_computeMLProdEigen, 12},
    {"_opera_RidgeCalibStep1", (DL_FUNC) &_opera_RidgeCalibStep1, 14},
    {"_opera_RidgeCalibStep2", (DL_FUNC) &_opera_RidgeCalibStep2, 5},
    {"_opera_computeRidgeCPP", (DL_FUNC) &_opera_computeRidgeCPP, 5},
    {"_opera_computeEWACalib", (DL_FUNC) &_opera_computeEWACalib, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_opera(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
