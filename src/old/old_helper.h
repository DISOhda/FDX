#include <RcppArmadillo.h>
using namespace Rcpp;

NumericVector stepfun(const NumericVector &x, const NumericVector &sfun);
NumericVector stepfun_leq(const NumericVector &x, const NumericVector &sfun);

NumericVector short_eff(const NumericVector &x, const double limit);

void colsortdec(NumericMatrix &mat);
void colsortasc(NumericMatrix &mat);
