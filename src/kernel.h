#include "dist_fun.h"
#include "helper.h"

////////// Discrete Lehmann-Romano

struct tau_m_results{
  double value;
  int index;
  std::vector<double> eval;
};

//' @rdname kernel
//' @export
// [[Rcpp::export]]
NumericVector kernel_DLR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const bool stepUp = false, const double zeta = 0.5, const NumericVector &support = NumericVector(), const Nullable<IntegerVector> &pCDFcounts = R_NilValue);

//' @rdname kernel
//' @export
// [[Rcpp::export]]
List kernel_DLR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5, const bool stepUp = false, const Nullable<IntegerVector> &pCDFcounts = R_NilValue);


////////// Weighted Lehmann-Romano

//' @rdname kernel
//' @export
// [[Rcpp::export]]
NumericVector kernel_wLR_fast(const NumericVector &qvalues, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false);


////////// Discrete Guo-Romano
  
//' @rdname kernel
//' @export
// [[Rcpp::export]]
List kernel_DGR_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const Nullable<IntegerVector> &pCDFcounts = R_NilValue);

//' @rdname kernel
//' @export
// [[Rcpp::export]]
List kernel_DGR_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5, const Nullable<IntegerVector> &pCDFcounts = R_NilValue);


////////// Weighted Guo-Romano

//' @rdname kernel
//' @export
// [[Rcpp::export]]
NumericVector kernel_wGR_fast(const NumericVector &qvalues, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false);


////////// Discrete Poisson-Binomial

//' @rdname kernel
//' @export
// [[Rcpp::export]]
NumericVector kernel_DPB_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const bool exact = true, const Nullable<IntegerVector> &pCDFcounts = R_NilValue);//, const Nullable<List> &pCDFindices = R_NilValue);

//' @rdname kernel
//' @export
// [[Rcpp::export]]
List kernel_DPB_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const bool adaptive = true, const double alpha = 0.05, const double zeta = 0.5, const bool exact = true, const Nullable<IntegerVector> &pCDFcounts = R_NilValue);//, const Nullable<List> &pCDFindices = R_NilValue);


////////// Weighted Poisson-Binomial

//' @rdname kernel
//' @export
// [[Rcpp::export]]
NumericVector kernel_wPB_fast(const NumericVector &qvalues, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false, bool exact = true);
