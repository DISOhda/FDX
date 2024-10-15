#include "dist_fun.h"
#include "helper.h"

//' @title Kernel functions
//' @name kernel
//' 
//' @keywords internal
//' 
//' @description
//' Kernel functions transform observed p-values or their support according to
//' \[HLR\], \[PB\] and \[HGR\]. The output is used by [`discrete.LR()`],
//' [`discrete.PB()`] and [`discrete.GR()`], respectively.
//' For each procedure, there is a kernel for fast computation and one for
//' calculation of critical values. Kernels followed by ".crit", e.g.
//' `kernel.DGR.crit`, compute and return these critical values, while
//' kernels ending in ".fast" only transform p-values and are therefore faster.
//' The end user should not use these functions directly.
//' 
//' **Note**: As of version 2.0, these functions are purely internal functions!
//' As a consequence, they have to be called directly via `:::`, e.g. 
//' `FDX:::kernel_DGR_fast()`. But users should **not** rely on them, as
//' parameters (including their names, order, etc.) may be changed without
//' notice!
//' 
//' @seealso
//' [`FDX`][`FDX-package`], [`discrete.LR()`]
//' [`discrete.GR()`], [`discrete.PB()`],
//' [`weighted.LR()`], [`weighted.GR()`],
//' [`discrete.PB()`]
//'
//' @templateVar pCDFlist TRUE
//' @templateVar adaptive TRUE
//' @templateVar alpha TRUE
//' @templateVar zeta TRUE
//' @templateVar exact TRUE
//' @templateVar weights FALSE
//' @template param
//' 
//' @param sorted_pv         numeric vector containing the raw p-values, sorted
//'                          in increasing order.
//' @param stepUp            single boolean specifying whether to conduct the
//'                          step-up (`TRUE`) or step-down (`FALSE`; the
//'                          default) version of the discrete Lehmann-Romano
//'                          procedure.
//' @param support           numeric vector, sorted in increasing order, that
//'                          contains the entirety of all observable values of
//'                          the p-value supports; for `kernel_DLR_fast()`, it
//'                          is ignored if `stepUp = FALSE`.
//' @param pCDFcounts        integer vector of counts that indicates to how many
//'                          p-values each **unique** p-value distribution
//'                          belongs.
//' @param sorted_w_pv       numeric vector containing the weighted p-values,
//'                          sorted in increasing order.
//' @param weights           numeric vector containing the **rescaled** weights,
//'                          sorted in **de**creasing order.
//' @param geom_weighting    a boolean specifying whether to conduct geometric
//'                          (`TRUE`) or arithmetic (`FALSE`)
//'                          weighting.
//'
//' @template example
//' @examples
//' 
//' alpha <- 0.05
//' 
//' # If not searching for critical constants, we use only the observed p-values
//' sorted.pvals <- sort(raw.pvalues)
//' y.DLR.fast <- FDX:::kernel_DLR_fast(pCDFlist, sorted.pvals, TRUE)
//' y.NDGR.fast <- FDX:::kernel_DGR_fast(pCDFlist, sorted.pvals, FALSE)$pval.transf
//' # transformed values
//' y.DLR.fast
//' y.NDGR.fast
//' 
//' # compute support
//' pv.list <- sort(unique(unlist(pCDFlist)))
//' y.DGR.crit <- FDX:::kernel_DGR_crit(pCDFlist, pv.list, sorted.pvals, TRUE)
//' y.NDPB.crit <- FDX:::kernel_DPB_crit(pCDFlist, pv.list, sorted.pvals, FALSE)
//' # critical constants
//' y.DGR.crit$crit.consts
//' y.NDPB.crit$crit.consts
//' # transformed values
//' y.DGR.crit$pval.transf
//' y.NDPB.crit$pval.transf
//' 
//' @return
//' For ".fast" kernels, a vector of transformed p-values is returned; ".crit"
//' kernels return a list object with critical constants (`$crit.consts`)
//' and transformed p-values (`$pval.transf`).
//' 

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
NumericVector kernel_wLR_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false);


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
NumericVector kernel_wGR_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false);


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
NumericVector kernel_wPB_fast(const NumericVector &sorted_w_pv, const NumericVector &weights, double alpha = 0.05, bool geom_weighting = false, bool exact = true);
