#' @name NDPB
#' 
#' @title
#' Wrapper Functions for the Non-Adaptive Discrete Guo-Romano Procedure
#' 
#' @description 
#' `NDPB()` is a wrapper function of [`discrete.PB()`] for computing 
#' non-adaptive \[DPB\]. It simply passes its arguments to [`discrete.PB()`]
#' with fixed `adaptive = FALSE`.
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar critical.values TRUE
#' @templateVar exact TRUE
#' @templateVar select.threshold TRUE
#' @templateVar pCDFlist.indices TRUE
#' @templateVar triple.dots TRUE
#' @templateVar weights FALSE
#' @template param
#' 
#' @template details_crit
#' 
#' @templateVar Adaptive TRUE
#' @templateVar Weighting FALSE
#' @template return
#' 
#' @seealso
#' [`discrete.PB()`], [`DPB()`], [`discrete.GR()`], [`DGR()`], [`NDGR()`],
#' [`discrete.LR()`], [`DLR()`], [`NDLR()`]
#' 
#' @references
#'  S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#'  Heterogeneous Tests.
#'  [arXiv:1912.04607v1](https://arxiv.org/abs/1912.04607v1).
#'  
#' @template example
#' @examples
#' 
#' # Non-adaptive DPB (exact) without critical values; using results object
#' NDPB.exact.fast <- NDPB(test.results)
#' summary(NDPB.exact.fast)
#' 
#' # Non-adaptive DPB (exact) with critical values; using extracted p-values and supports
#' NDPB.exact.crit <- NDPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(NDPB.exact.crit)
#' 
#' # Non-adaptive DPB (normal approx.) without critical values; using extracted p-values and supports
#' NDPB.norm.fast <- NDPB(raw.pvalues, pCDFlist, exact = FALSE)
#' summary(NDPB.norm.fast)
#' 
#' # Non-adaptive DPB (normal approx.) with critical values; using test results object
#' NDPB.norm.crit <- NDPB(test.results, critical.values = TRUE, exact = FALSE)
#' summary(NDPB.norm.crit)
#' 
#' @export
NDPB <- function(test.results, ...) UseMethod("NDPB")

#' @rdname NDPB
#' @export
NDPB.default <- function(
    test.results,
    pCDFlist,
    alpha            = 0.05,
    zeta             = 0.5,
    critical.values  = FALSE,
    exact            = TRUE,
    select.threshold = 1,
    pCDFlist.indices = NULL,
    ...
){
  out <- discrete.PB.default(
    test.results     = test.results,
    pCDFlist         = pCDFlist,
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = FALSE, 
    critical.values  = critical.values,
    exact            = exact,
    select.threshold = select.threshold,
    pCDFlist.indices = pCDFlist.indices,
    ...
  )
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(pCDFlist))
  )
  
  return(out)
}

#' @rdname NDPB
#' @export
NDPB.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    critical.values  = FALSE,
    exact            = TRUE,
    select.threshold = 1,
    ...
){
  out <- discrete.PB.DiscreteTestResults(
    test.results     = test.results,
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = FALSE, 
    critical.values  = critical.values,
    exact            = exact,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
