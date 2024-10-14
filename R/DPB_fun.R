#' @name DPB
#' 
#' @title
#' Wrapper Functions for the Discrete Guo-Romano Procedure
#' 
#' @description 
#' `DPB()` is a wrapper function of [`discrete.PB()`] for computing \[DPB\]. It
#' simply passes its arguments to [`discrete.PB()`] with fixed
#' `adaptive = TRUE`.
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
#' [`discrete.PB()`], [`NDPB()`], [`discrete.GR()`], [`DGR()`], [`NDGR()`],
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
#' # DPB (exact) without critical values; using results object
#' DPB.exact.fast <- DPB(test.results)
#' summary(DPB.exact.fast)
#' 
#' # DPB (exact) with critical values; using extracted p-values and supports
#' DPB.exact.crit <- DPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DPB.exact.crit)
#' 
#' # DPB (normal approximation) without critical values; using extracted p-values and supports
#' DPB.norm.fast <- DPB(raw.pvalues, pCDFlist, exact = FALSE)
#' summary(DPB.norm.fast)
#' 
#' # DPB (normal approximation) with critical values; using test results object
#' DPB.norm.crit <- DPB(test.results, critical.values = TRUE, exact = FALSE)
#' summary(DPB.norm.crit)
#' 
#' @export
DPB <- function(test.results, ...) UseMethod("DPB")

#' @rdname DPB
#' @export
DPB.default <- function(
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
    adaptive         = TRUE, 
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

#' @rdname discrete.PB
#' @export
DPB.DiscreteTestResults <- function(
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
    adaptive         = TRUE, 
    critical.values  = critical.values,
    exact            = exact,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
