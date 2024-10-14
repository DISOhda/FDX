#' @name NDLR
#' 
#' @title
#' Wrapper Functions for the Non-Adaptive Discrete Guo-Romano Procedure
#' 
#' @description 
#' `NDLR()` is a wrapper function of [`discrete.LR()`] for computing 
#' non-adaptive \[DLR\]. It simply passes its arguments to [`discrete.LR()`]
#' with fixed `adaptive = FALSE`.
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction TRUE
#' @templateVar critical.values TRUE
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
#' [`discrete.GR()`], [`DGR()`], [`discrete.LR()`], [`DLR()`], [`NDLR()`],
#' [`discrete.PB()`], [`DPB()`], [`NDPB()`]
#' 
#' @references
#'  S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#'  Heterogeneous Tests.
#'  [arXiv:1912.04607v1](https://arxiv.org/abs/1912.04607v1).
#'  
#' @template example
#' @examples
#' 
#' # Non-adaptive DLR without critical values; using results object
#' NDLR.sd.fast <- NDLR(test.results)
#' summary(NDLR.sd.fast)
#' 
#' # Non-adaptive DLR with critical values; using extracted p-values and supports
#' NDLR.sd.crit <- NDLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(NDLR.sd.crit)
#' 
#' # Non-adaptive DLR (step-up) without critical values; using extracted p-values and supports
#' NDLR.su.fast <- NDLR(raw.pvalues, pCDFlist, direction = "su")
#' summary(NDLR.su.fast)
#' 
#' # Non-adaptive DLR (step-up) with critical values; using test results object
#' NDLR.su.crit <- NDLR(test.results, direction = "su", critical.values = TRUE)
#' summary(NDLR.su.crit)
#' 
#' @export
NDLR <- function(test.results, ...) UseMethod("NDLR")

#' @rdname NDLR
#' @export
NDLR.default <- function(
    test.results,
    pCDFlist,
    alpha            = 0.05,
    zeta             = 0.5,
    direction        = "sd",
    critical.values  = FALSE,
    select.threshold = 1,
    pCDFlist.indices = NULL,
    ...
){
  out <- discrete.LR.default(
    test.results     = test.results,
    pCDFlist         = pCDFlist,
    alpha            = alpha,
    zeta             = zeta,
    direction        = direction,
    adaptive         = FALSE, 
    critical.values  = critical.values,
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

#' @rdname NDLR
#' @export
NDLR.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    direction        = "sd",
    critical.values  = FALSE,
    select.threshold = 1,
    ...
){
  out <- discrete.LR.DiscreteTestResults(
    test.results     = test.results,
    alpha            = alpha,
    zeta             = zeta,
    direction        = direction,
    adaptive         = FALSE, 
    critical.values  = critical.values,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
