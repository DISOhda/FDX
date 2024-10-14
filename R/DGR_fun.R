#' @name DGR
#' 
#' @title
#' Wrapper Functions for the Discrete Guo-Romano Procedure
#' 
#' @description 
#' `DGR()` is a wrapper function of [`discrete.GR()`] for computing \[DGR\]. It
#' simply passes its arguments to [`discrete.GR()`] with fixed
#' `adaptive = TRUE`.
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
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
#' [`discrete.GR()`], [`NDGR()`], [`discrete.LR()`], [`DLR()`], [`NDLR()`],
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
#' # DGR without critical values; using extracted p-values and supports
#' DGR.fast <- DGR(raw.pvalues, pCDFlist)
#' summary(DGR.fast)
#' 
#' # DGR with critical values; using test results object
#' DGR.crit <- DGR(test.results, critical.values = TRUE)
#' summary(DGR.crit)
#' 
#' @export
DGR <- function(test.results, ...) UseMethod("DGR")

#' @rdname DGR
#' @export
DGR.default <- function(
    test.results,
    pCDFlist,
    alpha            = 0.05,
    zeta             = 0.5,
    critical.values  = FALSE,
    select.threshold = 1,
    pCDFlist.indices = NULL,
    ...
){
  out <- discrete.GR.default(
    test.results     = test.results,
    pCDFlist         = pCDFlist,
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = TRUE,
    critical.values  = critical.values,
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

#' @rdname DGR
#' @export
DGR.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    critical.values  = FALSE,
    select.threshold = 1,
    ...
) {
  out <- discrete.GR.DiscreteTestResults(
    test.results     = test.results,
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = TRUE,
    critical.values  = critical.values,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
