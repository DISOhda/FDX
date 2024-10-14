#' @name DLR
#' 
#' @title
#' Wrapper Functions for the Discrete Guo-Romano Procedure
#' 
#' @description 
#' `DLR()` is a wrapper function of [`discrete.LR()`] for computing \[DLR\]. It
#' simply passes its arguments to [`discrete.LR()`] with fixed
#' `adaptive = TRUE`.
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
#' [`discrete.LR()`], [`NDLR()`], [`discrete.GR()`], [`DGR()`], [`NDGR()`],
#' [`discrete.PB()`], [`DPB()`], [`NDPB()`]
#' 
#' @references
#' Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
#'   Heterogeneous Tests. *Electronic Journal of Statistics*, *14*(2),
#'   pp. 4244-4272. \doi{10.1214/20-EJS1771}
#'  
#' @template example
#' @examples
#' 
#' # DLR without critical values; using results object
#' DLR.sd.fast <- DLR(test.results)
#' summary(DLR.sd.fast)
#' 
#' # DLR with critical values; using extracted p-values and supports
#' DLR.sd.crit <- DLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DLR.sd.crit)
#' 
#' # DLR (step-up) without critical values; using extracted p-values and supports
#' DLR.su.fast <- DLR(raw.pvalues, pCDFlist, direction = "su")
#' summary(DLR.su.fast)
#' 
#' # DLR (step-up) with critical values; using test results object
#' DLR.su.crit <- DLR(test.results, direction = "su", critical.values = TRUE)
#' summary(DLR.su.crit)
#' 
#' @export
DLR <- function(test.results, ...) UseMethod("DLR")

#' @rdname DLR
#' @export
DLR.default <- function(
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
    adaptive         = TRUE, 
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

#' @rdname DLR
#' @export
DLR.DiscreteTestResults <- function(
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
    adaptive         = TRUE, 
    critical.values  = critical.values,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
