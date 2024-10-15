#' @name continuous.LR
#' 
#' @title
#' Continuous Lehmann-Romano procedure
#' 
#' @description
#' Apply the usual (continuous) \[LR\] procedure, with or without computing the
#' critical values, to a set of p-values. A non-adaptive version is available as
#' well.
#' 
#' @templateVar test.results TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar adaptive TRUE
#' @templateVar critical.values TRUE
#' @templateVar select.threshold TRUE
#' @templateVar weights FALSE
#' @template param 
#' 
#' @details
#' `LR` and `NLR` are wrapper functions for `continuous.LR`. The
#' first one simply passes all its arguments to `continuous.LR` with
#' `adaptive = TRUE` and `NLR` does the same with
#' `adaptive = FALSE`.
#'  
#' @seealso
#' [`kernel()`], [`FDX`][FDX-package], [`continuous.GR()`],
#' [`discrete.LR()`], [`discrete.GR()`], 
#' [`discrete.PB()`], [`weighted.LR()`], 
#' [`weighted.GR()`], [`weighted.PB()`]
#' 
#' @references
#' Lehmann, E. L. & Romano, J. P. (2005). Generalizations of the familywise
#'   error rate. *The Annals of Statistics*, *33*(3), pp. 1138-1154.
#'   \doi{10.1214/009053605000000084}
#'   
#' @template example
#' @examples
#' 
#' # LR without critical values; using extracted p-values
#' LR.fast <- LR(raw.pvalues)
#' summary(LR.fast)
#' 
#' # LR with critical values; using test results object
#' LR.crit <- LR(test.results, critical.values = TRUE)
#' summary(LR.crit)
#' 
#' # Non-adaptive LR without critical values; using test results object
#' NLR.fast <- NLR(test.results)
#' summary(NLR.fast)
#' 
#' # Non-adaptive LR with critical values; using extracted p-values
#' NLR.crit <- NLR(raw.pvalues, critical.values = TRUE)
#' summary(NLR.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive TRUE
#' @templateVar Weighting FALSE
#' @template return
#' 
#' @importFrom checkmate assert check_numeric check_r6 qassert
#' @export
continuous.LR <- function(
    test.results,
    alpha = 0.05,
    zeta = 0.5,
    adaptive = TRUE,
    critical.values = FALSE,
    select.threshold = 1
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # test results (p-values)
  assert(
    check_numeric(
      x = test.results,
      lower = 0,
      upper = 1,
      any.missing = FALSE,
      min.len = 1
    ),
    check_r6(
      x = test.results,
      classes = "DiscreteTestResults",
      public = c("get_pvalues", "get_pvalue_supports", "get_support_indices")
    )
  )
  pvals <- if(is.numeric(test.results))
    test.results else
      test.results$get_pvalues()
  
  # FDP level
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # Exceedance probability
  qassert(x = zeta, rules = "N1[0, 1]")
  
  # adaptiveness
  qassert(adaptive, "B1")
  
  # compute and return critical values?
  qassert(critical.values, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- continuous.fdx.int(
    pvec        = pvals,
    method      = "LR",
    alpha       = alpha,
    zeta        = zeta,
    adaptive    = adaptive,
    crit.consts = critical.values,
    threshold   = select.threshold,
    data.name   = deparse(substitute(test.results))
  )
  
  return(output)
}

#' @rdname continuous.LR
#' @export
LR <- function(
    test.results,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- continuous.LR(test.results, alpha, zeta, TRUE, 
                       critical.values, select.threshold)
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}

#' @rdname continuous.LR
#' @export
NLR <- function(
    test.results,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- continuous.LR(test.results, alpha, zeta, FALSE, 
                       critical.values, select.threshold)
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}