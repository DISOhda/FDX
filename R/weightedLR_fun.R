#' @name weighted.LR
#' 
#' @title Weighted Lehmann-Romano Procedure
#' 
#' @description
#' Apply the weighted \[wLR\] procedure, with or without computing the
#' critical values, to a set of p-values. Both arithmetic and geometric
#' weighting are available.
#' 
#' @templateVar test.results TRUE
#' @templateVar weights TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar weighting.method TRUE
#' @templateVar critical.values TRUE
#' @templateVar select.threshold TRUE
#' @template param
#' 
#' @details
#' `wLR.AM` and `wLR.GM` are wrapper functions for `weighted.LR`.
#' The first one simply passes all its arguments to `weighted.LR` with
#' `weighting.method = "AM"` and `wLR.GM` does the same with
#' `weighting.method = "GM"`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][`FDX-package`], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.LR()`], 
#' [`discrete.GR()`], [`discrete.PB()`], 
#' [`weighted.GR()`], [`weighted.PB()`]
#' 
#' @references
#' Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
#'   Heterogeneous Tests. *Electronic Journal of Statistics*, *14*(2),
#'   pp. 4244-4272. \doi{10.1214/20-EJS1771}
#'   
#' @template exampleWeighted
#' @examples
#' 
#' # arithmetic-weighted Lehmann-Romano procedure without critical values
#' wLR.AM.fast <- wLR.AM(raw.pvalues.weighted, weights)
#' summary(wLR.AM.fast)
#' 
#' # arithmetic-weighted Lehmann-Romano procedure with critical values
#' wLR.AM.crit <- wLR.AM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wLR.AM.crit)
#' 
#' # geometric-weighted Lehmann-Romano procedure without critical values
#' wLR.GM.fast <- wLR.GM(raw.pvalues.weighted, weights)
#' summary(wLR.GM.fast)
#' 
#' # geometric-weighted Lehmann-Romano procedure with critical values
#' wLR.GM.crit <- wLR.GM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wLR.GM.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive FALSE
#' @templateVar Weighting TRUE
#' @template return
#' 
#' @importFrom checkmate assert assert_numeric check_numeric check_r6 qassert
#' @export
weighted.LR <- function(
    test.results,
    weights = NULL,
    alpha = 0.05,
    zeta = 0.5,
    weighting.method = c("AM", "GM"),
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
  n <- length(pvals)
  
  # weights
  assert_numeric(
    x = weights,
    lower = 0,
    finite = TRUE,
    len = n,
    all.missing = FALSE,
    null.ok = TRUE
  )
  if(is.null(weights)) weights <- rep(1, n)
  
  # FDP level
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # Exceedance probability
  qassert(x = zeta, rules = "N1[0, 1]")
  
  # Weighting method
  qassert(x = weighting.method, rules = "S1")
  match.arg(toupper(weighting.method), c("AM", "GM"))
  
  # compute and return critical values?
  qassert(critical.values, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- weighted.fdx.int(
    pvec        = pvals,
    weights     = weights,
    method      = "LR",
    weight.meth = weighting.method,
    alpha       = alpha,
    zeta        = zeta,
    crit.consts = critical.values,
    threshold   = select.threshold,
    data.name   = paste(
      deparse(substitute(test.results)),
      "and",
      deparse(substitute(weights))
    )
  )
  
  return(output)
}

#' @rdname weighted.LR
#' @export
wLR.AM <- function(test.results, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE, select.threshold = 1){
  out <- weighted.LR(test.results, weights, alpha, zeta, "AM", critical.values, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}

#' @rdname weighted.LR
#' @export
wLR.GM <- function(test.results, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE, select.threshold = 1){
  out <- weighted.LR(test.results, weights, alpha, zeta, "GM", critical.values, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}
