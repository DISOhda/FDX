#' @name weighted.PB
#' 
#' @title
#' Weighted Poisson-Binomial Procedure
#' 
#' @description
#' Apply the weighted \[wPB\] procedure, with or without computing the
#' critical values, to a set of p-values. Both arithmetic and geometric
#' weighting are available. Additionally, the user can choose between exact
#' computation of the Poisson-Binomial distribution or a refined normal
#' approximation.
#' 
#' @templateVar test.results TRUE
#' @templateVar weights TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar weighting.method TRUE
#' @templateVar critical.values TRUE
#' @templateVar exact TRUE
#' @templateVar select.threshold TRUE
#' @template param
#' 
#' @details
#' `wPB.AM` and `wPB.GM` are wrapper functions for `weighted.PB`.
#' The first one simply passes all its arguments to `weighted.PB` with
#' `weighting.method = "AM"` and `wPB.GM` does the same with
#' `weighting.method = "GM"`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][`FDX-package`], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.LR()`], 
#' [`discrete.GR()`], [`discrete.PB()`], 
#' [`weighted.LR()`], [`weighted.GR()`]
#' 
#' @references
#' Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
#'   Heterogeneous Tests. *Electronic Journal of Statistics*, *14*(2),
#'   pp. 4244-4272. \doi{10.1214/20-EJS1771}
#'   
#' @template exampleWeighted
#' @examples
#' 
#' # arithmetic-weighted Poisson-binomial procedure without critical values
#' wPB.AM.fast <- wPB.AM(raw.pvalues.weighted, weights)
#' summary(wPB.AM.fast)
#' 
#' # arithmetic-weighted Poisson-binomial procedure with critical values
#' wPB.AM.crit <- wPB.AM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wPB.AM.crit)
#' 
#' # geometric-weighted Poisson-binomial procedure without critical values
#' wPB.GM.fast <- wPB.GM(raw.pvalues.weighted, weights)
#' summary(wPB.GM.fast)
#' 
#' # geometric-weighted Poisson-binomial procedure with critical values
#' wPB.GM.crit <- wPB.GM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wPB.GM.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive FALSE
#' @templateVar Weighting TRUE
#' @template return
#' 
#' @importFrom checkmate assert assert_numeric check_numeric check_r6 qassert
#' @export
weighted.PB <- function(
    test.results,
    weights = NULL,
    alpha = 0.05,
    zeta = 0.5,
    weighting.method = c("AM", "GM"),
    critical.values = FALSE,
    exact = TRUE,
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
  
  # exact computation of Poisson-Binomial distribution or normal approximation?
  qassert(exact, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- weighted.fdx.int(
    pvec        = pvals,
    weights     = weights,
    method      = "PB",
    weight.meth = weighting.method,
    alpha       = alpha,
    zeta        = zeta,
    exact       = exact,
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

#' @rdname weighted.PB
#' @export
wPB.AM <- function(test.results, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE, exact = TRUE, select.threshold = 1){
  out <- weighted.PB(test.results, weights, alpha, zeta, "AM", critical.values, exact, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}

#' @rdname weighted.PB
#' @export
wPB.GM <- function(test.results, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE, exact = TRUE, select.threshold = 1){
  out <- weighted.PB(test.results, weights, alpha, zeta, "GM", critical.values, exact, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}
