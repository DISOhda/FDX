#' @name continuous.GR
#' 
#' @title
#' Continuous Guo-Romano procedure
#' 
#' @description
#' Apply the usual continuous \[GR\] procedure, with or without computing the
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
#' `GR` and `NGR` are wrapper functions for `continuous.GR`. The
#' first one simply passes all its arguments to `continuous.GR` with
#' `adaptive = TRUE` and `NGR` does the same with
#' `adaptive = FALSE`.
#' 
#' @seealso
#' [`kernel`], [`FDX-package`], [`continuous.LR()`],
#' [`discrete.LR()`], [`discrete.GR()`], 
#' [`discrete.PB()`], [`weighted.LR()`], 
#' [`weighted.GR()`], [`weighted.PB()`]
#'
#' @references
#' Guo, W. & Romano, J. P. (2007). A generalized Sidak-Holm procedure and
#'   control of generalized error rates under independence.
#'   *Statistical Applications in Genetics and Molecular Biology*, *6*(1),
#'   Art. 3, 35 pp. (eletronic). \doi{10.2202/1544-6115.1247}
#'   
#' @template example
#' @examples
#' 
#' # GR without critical values; using extracted p-values
#' GR.fast <- GR(raw.pvalues)
#' summary(GR.fast)
#' 
#' # LR with critical values; using test results object
#' GR.crit <- GR(test.results, critical.values = TRUE)
#' summary(GR.crit)
#' 
#' # Non-adaptive GR without critical values; using test results object
#' NGR.fast <- NGR(test.results)
#' summary(NGR.fast)
#' 
#' # Non-adaptive GR with critical values; using extracted p-values
#' NGR.crit <- NGR(raw.pvalues, critical.values = TRUE)
#' summary(NGR.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive TRUE
#' @templateVar Weighting FALSE
#' @template return
#' 
#' @importFrom checkmate assert check_numeric check_r6 qassert
#' @export
continuous.GR <- function(
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
    method      = "GR",
    alpha       = alpha,
    zeta        = zeta,
    adaptive    = adaptive,
    crit.consts = critical.values,
    threshold   = select.threshold,
    data.name   = deparse(substitute(test.results))
  )
  
  return(output)
}

#' @rdname continuous.GR
#' @export
GR <- function(
    test.results,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- continuous.GR(test.results, alpha, zeta, TRUE, 
                       critical.values, select.threshold)
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}

#' @rdname continuous.GR
#' @export
NGR <- function(
    test.results,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- continuous.GR(test.results, alpha, zeta, FALSE, 
                       critical.values, select.threshold)
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}