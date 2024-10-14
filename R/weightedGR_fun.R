#' @name weighted.GR
#' 
#' @title
#' Weighted Guo-Romano Procedure
#' 
#' @description
#' Apply the weighted \[wGR\] procedure, with or without computing the
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
#' `wGR.AM` and `wGR.GM` are wrapper functions for `weighted.GR`.
#' The first one simply passes all its arguments to `weighted.GR` with
#' `weighting.method = "AM"` and `wGR.GM` does the same with
#' `weighting.method = "GM"`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][`FDX-package`], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.LR()`], 
#' [`discrete.GR()`], [`discrete.PB()`], 
#' [`weighted.LR()`], [`weighted.PB()`]
#' 
#' @references
#' Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
#'   Heterogeneous Tests. *Electronic Journal of Statistics*, *14*(2),
#'   pp. 4244-4272. \doi{10.1214/20-EJS1771}
#' 
#' @template exampleWeighted
#' @examples
#' 
#' # arithmetic-weighted Guo-Romano procedure without critical values
#' wGR.AM.fast <- wGR.AM(raw.pvalues.weighted, weights)
#' summary(wGR.AM.fast)
#' 
#' # arithmetic-weighted Guo-Romano procedure with critical values
#' wGR.AM.crit <- wGR.AM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wGR.AM.crit)
#' 
#' # geometric-weighted Guo-Romano procedure without critical values
#' wGR.GM.fast <- wGR.GM(raw.pvalues.weighted, weights)
#' summary(wGR.GM.fast)
#' 
#' # geometric-weighted Guo-Romano procedure with critical values
#' wGR.GM.crit <- wGR.GM(raw.pvalues.weighted, weights, critical.values = TRUE)
#' summary(wGR.GM.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive FALSE
#' @templateVar Weighting TRUE
#' @template return
#' 
#' @importFrom checkmate assert assert_numeric check_numeric check_r6 qassert
#' @export
weighted.GR <- function(
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
    method      = "GR",
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

#' @importFrom pracma fzero
#' @importFrom stats qbeta
weighted.GR2 <- function(raw.pvalues, weights, alpha = 0.05, zeta = 0.5, weighting.method = "AM", critical.values = FALSE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(raw.pvalues) || !is.numeric(raw.pvalues) || all(is.na(raw.pvalues)) || any(raw.pvalues < 0 | raw.pvalues > 1))
    stop("All values of 'raw.pvalues' must be probabilities between 0 and 1!")
  
  m <- length(raw.pvalues)
  if(m != length(weights)) stop("The lengths of 'raw.pvalues' and 'weights' must be equal!")
  
  if(is.null(weights) || !is.numeric(weights) || all(is.na(weights)))
    stop("All values of 'weights' must be numeric!")
  
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  weighting.method <- match.arg(weighting.method, c("AM", "GM"))
  #--------------------------------------------
  #       remove NA's from raw p.values
  #--------------------------------------------
  non.missing <- !is.na(raw.pvalues) & !is.na(weights)
  raw.pvalues.ok <- raw.pvalues[non.missing]
  weights.ok <- weights[non.missing]
  #--------------------------------------------
  #       number of tests, m(l) and [alpha * m] + 1
  #--------------------------------------------
  m <- length(raw.pvalues.ok)
  a <- floor(alpha * 1:m)
  m.l <- m - 1:m + a + 1
  #--------------------------------------------
  #       Rescale weights
  #--------------------------------------------
  weights.rescaled <- weights.ok / mean(weights.ok)
  weights.decreasing <- sort(weights.rescaled, decreasing = TRUE)
  #--------------------------------------------
  #       Compute weighted p-values
  #--------------------------------------------
  switch(weighting.method,
         AM = {qvalues <- pmin(raw.pvalues.ok / weights.rescaled, 1)},
         GM = {qvalues <- geom_weight(raw.pvalues.ok, 1 / weights.rescaled)}
  )
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(qvalues)
  ro <- order(o)
  qvalues.sorted <- qvalues[o]
  #--------------------------------------------
  #        Compute [wLR-AM] or [wLR-GM] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  # Compute transformed weighted p-values
  y <- pmin(1, cummax(kernel_wGR_fast(qvalues.sorted, weights.decreasing, alpha, weighting.method == "GM")))
  # determine significant (transformed) p-values
  idx <- which(y > zeta)
  if(length(idx)){
    m.rej <- min(idx) - 1
    if (m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(qvalues <= qvalues.sorted[m.rej])
      pvec.rej <- raw.pvalues.ok[idx]
    }else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues.ok
  }
  # find critical values
  if(critical.values){
    if(weighting.method == "AM"){
      crit <- numeric(m)
      for(l in 1:m) crit[l] = fzero(function(x, w, b, z){ppbinom(b, pmin(w * x, 1), method = "GeoMeanCounter", lower.tail = FALSE) - z}, c(max(crit), 1), w = weights.decreasing[1:m.l[l]], b = a[l], z = zeta, tol = .Machine$double.neg.eps)$x
    }else{
      crit <- geom_weight(qbeta(zeta, a + 1, m.l - a), (1:m / cumsum(weights.decreasing))[m.l])
    }
  }
  
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej, Adjusted = y[ro], Weighted = qvalues)
  
  # add critical values to output list
  if(critical.values) output$Critical.values <- crit
  
  # include details of the used algorithm as strings
  alg <- "Weighted Guo-Romano procedure"
  output$Method <- if(weighting.method == "AM") paste(alg, "(arithmetic weighting)") else paste(alg, "(geometric weighting)")
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Weighting <- ifelse(weighting.method == "AM", "Arithmetic", "Geometric")
  
  # original test data
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$weights <- weights
  # object names of the data as strings
  output$Data$data.name <- deparse(substitute(raw.pvalues))
  
  class(output) <- "FDX"
  return(output)
}

#' @rdname weighted.GR
#' @export
wGR.AM <- function(
    test.results,
    weights,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- weighted.GR(test.results, weights, alpha, zeta, "AM",
                     critical.values, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}

#' @rdname weighted.GR
#' @export
wGR.GM <- function(
    test.results,
    weights,
    alpha = 0.05,
    zeta = 0.5,
    critical.values = FALSE,
    select.threshold = 1
) {
  out <- weighted.GR(test.results, weights, alpha, zeta, "GM",
                     critical.values, select.threshold)
  
  out$Data$Data.name <- paste(
    deparse(substitute(test.results)),
    "and",
    deparse(substitute(weights))
  )
  
  return(out)
}