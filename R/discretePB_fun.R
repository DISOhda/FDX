#' @name discrete.PB
#' 
#' @title
#' Discrete Poisson-Binomial procedure
#' 
#' @description
#' Apply the \[DPB\] procedure, with or without computing the critical values,
#' to a set of p-values and their discrete support. A non-adaptive version is
#' available as well. Additionally, the user can choose between exact
#' computation of the Poisson-Binomial distribution or a refined normal
#' approximation.
#' 
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar adaptive TRUE
#' @templateVar critical.values TRUE
#' @templateVar exact TRUE
#' @templateVar select.threshold TRUE
#' @templateVar pCDFlist.indices TRUE
#' @templateVar triple.dots TRUE
#' @templateVar weights FALSE
#' @template param
#'  
#' @details
#' `DPB` and `NDPB` are wrapper functions for `discrete.PB`.
#' The first one simply passes all its arguments to `discrete.PB` with
#' `adaptive = TRUE` and `NDPB` does the same with
#' `adaptive = FALSE`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][FDX-package], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.LR()`], 
#' [`discrete.GR()`], [`weighted.LR()`],
#' [`weighted.GR()`], [`weighted.PB()`]
#' 
#' @section References:
#'  S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#'  Heterogeneous Tests.
#'  [arXiv:1912.04607v1](https://arxiv.org/abs/1912.04607v1).
#' 
#' @template example
#' @examples
#' 
#' # DPB (exact) without critical values; using results object
#' DPB.exact.fast <- discrete.PB(test.results)
#' summary(DPB.exact.fast)
#' 
#' # DPB (exact) with critical values; using extracted p-values and supports
#' DPB.exact.crit <- discrete.PB(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DPB.exact.crit)
#' 
#' # DPB (normal approximation) without critical values; using extracted p-values and supports
#' DPB.norm.fast <- discrete.PB(raw.pvalues, pCDFlist, exact = FALSE)
#' summary(DPB.norm.fast)
#' 
#' # DPB (normal approximation) with critical values; using results object
#' DPB.norm.crit <- discrete.PB(test.results, critical.values = TRUE,
#'                              exact = FALSE)
#' summary(DPB.norm.crit)
#' 
#' # Non-adaptive DPB (exact) without critical values; using results object
#' NDPB.exact.fast <- discrete.PB(test.results, adaptive = FALSE)
#' summary(NDPB.exact.fast)
#' 
#' # Non-adaptive DPB (exact) with critical values; using extracted p-values and supports
#' NDPB.exact.crit <- discrete.PB(raw.pvalues, pCDFlist, adaptive = FALSE,
#'                                critical.values = TRUE)
#' summary(NDPB.exact.crit)
#' 
#' # Non-adaptive DPB (normal approx.) without critical values; using extracted p-values and supports
#' NDPB.norm.fast <- discrete.PB(raw.pvalues, pCDFlist, adaptive = FALSE,
#'                               exact = FALSE)
#' summary(NDPB.norm.fast)
#' 
#' # Non-adaptive DPB (normal approx.) with critical values; using results object
#' NDPB.norm.crit <- discrete.PB(test.results, adaptive = FALSE,
#'                               critical.values = TRUE, exact = FALSE)
#' summary(NDPB.norm.crit)
#' 
#' @templateVar Critical.values TRUE
#' @template return
#' 
#' @export
discrete.PB <- function(test.results, ...) UseMethod("discrete.PB")

#' @rdname discrete.PB
#' @importFrom checkmate assert_integerish assert_list assert_numeric qassert
#' @export
discrete.PB.default <- function(
    test.results,
    pCDFlist,
    alpha            = 0.05,
    zeta             = 0.5,
    adaptive         = TRUE,
    critical.values  = FALSE,
    exact            = TRUE,
    select.threshold = 1,
    pCDFlist.indices = NULL,
    ...
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # test results (p-values)
  qassert(x = test.results, rules = "N+[0, 1]")
  n <- length(test.results)
  
  # list structure of p-value distributions
  assert_list(
    x = pCDFlist,
    types = "numeric",
    any.missing = FALSE,
    min.len = 1,
    max.len = n
  )
  # individual p-value distributions
  for(i in seq_along(pCDFlist)){
    assert_numeric(
      x = pCDFlist[[i]],
      lower = 0,
      upper = 1,
      any.missing = FALSE,
      min.len = 1,
      sorted = TRUE
    )
    if(max(pCDFlist[[i]]) != 1)
      stop("Last value of each vector in 'pCDFlist' must be 1!")
  }
  m <- length(pCDFlist)
  
  # FDP level
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # Exceedance probability
  qassert(x = zeta, rules = "N1[0, 1]")
  
  # adaptiveness
  qassert(adaptive, "B1")
  
  # compute and return critical values?
  qassert(critical.values, "B1")
  
  # exact or approximate computation of Poisson-Binomial distribution
  qassert(exact, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  # list structure of indices
  assert_list(
    x = pCDFlist.indices,
    types = "numeric",
    any.missing = FALSE,
    len = m,
    unique = TRUE,
    null.ok = TRUE
  )
  # individual index vectors (if not NULL)
  if(is.null(pCDFlist.indices)){
    if(n != m){
      stop(
        paste(
          "If no indices for the p-value CDFs are provided, the lengths of",
          "'test.results' and 'pCDFlist' must be equal!"
        )
      )
    }
    pCDFlist.indices <- as.list(1:n)
    pCDFlist.counts <- rep(1, n)
  } else {
    set <- 1L:n
    for(i in seq_along(pCDFlist.indices)){
      pCDFlist.indices[[i]] <- assert_integerish(
        x = pCDFlist.indices[[i]],
        lower = 1,
        upper = n,
        any.missing = FALSE,
        min.len = 1,
        max.len = n,
        unique = TRUE,
        sorted = TRUE,
        coerce = TRUE
      )
      set <- setdiff(set, pCDFlist.indices[[i]])
    }
    if(length(set))
      stop("'pCDFlist.indices' must contain each p-value index exactly once!")
    pCDFlist.counts <- sapply(pCDFlist.indices, length)
  }
  
  #----------------------------------------------------
  #       check and prepare p-values for processing
  #----------------------------------------------------
  pvec <- match.pvals(test.results, pCDFlist, pCDFlist.indices)
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdx.int(
    pvec             = test.results,
    pCDFlist         = pCDFlist,
    pCDFlist.indices = pCDFlist.indices,
    method           = "PB",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    method.parameter = exact,
    crit.consts      = critical.values,
    threshold        = select.threshold,
    data.name        = paste(
                         deparse(substitute(test.results)),
                         "and",
                         deparse(substitute(pCDFlist))
                       )
  )
  
  return(output)
}

#' @rdname discrete.PB
#' @importFrom checkmate assert_r6 qassert
#' @export
discrete.PB.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    adaptive         = TRUE,
    critical.values  = FALSE,
    exact            = TRUE,
    select.threshold = 1,
    ...
) {
  #----------------------------------------------------
  #       check arguments
  #----------------------------------------------------
  # discrete test results object
  assert_r6(
    x = test.results,
    classes = "DiscreteTestResults",
    public = c("get_pvalues", "get_pvalue_supports", "get_support_indices")
  )
  
  # FDP level
  qassert(x = alpha, rules = "N1[0, 1]")
  
  # Exceedance probability
  qassert(x = zeta, rules = "N1[0, 1]")
  
  # adaptiveness
  qassert(adaptive, "B1")
  
  # compute and return critical values?
  qassert(critical.values, "B1")
  
  # exact or approximate computation of Poisson-Binomial distribution
  qassert(exact, "B1")
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdx.int(
    pvec             = test.results$get_pvalues(),
    pCDFlist         = test.results$get_pvalue_supports(unique = TRUE),
    pCDFlist.indices = test.results$get_support_indices(),
    method           = "PB",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    method.parameter = exact,
    crit.consts      = critical.values,
    threshold        = select.threshold,
    data.name        = deparse(substitute(test.results))
  )
  
  return(output)
}

discrete.PB2 <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE, exact = TRUE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  m <- as.integer(length(raw.pvalues))
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!") 
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(raw.pvalues, pCDFlist)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(critical.values){
    # compute transformed support
    y <- kernel_DPB_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, exact)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals > crit.constants)
  }else{
    # compute transformed sorted p-values
    y <- kernel_DPB_fast(pCDFlist, sorted.pvals, adaptive, alpha, exact)
    # determine significant (transformed) p-values
    idx <- which(y > zeta)
  }
  if(length(idx)){
    m.rej <- min(idx) - 1
    if(m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej])
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  
  if(critical.values){
    y <- y$pval.transf
  }
  # compute adjusted p-values
  pv.adj = cummax(pmin(y, 1))
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = pv.adj[ro]
  
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Poisson-Binomial procedure"
  alg <- if(exact) paste(alg, "(exact)") else paste(alg, "(normal approximation)")
  output$Method <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "FDX"
  return(output)
}
