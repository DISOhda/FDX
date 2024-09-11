#' @name discrete.GR
#' 
#' @title
#' Discrete Guo-Romano procedure
#' 
#' @description
#' Apply the \[DGR\] procedure, with or without computing the critical values, to
#' a set of p-values and their discrete support. A non-adaptive version is
#' available as well.
#' 
#' @details
#' `DGR` and `NDGR` are wrapper functions for `discrete.GR`.
#' The first one simply passes all its arguments to `discrete.GR` with
#' `adaptive = TRUE` and `NDGR` does the same with
#' `adaptive = FALSE`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][FDX-package], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.LR()`], 
#' [`discrete.PB()`], [`weighted.LR()`],
#' [`weighted.GR()`], [`weighted.PB()`]
#' 
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction FALSE
#' @templateVar adaptive TRUE
#' @templateVar critical.values TRUE
#' @templateVar exact FALSE
#' @templateVar pvalues FALSE
#' @templateVar sorted_pv FALSE
#' @templateVar stepUp FALSE
#' @templateVar support FALSE
#' @templateVar weights FALSE
#' @templateVar weighting.method FALSE
#' @template param
#'  
#' @section References:
#'  S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#'  Heterogeneous Tests.
#'  [arXiv:1912.04607v1](https://arxiv.org/abs/1912.04607v1).
#' 
#' @template example
#' @examples
#' 
#' DGR.fast <- DGR(raw.pvalues, pCDFlist)
#' summary(DGR.fast)
#' 
#' DGR.crit <- DGR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DGR.crit)
#' 
#' NDGR.fast <- NDGR(raw.pvalues, pCDFlist)
#' summary(NDGR.fast)
#' 
#' NDGR.crit <- NDGR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(NDGR.crit)
#' 
#' @templateVar Critical.values TRUE
#' @templateVar Adaptive TRUE
#' @templateVar Weighting FALSE
#' @template return
#'
#' @export
discrete.GR <- function(test.results, ...) UseMethod("discrete.GR")

#' @rdname discrete.GR
#' @importFrom checkmate assert_integerish assert_list assert_numeric qassert
#' @export
discrete.GR.default <- function(
    test.results,
    pCDFlist,
    alpha = 0.05,
    zeta = 0.5,
    adaptive = TRUE,
    critical.values = FALSE,
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
  pvec <- match.pvals(pCDFlist, test.results, pCDFlist.indices)
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdx.int(
    pvec             = test.results,
    pCDFlist         = pCDFlist,
    pCDFlist.indices = pCDFlist.indices,
    method           = "GR",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
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

#' @rdname discrete.GR
#' @importFrom checkmate assert_r6 qassert
#' @export
discrete.GR.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    adaptive         = TRUE,
    critical.values  = FALSE,
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
  
  # selection threshold
  qassert(x = select.threshold, rules = "N1(0, 1]")
  
  #----------------------------------------------------
  #       execute computations
  #----------------------------------------------------
  output <- discrete.fdx.int(
    pvec             = test.results$get_pvalues(),
    pCDFlist         = test.results$get_pvalue_supports(unique = TRUE),
    pCDFlist.indices = test.results$get_support_indices(),
    method           = "GR",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    crit.consts      = critical.values,
    threshold        = select.threshold,
    data.name        = deparse(substitute(test.results))
  )
  
  return(output)
}

#' @rdname discrete.GR
#' @export
discrete.GR2 <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE) {
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
  pvec <- match.pvals(pCDFlist, raw.pvalues)
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
    y <- kernel_DGR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals > crit.constants)
  }else{
    # compute transformed sorted p-values
    y <- kernel_DGR_fast(pCDFlist, sorted.pvals, adaptive, alpha)$pval.transf
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
  if(critical.values) output$Critical.values <- crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Guo-Romano procedure"
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

#' @rdname discrete.GR
#' @export
DGR <- function(test.results, ...) UseMethod("DGR")

#' @rdname discrete.GR
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

#' @rdname discrete.GR
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

#' @rdname discrete.GR
#' @export
NDGR <- function(test.results, ...) UseMethod("NDGR")

#' @rdname discrete.GR
#' @export
NDGR.default <- function(
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

#' @rdname discrete.GR
#' @export
NDGR.DiscreteTestResults <- function(
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
    adaptive         = FALSE,
    critical.values  = critical.values,
    select.threshold = select.threshold,
    ...
  )
  
  out$Data$Data.name <- deparse(substitute(test.results))
  
  return(out)
}
