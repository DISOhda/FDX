#' @name discrete.LR
#' 
#' @title
#' Discrete Lehmann-Romano procedure
#' 
#' @description
#' Apply the \[DLR\] procedure, with or without computing the critical values,
#' to a set of p-values and their discrete support. Both step-down and step-up
#' procedures can be computed and non-adaptive versions are available as well.
#' 
#' @details
#' `DLR` and `NDLR` are wrapper functions for `discrete.LR`.
#' The first one simply passes all its arguments to `discrete.LR` with
#' `adaptive = TRUE` and `NDLR` does the same with
#' `adaptive = FALSE`.
#' 
#' @seealso
#' [`kernel`], [`FDX`][FDX-package], [`continuous.LR()`],
#' [`continuous.GR()`], [`discrete.GR()`], 
#' [`discrete.PB()`], [`weighted.LR()`],
#' [`weighted.GR()`], [`weighted.PB()`]
#' 
#' @templateVar raw.pvalues TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction TRUE
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
#' DLR.sd.fast <- DLR(raw.pvalues, pCDFlist)
#' summary(DLR.sd.fast)
#' DLR.su.fast <- DLR(raw.pvalues, pCDFlist, direction = "su")
#' summary(DLR.su.fast)
#' 
#' DLR.sd.crit <- DLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DLR.sd.crit)
#' DLR.su.crit <- DLR(raw.pvalues, pCDFlist, direction = "su", critical.values = TRUE)
#' summary(DLR.su.crit)
#' 
#' NDLR.sd.fast <- NDLR(raw.pvalues, pCDFlist)
#' summary(NDLR.sd.fast)
#' NDLR.su.fast <- NDLR(raw.pvalues, pCDFlist, direction = "su")
#' summary(NDLR.su.fast)
#' 
#' NDLR.sd.crit <- NDLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(NDLR.sd.crit)
#' NDLR.su.crit <- NDLR(raw.pvalues, pCDFlist, direction = "su", critical.values = TRUE)
#' summary(NDLR.su.crit)
#' 
#' @templateVar Critical.values TRUE
#' @template return
#' 
#' @export
discrete.LR <- function(test.results, ...) UseMethod("discrete.LR")

#' @rdname discrete.LR
#' @importFrom checkmate assert_integerish assert_list assert_numeric qassert
#' @export
discrete.LR.default <- function(
    test.results,
    pCDFlist,
    alpha            = 0.05,
    zeta             = 0.5,
    direction        = "sd",
    adaptive         = TRUE,
    critical.values  = FALSE,
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
  
  # step-up/step-down direction
  direction <- match.arg(tolower(direction), c("su", "sd"))
  
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
    method           = "LR",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    method.parameter = (direction == "su"),
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

#' @rdname discrete.LR
#' @importFrom checkmate assert_r6 qassert
#' @export
discrete.LR.DiscreteTestResults <- function(
    test.results,
    alpha            = 0.05,
    zeta             = 0.5,
    direction        = "sd",
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
  
  # step-up/step-down direction
  direction <- match.arg(tolower(direction), c("su", "sd"))
  
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
    method           = "LR",
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    method.parameter = (direction == "su"),
    crit.consts      = critical.values,
    threshold        = select.threshold,
    data.name        = deparse(substitute(test.results))
  )
  
  return(output)
}

#' @rdname discrete.LR
#' @export
discrete.LR2 <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, direction = "sd", adaptive = TRUE, critical.values = FALSE){
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
  direction <- match.arg(direction, c("su", "sd"))
  if(direction == "su"){
    # SU case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, TRUE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      # compute transformed observed p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, TRUE, zeta, pv.list.all)
      # determine significant (transformed) p-values
      if(length(y)){
        idx <- which(y <= zeta * (floor(1:length(y) * alpha) + 1))
      }else{
        idx <- integer(0)
      }
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      m.rej <- 0
      idx <- integer(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    # SD case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, FALSE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      # compute transformed sorted p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, FALSE)
      # determine significant (transformed) p-values
      idx <- which(y > zeta * (floor(1:m * alpha) + 1)) 
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
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  if(direction == "sd"){
    if(critical.values){
      y <- y$pval.transf
    }
    # compute adjusted p-values
    pv.adj = cummax(pmin(y / (floor(alpha * 1:m) + 1), 1))
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted = pv.adj[ro]
  }
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Lehmann-Romano procedure"
  alg <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$Method <- paste(alg, switch(direction, su = "(step-up)", sd = "(step-down)"))
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

#' @rdname discrete.LR
#' @export
DLR <- function(test.results, ...) UseMethod("DLR")

#' @rdname discrete.LR
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

#' @rdname discrete.LR
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

#' @rdname discrete.LR
#' @export
NDLR <- function(test.results, ...) UseMethod("NDLR")

#' @rdname discrete.LR
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

#' @rdname discrete.LR
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
