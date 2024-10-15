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
#' @templateVar test.results TRUE
#' @templateVar pCDFlist TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
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
#' [`FDX`][FDX-package], [`discrete.GR()`], [`discrete.PB()`], 
#' [`continuous.LR()`], [`continuous.GR()`], [`weighted.LR()`],
#' [`weighted.GR()`], [`weighted.PB()`]
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
#' DLR.sd.fast <- discrete.LR(test.results)
#' summary(DLR.sd.fast)
#' 
#' # DLR with critical values; using extracted p-values and supports
#' DLR.sd.crit <- discrete.LR(raw.pvalues, pCDFlist, critical.values = TRUE)
#' summary(DLR.sd.crit)
#' 
#' # DLR (step-up) without critical values; using extracted p-values and supports
#' DLR.su.fast <- discrete.LR(raw.pvalues, pCDFlist, direction = "su")
#' summary(DLR.su.fast)
#' 
#' # DLR (step-up) with critical values; using results object
#' DLR.su.crit <- discrete.LR(test.results, direction = "su",
#'                            critical.values = TRUE)
#' summary(DLR.su.crit)
#' 
#' # Non-adaptive DLR without critical values; using results object
#' NDLR.sd.fast <- discrete.LR(test.results, adaptive = FALSE)
#' summary(NDLR.sd.fast)
#' 
#' # Non-adaptive DLR with critical values; using extracted p-values and supports
#' NDLR.sd.crit <- discrete.LR(raw.pvalues, pCDFlist, adaptive = FALSE,
#'                             critical.values = TRUE)
#' summary(NDLR.sd.crit)
#' 
#' # Non-adaptive DLR (step-up) without critical values; using extracted p-values and supports
#' NDLR.su.fast <- discrete.LR(raw.pvalues, pCDFlist, direction = "su",
#'                             adaptive = FALSE)
#' summary(NDLR.su.fast)
#' 
#' # Non-adaptive DLR (step-up) with critical values; using results object
#' NDLR.su.crit <- discrete.LR(test.results, direction = "su",
#'                             adaptive = FALSE, critical.values = TRUE)
#' summary(NDLR.su.crit)
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
  pvec <- match.pvals(test.results, pCDFlist, pCDFlist.indices)
  
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
