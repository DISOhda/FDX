#' @name fast.Discrete
#' 
#' @title
#' Fast application of discrete procedures
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' Applies the \[DLR\], \[DGR\] or \[DPB\] procedures, **without** computing the
#' critical values, to a data set of 2 x 2 contingency tables using Fisher's
#' exact test.
#' 
#' **Note**: These functions are deprecated and will be removed in a future
#' version. Please use [`direct.discrete.*()`][direct.discrete] with
#' `test.fun = DiscreteTests::fisher.test.pv` and (optional) 
#' `preprocess.fun = DiscreteDatasets::reconstruct_two` or 
#' `preprocess.fun = DiscreteDatasets::reconstruct_four` instead. Alternatively,
#' use a pipeline like\cr
#' `data |>`\cr
#' `  DiscreteDatasets::reconstruct_*(<args>) |>`\cr
#' `  DiscreteTests::*.test.pv(<args>) |>`\cr
#' `  discrete.*(<args>)`.
#'
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
#' @templateVar exact TRUE
#' @templateVar weights FALSE
#' @template param
#' 
#' @param counts        a data frame of 2 or 4 columns and any number of lines,
#'                      each line representing a 2 x 2 contingency table to
#'                      test. The number of columns and what they must contain
#'                      depend on the value of the `input` argument, see
#'                      Details of [`DiscreteFDR::fisher.pvalues.support()`].
#' @param alternative   same argument as in [`fisher.test()`]. The three
#'                      possible values are `"greater"` (default),
#'                      `"two.sided"` or `"less"`; may be abbreviated.
#' @param input         the format of the input data frame, see Details of
#'                      [`DiscreteFDR::fisher.pvalues.support()`]. The
#'                      three possible values are `"noassoc"` (default),
#'                      `"marginal"` or `"HG2011"`; may be 
#'                      abbreviated.
#' 
#' @examples
#' 
#' X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#' X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#' N1 <- rep(148, 9)
#' N2 <- rep(132, 9)
#' Y1 <- N1 - X1
#' Y2 <- N2 - X2
#' df <- data.frame(X1, Y1, X2, Y2)
#' df
#' 
#' # DLR
#' DLR.sd <- fast.Discrete.LR(counts = df, input = "noassoc")
#' summary(DLR.sd)
#' 
#' # DLR
#' DLR.su <- fast.Discrete.LR(counts = df, input = "noassoc", direction = "su")
#' summary(DLR.su)
#' 
#' # Non-adaptive DLR
#' NDLR.sd <- fast.Discrete.LR(counts = df, input = "noassoc", adaptive = FALSE)
#' summary(NDLR.sd)
#' 
#' # Non-adaptive DLR
#' NDLR.su <- fast.Discrete.LR(counts = df, input = "noassoc", direction = "su", adaptive = FALSE)
#' summary(NDLR.su)
#' 
#' # DGR
#' DGR <- fast.Discrete.GR(counts = df, input = "noassoc")
#' summary(DGR)
#' 
#' # Non-adaptive DGR
#' NDGR <- fast.Discrete.GR(counts = df, input = "noassoc", adaptive = FALSE)
#' summary(NDGR)
#' 
#' # DPB
#' DPB <- fast.Discrete.PB(counts = df, input = "noassoc")
#' summary(DPB)
#' 
#' # Non-adaptive DPB
#' NDPB <- fast.Discrete.PB(counts = df, input = "noassoc", adaptive = FALSE)
#' summary(NDPB)
#' 
#' @templateVar Critical.values FALSE
#' @templateVar Adaptive TRUE
#' @templateVar Weighting FALSE
#' @template return
#' 
#' @importFrom DiscreteFDR fisher.pvalues.support
#' @importFrom lifecycle deprecate_soft
#' @export
fast.Discrete.LR <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, direction = "sd", adaptive = TRUE){
  deprecate_soft("2.0.0", "fast.Discrete.LR()", "direct.discrete.LR()")
  
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.LR(raw.pvalues, pCDFlist, alpha, zeta, direction, adaptive, FALSE)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}

#' @rdname fast.Discrete
#' @export
fast.Discrete.GR <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, adaptive = TRUE){
  deprecate_soft("2.0.0", "fast.Discrete.GR()", "direct.discrete.GR()")
  
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.GR(raw.pvalues, pCDFlist, alpha, zeta, adaptive, FALSE)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}

#' @rdname fast.Discrete
#' @export
fast.Discrete.PB <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, adaptive = TRUE, exact = FALSE){
  deprecate_soft("2.0.0", "fast.Discrete.PB()", "direct.discrete.PB()")
  
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.PB(raw.pvalues, pCDFlist, alpha, zeta, adaptive, FALSE, exact)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}

#' @name direct.discrete
#' 
#' @title 
#' Direct Application of Multiple Testing Procedures to Dataset
#' 
#' @description
#' Apply the \[DLR\], \[NDLR\], \[DGR\], \[NDGR\], \[PB\] or \[NPB\] procedure,
#' with or without computing the critical constants, to a data set of 2x2
#' contingency tables using a hypothesis test function from package
#' [DiscreteTests][DiscreteTests::DiscreteTests-package].
#' 
#' @templateVar dat TRUE
#' @templateVar test.fun TRUE
#' @templateVar test.args TRUE
#' @templateVar alpha TRUE
#' @templateVar zeta TRUE
#' @templateVar direction TRUE
#' @templateVar adaptive TRUE
#' @templateVar critical.values TRUE
#' @templateVar exact TRUE
#' @templateVar select.threshold TRUE
#' @templateVar preprocess.fun TRUE
#' @templateVar preprocess.args TRUE
#' @templateVar weights FALSE
#' @template param
#' 
#' @param dat               input data; must be suitable for the first parameter
#'                          of the provided `preprocess.fun` function or, if
#'                          `preprocess.fun = NULL`, for the first parameter of
#'                          the `test.fun` function.
#' @param test.fun          function **from package
#'                          [`DiscreteTests`][DiscreteTests::DiscreteTests-package]**,
#'                          i.e. one whose name ends with `*_test_pv` and which
#'                          performs hypothesis tests and provides an object
#'                          with p-values and their support sets; can be
#'                          specified by a single character string (which is
#'                          automatically checked for being a suitable function
#'                          **from that package** and may be abbreviated) or a
#'                          single function object.
#' @param test.args         optional named list with arguments for `test.fun`;
#'                          the names of the list fields must match the test
#'                          function's parameter names. The first parameter of
#'                          the test function (i.e. the data) **MUST NOT** be
#'                          included!
#' @param preprocess.fun    optional function for pre-processing the input
#'                          `data`; its result must be suitable for the first
#'                          parameter of the `test.fun` function.
#' @param preprocess.args   optional named list with arguments for
#'                          `preprocess.fun`; the names of the list fields must
#'                          match the pre-processing function's parameter names.
#'                          The first parameter of the test function (i.e. the
#'                          data) **MUST NOT** be included!
#' 
#' @template example
#' @examples
#' 
#' # DLR
#' DLR.sd <- direct.discrete.LR(df, "fisher")
#' summary(DLR.sd)
#' 
#' # Non-adaptive DLR (step-up variant; adjusted p-values do not exist here!)
#' NDLR.su <- direct.discrete.LR(df, "fisher", direction = "su", adaptive = FALSE)
#' summary(NDLR.su)
#' 
#' # DGR
#' DGR <- direct.discrete.GR(df, "fisher")
#' summary(DGR)
#' 
#' # Non-adaptive DGR
#' NDGR <- direct.discrete.GR(df, "fisher", adaptive = FALSE)
#' summary(NDGR)
#' 
#' # DPB (normal approximation)
#' PB.approx <- direct.discrete.PB(df, "fisher", exact = FALSE)
#' summary(DGR)
#' 
#' # Non-adaptive DPB
#' NPB.exact <- direct.discrete.GR(df, "fisher", adaptive = FALSE)
#' summary(NDGR)
#' 
#' @importFrom DiscreteFDR generate.pvalues
#' @export
direct.discrete.LR <- function(
    dat,
    test.fun, 
    test.args = NULL,
    alpha = 0.05,
    zeta = 0.5,
    direction = "su",
    adaptive = FALSE,
    critical.values = FALSE,
    select.threshold = 1,
    preprocess.fun = NULL, 
    preprocess.args = NULL
) {
  out <- discrete.LR.DiscreteTestResults(
    test.results = generate.pvalues(
      dat             = dat,
      test.fun        = test.fun,
      test.args       = test.args,
      preprocess.fun  = preprocess.fun,
      preprocess.args = preprocess.args
    ),
    alpha            = alpha,
    zeta             = zeta,
    direction        = direction,
    adaptive         = adaptive,
    critical.values  = critical.values,
    select.threshold = select.threshold
  )
  
  out$Data$Data.name <- deparse(substitute(dat))
  
  return(out)
}

#' @name direct.discrete
#' @importFrom DiscreteFDR generate.pvalues
#' @export
direct.discrete.GR <- function(
    dat,
    test.fun, 
    test.args = NULL,
    alpha = 0.05,
    zeta = 0.5,
    adaptive = FALSE,
    critical.values = FALSE,
    select.threshold = 1,
    preprocess.fun = NULL, 
    preprocess.args = NULL
) {
  out <- discrete.GR.DiscreteTestResults(
    test.results = generate.pvalues(
      dat             = dat,
      test.fun        = test.fun,
      test.args       = test.args,
      preprocess.fun  = preprocess.fun,
      preprocess.args = preprocess.args
    ),
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    critical.values  = critical.values,
    select.threshold = select.threshold
  )
  
  out$Data$Data.name <- deparse(substitute(dat))
  
  return(out)
}

#' @name direct.discrete
#' @importFrom DiscreteFDR generate.pvalues
#' @export
direct.discrete.PB <- function(
    dat,
    test.fun, 
    test.args        = NULL,
    alpha            = 0.05,
    zeta             = 0.5,
    adaptive         = FALSE,
    critical.values  = FALSE,
    exact            = TRUE,
    select.threshold = 1,
    preprocess.fun   = NULL, 
    preprocess.args  = NULL
) {
  out <- discrete.PB.DiscreteTestResults(
    test.results = generate.pvalues(
      dat             = dat,
      test.fun        = test.fun,
      test.args       = test.args,
      preprocess.fun  = preprocess.fun,
      preprocess.args = preprocess.args
    ),
    alpha            = alpha,
    zeta             = zeta,
    adaptive         = adaptive,
    critical.values  = critical.values,
    exact            = exact,
    select.threshold = select.threshold
  )
  
  out$Data$Data.name <- deparse(substitute(dat))
  
  return(out)
}