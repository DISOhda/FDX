#' @name FDX-package
#' 
#' @title
#' False Discovery Exceedance (FDX) Control for Heterogeneous and Discrete Tests
#' 
#' @import Rcpp
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib FDX
#' @description
#' This package implements the \[HLR\], \[HGR\] and \[HPB\] procedures for both
#' heterogeneous and discrete tests (see Reference). 
#' 
#' @details
#' The functions are reorganized from the reference paper in the following way.
#' [`discrete.LR()`] (for Discrete Lehmann-Romano) implements \[DLR\],
#' [`discrete.GR()`] (for Discrete Guo-Romano) implements \[DGR\] and
#' [`discrete.PB()`] (for Discrete Poisson-Binomial) implements \[DPB\].
#' [`DLR()`] and [`NDLR()`] are wrappers for [`discrete.LR()`] to access
#' \[DLR\] and its non-adaptive version directly. Likewise, [`DGR()`],
#' [`NDGR()`], [`DPB()`] and [`NDPB()`] are wrappers for
#' [`discrete.GR()`] and [`discrete.PB()`], respectively. Their main
#' parameters are a vector of raw observed p-values and a list of the same
#' length, whose elements are the discrete supports of the CDFs of the p-values.
#' 
#' In the same fashion, [`weighted.LR()`] (for Weighted Lehmann-Romano),
#' [`weighted.GR()`] (for Weighted Guo-Romano) and [`weighted.PB()`]
#' (for Weighted Poisson-Binomial) implement \[wLR\], \[wGR\] and \[wGR\],
#' respectively. They also possess wrapper functions, namely [`wLR.AM()`],
#' [`wGR.AM()`] and [`wPB.AM()`] for arithmetic weighting, and [`wLR.GM()`],
#' [`wPB.GM()`] and [`wPB.GM()`] for geometric weighting.
#'  
#' The functions [`fast.Discrete.LR()`], [`fast.Discrete.GR()`]
#' and [`fast.Discrete.PB()`] are wrappers for
#' [`DiscreteFDR::fisher.pvalues.support()`] and [`discrete.LR()`],
#' [`discrete.GR()`] and [`discrete.PB()`], respectively, which allow to apply
#' discrete procedures directly to a data set of contingency tables.
#'  
#' @section References:
#'  S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#'  Heterogeneous Tests.
#'  [arXiv:1912.04607v1](https://arxiv.org/abs/1912.04607v1).
"_PACKAGE"
