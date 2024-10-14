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
#' @references
#' Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
#'   Heterogeneous Tests. *Electronic Journal of Statistics*, *14*(2),
#'   pp. 4244-4272. \doi{10.1214/20-EJS1771}
#'   
#' Lehmann, E. L. & Romano, J. P. (2005). Generalizations of the familywise
#'   error rate. *The Annals of Statistics*, *33*(3), pp. 1138-1154.
#'   \doi{10.1214/009053605000000084}
#'   
#' Guo, W. & Romano, J. P. (2007). A generalized Sidak-Holm procedure and
#'   control of generalized error rates under independence.
#'   *Statistical Applications in Genetics and Molecular Biology*, *6*(1),
#'   Art. 3, 35 pp. (eletronic). \doi{10.2202/1544-6115.1247}
"_PACKAGE"
