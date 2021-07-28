#'@return
#'A \code{FDX} S3 class object whose elements are:
#'\item{Rejected}{Rejected raw p-values.}
#'\item{Indices}{Indices of rejected hypotheses.}
#'\item{Num.rejected}{Number of rejections.}
#'\item{Adjusted}{Adjusted p-values (only for step-down direction).}
#'<%=ifelse(exists("Weighting") && Weighting, "\\item{Weighted}{Weighted p-values.}","") %>
#'<%=ifelse(exists("Critical.values") && Critical.values, "\\item{Critical.values}{Critical values (if requested).}","") %>
#'\item{Method}{A character string describing the used algorithm, e.g. 'Discrete Lehmann-Romano procedure (step-up)'.}
#'\item{FDP.threshold}{FDP threshold \code{alpha}.}
#'\item{Exceedance.probability}{Probability \code{zeta} of FDP exceeding \code{alpha}; thus, FDP is being controlled at level \code{alpha} with confidence \code{1 - zeta}.}
#'<%=ifelse(exists("Adaptive") && Adaptive, "\\item{Adaptive}{A boolean specifying whether an adaptive procedure was conducted or not.}","") %>
#'<%=ifelse(exists("Weighting") && Weighting, "\\item{Weighting}{A character string describing the weighting method.}","") %>
#'\item{Data$raw.pvalues}{The values of \code{raw.pvalues}.}
#'<%=ifelse(exists("pCDFlist") && pCDFlist, "\\item{Data$pCDFlist}{The values of \\code{pCDFlist}.}","") %>
#'<%=ifelse(exists("weights") && weights, "\\item{Data$weights}{The values of \\code{weights}.}","") %>
#'\item{Data$data.name}{The respective variable names of \code{raw.pvalues} and \code{pCDFlist}.}
