#' @return
#' A \code{FDX} S3 class object whose elements are:
#' \item{Rejected}{rejected raw \eqn{p}-values.}
#' \item{Indices}{indices of rejected \eqn{p}-values.}
#' \item{Num.rejected}{number of rejections.}
#' <%=ifelse(exists("Weighting") && Weighting, "\\item{Weighted}{weighted \\eqn{p}-values.}","") %>
#' \item{Adjusted}{adjusted \eqn{p}-values<%=ifelse(exists("direction") && direction, " (only for step-down direction)", "")%>.}
#' \item{Critical.values}{critical values (only exists if computations where performed with `critical.values = TRUE`).}
#' \item{Select}{list with data related to \eqn{p}-value selection; only exists if `threshold < 1`.}
#' \item{Select$Threshold}{\eqn{p}-value selection `threshold`.}
#' \item{Select$Effective.Thresholds}{results of each \eqn{p}-value CDF evaluated at the selection threshold.}
#' \item{Select$Pvalues}{selected \eqn{p}-values that are \eqn{\leq} selection `threshold`.}
#' \item{Select$Indices}{indices of \eqn{p}-values \eqn{\leq} selection `threshold`.}
#' \item{Select$Scaled}{scaled selected \eqn{p}-values.}
#' \item{Select$Number}{number of selected \eqn{p}-values \eqn{\leq} `threshold`.}
#' \item{Data}{list with input data.}
#' \item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Lehmann-Romano procedure (step-up)'.}
#' \item{Data$Raw.pvalues}{all observed raw \eqn{p}-values.}
#' <%=ifelse(exists("weights") && weights, "\\item{Data$Weights}{the weights for the raw \\eqn{p}-values.}","") %>
#' <%=ifelse(exists("pCDFlist") && pCDFlist, "\\item{Data$pCDFlist}{list of the \\eqn{p}-value supports.}","") %>
#' \item{Data$FDP.threshold}{FDP threshold \code{alpha}.}
#' \item{Data$Exceedance.probability}{probability `zeta` of FDP exceeding `alpha`; thus, FDP is being controlled at level `alpha` with confidence 1 - `zeta`.}
#' <%=ifelse(exists("Adaptive") && Adaptive, "\\item{Data$Adaptive}{boolean indicating whether an adaptive procedure was conducted or not.}","") %>
#' <%=ifelse(exists("Weighting") && Weighting, "\\item{Data$Weighting}{character string describing the weighting method.}","") %>
#' \item{Data$Data.name}{the respective variable name(s) of the input data.}
