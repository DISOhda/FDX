% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/continuousGR.R
\name{continuous.GR}
\alias{continuous.GR}
\alias{GR}
\alias{NGR}
\title{Continuous Guo-Romano procedure}
\usage{
continuous.GR(
  test.results,
  alpha = 0.05,
  zeta = 0.5,
  adaptive = TRUE,
  critical.values = FALSE,
  select.threshold = 1
)

GR(
  test.results,
  alpha = 0.05,
  zeta = 0.5,
  critical.values = FALSE,
  select.threshold = 1
)

NGR(
  test.results,
  alpha = 0.05,
  zeta = 0.5,
  critical.values = FALSE,
  select.threshold = 1
)
}
\arguments{
\item{test.results}{either a numeric vector with p-values or an R6 object of class \code{\link[DiscreteTests]{DiscreteTestResults}} from package \code{\link[DiscreteTests:DiscreteTests-package]{DiscreteTests}} for which a discrete FDR procedure is to be performed.}

\item{alpha}{single real number strictly between 0 and 1 specifying the target FDP.}

\item{zeta}{single real number strictly between 0 and 1 specifying the target probability of not exceeding the desired FDP. If \code{zeta = NULL} (the default), then \code{zeta} is chosen equal to \code{alpha}.}

\item{adaptive}{single boolean indicating whether to conduct an adaptive procedure or not.}

\item{critical.values}{single boolean indicating whether critical constants are to be computed.}

\item{select.threshold}{single real number strictly between 0 and 1 indicating the largest raw p-value to be considered, i.e. only p-values below this threshold are considered and the procedures are adjusted in order to take this selection effect into account; if \code{threshold = 1} (the default), all raw p-values are selected.}
}
\value{
A \code{FDX} S3 class object whose elements are:
\item{Rejected}{rejected raw \eqn{p}-values.}
\item{Indices}{indices of rejected \eqn{p}-values.}
\item{Num.rejected}{number of rejections.}

\item{Adjusted}{adjusted \eqn{p}-values.}
\item{Critical.values}{critical values (only exists if computations where performed with \code{critical.values = TRUE}).}
\item{Select}{list with data related to \eqn{p}-value selection; only exists if \code{select.threshold < 1}.}
\item{Select$Threshold}{\eqn{p}-value selection threshold.}
\item{Select$Effective.Thresholds}{results of each \eqn{p}-value CDF evaluated at the selection threshold.}
\item{Select$Pvalues}{selected \eqn{p}-values that are \eqn{\leq} selection threshold.}
\item{Select$Indices}{indices of \eqn{p}-values \eqn{\leq} selection threshold.}
\item{Select$Scaled}{scaled selected \eqn{p}-values.}
\item{Select$Number}{number of selected \eqn{p}-values \eqn{\leq} selection threshold.}
\item{Data}{list with input data.}
\item{Data$Method}{character string describing the used algorithm, e.g. 'Discrete Lehmann-Romano procedure (step-up)'.}
\item{Data$Raw.pvalues}{all observed raw \eqn{p}-values.}

\item{Data$FDP.threshold}{FDP threshold \code{alpha}.}
\item{Data$Exceedance.probability}{probability \code{zeta} of FDP exceeding \code{alpha}; thus, FDP is being controlled at level \code{alpha} with confidence 1 - \code{zeta}.}
\item{Data$Adaptive}{boolean indicating whether an adaptive procedure was conducted or not.}

\item{Data$Data.name}{the respective variable name(s) of the input data.}
}
\description{
Apply the usual continuous [GR] procedure, with or without computing the
critical values, to a set of p-values. A non-adaptive version is available as
well.
}
\details{
\code{GR} and \code{NGR} are wrapper functions for \code{continuous.GR}. The
first one simply passes all its arguments to \code{continuous.GR} with
\code{adaptive = TRUE} and \code{NGR} does the same with
\code{adaptive = FALSE}.
}
\examples{
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
df

# Construction of the p-values and their supports with Fisher's exact test
library(DiscreteTests)  # for Fisher's exact test
test.results <- fisher_test_pv(df)
raw.pvalues <- test.results$get_pvalues()
pCDFlist <- test.results$get_pvalue_supports()

# GR without critical values; using extracted p-values
GR.fast <- GR(raw.pvalues)
summary(GR.fast)

# LR with critical values; using test results object
GR.crit <- GR(test.results, critical.values = TRUE)
summary(GR.crit)

# Non-adaptive GR without critical values; using test results object
NGR.fast <- NGR(test.results)
summary(NGR.fast)

# Non-adaptive GR with critical values; using extracted p-values
NGR.crit <- NGR(raw.pvalues, critical.values = TRUE)
summary(NGR.crit)

}
\references{
Guo, W. & Romano, J. P. (2007). A generalized Sidak-Holm procedure and
control of generalized error rates under independence.
\emph{Statistical Applications in Genetics and Molecular Biology}, \emph{6}(1),
Art. 3, 35 pp. (electronic). \doi{10.2202/1544-6115.1247}
}
\seealso{
\code{\link{kernel}}, \code{\link{FDX-package}}, \code{\link[=continuous.LR]{continuous.LR()}},
\code{\link[=discrete.LR]{discrete.LR()}}, \code{\link[=discrete.GR]{discrete.GR()}},
\code{\link[=discrete.PB]{discrete.PB()}}, \code{\link[=weighted.LR]{weighted.LR()}},
\code{\link[=weighted.GR]{weighted.GR()}}, \code{\link[=weighted.PB]{weighted.PB()}}
}
