% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DPB_fun.R, R/discretePB_fun.R
\name{DPB.DiscreteTestResults}
\alias{DPB.DiscreteTestResults}
\alias{discrete.PB}
\alias{discrete.PB.default}
\alias{discrete.PB.DiscreteTestResults}
\title{Discrete Poisson-Binomial procedure}
\usage{
\method{DPB}{DiscreteTestResults}(
  test.results,
  alpha = 0.05,
  zeta = 0.5,
  critical.values = FALSE,
  exact = TRUE,
  select.threshold = 1,
  ...
)

discrete.PB(test.results, ...)

\method{discrete.PB}{default}(
  test.results,
  pCDFlist,
  alpha = 0.05,
  zeta = 0.5,
  adaptive = TRUE,
  critical.values = FALSE,
  exact = TRUE,
  select.threshold = 1,
  pCDFlist.indices = NULL,
  ...
)

\method{discrete.PB}{DiscreteTestResults}(
  test.results,
  alpha = 0.05,
  zeta = 0.5,
  adaptive = TRUE,
  critical.values = FALSE,
  exact = TRUE,
  select.threshold = 1,
  ...
)
}
\arguments{
\item{test.results}{either a numeric vector with p-values or an R6 object of class \code{\link[DiscreteTests]{DiscreteTestResults}} from package \code{\link[DiscreteTests:DiscreteTests-package]{DiscreteTests}} for which a discrete FDR procedure is to be performed.}

\item{alpha}{single real number strictly between 0 and 1 specifying the target FDP.}

\item{zeta}{single real number strictly between 0 and 1 specifying the target probability of not exceeding the desired FDP. If \code{zeta = NULL} (the default), then \code{zeta} is chosen equal to \code{alpha}.}

\item{critical.values}{single boolean indicating whether critical constants are to be computed.}

\item{exact}{single boolean indicating whether to compute the Poisson-Binomial distribution exactly or by normal approximation.}

\item{select.threshold}{single real number strictly between 0 and 1 indicating the largest raw p-value to be considered, i.e. only p-values below this threshold are considered and the procedures are adjusted in order to take this selection effect into account; if \code{threshold = 1} (the default), all raw p-values are selected.}

\item{...}{further arguments to be passed to or from other methods. They are ignored here.}

\item{pCDFlist}{list of the supports of the CDFs of the p-values; each list item must be a numeric vector, which is sorted in increasing order and whose last element equals 1.}

\item{adaptive}{single boolean indicating whether to conduct an adaptive procedure or not.}

\item{pCDFlist.indices}{list of numeric vectors containing the test indices that indicate to which raw p-value each \strong{unique} support in \code{pCDFlist} belongs; ignored if the lengths of \code{test.results} and \code{pCDFlist} are equal.}
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

\item{Data$pCDFlist}{list of the \eqn{p}-value supports.}
\item{Data$FDP.threshold}{FDP threshold \code{alpha}.}
\item{Data$Exceedance.probability}{probability \code{zeta} of FDP exceeding \code{alpha}; thus, FDP is being controlled at level \code{alpha} with confidence 1 - \code{zeta}.}

\item{Data$Data.name}{the respective variable name(s) of the input data.}
}
\description{
Apply the [DPB] procedure, with or without computing the critical values,
to a set of p-values and their discrete support. A non-adaptive version is
available as well. Additionally, the user can choose between exact
computation of the Poisson-Binomial distribution or a refined normal
approximation.
}
\details{
\code{DPB} and \code{NDPB} are wrapper functions for \code{discrete.PB}.
The first one simply passes all its arguments to \code{discrete.PB} with
\code{adaptive = TRUE} and \code{NDPB} does the same with
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

# DPB (exact) without critical values; using results object
DPB.exact.fast <- discrete.PB(test.results)
summary(DPB.exact.fast)

# DPB (exact) with critical values; using extracted p-values and supports
DPB.exact.crit <- discrete.PB(raw.pvalues, pCDFlist, critical.values = TRUE)
summary(DPB.exact.crit)

# DPB (normal approximation) without critical values; using extracted p-values and supports
DPB.norm.fast <- discrete.PB(raw.pvalues, pCDFlist, exact = FALSE)
summary(DPB.norm.fast)

# DPB (normal approximation) with critical values; using results object
DPB.norm.crit <- discrete.PB(test.results, critical.values = TRUE,
                             exact = FALSE)
summary(DPB.norm.crit)

# Non-adaptive DPB (exact) without critical values; using results object
NDPB.exact.fast <- discrete.PB(test.results, adaptive = FALSE)
summary(NDPB.exact.fast)

# Non-adaptive DPB (exact) with critical values; using extracted p-values and supports
NDPB.exact.crit <- discrete.PB(raw.pvalues, pCDFlist, adaptive = FALSE,
                               critical.values = TRUE)
summary(NDPB.exact.crit)

# Non-adaptive DPB (normal approx.) without critical values; using extracted p-values and supports
NDPB.norm.fast <- discrete.PB(raw.pvalues, pCDFlist, adaptive = FALSE,
                              exact = FALSE)
summary(NDPB.norm.fast)

# Non-adaptive DPB (normal approx.) with critical values; using results object
NDPB.norm.crit <- discrete.PB(test.results, adaptive = FALSE,
                              critical.values = TRUE, exact = FALSE)
summary(NDPB.norm.crit)

}
\references{
Döhler, S. & Roquain, E. (2020). Controlling False Discovery Exceedance for
Heterogeneous Tests. \emph{Electronic Journal of Statistics}, \emph{14}(2),
pp. 4244-4272. \doi{10.1214/20-EJS1771}
}
\seealso{
\code{\link{kernel}}, \code{\link[=FDX-package]{FDX}}, \code{\link[=continuous.LR]{continuous.LR()}},
\code{\link[=continuous.GR]{continuous.GR()}}, \code{\link[=discrete.LR]{discrete.LR()}},
\code{\link[=discrete.GR]{discrete.GR()}}, \code{\link[=weighted.LR]{weighted.LR()}},
\code{\link[=weighted.GR]{weighted.GR()}}, \code{\link[=weighted.PB]{weighted.PB()}}
}
