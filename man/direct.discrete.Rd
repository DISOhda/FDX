% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_discrete.R
\name{direct.discrete}
\alias{direct.discrete}
\alias{direct.discrete.LR}
\alias{direct.discrete.GR}
\alias{direct.discrete.PB}
\title{Direct Application of Multiple Testing Procedures to Dataset}
\usage{
direct.discrete.LR(
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
)

direct.discrete.GR(
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
)

direct.discrete.PB(
  dat,
  test.fun,
  test.args = NULL,
  alpha = 0.05,
  zeta = 0.5,
  adaptive = FALSE,
  critical.values = FALSE,
  exact = TRUE,
  select.threshold = 1,
  preprocess.fun = NULL,
  preprocess.args = NULL
)
}
\arguments{
\item{dat}{input data; must be suitable for the first parameter
of the provided \code{preprocess.fun} function or, if
\code{preprocess.fun = NULL}, for the first parameter of
the \code{test.fun} function.}

\item{test.fun}{function \strong{from package
\code{\link[DiscreteTests:DiscreteTests-package]{DiscreteTests}}},
i.e. one whose name ends with \verb{*_test_pv} and which
performs hypothesis tests and provides an object
with p-values and their support sets; can be
specified by a single character string (which is
automatically checked for being a suitable function
\strong{from that package} and may be abbreviated) or a
single function object.}

\item{test.args}{optional named list with arguments for \code{test.fun};
the names of the list fields must match the test
function's parameter names. The first parameter of
the test function (i.e. the data) \strong{MUST NOT} be
included!}

\item{alpha}{single real number strictly between 0 and 1 specifying the target FDP.}

\item{zeta}{single real number strictly between 0 and 1 specifying the target probability of not exceeding the desired FDP. If \code{zeta = NULL} (the default), then \code{zeta} is chosen equal to \code{alpha}.}

\item{direction}{single character string specifying whether to perform the step-up (\verb{"su}) or step-down (\code{"sd"}; the default) version of the Lehmann-Romano procedure.}

\item{adaptive}{single boolean indicating whether to conduct an adaptive procedure or not.}

\item{critical.values}{single boolean indicating whether critical constants are to be computed.}

\item{select.threshold}{single real number strictly between 0 and 1 indicating the largest raw p-value to be considered, i.e. only p-values below this threshold are considered and the procedures are adjusted in order to take this selection effect into account; if \code{threshold = 1} (the default), all raw p-values are selected.}

\item{preprocess.fun}{optional function for pre-processing the input
\code{data}; its result must be suitable for the first
parameter of the \code{test.fun} function.}

\item{preprocess.args}{optional named list with arguments for
\code{preprocess.fun}; the names of the list fields must
match the pre-processing function's parameter names.
The first parameter of the test function (i.e. the
data) \strong{MUST NOT} be included!}

\item{exact}{single boolean indicating whether to compute the Poisson-Binomial distribution exactly or by normal approximation.}
}
\description{
Apply the [DLR], [NDLR], [DGR], [NDGR], [PB] or [NPB] procedure,
with or without computing the critical constants, to a data set of 2x2
contingency tables using a hypothesis test function from package
\link[DiscreteTests:DiscreteTests-package]{DiscreteTests}.
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

# DLR
DLR.sd <- direct.discrete.LR(df, "fisher")
summary(DLR.sd)

# Non-adaptive DLR (step-up variant; adjusted p-values do not exist here!)
NDLR.su <- direct.discrete.LR(df, "fisher", direction = "su", adaptive = FALSE)
summary(NDLR.su)

# DGR
DGR <- direct.discrete.GR(df, "fisher")
summary(DGR)

# Non-adaptive DGR
NDGR <- direct.discrete.GR(df, "fisher", adaptive = FALSE)
summary(NDGR)

# DPB (normal approximation)
PB.approx <- direct.discrete.PB(df, "fisher", exact = FALSE)
summary(DGR)

# Non-adaptive DPB
NPB.exact <- direct.discrete.GR(df, "fisher", adaptive = FALSE)
summary(NDGR)

}
