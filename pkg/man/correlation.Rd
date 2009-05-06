\name{correlation}
\alias{cor.rlogit}
\alias{cov.rlogit}
\alias{sd.rlogit}
\title{Correlation structure of the random parameters}
\description{
  Functions that extract the correlation structure of a rlogit object
}
\usage{
cor.rlogit(x)
cov.rlogit(x)
sd.rlogit(x)
}

\arguments{
  \item{x}{an \code{rlogit} object with \code{correlation=TRUE}.}
}
\value{
  A matrix for \code{cor.rlogit} and \code{cov.rlogit}, a vector for \code{sd.rlogit}
}
\details{
  These functions extract the standard deviations (\code{sd.rlogit}), the
  covariance (\code{cov.rlogit}) or the correlations (\code{cor.rlogit}).
}

\author{Yves Croissant}

\keyword{regression}
