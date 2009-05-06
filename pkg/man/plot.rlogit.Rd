\name{plot.rlogit}
\alias{plot.rlogit}
\alias{plot.rpar}
\title{Plot of the distribution of estimated random parameters}
\description{
  Methods for \code{rpar} and \code{rlogit} objects which provide a
  plot of the distribution of one or all of the estimated random parameters
}
\usage{
\method{plot}{rlogit}(x, ...)
\method{plot}{rpar}(x, type = c("density","probability"), norm = TRUE, ...)
}
\arguments{
  \item{x}{a \code{rlogit} or a \code{rpar} object,}
  \item{type}{the function to be plotted, whether the density or the
    probability density function,}
  \item{norm}{if true, the random parameter is divided by the
    \code{norm} variable if exists.}
  \item{...}{further arguments, passed to \code{plot.rpar} for the
    \code{rlogit} method and to \code{plot} for the \code{rpar} method.}
}
\details{
  For the \code{rpar} method, one plot is drawn. For the \code{rlogit}
  method, one plot for each random parameter is drawn.
}

\seealso{
  \code{\link{rlogit}} for the estimation of random parameters logit models
  and \code{\link{rpar}} for the description of \code{rpar} objects and
  \code{\link{distribution}} for functions which return
  informations about the distribution of random parameters.
}

\author{Yves Croissant}
\keyword{regression}