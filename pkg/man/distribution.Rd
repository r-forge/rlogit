\name{distribution}
\alias{distribution}
\alias{med}
\alias{stdev}
\alias{mean.rpar}
\alias{med.rpar}
\alias{stdev.rpar}
\alias{qrlogit}
\alias{prlogit}
\alias{drlogit}
\alias{qrlogit.rpar}
\alias{prlogit.rpar}
\alias{drlogit.rpar}
\alias{qrlogit.rlogit}
\alias{prlogit.rlogit}
\alias{drlogit.rlogit}



\title{Functions used to describe the characteristics of estimated
  random parameters}
\description{
  \code{rpar} objects contain all the relevant information about
  the distribution of random parameters. These functions enables to
  obtain easily descriptive statistics, density, probability and
  quantiles of the distribution.
}
\usage{
med(x, ...)
stdev(x, ...)
qrlogit(x, ...)
prlogit(x, ...)
drlogit(x, ...)
\method{mean}{rpar}(x, norm = TRUE, ...)
\method{med}{rpar}(x, norm = TRUE, ...)
\method{stdev}{rpar}(x, norm = TRUE, ...)
\method{qrlogit}{rpar}(x, norm = TRUE, ...)
\method{prlogit}{rpar}(x, norm = TRUE, ...)
\method{drlogit}{rpar}(x, norm = TRUE, ...)
\method{qrlogit}{rlogit}(x, par = 1, y = NULL, norm = TRUE, ...)
\method{prlogit}{rlogit}(x, par = 1, y = NULL, norm = TRUE, ...)
\method{drlogit}{rlogit}(x, par = 1, y = NULL, norm = TRUE, ...)

}
\arguments{
  \item{x}{ a \code{rlogit} or a \code{rpar} object,}
  \item{norm}{the variable used for normalization if any,}
  \item{par}{the required parameter for the \code{rlogit} methods (either
    the name or the position of the parameter,}
  \item{y}{values for which the function has to be evaluated,}
  \item{...}{further arguments.}
}
\value{
  a numeric vector for \code{qrlogit}, \code{drlogit} and \code{prlogit}, a
  numeric value for \code{mean}, \code{stdev} and \code{med}.
}
\details{
  \code{mean}, \code{med} and \code{stdev} compute respectively the
  mean, the median and the standard deviation of the random
  parameter. \code{qrlogit}, \code{prlogit}, \code{drlogit} return functions that
  compute the quantiles, the probability and the density of the random
  parameters.
}
\seealso{
  \code{\link{rlogit}} for the estimation of random parameters logit models
  and \code{\link{rpar}} for the description of \code{rpar} objects.
}
\examples{
}
\author{Yves Croissant}
\keyword{regression}