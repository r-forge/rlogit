\name{rpar}
\alias{rpar}
\alias{print.rpar}
\alias{summary.rpar}

\title{random parameter objects}
\description{\code{rpar} objects contain the relevant information about
  estimated random parameters. The homonymous function extract on
  \code{rpar} object from a \code{rlogit} object.
}
\usage{
rpar(x, par, ...)
}
\arguments{
  \item{x}{a \code{rlogit} object,}
  \item{par}{the name of the parameter to be extracted,}
  \item{...}{further arguments.}
}
\value{
  a \code{rpar} object, which contain :
  \item{dist}{the name of the distribution,}
  \item{mean}{the first parameter of the distribution,}
  \item{sigma}{the second parameter of the distribution,}
  \item{name}{the name of the parameter,}
  \item{norm}{the name of the coefficient which is used for
    normalisation, if any.}
}

\details{
  \code{rlogit} objects contain an element called rpar which contain a list
  of \code{rpar} objects, one for each estimated random parameter. The
  \code{print} method prints the name of the distribution and the
  parameter, the \code{summary} behave like the one for numeric vectors.
}
\seealso{
  \code{\link{rlogit}} for the estimation of a random parameters logit model.
}

\author{Yves Croissant}
\keyword{regression}