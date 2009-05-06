\name{rlogit}
\alias{rlogit}
\alias{summary.rlogit}
\alias{print.summary.rlogit}
\title{Random parameter logit model}
\description{
  Estimation by simulated likelihood of the random parameter (or mixed)
  logit model
}
\usage{
rlogit(formula, data, rpar = NULL, correlation = FALSE, weights = NULL,
    id = NULL, start = NULL, R = 50, eta = NULL, halton = NULL,
    fixed = FALSE, norm = NULL, ...)
\method{summary}{rlogit}(object, ...)
\method{print}{summary.rlogit}(x, digits = max(3, getOption("digits") - 2),
    width = getOption("width"), ...)
}

\arguments{
  \item{x, object}{an object of class rlogit}
  \item{formula}{a symbolic description for the model to be estimated,}
  \item{data}{the data,}
  \item{rpar}{a vector indicating which coefficients are random and
    which distribution they are assumed to follow,}
  \item{correlation}{a boolean indicating whether random parameters are assumed
    to be corelated,}
  \item{weights}{a vector of weights,}
  \item{id}{the individual index in case of panel data,}
  \item{start}{a vector of starting values,}
  \item{R}{the number of draws for the random parameters,}
  \item{eta}{a matrix of random numbers,}
  \item{halton}{a boolean indivating whether Halton sequences should be
    used,}
  \item{fixed}{a vector of boolean indicating which coefficients should
    be fixed,}
  \item{norm}{the name of the variable used for normalisation,}
  \item{digits}{the number of digits,}
  \item{width}{the width of the printing,}
  \item{...}{further arguments.}
}

\value{
An object of class \code{"rlogit"}, a list with elements
\code{coefficients}, \code{vcov}, \code{X}, \code{hessian},
\code{gradient}, \code{call}, \code{model}, \code{df.residual},
\code{logLik} and \code{type}.
}
\details{
  \code{rlogit} estimate a random parameter logit model using simulations.
  
}

\examples{
\dontrun{
data("Mode",package="Ecdat")
Mode <- mlogit.data(Mode,choice="choice",shape="wide",varying=2:9)
za <- rlogit(choice~cost+time,Mode,rpar=c(time="n"))
summary(za)

data("Heating",package="Ecdat")
Hea <- mlogit.data(Heating,choice="depvar",varying=c(3:12,17:21),shape="wide")
tt <- rlogit(depvar~ic+oc+(income+agehed):chdim,data=Hea)


}
}

\author{Yves Croissant}

\references{
  Train, K. (2004) \emph{Discrete Choice Modelling, whith Simulations},
  Cambridge University Press.
}
\keyword{regression}
