\name{halton}
\alias{halton}

\title{Halton sequences}

\description{
  This function enables the creation of a halton sequence with a user
  supplied length and prime number
}

\usage{
halton(prime = 3, length = 100, drop = 10)
}
\arguments{
  \item{prime}{a prime number,}
  \item{length}{the length of the sequence,}
  \item{drop}{the number of numbers droped at the begining of the
    sequence.}

}

\value{
  a numeric vector
}

\author{Yves Croissant}
\keyword{regression}