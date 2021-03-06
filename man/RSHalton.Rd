% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HaltonFunctions.r
\name{RSHalton}
\alias{RSHalton}
\title{Generate numbers from a Halton Sequence with a random start}
\usage{
RSHalton(n = 10, seeds = c(0, 0), bases = c(2, 3), boxes = 0, J = c(0, 0))
}
\arguments{
\item{n}{Number of points required}

\item{seeds}{Random starting point in each dimension}

\item{bases}{Co-prime base for the Halton Sequence}

\item{boxes}{Halton boxes that points are required to be generated in}

\item{J}{Defines the Halton frame, and relates to the number of boxes used.}
}
\value{
Matrix with the columns, order of point, x in [0,1) and y in [0,1)
}
\description{
For efficiency, this function can generate points along a random start Halton Sequence for
predefined Halton
}
\examples{
\dontrun{
# First 10 points in the Halton Sequence for base 2,3
pts <- RSHalton(n = 10)
# First 10 points in the Halton Sequence for base 2,3 with starting point at the 15th and 22nd index.
pts <- RSHalton(n = 10, seeds = c(14, 21))
}
}
