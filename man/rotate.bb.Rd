% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.r
\name{rotate.bb}
\alias{rotate.bb}
\title{Rotate Bounding box by theta radians}
\usage{
rotate.bb(shp, theta)
}
\arguments{
\item{shp}{A spatial file with the spatial boundary of the sample.}

\item{theta}{Radians of rotation. Positive to the right of pi/2, negative to the left.}
}
\value{
Matrix
}
\description{
Given some shp defined as the boundary of interest, rotate it around the centroid
and return the rotation and the centroid as attributes. This is used for
defining a Master Sample bounding box that has random rotation while ensuring that
the new rotated bounding box fits the shp.
}
\examples{
\dontrun{
data(NS_bioregion)
bb.new <- rotate.bb(NS_bioregion, -pi/3)
}
}
