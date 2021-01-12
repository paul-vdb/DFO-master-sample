## Main package-level documentation.

#' BASMasterSample: A package to use Balanced Acceptance Sampling (BAS) as a 
#' Master Sample for environmental monitoring.
#'
#'
#' The \code{BASMasterSample} package uses a fixed bounding box, random seed
#' and random rotation to select BAS sample sites from a Master Sample.
#' The default is to select a Master Sample from Canada's Western Marine Master Sample.
#' However, the package is general for any master sample passed to it in the correct format and
#' can generate a new master sample given the boundary shape of interest. To date it can 
#' select random sites, or by stratification.
#'
#' At this stage the package is just for continuous resources but in the future we 
#' will add additional functions for selecting master samples for linear and point features.
#'
#' @section Data Structure:
#'
#' Add information for how to set up your spatial files to be used by this package.
#'
#' @section Selecting Sites:
#'
#' How to use the main function to select sites
#' 
#' @section Stratification:
#'
#' How to use the function for stratification.
#'
#' @section New Master Sample:
#' 
#' How to make a new master sample.
#'
#' @references Robertson, B. L., Brown, J. A., McDonald, T., and Jaksons, P. (2013). 
#'	BAS: Balanced acceptance sampling of natural resources. 
#'	\emph{Biometrics}, 69(3), 776--784.
#' @references Robertson, B. L., McDonald, T., Price, C. J., and Brown, J. A. (2017).
#' A modification of balanced acceptance sampling. 
#' \emph{Statistics & Probability Letters}, 129, 107--112.
#' @references van Dam-Bates, P., Gansell, O., and Robertson, B. (2018). 
#'	Using balanced acceptance sampling as a master sample for environmental surveys. 
#' \emph{Methods in Ecology and Evolution}, 9(7), 1718--1726.
#'
#' @docType package
#' @name ascr
NULL
