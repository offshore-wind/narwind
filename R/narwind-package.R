#' Population consequences of offshore wind
#'
#' The \code{narwind} package provides methods for fitting and selecting among behavioural dose-response models using Bayesian reversible jump Markov Chain Monte Carlo (rjMCMC). 
#' 
#'
#' @author Phil J. Bouchet
#' @docType package
#' @name narwind
#' @useDynLib narwind
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "4.0.0")  utils::globalVariables(c("."))

#' North Atlantic right whale density
#'
#' Data from .
#'
#' @references 
#'
#' Smith J (2006) Title, 41 p.
#'
#' @format A list of twelve \code{SpatialGridDataFrame}
#' \describe{
#'   \item{Jan}{Density surface for January}
#'   \item{Feb}{Density surface for February}
#'   \item{Mar}{Density surface for March}
#'   \item{Apr}{Density surface for April}
#'   \item{May}{Density surface for May}
#'   \item{Jun}{Density surface for June}
#'   \item{Jul}{Density surface for July}
#'   \item{Aug}{Density surface for August}
#'   \item{Sep}{Density surface for September}
#'   \item{Oct}{Density surface for October}
#'   \item{Nov}{Density surface for November}
#'   \item{Dec}{Density surface for December}
#' }
#'
#' @name density_narw
#' @docType data
#' @source Data provided by Jason Roberts (Marine Geospatial Ecology Lab, Duke University).
#' @keywords datasets
NULL