#' Assessment of offshore wind impacts on North Atlantic right whales
#'
#' The \code{narwind} package integrates a spatially explicit individual-based model and a stochastic population model to explore how the construction of wind farms along the U.S. east coast may affect the future abundance of the critically endangered population of North Atlantic right whales (Eubalaena glacialis). The software provides functionality for crafting custom, testable offshore wind development scenarios reflecting real-world offshore development activities.   
#'
#' @author Phil J. Bouchet
#' @docType package
#' @name narwind
#' @useDynLib narwind
NULL

# Quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "4.0.0")  utils::globalVariables(c("."))