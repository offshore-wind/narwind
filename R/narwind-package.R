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