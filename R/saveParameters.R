#' Update scenario object when app stops and return to R environment
#'
#' @param obj An object of class \code{narwscenario}, as returned by \code{\link{scenario}}.
#'
#' @return An updated \code{narwscenario} object.
#' @export
#' @author Rob Schick

saveParameters <- function(obj){
  
  assign("scenario_01", obj, envir = .GlobalEnv)
  
}