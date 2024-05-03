#' Update scenario object when app stops and return to R environment
#'
#' @param obj An object of class \code{narwscenario}, as returned by \code{\link{scenario}}.
#'
#' @return An updated \code{narwscenario} object.
#' @export
#' @author Rob Schick

saveParameters <- function(obj){
  cat("saveParameters is being called.\n")
  # print(str(obj))
  assign("scenario_custom", obj, envir = .GlobalEnv)
  
}