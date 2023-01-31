#' Overview of model results
#'
#' @param obj Object returned by run_model
#' @param n.rows Number of rows to print
#'
#' @export
#' 
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' animals <- run_model(10)
#' plot(animals)
#' }

print.narwsim <- function(obj, n.rows = 4) {
  
  for(k in seq_along(obj$param$cohort.id)){
    
    if(!obj$param$cohort.id[k] == 5) obj[["sim"]][[k]] <- list(obj[["sim"]][[k]])
    
    cat("\n++++++++++++++++++++++++++++++++++", obj$param$cohort.names[k], "++++++++++++++++++++++++++++++++++\n\n")
    
    purrr::walk(
      .x = seq_along(obj[["sim"]][[k]]),
      .f = ~ {
        if(obj$param$cohort.id[k] == 5){
        cat("--------------------------\n")
        cat(c("Adults", "Calves")[.x], "\n")
        cat("--------------------------\n")
        }
        print(obj[["sim"]][[k]][[.x]][1:n.rows, , 1])
        cat("\n")
      }
    )
  }
}