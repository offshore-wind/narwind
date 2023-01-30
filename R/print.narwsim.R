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
  purrr::walk(
    .x = names(obj[["attrib"]]),
    .f = ~ {
      cat("--------------------------\n")
      cat(.x, "\n")
      cat("--------------------------\n")
      print(obj[["attrib"]][[.x]][1:n.rows, , 1])
      cat("\n")
    }
  )
}