#' Update scenario object
#'
#' @param obj An object of class \code{narwscenario}, as returned by \code{\link{scenario}}.
#' @param vessel.tbl A \code{data.frame} object defining vessel traffic parameters from custom user input to the Shiny app. 
#'
#' @return An updated \code{narwscenario} object.
#' @export
#' @author Phil J. Bouchet
update.narwscenario <- function(obj, 
                                vessel.tbl){
  
  # Convert to data.table and match column names
  vessel.tbl <- data.table::as.data.table(vessel.tbl) |> 
    janitor::clean_names(parsing_option = 0) |> 
    dplyr::rename(Nvessels = nvessels,
                  speed_knt = speed.knt,
                  roundtrips_foundation = roundtrips.foundation)
  
  # Update vessel transit data
  obj$vessels[vessel.tbl, on = "vesselclass", Nvessels := i.Nvessels]
  obj$vessels[vessel.tbl, on = "vesselclass", speed_knt := i.speed_knt]
  obj$vessels[vessel.tbl, on = "vesselclass", roundtrips_foundation := i.roundtrips_foundation]
  
  return(obj)
}