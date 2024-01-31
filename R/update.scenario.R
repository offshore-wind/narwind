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
  
  # vessel.tbl <- data.frame(
  #   VesselClass = c("Cable Lay", "Construction/Crane", "Crew Transfer", "Heavy Cargo", "Support Vessels", "Survey", "Tugs"),
  #   NumberOfVessels = 5,
  #   SpeedOfVessels = c(rep(99,7)),
  #   NumberOfTrips = 10
  # )
  
  # Convert to data.table and match column names
  vessel.tbl <- data.table::as.data.table(vessel.tbl) |> 
    janitor::clean_names(parsing_option = 0) |> 
    dplyr::rename(Nvessels = numberofvessels,
                  speed_knt = speedofvessels,
                  roundtrips_foundation = numberoftrips)
  
  # Update vessel transit data
  obj$vessels[vessel.tbl, on = "vesselclass", Nvessels := i.Nvessels]
  obj$vessels[vessel.tbl, on = "vesselclass", speed_knt := i.speed_knt]
  obj$vessels[vessel.tbl, on = "vesselclass", roundtrips_foundation := i.roundtrips_foundation]
  
  return(obj)
}