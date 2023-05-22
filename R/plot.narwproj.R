#' Plot
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @import ggplot2
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

plot.narwproj <- function(obj, cohort = FALSE, scales = "free"){
  
  start.year <- lubridate::year(lubridate::now())
  end.year <- start.year + obj$param$yrs
  
  df <- obj$prj$mean
  
  if(cohort){
    df <- df[cohort != "North Atlantic right whales"]
  } else {
    df <- df[cohort == "North Atlantic right whales"]
  }
  
  ggplot2::ggplot() +
    ggplot2::geom_path(data = df, aes(x = year, y = mean, group = "cohort"), colour = "#1565C0") +
    ggplot2::geom_ribbon(data = df, aes(x = year, y = mean, ymin = lwr, ymax = uppr), alpha = 0.25, fill = "#1565C0") +
    ggplot2::facet_wrap(~cohort, scales = scales) + 
    ggplot2::scale_x_continuous(breaks = pretty(c(start.year, end.year))) +
    xlab("") + ylab("Abundance") +
    theme_narw()
}
