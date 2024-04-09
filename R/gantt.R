#' Schedule of piling activities
#'
#' Plots a gantt chart showing the  
#' 
#' @param obj An object of class \code{narwscenario}.
#' @export
#' @author Phil J. Bouchet
#' @return A ggplot2 object
gantt <- function(obj){
  
  routes <- obj$routes |> do.call(what = rbind)
  turbines <- obj$locs
  start.day <- obj$start.day
  start.month <- obj$start.month
  piles.per.day <- obj$piles.per.day
  nfarms <- length(unique(routes$windfarm))
  
  # Determine if dates have been supplied
  has.dates <- "date" %in% names(turbines)
  
  # Number of turbines per site
  n.turbines <- if(nrow(turbines) == 0) 0 else turbines[, .N, windfarm]$N
  
  if(!has.dates){
    
    # If dates are not provided, need to compute them
    if(!nrow(turbines)==0){
      
      # Format piling dates
      start.day <- stringr::str_pad(start.day, width = 2, pad = "0")
      
      # Duration of piling activities
      piling.durations <- ceiling(n.turbines / piles.per.day)
      
      start.date <- paste0(lubridate::year(lubridate::now())+1,
                           stringr::str_pad(start.month, width = 2, pad = "0"), start.day)
      
      piling.dates <- purrr::map2(
        .x = start.date,
        .y = piling.durations,
        .f = ~ {
          out <- seq(
            from = lubridate::ymd(.x),
            by = "day",
            length.out = .y
          )
          if(any(grepl(pattern = "02-29", x = out))){
            out <- out[which(!grepl(pattern = "02-29", x = out))]
            out <- c(out, out[length(out)] + 1)
          }
          out
        })
      
    } else {
      
      piling.durations <- rep(365, nfarms)
      piling.dates <- purrr::map(.x = seq_len(nfarms), .f = ~get_dates(gantt = TRUE))
      
    }
    
  } else {
    
    piling.durations <- turbines[,list(duration = length(unique(date))),windfarm]$duration
    piling.dates <- lapply(X = seq_len(nfarms), FUN = function(a) turbines[windfarm == a,]$date)
    
  }
  
  wind.sites <- routes |> 
    data.table::as.data.table() |> 
    dplyr::group_by(site, windfarm) |>
    dplyr::summarise(n = 1, .groups = "keep")
  
  # Create a gantt chart
  gantt <- data.table::data.table(
    wp = wind.sites$site,
    activity = paste0("Wind farm ", wind.sites$windfarm),
    start_date = do.call(c, lapply(X = piling.dates, FUN = function(x) min(x))),
    end_date = do.call(c, lapply(X = piling.dates, FUN = function(x) max(x)))
  )
  
  gg.gantt <- create_gantt(project = gantt,
                           project_start_date = "2024-01",
                           font_family = "sans",
                           by_date = TRUE,
                           exact_date = TRUE,
                           mark_quarters = TRUE,
                           mark_years = F,
                           size_wp = 3,
                           size_activity = 2,
                           line_end = "butt", # "butt", "square"
                           month_breaks = 1,
                           hide_wp = FALSE,
                           show_vertical_lines = TRUE,
                           alpha_wp = 1,
                           size_text = 1.3,
                           alpha_activity = 0.35, 
                           colour_palette = rep("black", 5),
                           axis_text_align = "right")
  
  print(gg.gantt)
  return(gg.gantt)
}

