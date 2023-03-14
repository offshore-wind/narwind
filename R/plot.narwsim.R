#' Plot
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @param n.ind Number of animals to plot.
#' @param id ID of animals to plot
#' @param animate Logical. If TRUE, creates a .gif animation of the tracks.
#' @param web If TRUE, returns a plotly graph
#' @import ggplot2
#' @import patchwork
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

plot.narwsim <- function(obj,
                         whale.id = NULL,
                         animate = FALSE,
                         web = FALSE,
                         nmax = 100
){
  
  # obj = m
  # whale.id = NULL
  # animate = FALSE
  # web = FALSE
  # nmax = 100
  
  if(!"narwsim" %in% class(obj)) stop("Input must be an object of class <narwsim>")
  
  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  cohort.id <- obj$param$cohort.id
  n.ind <- obj$param$nsim
  if(n.ind > nmax) {warning("Plotting only the first ", nmax, " tracks"); n.ind <- nmax}
  if(is.null(whale.id)) whale.id <- seq_len(n.ind)
  
  if(any(cohort.id == 5)) cohort.names[which(cohort.id == 5)] <- paste0(cohort.names[which(cohort.id == 5)], " + calves")

  sim <- obj$sim
  locs.dead <- obj$mrt$locs
  locations <- purrr::set_names(x = cohort.ab) |>
    purrr::map(.f = ~sim[[.x]][, list(date, whale, easting, northing, region, cohort_name)])


  
  # ....................................
  # Load required objects
  # ....................................
  
  support_grd <- support_poly |> sf::st_as_sf()
  
  # ....................................
  # Plotting options
  # ....................................
  
  gg.opts <- ggplot2::ggplot() + 
    ggplot2::geom_sf(data = support_grd, fill = "orange", linewidth = 0, na.rm = TRUE, alpha = 0.5) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5),
                   axis.text = ggplot2::element_text(colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
                   axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                   axis.title = element_text(size = 12),
                   axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                   axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                   strip.background = element_rect(fill = "grey20"),
                   strip.text = element_text(colour = 'white', size = 12))
  
  # ....................................
  # Generate plots
  # ....................................
  
  plot.list <- purrr::set_names(cohort.ab) |> 
    purrr::map2(.y = cohort.names, .f = ~{
    
    tracks <- lapply(seq_len(n.ind), function(x) data.frame(locations[[.x]][whale == x,])) |>
      do.call(what = rbind) |> 
      dplyr::filter(whale %in% whale.id)
    
    locs.dead <- purrr::map(.x = locs.dead, .f = ~ .x |> 
                              dplyr::filter(whale %in% whale.id) |> 
                              dplyr::mutate(cause_death = ifelse(strike == 1 & bc > 0.05, "strike", "starve")))
    
    # ....................................
    # Plot base map (land)
    # ....................................
    
    base_p <- gg.opts + ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", linewidth = 0.25)
    
    # ....................................
    # Plot tracks
    # ....................................
    
    track_p <- base_p +
      
      ggplot2::geom_path(
        data = tracks,
        mapping = ggplot2::aes(x = easting, y = northing, group = whale), alpha = 0.7, linewidth = 0.2) +
      
      {if(nrow(locs.dead[[.x]]) > 0) 
        ggplot2::geom_point(data = locs.dead[[.x]], aes(x = easting, y = northing, colour = factor(cause_death)))} +
      
      ggplot2::scale_color_manual(values = c("firebrick", "royalblue4")) +
      ggplot2::theme(legend.position = "bottom",
                     legend.text = element_text(size = 10),
                     legend.key = element_rect(fill = "transparent")) +
      labs(color = NULL) +
      ggplot2::coord_sf(expand = FALSE) +
      # ggtitle(.y) +
      ggplot2::facet_wrap(~cohort_name)

    
    track_p
    
  })
  
  if(web){
    web_p <- purrr::map(.x = plot.list, .f = ~plotly::ggplotly(.x))
    purrr::walk(web_p, print)
  } else {
    print(patchwork::wrap_plots(plot.list))
  }
  
  if(animate){
    
    # ....................................
    # If multiple individuals 
    # ....................................
    
    # Create all combinations of data
    # track.df <- tidyr::expand_grid(
    #   id = unique(locs$trackID),
    #   date = seq.Date(from = min(locs$date), to = max(locs$date), by = 1))
    
    # # Create complete dataset
    # df_all <- dplyr::left_join(track.df, locs)
    
    # mov <- 
    #   
    #   # Base map
    #   base_p + 
    #   
    #   # Simulated track
    #   ggplot2::geom_path(
    #     data = locs,
    #     mapping = ggplot2::aes(x = easting, y = northing, group = trackID, fill = trackID), alpha = 0.3) +
    #   
    #   geom_point(data = locs,
    #              aes(x = easting, y = northing, group = trackID, fill = trackID), alpha = 0.7, shape = 21, size = 2)
    
    # lines and points
    # geom_path(data = df_all, 
    #           aes(x=lon,y=lat,group=id,color=spd), 
    #           alpha = 0.3)+
    # geom_point(data = df_all, 
    #            aes(x=lon,y=lat,group=id,fill=spd),
    #            alpha = 0.7, shape=21, size = 2)+
    
    # anim <- mov + 
    #   
    #   # Reveal data along given dimension - here, allows data to gradually appear through time.
    #   gganimate::transition_reveal(along = date) +
    #   
    #   # Control easing of aesthetics – easing defines how a value changes to another during tweening
    #   gganimate::ease_aes('linear') +
    #   
    #   # Use date as title
    #   ggplot2::ggtitle("Date: {frame_along}")
    # 
    # anim.out <- gganimate::animate(anim, 
    #                                nframes = 5, 
    #                                fps = 10,
    #                                height = 2, 
    #                                width = 3, 
    #                                units = "in", 
    #                                res = 150)
    
    
    # movevis
    
    locs$time <- as.POSIXct(locs$date)
    locs$date <- lubridate::parse_date_time(locs$date, "Ymd")
    locs.m <- moveVis::df2move(df = locs, proj = narw_crs(), x = "easting", y = "northing", time = "date") 
    m <- moveVis::align_move(locs.m, res = 1, digit = 0, unit = "days", spaceMethod = "greatcircle")
    frames <- moveVis::frames_spatial(m,  # move object
                                      map_service = "carto",
                                      map_type = "light",  # base map
                                      # map_token = "pk.eyJ1IjoicGIyODIiLCJhIjoiY2xhc2I3ZmNzMW50azNvbnYxZ3M4NDZyeCJ9.0ECTFIih0onl92DP4-ZosA",
                                      path_size = 2, 
                                      path_colours = c("orange"), 
                                      alpha = 0.5) |>   # path
      moveVis::add_labels(x = "Longitude", y = "Latitude") |>  # add some customizations
      moveVis::add_northarrow(colour = "black", position = "bottomright") |>  
      moveVis::add_scalebar(colour = "black", position = "bottomleft") |>  
      moveVis::add_timestamps(m, type = "label") |> 
      moveVis::add_progress(size = 2)
    
    # world_sf <- get_world() |> sf::st_as_sf()
    # frames.gg <- moveVis::add_gg(frames, gg = ggplot2::expr(ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", size = 0.25)))
    # moveVis::animate_frames(frames.gg, out_file = "animation2.gif", overwrite = TRUE, display = FALSE)
    
    stamp <- gsub(pattern = "-", replacement = "", x = Sys.time()) |> 
      gsub(pattern = " ", replacement = "_") |> 
      gsub(pattern = ":", replacement = "")
    moveVis::animate_frames(frames, out_file = paste0("NARW_animation_", stamp, ".gif"), overwrite = TRUE, display = FALSE)
  }
  
}