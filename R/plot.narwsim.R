#' Plot
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @param n Number of animals to plot.
#' @param id ID of animals to plot
#' @param animate Logical. If TRUE, creates a .gif animation of the tracks.
#' @param color.by Whether to color tracks by animal ID.
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
                         id = NULL,
                         animate = FALSE,
                         color.by = FALSE,
                         web = FALSE,
                         world.fill = "lightgrey",
                         nmax = 100
){
  
  # obj = m
  # id = NULL
  # animate = FALSE
  # color.by = FALSE
  # web = FALSE
  # world.fill = "lightgrey"
  # nmax = 100
  
  if(!"narwsim" %in% class(obj)) stop("Input must be an object of class <narwsim>")

  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  cohort.id <- obj$param$cohort.id
  if(any(cohort.id == 5)) cohort.names[which(cohort.id == 5)] <- paste0(cohort.names[which(cohort.id == 5)], " + calves")

  locations <- obj$locs
  
  n <- dim(locations[[1]])[3]
  if(n > nmax) {warning("Plotting only the first ", nmax, " tracks"); n <- nmax}
  
  # Provide id (integer) to pick an animal
  # Set id = -1 for a random animal
  
  # ....................................
  # Load required objects
  # ....................................
  
  support_grd <- targets::tar_read(support_poly) |> sf::st_as_sf()
  wrld <- targets::tar_read(world)
  
  # ....................................
  # Plotting options
  # ....................................
  
  gg.opts <- ggplot() + 
    ggplot2::geom_sf(data = support_grd, fill = "orange", linewidth = 0, na.rm = TRUE, alpha = 0.5) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 0.5),
                   axis.text = ggplot2::element_text(colour = "black"),
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5))
  
  # ....................................
  # Generate plots
  # ....................................
  
  plot.list <- purrr::set_names(cohort.ab) |> 
    purrr::map2(.y = cohort.names, .f = ~{
    
    tracks <- lapply(seq_len(n), function(x) dplyr::mutate(data.frame(locations[[.x]][,,x]), trackID = dimnames(locations[[.x]])[[3]][x]))
    
    # Add date as column
    locs <- purrr::map(.x = tracks, .f = ~{
      tibble::rownames_to_column(.x, var = "date") |> 
        dplyr::mutate(date = lubridate::as_date(date))}) |> 
      do.call(what = rbind)
    
    # Select tracks
    if(is.null(id)) animal_id <- 1:n else animal_id <- unique(locs$trackID)[id]
    locs <- locs |> dplyr::filter(trackID %in% paste0("whale.", animal_id))
    
    # ....................................
    # Plot base map (land)
    # ....................................
    
    base_p <- gg.opts + ggplot2::geom_sf(data = sf::st_as_sf(wrld), fill = world.fill, color = "black", linewidth = 0.25)
    
    # ....................................
    # Plot tracks
    # ....................................
    
    track_p <- base_p +
      ggplot2::geom_path(
        data = locs,
        mapping = ggplot2::aes(x = easting, y = northing, colour = trackID), alpha = 0.7, linewidth = 0.2) +
      {if(color.by) ggplot2::scale_color_manual(name = "Animal ID", values = pals::glasbey(length(animal_id)))} +
      {if(!color.by) ggplot2::scale_color_manual(name = "Animal ID", values = rep("black", length(animal_id)))} +
      ggplot2::theme(legend.position = "none") + 
      ggplot2::coord_sf(expand = FALSE) +
      ggtitle(.y)
    
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
    # frames.gg <- moveVis::add_gg(frames, gg = ggplot2::expr(ggplot2::geom_sf(data = sf::st_as_sf(wrld), fill = "lightgrey", color = "black", size = 0.25)))
    # moveVis::animate_frames(frames.gg, out_file = "animation2.gif", overwrite = TRUE, display = FALSE)
    
    stamp <- gsub(pattern = "-", replacement = "", x = Sys.time()) |> 
      gsub(pattern = " ", replacement = "_") |> 
      gsub(pattern = ":", replacement = "")
    moveVis::animate_frames(frames, out_file = paste0("NARW_animation_", stamp, ".gif"), overwrite = TRUE, display = FALSE)
  }
  
}