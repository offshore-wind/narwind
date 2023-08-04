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
#' @import ggnewscale
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
                         what = "map",
                         whaleID = NULL,
                         animate = FALSE,
                         web = FALSE,
                         ...
){
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default optional arguments
  nmax <- 100
  lwd <- 0.2
  bymonth <- FALSE
  bywhale <- FALSE
  cohortID <- obj$param$cohortID
  if(is.null(whaleID)) whaleID <- 1:obj$param$nsim
  
  if(bymonth & bywhale) stop("<bymonth> and <bywhale> cannot both be set to TRUE")
  if(!what %in% c("inits", "map", "gam", "feed", "migrate")) stop("Unrecognized input to <what> argument")
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")
  
  # Default values
  if(length(args) > 0){
    if("nmax" %in% names(args)) nmax <- args[["nmax"]]
    if("lwd" %in% names(args)) lwd <- args[["lwd"]]
    if("bymonth" %in% names(args)) bymonth <- args[["bymonth"]]
    if("bywhale" %in% names(args)) bywhale <- args[["bywhale"]]
    if("cohortID" %in% names(args)) cohortID <- args[["cohortID"]]
    if("whaleID" %in% names(args)) whaleID <- args[["whaleID"]]
  }
  
  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohortID, abb]
  cohort.names <- cohorts[id %in% cohortID, name]
  n.ind <- obj$param$nsim

  #' -------------------------------------------------------
  # INITIAL LOCATIONS ---- 
  #' -------------------------------------------------------
  
  if(what == "inits"){
    
    if (!is.null(whaleID)) {
      for (k in 1:length(obj$init$xy)) {
        if (!"xyinits" %in% class(obj$init$xy)) stop("Input not recognized")
        par(mfrow = c(3, 4))
        for (h in 1:12) {
          sp:::plot.SpatialPolygons(world, col = "grey", axes = TRUE, main = month.name[h])
          xpos <- sapply(obj$init$xy[[k]][h, ], FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
          ypos <- sapply(obj$init$xy[[k]][h, ], FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
          points(cbind(xpos, ypos), pch = 16, col = "orange")
        }
      }
      par(mfrow = c(1, 1))
    } else {
      coords <- sp::coordinates(density_narw[[1]])
      colnames(coords) <- c("x", "y")
      month.colours <- c("#962938", "#be3447", "#296f96", "#348cbe", "#7db9db", "#b77700", "#ea9800", "#ffb01e", "#ffd104", "#41cbb8", "#d05566", "#db7d8a")
      for (k in 1:length(obj$init$xy)) {
        xinit <- matrix(data = coords[obj$init$xy[[k]], "x"], nrow = nrow(obj$init$xy[[k]]))
        yinit <- matrix(data = coords[obj$init$xy[[k]], "y"], nrow = nrow(obj$init$xy[[k]]))
        sp::plot(regions, main = paste0(cohort.names[k], " -- Individual # ", whaleID))
        sp::plot(world, col = "grey", add = TRUE)
        points(xinit[, whaleID], yinit[, whaleID], pch = 16, col = month.colours)
        legend("bottomright", legend = month.abb, fill = month.colours)
      }
    }
    
    #' -------------------------------------------------------
    # FITTED GAMs ----
    #' -------------------------------------------------------  
    
    } else if (what == "gam") {

      plot(obj$gam$fit$surv, scheme = 1, scale = 0, pages = 1)
      plot(obj$gam$fit$bc, scheme = 1, scale = 0, pages = 1)
    
  # } else if (what == "prob") {
  #   
  #   cohorts <- obj$param$cohorts
  #   
  #   x_axis <- seq(0,1,0.01)
  #   
  #   bc_preds <- obj$gam$pred$bc
  #   surv_preds <- obj$gam$pred$surv
  #   
  #   if(5 %in% cohortID) which.cohorts <- c(0, cohortID) else which.cohorts <- cohortID
  #   
  #   pred.df <- tibble::tibble(cohort = factor(do.call(c, purrr::map(.x = which.cohorts, .f = ~rep(.x, each = length(x_axis))))),
  #                             start_bc = rep(x_axis, length(which.cohorts)))
  #   
  #   pred.df$surv <- apply(X = pred.df, MARGIN = 1, FUN = function(x) surv_preds[[as.character(x[1])]](x[2]))
  #   pred.df$bc <- apply(X = pred.df, MARGIN = 1, FUN = function(x) bc_preds[[as.character(x[1])]](x[2]))
  #   
  #   pred.df <- pred.df |> tidyr::pivot_longer(!c(cohort, start_bc), names_to = "variable", values_to = "value")
  #   
  # 
  #   to_string <- ggplot2::as_labeller(c('bc' = "Body condition", 'surv' = "Survival"))
  #   linecol <- c("black", "#f46a9b", "#ef9b20", "#edbf33", "#87bc45", "#27aeef", "#b33dc6")
  #   
  #   
  #   ggplot2::ggplot(data = pred.df) +
  #     ggplot2::geom_line(aes(x = start_bc, y = value, col = cohort), linewidth = 0.65) +
  #     theme_narw() +
  #     xlab("Body condition") +
  #     ylab("Probability") +
  #     ggplot2::scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  #     ggplot2::scale_color_manual(values = linecol, labels = cohorts$name) +
  #     ggplot2::facet_wrap(~variable, scales = "free", labeller = to_string) +
  #     ggplot2::theme(legend.position = "right",
  #                    legend.title = ggplot2::element_blank())
  #   
    #' -------------------------------------------------------
    # MOVEMENT TRACKS ----
    #' -------------------------------------------------------  
    
  } else if (what %in% c("map", "feed", "migrate")) {
  
  if(n.ind > nmax) {warning("Plotting only the first ", nmax, " tracks"); n.ind <- nmax}
  if(is.null(whaleID)) whaleID <- seq_len(n.ind)

  sim <- obj$sim
  locs.dead <- obj$dead
  locs.birth <- obj$birth
  
  locations <- purrr::set_names(x = cohort.ab) |>
    purrr::map(.f = ~sim[[.x]][day > 0, list(date, month, whale, easting, northing, region, cohort_name, feed, north)])
  
  tracks <- data.table::rbindlist(locations)

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
    ggplot2::ylab("")
  
  # ....................................
  # Generate plots
  # ....................................
  
  COLORS <- c('starve' = 'firebrick', 'strike' = 'deepskyblue4', 'birth' = "#15A471", 'other' = 'darkorange')
  SHAPES <- c('Calves' = 1, 'Juveniles' = 16, 'Adults' = 16)
  
  plot.list <- purrr::set_names(cohortID) |>
  purrr::map(.f = ~ {

    # Extract movement tracks
    tracks.cohort <- tracks[whale %in% whaleID & cohort_name == cohorts[id == .x, name]]
    
    if(.x == 5){
      new.name <- paste0(unique(tracks.cohort$cohort_name), " + calves")
      # tracks[cohort_name == .y, cohort_name:= new.name]
      tracks.cohort$cohort_name <- new.name
      locs.dead[cohort %in% c(0,5), cohort_name:= new.name]
      # if("data.table" %in% class(locs.birth)) 
      locs.birth[cohort %in% c(0,5), cohort_name:= new.name]
    }

    # Plot base map (land)
    base_p <- gg.opts + 
      ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", linewidth = 0.25) +
      theme_narw(vertical = TRUE)

    # Produce required plot(s)
    if(bymonth) {
      
      track_p <- 
        base_p +
        ggplot2::geom_path(data = tracks.cohort, 
                           mapping = ggplot2::aes(x = easting, y = northing, group = whale, colour = factor(month)), 
                                               alpha = 0.7, linewidth = lwd) +
        ggplot2::scale_color_manual(values = pals::viridis(12), labels = month.abb) +
        labs(color = NULL) +
        ggplot2::coord_sf(expand = FALSE) +
        ggplot2::facet_wrap(~cohort_name) +
        ggplot2::theme(
          legend.title = element_blank(),
          legend.position = c(1000, -500),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = "transparent")) 
      
    } else {
      
      track_p <- 
        
        base_p +
        
        # ............................................................
        # Color by individual
        # ............................................................
        
        {if(bywhale) ggplot2::geom_path(data = tracks.cohort,
                mapping = ggplot2::aes(x = easting, y = northing, group = whale, 
                                       col = factor(whale)), alpha = 0.7, linewidth = lwd) } +
        
        {if(bywhale & length(whaleID) <= 10) ggplot2::scale_color_manual(values = whaleID)} +
        
        {if(bywhale & length(whaleID) <= 10) ggnewscale::new_scale_colour() }
      
      if(what == "feed"){
        
        track_p <- track_p +
          ggplot2::geom_point(data = tracks.cohort,
                              mapping = ggplot2::aes(x = easting, y = northing, group = whale, colour = factor(feed))) 
        
      } else if(what == "migrate"){
        
        track_p <- track_p +
          ggplot2::geom_path(data = tracks.cohort,
            mapping = ggplot2::aes(x = easting, y = northing, group = whale, colour = factor(north))) +
          
          {if(.x == 5 & "data.table" %in% class(locs.birth))
            ggplot2::geom_point(data = locs.dead[cohort == 0 & whale %in% whaleID & abb %in% cohorts[id == 0, abb]],
                                aes(x = easting, y = northing, colour = factor(cause_death), shape = factor(class)))}
        
      } else {
        
        track_p <- track_p +
          
          # ............................................................
        # Same color
        # ............................................................
        
        {if(!bywhale) ggplot2::geom_path(data = tracks.cohort,
             mapping = ggplot2::aes(x = easting, y = northing, group = whale), alpha = 0.7, linewidth = lwd)} +
        
        # ............................................................
        # Mortality - adults
        # ............................................................
        
        # Adults / juveniles
        ggplot2::geom_point(data = locs.dead[cohort > 0 & whale %in% whaleID & abb %in% cohorts[id == .x, abb]],
                            aes(x = easting, y = northing, colour = factor(cause_death), shape = factor(class))) +
        
        # ............................................................
        # Calving events
        # ............................................................
        
        {if(.x == 5 & "data.table" %in% class(locs.birth))
            ggplot2::geom_point(data = locs.dead[cohort == 0 & whale %in% whaleID & abb %in% cohorts[id == 0, abb]],
                                aes(x = easting, y = northing, colour = factor(cause_death), shape = factor(class)))} +
        
        {if(.x == 5 & "data.table" %in% class(locs.birth))
            ggplot2::geom_point(data = locs.birth[whale %in% whaleID], aes(x = easting, y = northing, colour = factor(event)))
        } +
        
        # ............................................................
        # Color and theme options
        # ............................................................
        
        ggplot2::scale_color_manual(values = COLORS) +
        ggplot2::scale_shape_manual(values = SHAPES)

    }
     
    track_p <- track_p +
      ggplot2::theme(
          legend.title = element_blank(),
          legend.position = c(1000, -500),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = "transparent"))  +
        labs(color = NULL) +
        ggplot2::facet_wrap(~cohort_name) +
        ggplot2::coord_sf(expand = FALSE) +
        ggplot2::theme(plot.margin = unit(c(-0.30,0,0,0), "null"))
    }

    print(track_p)
  })
  
  if(web){
    web_p <- purrr::map(.x = plot.list, .f = ~plotly::ggplotly(.x))
    purrr::walk(web_p, print)
  } else {
    print(patchwork::wrap_plots(plot.list) + patchwork::plot_layout(guides = "collect"))
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
}