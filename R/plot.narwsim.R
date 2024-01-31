#' Model plots
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @param what 
#' @param animate Logical. If TRUE, creates a .gif animation of the tracks.
#' @param web If TRUE, returns a plotly graph
#' @import ggplot2
#' @import patchwork
#' @import ggnewscale
#'
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(10000)
#' plot(m)
#' }

plot.narwsim <- function(obj,
                         what = "map",
                         animate = FALSE,
                         web = FALSE,
                         ...
){
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default optional arguments
  if(!what %in% c("inits", "map", "pred", "feed", "migrate")) stop("Unrecognized input to <what> argument")
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")
  
  # Default values
  if("nL" %in% names(args)) nL <- args[["nL"]] else nL <- 100
  if("lwd" %in% names(args)) lwd <- args[["lwd"]] else lwd <- 0.2
  if("bymonth" %in% names(args)) bymonth <- args[["bymonth"]] else bymonth <- FALSE
  if("bywhale" %in% names(args)) bywhale <- args[["bywhale"]] else bywhale <- FALSE
  if("cohort" %in% names(args)) cohort <- args[["cohort"]] else  cohort <- obj$param$cohort
  if("whale" %in% names(args)) whaleID <- args[["whale"]] else whaleID <- 1:obj$param$nsim
  if("vid" %in% names(args)) vid <- args[["vid"]] else vid <- "gif"
  
  if(bymonth & bywhale) stop("<bymonth> and <bywhale> cannot both be set to TRUE")
  
  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohort, abb]
  cohort.names <- cohorts[id %in% cohort, name]
  n.ind <- length(whaleID)

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
        sp::plot(regions, main = paste0(cohort.names[k], " -- Individual # ", whale))
        sp::plot(world, col = "grey", add = TRUE)
        points(xinit[, whaleID], yinit[, whaleID], pch = 16, col = month.colours)
        legend("bottomright", legend = month.abb, fill = month.colours)
      }
    }
    
    #' -------------------------------------------------------
    # FITTED GAMs ----
    #' -------------------------------------------------------  
    
    } else if (what == "pred") {

      if(!"post" %in% names(obj)){
       
      plot(obj$gam$fit$surv, scheme = 1, scale = 0, pages = 1, ylab = "s(start_bc)", main = "Survival")
      plot(obj$gam$fit$bc, scheme = 1, scale = 0, pages = 1, ylab = "s(start_bc)", main = "Body condition")
      
      } else {
        
        # Extract objects
        modlist <- obj$gam$fit
        pred.x <- obj$post$pred.x
        bc.range <- obj$post$bc.range
        nsamples <- obj$post$nsamples
        mbc.x <- obj$post$mbc.x
        preds.surv <- obj$post$preds$survbc$surv
        preds.bc <- obj$post$preds$survbc$bc
        preds.mbc <- obj$post$preds$mbc
        
        # Create prediction data.frame
        pred.df <- expand.grid(start_bc = seq(bc.range[1], bc.range[2], length.out = 100), 
                               cohort = cohorts$id, 
                               mcmc = seq_len(nsamples))
        
        pred.df$pred_surv <- as.numeric(t(preds.surv))
        pred.df$pred_bc <- as.numeric(t(preds.bc))
        
        linewidth <- 0.65
        
        sample.of.five <- purrr::map(.x = obj$post$samples, .f = ~sample(1:ncol(.x), size = 5, replace = FALSE))
        par(mfrow = c(5,2))
        plot(coda::as.mcmc(obj$post$samples$surv[, sample.of.five$surv]), auto.layout = FALSE, main = "Survival")
        plot(coda::as.mcmc(obj$post$samples$bc[, sample.of.five$bc]), auto.layout = FALSE, main = "Body condition")
        plot(coda::as.mcmc(obj$post$samples$mbc[, sample.of.five$mbc]), auto.layout = FALSE, main = "Gestation")
        par(mfrow = c(1,1))
        
        # plot(1:1, xaxt='n', yaxt = 'n', type = 'n', xlab = '', ylab = '', frame.plot=F)
        # text(1,1, "Survival model", cex = 2)
        # par(mfrow = c(5,3))
        # for(pp in sample.of.nine$surv){plot(obj$post$samples$surv[, pp], type = "l", xlab = "", ylab = "")}
        # par(mfrow = c(1,1))
        # 
        # plot(1:1, xaxt='n', yaxt = 'n', type = 'n', xlab = '', ylab = '', frame.plot=F)
        # text(1,1, "Health model", cex = 2)
        # par(mfrow = c(3,3))
        # for(pp in sample.of.nine$bc){plot(obj$post$samples$bc[, pp], type = "l", xlab = "", ylab = "")}
        # par(mfrow = c(1,1))
        # 
        # plot(1:1, xaxt='n', yaxt = 'n', type = 'n', xlab = '', ylab = '', frame.plot=F)
        # text(1,1, "Gestation model", cex = 2)
        # par(mfrow = c(3,3))
        # for(pp in sample.of.nine$mbc){plot(obj$post$samples$mbc[, pp], type = "l", xlab = "", ylab = "")}
        # par(mfrow = c(1,1))
        
        # Create plots for survival and body condition
        p <- purrr::map(
          .x = names(modlist),
          .f = ~ {
            
            # Add predictions to data.frame
            mean_pred <- cbind(pred.x, pred = predict(modlist[[.x]], pred.x, "response"))
            
            # Sample a subset of predictions
            sample.df <- pred.df[pred.df$mcmc %in% sample(nsamples, nL, replace = FALSE), ] |> 
              dplyr::select(start_bc, cohort, mcmc, paste0("pred_", .x))
            
            # Rename columns (remove pred_)
            lookup <- c(pred = paste0("pred_", .x))
            sample.df <- dplyr::rename(sample.df, all_of(lookup))
            
            # Define facet names
            facet.names <- cohorts$name
            names(facet.names) <- cohorts$id
            
            ggplot2::ggplot(data = sample.df, aes(x = start_bc, y = pred)) +
              ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
              ggplot2::geom_line(data = mean_pred, linewidth = linewidth) +
              ggplot2::facet_wrap(vars(cohort), scales = "free_y", 
                                  labeller = ggplot2::labeller(cohort = facet.names)) +
              # ggplot2::scale_y_continuous(limits = c(0,1)) +
              theme_narw() +
              xlab("Starting body condition (%)") +
              ylab(dplyr::case_when(
                .x == "surv" ~ "p(survival)",
                .x == "bc" ~ "Final body condition (%)",
                .default = ""
              ))
          }
        ) |> purrr::set_names(nm = names(modlist))
        
        purrr::walk(.x = p, .f = ~print(.x))
        
        # Create plot for minimum BC required for gestation
        mbc.df <- expand.grid(mass = mbc.x, mcmc = seq_len(nsamples))
        mbc.df$min_bc <- as.numeric(t(preds.mbc))
        
        # Sample a subset of predictions 
        mbc.df <- mbc.df[mbc.df$mcmc %in% sample(nsamples, nL, replace = FALSE), ]
        
        # Add mean response
        pdf <- data.frame(mass = mbc.x, min_bc = predict(gam_gest, data.frame(mass = mbc.x), "response"))
        
        # Create plot
        gp <- ggplot2::ggplot(data = mbc.df, aes(x = mass, y = min_bc)) +
          ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
          ggplot2::geom_line(data = pdf, linewidth = linewidth) +
          theme_narw() +
          scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
          xlab("Total mass (kg)") +
          ylab("Min body condition (%)")
        
        print(gp)

      }
  # } else if (what == "prob") {
  #   
  #   cohorts <- obj$param$cohorts
  #   
  #   x_axis <- seq(0,1,0.01)
  #   
  #   bc_preds <- obj$gam$pred$bc
  #   surv_preds <- obj$gam$pred$surv
  #   
  #   if(5 %in% cohort) which.cohorts <- c(0, cohort) else which.cohorts <- cohort
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
  
  if(n.ind > nL & length(whaleID) > nL) {warning("Plotting only the first ", nL, " tracks"); n.ind <- nL}

  sim <- obj$sim
  locs.dead <- obj$dead
  locs.birth <- obj$birth
  
  locations <- purrr::set_names(x = cohort.ab) |>
    purrr::map(.f = ~sim[[.x]][day > 0, list(date, month, whale, easting, northing, region, cohort_name, feed, north)])
  
  tracks <- data.table::rbindlist(locations)
  tracks <- tracks[whale %in% whaleID]

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
  
  plot.list <- purrr::set_names(cohort) |>
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
        
        {if(bywhale & length(whaleID) <= 10) ggplot2::scale_color_manual(values = whale)} +
        
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

    # Prepare data for animation
    dates <- data.frame(date = get_dates(strip = TRUE), fulldate = get_dates(strip = FALSE))
    animxy <- dplyr::left_join(x = tracks, y = dates, by = "date")
    animxy$fulldate <- lubridate::parse_date_time(animxy$fulldate, "Ymd")
    animxy$individual <- paste0(animxy$cohort_name, "_", animxy$whale)
    animxy <- add_latlon(animxy)
    
    # Convert data into moveStack object
    animxy.mstack <- suppressWarnings(moveVis::df2move(df = animxy,
                                 proj = narw_crs(latlon = TRUE), 
                                 x = "long", 
                                 y = "lat",
                                 time = "fulldate", 
                                 track_id = "individual"))
    
    # Align movement data to a uniform time scale with a uniform temporal resolution throughout 
    # the complete movement sequence -- to prepare the movement data to be interpretable by famres_sp
    animxy.align <- suppressWarnings(alignmove(animxy.mstack))
    
    # Correct time stamps
    animxy.align$time <- animxy.align$time + lubridate::years(lubridate::year(lubridate::now())-2021)    

    # Define animat colours
    coh <- gsub("_", "", unique(gsub("[0-9]+", "", animxy.align@trackId)))
    path_colours <- data.frame(name = coh, cohort = gsub("\\.", ")", gsub("\\.\\.", " (", coh)))
    path_colours$cohort <- gsub(pattern = "female \\(", replacement = "female, ", x = path_colours$cohort)
    path_colours <- dplyr::left_join(path_colours, cohorts[, c("name", "colour")], by = dplyr::join_by(cohort == name))
    animxy.align$colour <- sapply(X = gsub("_", "", gsub("[0-9]+", "", animxy.align@trackId)), FUN = function(x) path_colours[path_colours$name == x, "colour"])
    
    # Code to download basemap tiles
    # basemap <- basemaps::basemap_png(ext = sf::st_bbox(sf::st_as_sf(world)),
    #                                  map_service = "esri",
    #                                  map_type = "world_ocean_base",
    #                                  map_res = 1,
    #                                  map_dir = getwd())
    
    # Create frames of spatial movement maps for animation
    frames <- suppressWarnings(frames_sp(mobj = animxy.align,
                        map_service = "esri",
                        map_type = "world_ocean_base",
                        path_size = 2, 
                        path_legend = FALSE,
                        alpha = 1, 
                        maxpixels = 5e8)) # maxpixels controls basemap resolution
    
    path_colours$longitude <- 0
    path_colours$latitude <- 0
    path_colours$cohort <- factor(path_colours$cohort)
    
    fcol <- sapply(X = levels(path_colours$cohort), FUN = function(x) path_colours$colour[which(path_colours$cohort == x)])

    # print(path_colours)
    
    frames.out <- frames |>
      moveVis::add_labels(x = "", y = "") |>  # add some customizations
      addtimestamps(animxy.align, type = "label", size = 4, alpha = 0.5) |> 
      moveVis::add_progress(size = 4, colour = "#104E8B")
    
    frames.out <- moveVis::add_gg(frames.out,
      gg = expr(geom_point(aes(x = longitude, y = latitude, colour = cohort), data = path_colours)), data = path_colours)

    frames.out <- moveVis::add_gg(frames.out, gg = expr(scale_color_manual(values = fcol, name = "")), data = path_colours)

    frames.out <- moveVis::add_gg(frames.out, gg = expr(theme(
      legend.position = c(0.75, 0.2),
      legend.text = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = "transparent")
    )), data = path_colours)
    
    stamp <- gsub(pattern = "-", replacement = "", x = Sys.time()) |> 
      gsub(pattern = " ", replacement = "_") |> 
      gsub(pattern = ":", replacement = "")
    
    # Export animation as GIF
    moveVis::animate_frames(frames.out, 
                            out_file = paste0("NARW_animation_", stamp, ".", vid), 
                            overwrite = TRUE, 
                            display = FALSE)
  }
  }
}