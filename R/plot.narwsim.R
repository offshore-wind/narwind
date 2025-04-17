#' Model plots
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @param what 
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
  if("alpha" %in% names(args)) alpha <- args[["alpha"]] else alpha <- 0.7
  if("bymonth" %in% names(args)) bymonth <- args[["bymonth"]] else bymonth <- FALSE
  if("bywhale" %in% names(args)) bywhale <- args[["bywhale"]] else bywhale <- FALSE
  if("cohort" %in% names(args)) cohortID <- args[["cohort"]] else cohortID <- obj$param$cohort
  if("whale" %in% names(args)) whaleID <- args[["whale"]] else whaleID <- 1:obj$param$nsim
  if("vid" %in% names(args)) vid <- args[["vid"]] else vid <- "gif"
  if("bbox" %in% names(args)) bbox <- args[["bbox"]] else bbox <- TRUE
  
  if(bymonth & bywhale) stop("<bymonth> and <bywhale> cannot both be set to TRUE")
  
  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohortID, abb]
  cohort.names <- cohorts[id %in% cohortID, name]
  n.ind <- length(whaleID)
  
  #' -------------------------------------------------------
  # INITIAL LOCATIONS ---- 
  #' -------------------------------------------------------
  
  if(what == "inits"){
    
    if (!is.null(whaleID)) {
      for (k in 1:length(obj$init$xy)) {
        if (!"xyinits" %in% class(obj$init$xy)) stop("Input not recognized")
        par(mfrow = c(4, 4))
        for (h in 1:15) {
          sp:::plot.SpatialPolygons(world, col = "grey", axes = TRUE, main = c(month.name, month.name[10:12])[h])
          xpos <- sapply(obj$init$xy[[k]][h, whaleID], 
                         FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
          ypos <- sapply(obj$init$xy[[k]][h, whaleID], 
                         FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
          points(cbind(xpos, ypos), pch = 16, col = "orange")
        }
      }
      par(mfrow = c(1, 1))
    } else {
      coords <- sp::coordinates(density_narw[[1]])
      colnames(coords) <- c("x", "y")
      month.colours <- c("#962938", "#be3447", "#296f96", "#348cbe", "#7db9db", "#b77700", "#ea9800", "#ffb01e", "#ffd104", "#41cbb8", "#d05566", "pink", "royalblue4", "black", "purple")
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
      
      plot(obj$gam$fit$surv, scheme = 1, scale = 0, pages = 1, ylab = "s(Start BC)", xlab = "Start BC", main = "Survival")
      plot(obj$gam$fit$bc, scheme = 1, scale = 0, pages = 1, ylab = "s(Start BC)", xlab = "Start BC",main = "Body condition")
      
      plot(obj$gam$fit_calf$surv, scheme = 1, scale = 0, ylab = "s(Start BC)", xlab = "Start BC", main = "Survival")
      plot(obj$gam$fit_calf$bc, scheme = 1, scale = 0, ylab = "s(Start BC)", xlab = "Start BC", main = "Body condition")
      
    } else {
      
      # Extract objects
      modlist <- obj$gam$fit # Fitted GAMs
      modlist_calf <- obj$gam$fit_calf
      
      pred.x <- obj$post$pred.x$all # Prediction points -- used for mean curve
      pred.x_calf <- obj$post$pred.x$calf
      
      bc.range <- obj$post$bc.range # Allowable range of BC
      nsamples <- obj$post$nsamples # Number of posterior draws
      
      # mbc.x <- obj$post$mbc.x
      preds.surv <- obj$post$preds$survbc$surv # Posterior samples from survival model
      preds.bc <- obj$post$preds$survbc$bc # Posterior samples from body condition model
      # preds.mbc <- obj$post$preds$mbc
      
      preds.surv_calf <- obj$post$preds$survbc_calf$surv # Posterior samples from survival model
      preds.bc_calf <- obj$post$preds$survbc_calf$bc # Posterior samples from body condition model
      
      # Create prediction data.frame
      pred.df <- expand.grid(start_bc = seq(bc.range[1], bc.range[2], length.out = 100), 
                             cohort = cohorts$id, 
                             mcmc = seq_len(nsamples))
      
      pred.df$pred_surv <- as.numeric(t(preds.surv))
      pred.df$pred_bc <- as.numeric(t(preds.bc))
      
      pred.df$pred_surv_calf <- as.numeric(t(preds.surv_calf))
      pred.df$pred_bc_calf <- as.numeric(t(preds.bc_calf))
      
      # Thickness of lines on plot
      linewidth <- 0.65
      
      # Randomly picks 5 spline parameters and plots their trace
      # sample.of.five <- purrr::map(.x = obj$post$samples, .f = ~sample(1:ncol(.x), size = 5, replace = FALSE))
      # par(mfrow = c(5,2))
      # plot(coda::as.mcmc(obj$post$samples$surv[, sample.of.five$surv]), auto.layout = FALSE, main = "Survival")
      # plot(coda::as.mcmc(obj$post$samples$bc[, sample.of.five$bc]), auto.layout = FALSE, main = "Body condition")
      # # plot(coda::as.mcmc(obj$post$samples$mbc[, sample.of.five$mbc]), auto.layout = FALSE, main = "Gestation")
      # par(mfrow = c(1,1))
      
      # Create plots for survival and body condition
      p <- purrr::map(
        .x = names(modlist),
        .f = ~ {
          
          # Add mean predictions to data.frame
          mean_pred <- cbind(pred.x, pred = predict(modlist[[.x]], pred.x, "response"))
          if(.x == "bc") mean_pred$pred <- clamp(mean_pred$pred, find_maxBC())
          
          # Sample a subset of predictions
          sample.df <- pred.df[pred.df$mcmc %in% sample(nsamples, nL, replace = FALSE), ] |> 
            dplyr::select(start_bc, cohort, mcmc, paste0("pred_", .x))
          
          if(.x == "bc") sample.df$pred_bc <- clamp(sample.df$pred_bc, find_maxBC())
          
          # Rename columns (remove pred_)
          lookup <- c(pred = paste0("pred_", .x))
          sample.df <- dplyr::rename(sample.df, all_of(lookup))
          
          # Define facet names
          facet.names <- cohorts$name
          names(facet.names) <- cohorts$id
          
          if(.x == "surv"){
            nocalf <- obj$gam$dat[cohort > 0]
            calfdat <- obj$gam$dat[cohort== 0 & event == "birth"]
            data.pts <- rbind(nocalf, calfdat)
          } else {
            alive.calf <- unique(obj$gam$dat[cohort == 0 & alive == 1, whale])
            adat <- obj$gam$dat[alive == 1,]
            nolac <- adat[!cohort == 5]
            lacdat <- adat[cohort == 5 & whale  %in% alive.calf]
            data.pts <- rbind(nolac, lacdat)
          }

          names(data.pts) <- gsub(pattern = ifelse(.x == "bc", "end_bc", "alive"), replacement = "pred", names(data.pts))
          
          ggplot2::ggplot(data = sample.df, aes(x = start_bc, y = pred)) +
            ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
            {if(.x == "bc") ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") } +
            ggplot2::geom_line(data = mean_pred, linewidth = linewidth) +
            ggplot2::geom_rug(data = data.pts, sides = "b") +
            ggplot2::facet_wrap(vars(cohort), scales = "free_y", 
                                labeller = ggplot2::labeller(cohort = facet.names)) +
            {if(.x == "surv") ggplot2::scale_y_continuous(limits = c(0,1))} +
            {if(.x == "bc") ggplot2::scale_y_continuous(limits = c(0, find_maxBC()))} +
            theme_narw(bbox = bbox) +
            xlab("Starting body condition") +
            ylab(dplyr::case_when(
              .x == "surv" ~ "p(survival)",
              .x == "bc" ~ "Final body condition",
              .default = ""
            ))
        }
      ) |> purrr::set_names(nm = names(modlist))
      
      
      p_calf <- purrr::map(
        .x = names(modlist_calf),
        .f = ~ {
          
          # Add mean predictions to data.frame
          mean_pred_calf <- cbind(pred.x_calf, pred = predict(modlist_calf[[.x]], pred.x_calf, "response"))
          if(.x == "bc") mean_pred_calf$pred <- clamp(mean_pred_calf$pred, find_maxBC())
          
          # Sample a subset of predictions
          sample.df <- pred.df[pred.df$mcmc %in% sample(nsamples, nL, replace = FALSE) & pred.df$cohort == 0, ] |> 
            dplyr::select(start_bc, cohort, mcmc, paste0("pred_", .x, "_calf"))
          
          if(.x == "bc") sample.df$pred_bc_calf <- clamp(sample.df$pred_bc_calf, find_maxBC())
          
          # Rename columns (remove pred_)
          lookup <- c(pred = paste0("pred_", .x, "_calf"))
          sample.df <- dplyr::rename(sample.df, all_of(lookup))
          
          # Define facet names
          facet.names <- cohorts$name
          names(facet.names) <- cohorts$id
          
          data.pts <- dplyr::left_join(x = obj$gam$dat[cohort == 0 & event == "birth"],
                                       y = obj$gam$dat[cohort == 5], by = "whale") |> 
            dplyr::select(alive.x, start_bc.y, end_bc.x)
          
          if(.x == "surv"){
            
            data.pts <- data.pts |> dplyr::rename(pred = alive.x, start_bc = start_bc.y, end_bc = end_bc.x)
            
          } else {
           
            data.pts <- data.pts |> dplyr::rename(alive = alive.x, start_bc = start_bc.y, pred = end_bc.x) |> 
              dplyr::filter(alive == 1)
            
          }
          
          ggplot2::ggplot(data = sample.df, aes(x = start_bc, y = pred)) +
            ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
            {if(.x == "bc") ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") } +
            ggplot2::geom_line(data = mean_pred_calf, linewidth = linewidth) +
            ggplot2::geom_rug(data = data.pts, sides = "b") +
            {if(.x == "surv") ggplot2::scale_y_continuous(limits = c(0,1))} +
            {if(.x == "bc") ggplot2::scale_y_continuous(limits = c(0, find_maxBC()))} +
            theme_narw(bbox = bbox) +
            xlab("Female body condition") +
            ylab(dplyr::case_when(
              .x == "surv" ~ "p(survival)",
              .x == "bc" ~ "Calf body condition",
              .default = ""
            ))
        }
      ) |> purrr::set_names(nm = names(modlist))
      
      purrr::walk(.x = p, .f = ~print(.x))
      purrr::walk(.x = p_calf, .f = ~print(.x))
      
    }
    
    #' -------------------------------------------------------
    # MOVEMENT TRACKS ----
    #' -------------------------------------------------------  
    
  } else if (what %in% c("map", "feed", "migrate")) {
    
    if(n.ind > nL & length(whaleID) > nL) {warning("Plotting only the first ", nL, " tracks"); n.ind <- nL}
    
    sim <- obj$sim
    locs.dead <- obj$dead
    locs.birth <- obj$birth
    
    locations <- purrr::set_names(x = cohort.ab) |>
      purrr::map(.f = ~sim[[.x]][day > 0, list(date, day, month, whale, easting, northing, region, cohort, cohort_name, feed, north)])
    
    tracks <- data.table::rbindlist(locations)
    
    # Only retain animals of interest
    tracks <- tracks[whale %in% whaleID]
    
    # ....................................
    # Load required objects
    # ....................................
    
    support_grd <- support_poly |> sf::st_as_sf()
    
    # ....................................
    # Plotting options
    # ....................................
    
    gg.opts <- ggplot2::ggplot() + 
      ggplot2::geom_sf(data = support_grd, fill = "orange", linewidth = 0, na.rm = TRUE, alpha = 0.6) +
      ggplot2::xlab("") +
      ggplot2::ylab("")
    
    # ....................................
    # Generate plots
    # ....................................
    
    # Define legend colors and shapes
    COLORS <- c('Mortality (starve)' = 'firebrick',
                "Mortality (natural)" = "darkorange",
                'Mortality (strike)' = 'deepskyblue4',
                'Birth' = "#15A471")
    
    SHAPES <- c('Calves' = 1, 'Adults/Juveniles' = 16)
    
    locs.dead$cause_death <- sapply(locs.dead$cause_death, switch,
                                    "strike (female)" = "Mortality (strike)",
                                    "strike" = "Mortality (strike)",
                                    "starve" = "Mortality (starve)",
                                    "starve (female)" = "Mortality (starve)",
                                    "natural" = "Mortality (natural)",
                                    "natural (female)" = "Mortality (natural)")
    
    locs.dead$class <- sapply(locs.dead$class, switch,
                              "Adults" = "Adults/Juveniles",
                              "Juveniles" = "Adults/Juveniles",
                              "Calves" = "Calves")
    
    locs.birth$event <- "Birth"
    
    # Extract movement tracks
    tracks.cohort <- tracks[cohort %in% cohortID & whale %in% whaleID]
    month.names <- c(month.abb[1:9], paste0(month.abb[10:12], " (burn in)"), month.abb[10:12])
    tracks.cohort <- tracks.cohort |> 
      dplyr::mutate(month = month.names[month])|> 
      dplyr::mutate(month = factor(month, levels = c(month.names[10:12], month.names[1:9],
                                                     month.names[13:15])))
    
    # Rename lactating cohort so that calves are plotted together with lactating females
    if(5 %in% cohortID){
      new.name <- paste0(unique(tracks.cohort[cohort == 5,]$cohort_name), " + calves")
      tracks.cohort[cohort == 5,]$cohort_name <- new.name
      locs.dead[cohort %in% c(0,5), cohort_name:= new.name]
      locs.birth[cohort %in% c(0,5), cohort_name:= new.name]
    }
    
    # List format for leaflet plotting
    if(web) tracks.cohort <- split(tracks.cohort, f = tracks.cohort$cohort_name)
    
    # Basemap (land)
    base_p <- gg.opts + 
      ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", linewidth = 0.25) +
      theme_narw(vertical = TRUE, bbox = bbox)
    
    create_map <- function(input.tracks){
      
      # Color by month 
      if(bymonth) {
        
        out_p <- base_p +
          ggplot2::geom_path(data = input.tracks, 
                             mapping = ggplot2::aes(x = easting, y = northing, group = whale), 
                             linewidth = 1.5) +
          ggplot2::facet_wrap(~month) +
          labs(color = NULL)

        
      } else {
        
        # Color by individual
        if(bywhale){
          
          out_p <- base_p +
            ggplot2::geom_path(
              data = input.tracks,
              mapping = ggplot2::aes(
                x = easting, y = northing, group = whale,
                col = factor(whale)), alpha = alpha, linewidth = lwd) +
            {if (length(whaleID) <= 10) ggplot2::scale_color_manual(values = whale)} +
            {if (length(whaleID) <= 10) ggnewscale::new_scale_colour()}
          
        } else {
          
          if(what == "feed"){
            
            out_p <- base_p +
              ggplot2::geom_point(data = input.tracks,
                                  mapping = ggplot2::aes(x = easting, 
                                                         y = northing, 
                                                         group = whale, 
                                                         colour = factor(feed))) 
            
          } else if(what == "migrate"){
            
            out_p <- base_p +
              ggplot2::geom_path(data = input.tracks,
                                 mapping = ggplot2::aes(x = easting, 
                                                        y = northing, 
                                                        group = whale, 
                                                        colour = factor(north))) 
            
          } else {
            
            out_p <- base_p +
              ggplot2::geom_path(
                data = input.tracks,
                mapping = ggplot2::aes(x = easting, 
                                       y = northing, 
                                       group = whale), 
                alpha = alpha, 
                linewidth = lwd) +
              
              # Mortality - adults
              ggplot2::geom_point(
                data = locs.dead[cohort > 0 & whale %in% whaleID & cohort %in% unique(input.tracks$cohort)],
                aes(x = easting, 
                    y = northing, 
                    colour = factor(cause_death), 
                    shape = factor(class))) +
              
              # Calf mortality
              {if (5 %in% unique(input.tracks$cohort) & "data.table" %in% class(locs.birth)) {
                ggplot2::geom_point(
                  data = locs.dead[cohort == 0 & whale %in% whaleID],
                  aes(x = easting, 
                      y = northing, 
                      colour = factor(cause_death), 
                      shape = factor(class)))}} +
              
              # Calving events
              {if (5 %in% unique(input.tracks$cohort) & "data.table" %in% class(locs.birth)) {
                ggplot2::geom_point(data = locs.birth[whale %in% whaleID], 
                                    aes(x = easting, 
                                        y = northing, 
                                        colour = factor(event)))}} +
              
              ggplot2::scale_color_manual(values = COLORS) +
              ggplot2::scale_shape_manual(values = SHAPES) +
              ggplot2::labs(shape = "Age class", colour = "Event")
          }
        }
      }
      
        out_p +
          ggplot2::theme(
            legend.background = ggplot2::element_rect(fill = "white", colour = NA),
            legend.text = ggplot2::element_text(size = 10),
            legend.key = ggplot2::element_rect(fill = "transparent"))  +
          ggplot2::theme(legend.title = element_text(face = "bold")) +
          {if(!bymonth) ggplot2::facet_wrap(~cohort_name)} +
          ggplot2::coord_sf(expand = FALSE) +
          {if(!5 %in% unique(input.tracks$cohort)) ggplot2::guides(shape = "none")}
      
    } # End create_map
    
    if(web){
      
      tracks_p <- purrr::map(.x = tracks.cohort, .f = ~create_map(.x))
      web_p <- purrr::map(.x = tracks_p, .f = ~plotly::ggplotly(.x))
      purrr::walk(web_p, print)
      
    } else {
      
      tracks_p <- create_map(tracks.cohort)
      print(tracks_p)
      
    }
  }
}