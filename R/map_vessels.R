#' Vessel strikes
#'
#' Generates monthly rasters describing the risk of vessel strikes from wind farm development activities, as defined in a scenario object. 
#' 
#' @param obj An object of class \code{narwscenario}.
#' @param z Character. Output variable; can be one of "risk" (relative strike risk), "density" (transit density), "transits" (total number of transits per cell), or "dist" (total distance traveled per cell).
#' @param strike_scalar Numeric. Value of the factor used to scale the resulting strike surface(s).
#' @param which.month Integer, or vector of integers. Months of the year for which strike rasters are returned; 1 corresponds to January, 2 to February etc.
#' @param vessel.speed Numeric vector. Nominal transit speeds for each vessel class (Cable Lay, Construction/Crane, Crew Transfer,Heavy Cargo, Support Vessels <100 m, Survey, Tugs).
#' @param speed.limit Numeric vector of length 12, indicating the maximum speeds imposed during slow-downs for each month of the year. Months with \code{NA} values are subject to no speed restrictions.
#' @param baseline Logical. Defaults to \code{FALSE}. If \code{TRUE}, baseline conditions (excluding wind farm development activities) are included in the output.
#' @param do.plot Logical. Defaults to \code{FALSE}. If \code{TRUE}, strike maps are generated.
#' @param spgdf Defaults to \code{FALSE}. If \code{TRUE}, output surfaces are returned as \code{SpatialGridDataFrame} objects.
#' @export
#' @author Phil J. Bouchet
#' @return A list containing strike risk rasters for each month, in \code{SpatialGridDataFrame} or \code{raster} format.
#' @examples
#' \dontrun{
#' library(narwind)
#' vesselmaps <- map_vessels(obj = scenario_01,
#'                           z = "dist",
#'                           strike_scalar = 1e-7,
#'                           which.month = 1:12,
#'                           vessel.speed = NULL,
#'                           speed.limit = c(10,10,10,10,rep(NA,8)),
#'                           baseline = TRUE,
#'                           gantt = TRUE,
#'                           do.plot = TRUE,
#'                           spgdf = FALSE)
#' }

map_vessels <- function(obj = NULL,
                        z = "dist",
                        strike_scalar = 1e-7,
                        which.month = 1:12,
                        vessel.speed = NULL,
                        speed.limit = NULL,
                        baseline = TRUE,
                        do.plot = FALSE,
                        spgdf = FALSE){
  
  # Requires the following packages: 
  # sf, data.table, tidyverse, terra (raster), purrr, rgeoda, lubridate, fasterize, grid
  
  ##'...............................
  ## DATA IMPORT
  ##'...............................
  
  if(!inherits(obj, "narwscenario")) stop("Input object must be of class <narwscenario>.")
  if(!is.null(speed.limit)){
    if(!length(speed.limit) == 12) stop("<speed.limit> must be a numeric vector of length 12.")
    if(all(is.na(speed.limit))) speed.limit <- NULL
  }
  
  if(is.numeric(vessel.speed)){
    if(!length(speed.limit) == 7) stop("<vessel.speed> must be a numeric vector of length 7")
  }
  
  if(is.null(speed.limit)) speed.limit <- rep(NA, 12)
  
  turbines <- obj$locs
  if(is.list(obj$routes)) routes <- obj$routes |> do.call(what = rbind)
  if(!raster::compareCRS(routes, narw_crs())) routes <- sp::spTransform(routes, CRSobj = narw_crs())
  routes <- routes |> sf::st_as_sf()

  vparams <- obj$vessels
  start.month <- obj$start.month
  start.day <- obj$start.day
  piles.per.day <- obj$piles.per.day
  
  # Determine if dates have been supplied
  has.dates <- "date" %in% names(turbines)
  
  # Add route ID
  routes$routeID <- paste0(routes$windfarm, "-", routes$routeID)
  
  # Adjust speeds if supplied
  if(!is.null(vessel.speed)) vparams$speed_knt <- do.call(c, purrr::map(.x = vessel.speed, .f = ~rep(.x, each = 5)))
  
  ##'...............................
  ## PILING
  ##'...............................
  
  # Number of wind farm sites
  nfarms <- length(unique(routes$windfarm))
  
  # Number of turbines per site
  n.turbines <- if(nrow(turbines) == 0) 0 else turbines[, .N, windfarm]$N
  
  ### Function checks ====
  if(!z %in% c("risk", "density", "transits", "dist")) stop("Unrecognised input in <var>.")
  if(!nrow(turbines)==0 & !all(c("windfarm", "longitude", "latitude") %in% names(turbines))) stop("Missing columns in <turbine.csv>")
  if(!all(c("routeID", "speed_knt", "Nvessels") %in% names(vparams))) stop("Missing columns in <ship.params>")
  
  ### Piling activities ====
  
  # Add projected coordinates if missing
  if(nrow(turbines) > 0 & !all(c("x", "y") %in% names(turbines))) turb <- add_xy(turbines)
  
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
  
  
  # Determine the number of days during which piling occurs in each month
  if (!nrow(turbines)==0) {
    
    piling.mdays <- lapply(X = piling.dates, FUN = function(x) {
      tb <- table(lubridate::month(x))
      out <- rep(0, 12)
      out[as.numeric(names(tb))] <- as.numeric(tb)
      out
    })
    
    if (!all.equal(sapply(X = piling.dates, length), piling.durations)) warning("Dates do not match piling durations")

    # Calculate the proportions of piling days in each month of piling
    piling.prop <- lapply(X = piling.mdays, FUN = function(x) x / sum(x))
    
  } else {
    
    # During the operations phase
    piling.prop <- lapply(1:nfarms, FUN = function(x) rep(1/12,12))
  }

  ##'...............................
  ## VESSEL TRAFFIC
  ##'...............................
  
  # Check CRS
  if(!identical(narw_crs()@projargs, sf::st_crs(routes)$input)) routes <- sf::st_transform(routes, crs = narw_crs())
  
  # Total lengths of routes (km)
  routes$route_km <- as.numeric(sf::st_length(routes))
  vparams[sf::st_drop_geometry(routes), on = "routeID", route_km := i.route_km]
  
  if("roundtrips_foundation" %in% names(vparams)){
    
    # Retrieve number of foundations
    vparams[, nturbines:= n.turbines[windfarm]]
    vparams[, roundtrips := Nvessels * round(nturbines * roundtrips_foundation,0)]
    
  } else {
    
    # Determine the time needed to travel each route (both ways) based on vessel speed (1 knot is 1.852 km/hr)
    vparams[, roundtrip_days := ceiling((2 * route_km / (1.852 * speed_knt))/24)]
    
    # Retrieve piling durations and determine maximum number of vessel voyages per month
    vparams[, piling_days := piling.durations[windfarm]]
    vparams[, roundtrips := Nvessels * floor(piling_days / roundtrip_days)]
    
  }
  
  # Intersect with grid
  routes <- suppressWarnings(sf::st_intersection(routes, vessel_grid$sf))
  
  # Calculate segment lengths
  routes$dist_km <- as.numeric(sf::st_length(routes))
  
  # Generate rasters
  strike_r <- lapply(X = which.month, FUN = function(mo){
    
    # Make a copy of the vparams object so that speed restrictions can be applied
    v.params <- vparams
    
    # Apply speed restriction to all vessels
    if(!is.na(speed.limit[mo])){
      v.params[, speed_knt := speed.limit[mo]]
      warning("Speed restrictions applied")
    }
    
    # Determine number of single transits
    v.params[, ntrans := 2 * floor(roundtrips * sapply(X = windfarm, FUN = function(x) piling.prop[[x]][mo]))]
    
    # Calculate p(lethal) strikes based on vessel speed
    v.params[, PLETH := 1 / (1 + exp(-(-1.905 + 0.217 * speed_knt)))]
    
    # Tally up
    vp <- v.params[, list(PLETH = sum(ntrans * PLETH)), routeID]
    vn <- v.params[, list(ntrans = sum(ntrans)), routeID]
    
    # Convert routes to data.table and add relevant attributes
    ro <- data.table::as.data.table(routes)
    
    # Join attributes
    ro[vp, on = "routeID", PLETH := i.PLETH]
    ro[vn, on = "routeID", ntrans := i.ntrans]
    
    # Weigh by distance of voyage segment
    ro[, PLETH_w := PLETH * dist_km]
    
    # Sum by cell
    ro_sum <- ro[, list(area_km2 = mean(area_km2),
                        tot_PLETH = sum(PLETH_w),
                        tot_dist_km = sum(dist_km * ntrans),
                        ntrans = sum(ntrans)), "cellID"] 
    
    # Add to baseline
    grid_out <- data.table::rbindlist(list(ro_sum, strike_layer[month == mo, names(ro_sum), with = FALSE]))
    
    # Re-calculate totals per cell
    grid_out <- grid_out[, list(area_km2 = mean(area_km2),
                                tot_PLETH = sum(tot_PLETH),
                                tot_dist_km = sum(tot_dist_km),
                                ntrans = sum(ntrans)), "cellID"]
    
    # Apply estimated scalar to calculate risk [TO UPDATE, set to 1 for now]
    # and derive transit densities
    if(!is.null(strike_scalar)){
     grid_out[, risk:= strike_scalar * tot_PLETH]
     grid_out[, risk:= ifelse(risk > 1, 1, risk)]
    } else {
      grid_out[, risk:= tot_PLETH]
    }
    grid_out[, tdens := ntrans / area_km2]
    
    # Baseline
    baseline_out <-
      dplyr::left_join(
        x = vessel_grid$sf,
        y = strike_layer[month == mo, names(ro_sum)[!names(ro_sum) == "area_km2"], with = FALSE], by = "cellID") |>
      dplyr::mutate(tdens = ntrans / area_km2) |> 
      {\(.) if(!is.null(strike_scalar))
        dplyr::mutate(., risk = strike_scalar * tot_PLETH)
        else dplyr::mutate(., risk = tot_PLETH)}() |> 
      {\(.) if(!is.null(strike_scalar))
        dplyr::mutate(., risk = ifelse(risk > 1, 1, risk))
        else . }() |> 
      dplyr::mutate(
        tot_PLETH = dplyr::coalesce(tot_PLETH, 0L),
        ntrans = dplyr::coalesce(ntrans, 0L),
        tot_dist_km = dplyr::coalesce(tot_dist_km, 0L),
        tdens = dplyr::coalesce(tdens, 0L)) |>
      dplyr::rename(transits = ntrans, density = tdens, dist = tot_dist_km)
    
    # With added traffic
    strike_out <- dplyr::left_join(x = vessel_grid$sf, y = grid_out[, !c("area_km2"), with = FALSE], by = "cellID") |> 
      dplyr::mutate(tot_PLETH = dplyr::coalesce(tot_PLETH, 0L),
                    ntrans = dplyr::coalesce(ntrans, 0L),
                    tot_dist_km = dplyr::coalesce(tot_dist_km, 0L),
                    tdens = dplyr::coalesce(tdens, 0L)) |> 
      dplyr::rename(transits = ntrans, density = tdens, dist = tot_dist_km)
    
    # With wind farm traffic
    out <- fasterize::fasterize(
      sf = strike_out[z],
      raster = raster::subset(vessel_grid$raster, "cellID"),
      field = z)
    names(out) <- z
    
    if(baseline){
      out <- list(scenario = out)
      names(out[["scenario"]]) <- z
    } 
    
    # Without wind farm traffic
    if(baseline){
      out[["baseline"]] <- 
        fasterize::fasterize(
          sf = baseline_out[z],
          raster = raster::subset(vessel_grid$raster, "cellID"),
          field = z)
      names(out[["baseline"]]) <- z
    }
    
    out
    
  })
  
  names(strike_r) <- month.abb[which.month]
  
  if(do.plot){
    
    # Generate maps of vessel strike risk for each month of the year
    # Accuracy of natural breaks
    acc <- switch(z,
                  "risk" = ifelse(!is.null(strike_scalar), 1e-5, 10),
                  "density" = 0.05,
                  "transits" = 10,
                  "dist" = 100)
    
    z.title <- switch(z,
                      "risk" = "Relative strike risk",
                      "density" = expression(Transits~km^{-2}),
                      "transits" = "No. transits",
                      "dist" = "Distance (km)")
    
    for(mo in which.month){
      
      # Breaks for colour legend
      z_brks <- raster::as.data.frame(strike_r[[month.abb[mo]]][[1]], na.rm = TRUE)[z]
      if(baseline) z_brks <- rbind(z_brks, raster::as.data.frame(strike_r[[month.abb[mo]]][[2]], na.rm = TRUE)[z])
      z_brks <-  unique(roundany(
        c(0, rgeoda::natural_breaks(k = 20, z_brks[z]
        ), max(z_brks[z])),
        accuracy = acc, f = ceiling
      ))
      
      if(baseline){
        
        par(mfrow = c(1,2))
        
        # Baseline risk
        terra::plot(terra::rast(strike_r[[month.abb[mo]]][[2]]),
                    pax = list(cex.axis = 0.85),
                    plg = list(title = z.title, title.adj = 0.25),
                    legend = "bottomright",
                    main = paste0("Baseline", " [", month.name[mo], "]"),
                    col = pals::viridis(20),
                    breaks = z_brks)
        
        terra::plot(terra::vect(world), col = "grey", add = TRUE)
        
      }
      
      # Baseline + wind farm risk
      terra::plot(terra::rast(strike_r[[month.abb[mo]]][[1]]),
                  pax = list(cex.axis = 0.85),
                  plg = list(title = z.title, title.adj = 0.25),
                  legend = "bottomright",
                  main = paste0("Wind scenario", " [", month.name[mo], "]"),
                  col = pals::viridis(20),
                  breaks = z_brks)
      
      terra::plot(terra::vect(world), col = "grey", add = TRUE)
      
      par(mfrow = c(1,1))
      
    }
  }
  
  if(spgdf){
    if(baseline) strike_r <- purrr::map_depth(.x = strike_r, .depth = 2, .f = ~ as(.x, "SpatialGridDataFrame")) else strike_r <- purrr::map(.x = strike_r, .f = ~ as(.x, "SpatialGridDataFrame"))
  }
  return(strike_r)
  
} # End function
