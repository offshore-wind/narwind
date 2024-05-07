#' Individual-based model of North Atlantic right whales
#'
#' Simulates the movements, behavior, and energy budgets of North Atlantic right whales (Eubalaena glacialis) exposed to offshore wind development activities within U.S. waters. The model operates in daily time steps and runs for a full calendar year.
#' 
#' @param nsim Integer. Number of simulated animals.
#' @param scenario An object of class \code{narwscenario}, as returned by \code{\link{scenario}}.
#' @param label Character. Text label assigned to the simulation object; used for plotting.
#' @param piling.hrs Numeric. Length of time (hours) during which whales cease foraging following a response to pile-driving noise exposure.
#' @param n.cores Integer. Number of cores to use for parallel processing. By default, each population cohort is assigned to a separate worker if enough cores are available. If not, N-2 cores are used, where N is the number of cores available.
#' @param progress Logical. Whether a progress bar should be printed to the R console during execution. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to the simulator.
#' @import data.table
#' @import sf
#' @importFrom foreach `%dopar%`
#' @return A list object of class \code{narwsim}.
#' @export
#' @author Phil J. Bouchet
#' @seealso \code{\link{initialize}}
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' }

narw <- function(nsim = 1e3,
                 scenario = NULL,
                 label = "",
                 piling.hrs = 4,
                 n.cores = NULL,
                 progress = TRUE,
                 ...){
  
  
  # pair: An object of class \code{narwsim}, to which the current simulation must be matched. With the exception of wind farm parameters, simulation conditions between paired runs are identical; pairing is therefore useful for comparative assessments of competing offshore wind scenarios.
  # pair = NULL
  
  # Match random seed if pairing scenario is provided
  # if(!is.null(attr(pair,"seed"))){
  #   assign(".Random.seed",  attr(pair, "seed"), .GlobalEnv)
  # }
  
  # if(!is.null(attr(pair,"seed"))){
  #   rdseed <- attr(pair,"seed")
  # } else {
  #   rdseed <- sample(1:10000, size = 1)
  # }
  # 
  # set.seed(rdseed)
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Start stopwatch
  start.time <- Sys.time()
  
  ##'...............................
  ## Function checks ----
  ##''...............................
  
  # Scenario must be an object of class <narwscenario>, i.e., a list with the following items:
  # 
  # ** VESSEL TRAFFIC **
  # - Turbine locations and piling dates (data.frame) -- compulsory fields: windfarm, x, y, date
  # - Vessel routes (shapefile) -- compulsory fields: routeID, windfarm
  # - Vessel traffic parameters (data.frame) -- compulsory fields: category, windfarm, routeID, speed_knt, Nvessels
  # 
  # ** PILING NOISE **
  # - Ambient noise level
  # - Source level
  # - Noise attenuation
  # - Spherical/cylindrical spreading
  # 
  # ** DISTURBANCE REGIME **
  # - Time horizon for projections
  # - Schedule of activities (yearly)
  
  if(!is.null(scenario)){
    if(!inherits(scenario, "narwscenario")){
      if(is.numeric(scenario)){
        if(!scenario %in% 1:3) stop("Unrecognized scenario.")
      } else {
        stop("Erroneous input to <scenario>")
      }
    }
  }
  
  # if(!is.null(pair) & !inherits(pair, "narwsim")) stop("Input to <pair> must be an object of class <narwsim>")
  
  if(piling.hrs < 0 | piling.hrs > 24) stop("<piling.hrs> must be between 0 and 24.")
  
  ##'...............................
  ## Initialization & DEFAULT VALUES ----
  ##''...............................
  
  # Scenario object
  if(!is.null(scenario) & is.numeric(scenario)) scenario <- get(paste0("scenario_0", scenario))
  
  # Start simulation in October
  init.month <- 10 
  
  # Ambient noise levels
  if(!is.null(scenario)) ambient.dB <- scenario$ambient else ambient.dB <- 80
  
  # Default values
  # Starvation threshold expressed as relative blubber mass
  # As per Pirotta et al. (2018) - to test extreme conditions of leanness
  if("cohort" %in% names(args)) cohort <- args[["cohort"]] else cohort <- 1:6
  if("cease.nursing" %in% names(args)) cease.nursing <- args[["cease.nursing"]] else cease.nursing <- c(FALSE, piling.hrs)
  
  # For package development purposes
  if("stressors" %in% names(args)) stressors <- args[["stressors"]] else stressors <- TRUE
  if("growth" %in% names(args)) growth <- args[["growth"]] else growth <- TRUE
  if("pkgdev" %in% names(args)) pkgdev <- args[["pkgdev"]] else pkgdev <- FALSE
  
  # Vessel strike risk scalar.
  # Factor by which vessel strike risk rasters are multiplied to scale to strike probabilities.
  strike.scalar <- 1e-07
  
  # Scalar for prey surfaces
  prey.scalar <- 29
  
  # Mortality parameters
  starvation.death <- 0.005
  starvation.onset <- 0.05
  
  # Step size for availability radius movement model
  # Need to divide by 2 as there are two half steps (Michelot, 2020)
  # UNUSED - as using Normal kernel
  # step.size <- 120/2
  
  # Number of proposals when sampling new location
  # In line with sample sizes used in case studies from Michelot (2020)
  n.prop <- 50
  
  cohort <- cohort[cohort > 0]
  if(any(!cohort %in% 1:6)) stop("Unrecognized cohort")
  
  # ............................................................
  # Define cohorts and their parameters
  # ............................................................
  
  cohorts <- data.table::data.table(id = 0:6,
                                    name = c("Calves (male, female)",
                                             "Juveniles (male)",
                                             "Juveniles (female)",
                                             "Adults (male)",
                                             "Adults (female, pregnant)",
                                             "Adults (female, lactating)",
                                             "Adults (female, resting)"),
                                    class = c("Calves", "Juveniles", "Juveniles", "Adults", "Adults", "Adults", "Adults"))
  cohorts[, abb:= abbreviate(tolower(name), minlength = 6)]
  cohorts[, abb:= gsub(pattern = "j\\(", replacement = "jv\\(", x = abb)]
  cohorts[, abb:= gsub(pattern = "a\\(", replacement = "ad\\(", x = abb)]
  cohorts[, colour := c("black","#104E8B","#F69554","#22BA9C","#84375A","#EEB422","#942F33")]
  
  cat("--------------------------------------------------------------------------------\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("\n")
  cat("                 NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("                           *** BIOENERGETIC MODEL ***\n")
  cat("\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("--------------------------------------------------------------------------------\n\n")
  
  date_time()
  
  cat("––– Parameters:\n\n")
  
  cat("+ Animats: N =", formatC(nsim, big.mark = ","), "[individual(s) / cohort]\n")
  cat("+ Label:", ifelse(!label == "", label, "None"), "\n")
  # cat("+ Pairing:", ifelse(is.null(pair), "None", ifelse(pair$param$label == "", deparse(substitute(pair)), pair$param$label)), "\n\n")
  
  cat("+ Cohorts: \n\n")
  for (h in cohort) cat(cohorts[id==h, abb], ": ", cohorts[id==h, name], "\n", sep = "")
  cat("\n")
  
  cat("––– Scenario:\n\n")
  cat("+ Phase:", ifelse(is.null(scenario), "Baseline", switch(scenario$phase,
                                                               '1' = "Construction",
                                                               '2' = "Operation & Maintenance")), "\n")
  if(!is.null(scenario)){
    
    cat("+ Windfarms: N =", length(unique(scenario$locs$windfarm)), "\n")
    
    if(scenario$phase < 2){
    cat("\n")
    cat("––– Piling activities:\n\n")
    purrr::walk(.x = unique(scenario$locs$windfarm),
                .f = ~{
                  cat("+ Wind farm ", .x, ": ", 
                      paste0(as.character(format(as.Date(min(scenario$locs[windfarm == .x, ]$date)), "%b-%d")),
                      " (start) ----> ", 
                      as.character(format(as.Date(max(scenario$locs[windfarm == .x, ]$date)), "%b-%d")),
                      " (end) @", scenario$piles.per.day, " [pile(s) / day]\n"), sep = "")
                })
    cat("\n")
    cat("––– Piling noise:\n\n")
    cat("+ Ambient noise level:", scenario$ambient, "[dB]\n")
    cat("+ Source level:", scenario$sourceLvL, "[dB]\n")
    cat("+ Noise mitigation:", scenario$lowerdB, "[dB]\n")
    cat("+ Log range coefficient:", scenario$logrange, "\n")
    cat("+ Sound absorption coefficient:", scenario$absorb, "[dB/km]\n")
    cat("+ Duration of foraging cessation:", piling.hrs, "[hrs]\n")
    }
    
    cat("\n")
    cat("––– Vessel traffic:\n\n")
    cat("+ Routes: N =", length(unique(scenario$vessels$routeID)))
    vs <- scenario$vessels[, list(mean_speed = mean(speed_knt),
                                  min_speed = min(speed_knt),
                                  max_speed = max(speed_knt)), vesselclass]
    vs[, speed := paste0(mean_speed, " [", min_speed, "-", max_speed, "]")]
    vs <- vs[, c(1,5)]
    
    vn <- scenario$vessels[, list(median_fleet = median(Nvessels),
                                  min_fleet = min(Nvessels),
                                  max_fleet = max(Nvessels)), vesselclass]
    vn[, N := paste0(median_fleet, " [", min_fleet, "-", max_fleet, "]")]
    vn <- vn[, 5]
    
    vt <- scenario$vessels[, list(mean_trips = round(mean(roundtrips_foundation), 2),
                                  min_trips = round(min(roundtrips_foundation), 2),
                                  max_trips = round(max(roundtrips_foundation), 2)), vesselclass]
    vt[, trips := paste0(mean_trips, " [", min_trips, "-", max_trips, "]")]
    vt <- vt[, 5]
    vout <- cbind(vs, vn, vt) |> dplyr::rename(class = vesselclass)
    
    print(knitr::kable(vout, format = "simple"))
    cat("\n")
    
  }
  cat("\n––– Execution:\n\n")
  
  console(msg = "Starting up")
  
  # ............................................................
  # Set up timeline
  # ............................................................
  
  # Define time steps and sequence of density maps for each day
  # Remove leap days if necessary, as vessel strike layer built from data obtained for 2019
  date_seq <- get_dates()
  map_seq <- as.numeric(stringr::str_sub(date_seq,-2,-1))
  monthIDs <- c(map_seq[1:366], data.table::rleid(map_seq)[367:length(map_seq)])
  
  # ............................................................
  # Starvation mortality
  # ............................................................
  
  # Get coefficients of exponential decay curve
  mortality_coefs <- mortality_curve(starvation.death, starvation.onset)
  
  # ............................................................
  # Dose-response
  # ............................................................
  
  # Fit spline to dose-response data
  spline_dr <- splinefun(x = doseresponse[,2], y = doseresponse[,1])
  
  # Use inverse transform sampling to generate thresholds of response
  n.inv <- max(10000, nsim)
  median_doseresponse <- spline_dr(runif(n.inv))
  
  # doseresponse_IDs <- sample(1:10000, size = nsim * length(date_seq), replace = TRUE) |> 
  #   as.matrix(nrow = length(date_seq), ncol = nsim)
  
  # Checking that we get the dose-response back
  # y <- sapply(X = doseresponse[,1], FUN = function(x) sum(dose_n <= x) / 10000000)
  # plot(doseresponse[,1], doseresponse[,2], type = 'l')
  # lines(doseresponse[,1], y, col = "orange", lty = 2)
  
  # Generate random seeds for behavioral responses
  db.seed <- lapply(X = 1:6, FUN = function(x) sample.int(n.inv, nsim)) |> do.call(what = cbind)
  
  # ............................................................
  # Load spatial support and regions
  # ............................................................
  
  # Convert from logical to numeric; also need to flip when converting from SpatialGridDataFrame to matrix
  geomap <- t(1*raster::as.matrix(density_support))
  
  # ............................................................
  # Load density surfaces
  # ............................................................
  
  # Creates a sequence of month IDs so that the correct monthly surfaces are retrieved in the simulator
  layerIDs <- rle(map_seq)$values
  
  maps <- lapply(density_narw[layerIDs], raster::as.matrix)
  if(!all(sapply(maps, class)[1,] %in% "matrix")) stop("Cannot find input of class <matrix>.")
  
  maps.weighted.seus <- lapply(density_weighted_seus[layerIDs], raster::as.matrix)
  maps.weighted.gsl <- lapply(density_weighted_gsl[layerIDs], raster::as.matrix)
  maps.weighted.gom <- lapply(density_weighted_gom[layerIDs], raster::as.matrix)
  
  # ............................................................
  # Load prey, daylight, and stressor surfaces
  # ............................................................
  
  coords <- raster::coordinates(density_narw[[1]])
  colnames(coords) <- c("x", "y")
  map_limits <- get_limits(density_narw[[1]])
  map_resolution <- density_narw[[1]]@grid@cellsize
  
  # Daylight hours
  
  # daylightmaps <- lapply(daylight, raster::as.matrix) # Monthly
  # map_limits_daylight <- get_limits(daylight[[1]])
  # map_resolution_daylight <- daylight[[1]]@grid@cellsize
  
  # Regions
  
  regionsmap <- raster::as.matrix(regions_m) # Monthly
  map_limits_regions <- get_limits(regions_m)
  map_resolution_regions <- regions_m@grid@cellsize
  
  # Prey
  
  preymaps <- lapply(prey_layer[layerIDs], raster::as.matrix) # Monthly
  map_limits_prey <- get_limits(prey_layer[[1]])
  map_resolution_prey <- prey_layer[[1]]@grid@cellsize
  
  # Entanglement risk
  
  fishingmaps <- lapply(fishing_layer[layerIDs], raster::as.matrix) # Monthly
  map_limits_fishing <- get_limits(fishing_layer[[1]])
  map_resolution_fishing <- fishing_layer[[1]]@grid@cellsize
  
  # Vessel strike risk
  
  if(!is.null(scenario)){
    
    vesselmaps <- map_vessels(obj = scenario, z = "risk", baseline = FALSE, spgdf = TRUE, strike_scalar = strike.scalar)
    
  } else { # Baseline
    
    vesselmaps <- lapply(X = 1:12, FUN = function(x){
      tmp <- dplyr::left_join(x = vessel_grid$sf[, "cellID"], y = strike_layer[month == x], by = "cellID")
      r <- fasterize::fasterize(
        sf = tmp["tot_PLETH"],
        raster = raster::subset(vessel_grid$raster, "cellID"),
        field = "tot_PLETH")
      r <- r * strike.scalar
      r[r>1] <- 1
      as(r, "SpatialGridDataFrame")
    }) |> purrr::set_names(nm = month.abb) # Monthly
  }
  
  vesselmaps.spdf <- vesselmaps
  
  map_limits_vessels <- get_limits(vesselmaps[[1]])
  map_resolution_vessels <- vesselmaps[[1]]@grid@cellsize
  vesselmaps <- lapply(vesselmaps[layerIDs], raster::as.matrix)
  
  # Piling noise - daily rasters
  
  console(msg = "Starting up", suffix = tickmark())
  
  if(!is.null(scenario)){
    
    noisemaps <- map_noise(obj = scenario)
    
  } else {
    
    ambient.r <- as(dummy_raster(value = ambient.dB), "SpatialGridDataFrame")
    noisemaps <- lapply(date_seq, FUN = function(x) ambient.r)
    names(noisemaps) <- date_seq
  }
  
  noisemaps.spdf <- noisemaps
  
  map_limits_noise <- map_limits
  map_resolution_noise <- map_resolution
  noisemaps <- lapply(noisemaps, raster::as.matrix)
  
  # Test coarser resolution
  # preymaps <- lapply(dummy_prey, FUN = function(x){
  #   raster::as.matrix(as(raster::aggregate(raster::raster(x), fact = 2), "SpatialGridDataFrame"))
  # })
  
  # mo <- 1
  # evalEnvironment(maps[[mo]],
  #                 maps.weighted.seus[[mo]],
  #                 maps.weighted.gsl[[mo]],
  #                 preymaps[[mo]],
  #                 fishingmaps[[mo]],
  #                 vesselmaps[[mo]],
  #                 noisemaps[[mo]],
  #                 daylightmaps[[mo]],
  #                 regionsmap,
  #                 map_limits,
  #                 map_resolution,
  #                 -99.0587, -136.079, 'F')
  # 
  # raster::extract(raster::raster(density_narw[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(density_weighted_seus[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(density_weighted_gsl[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(fishing_layer[[mo]], y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(dummy_vessels[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(density_narw[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(density_narw[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  # raster::extract(raster::raster(density_narw[[mo]]), y = data.frame(x = -99.0587, y = -136.079))
  
  
  # ............................................................
  # Migratory behavior
  # ............................................................
  
  migration.df <- purrr::map(
    .x = cohort,
    .f = ~ {
      data.frame(
        seus = prob_migration(nsim, "SEUS", .x),
        gsl = prob_migration(nsim, "GSL", .x)
      )
    }
  ) |> purrr::set_names(nm = cohorts[id %in% cohort, abb])
  
  # ............................................................
  # Initial locations
  # ............................................................
  
  init_maps <- purrr::map(.x = maps, .f = ~{
    p <- as.numeric(.x)
    p[!is.finite(p)] <- 0
    p})
  
  # Simulate initial locations for latent animals
  # The warp function re-scale density surfaces within pre-defined bounds.
  # Return cell ID for each animal (columns) x month (rows)
  
  maps.seus <- purrr::map(.x = init_maps,
                          .f = ~ warp(.x, coords, lower = -5000, upper = -300, left = -1000, right = 3000))
  maps.gsl <- purrr::map(.x = init_maps, 
                         .f = ~ warp(.x, coords, lower = 1475, upper = 5000, left = -1000, right = 1300))
  maps.gom <- purrr::map(.x = init_maps, 
                         .f = ~ warp(.x, coords, lower = 900, upper = 1500, left = -1000, right = 1500))
  maps.sne <- purrr::map(.x = init_maps, 
                         .f = ~ warp(.x, coords, lower = 700, upper = 850, left = 500, right = 750))
  
  init.inds <- purrr::map(
    .x = cohorts[id %in% cohort, abb],
    .f = ~ initiate_xy(
      maplist = list("sne" = maps.sne,
                     "seus" = maps.seus,
                     "gsl" = maps.gsl,
                     "gom" = maps.gom,
                     "main" = init_maps),
      coords = coords,
      migrate = migration.df[[.x]],
      init.month = init.month,
      nsim = nsim
    )
  ) |> purrr::set_names(nm = cohorts[id %in% cohort, abb])
  
  # plot(world, axes = TRUE, col = "grey")
  # xy1 <- purrr::map(.x = init.inds[[1]][,1], .f = ~raster::xyFromCell(raster::raster(density_narw[[1]]), .x)) |> 
  #   do.call(what = rbind)
  # points(xy1)
  
  # ............................................................
  # Run simulation
  # ............................................................
  
  if(length(cohort) > 1){
    
    Ncores <- parallelly::availableCores()
    
    # Set number of cores to use
    if(is.null(n.cores)){
      if(length(cohort) <= (Ncores - 2)) n.cores <- length(cohort) else n.cores <- Ncores - 2
    } else {
      if(n.cores > Ncores)
        stop("Insufficient cores.")
    }
    
    console(msg = paste0("Initializing parallel computing (", n.cores, " cores)"))
    
    # Initializing workers
    cl <- snow::makeSOCKcluster(n.cores)
    doSNOW::registerDoSNOW(cl)
    on.exit(snow::stopCluster(cl))
    # parallel::clusterSetRNGStream(cl, rdseed)
    
    console(msg = paste0("Initializing parallel computing (", n.cores, " cores)"), suffix = tickmark())
    
    console(msg = "Running simulations", suffix = "")
    
    # Compute estimates
    out <- foreach::foreach(i = seq_along(cohort), 
                            .packages = c("Rcpp", "narwind")) %dopar% {
                              
                              # Set environments
                              .GlobalEnv$coords <- coords
                              
                              # # Match random seed if baseline scenario is provided
                              # if(!is.null(attr(pair,"seed"))){
                              #   assign(".Random.seed",  attr(pair, "seed"))
                              # }
                              
                              # Only used for package development purposes
                              if(pkgdev) Rcpp::sourceCpp("src/simtools.cpp")
                              
                              NARW_simulator(
                                cohortID = cohort[i],
                                seus = migration.df[[i]][,1],
                                gsl = migration.df[[i]][,2],
                                densities = maps,
                                densities_seus = maps.weighted.seus,
                                densities_gsl = maps.weighted.gsl,
                                densities_gom = maps.weighted.gom,
                                densitySeq = map_seq,
                                prey = preymaps,
                                fishing = fishingmaps,
                                vessels = vesselmaps,
                                noise = noisemaps,
                                doseresp_seed = db.seed,
                                doseresp = median_doseresponse,
                                regions = regionsmap,
                                M = n.prop, 
                                # stepsize = step.size,
                                xinit = matrix(data = coords[init.inds[[i]], 'x'], nrow = nrow(init.inds[[i]])),
                                yinit = matrix(data = coords[init.inds[[i]], 'y'], nrow = nrow(init.inds[[i]])),
                                support = geomap,
                                limits = map_limits,
                                limits_regions = map_limits_regions,
                                limits_prey = map_limits_prey,
                                limits_fishing = map_limits_fishing,
                                limits_vessels = map_limits_vessels,
                                limits_noise = map_limits_noise,
                                resolution = map_resolution,
                                resolution_regions = map_resolution_regions,
                                resolution_prey = map_resolution_prey,
                                resolution_fishing = map_resolution_fishing,
                                resolution_vessels = map_resolution_vessels,
                                resolution_noise = map_resolution_noise,
                                stressors = stressors,
                                growth = growth,
                                prey_scale = prey.scalar,
                                starvation_coefs = mortality_coefs,
                                starvation_death = starvation.death,
                                starvation_onset = starvation.onset,
                                nursing_cessation = cease.nursing,
                                piling_hrs = piling.hrs,
                                progress = FALSE
                              )
                            } # End foreach loop
    
    names(out) <- cohorts[id %in% cohort, abb]
    # close(pb) # Close progress bar
    
  } else {
    
    # d <- maps
    # if(cohort == 5) d <- maps.weighted else d <- maps
    # if(cohort > 0) d <- maps.weighted else d <- maps
    
    # Simulated animals begin at the position of their corresponding latent animal for the starting month
    # This can be checked by running the below code
    # xpos <- apply(init.inds[[1]], 2, FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
    # ypos <- apply(init.inds[[1]], 2, FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
    
    out <- list(NARW_simulator(
      cohortID = cohort,
      seus = migration.df[[1]][, 1],
      gsl = migration.df[[1]][, 2],
      densities = maps,
      densities_seus = maps.weighted.seus,
      densities_gsl = maps.weighted.gsl,
      densities_gom = maps.weighted.gom,
      densitySeq = monthIDs,
      prey = preymaps,
      fishing = fishingmaps,
      vessels = vesselmaps,
      noise = noisemaps,
      doseresp_seed = db.seed,
      doseresp = median_doseresponse,
      regions = regionsmap,
      M = n.prop,
      # stepsize = step.size,
      xinit = matrix(data = coords[init.inds[[1]], "x"], nrow = nrow(init.inds[[1]])),
      yinit = matrix(data = coords[init.inds[[1]], "y"], nrow = nrow(init.inds[[1]])),
      support = geomap,
      limits = map_limits,
      limits_regions = map_limits_regions,
      limits_prey = map_limits_prey,
      limits_fishing = map_limits_fishing,
      limits_vessels = map_limits_vessels,
      limits_noise = map_limits_noise,
      resolution = map_resolution,
      resolution_regions = map_resolution_regions,
      resolution_prey = map_resolution_prey,
      resolution_fishing = map_resolution_fishing,
      resolution_vessels = map_resolution_vessels,
      resolution_noise = map_resolution_noise,
      stressors = stressors,
      growth = growth,
      prey_scale = prey.scalar,
      starvation_coefs = mortality_coefs,
      starvation_death = starvation.death,
      starvation_onset = starvation.onset,
      nursing_cessation = cease.nursing,
      piling_hrs = piling.hrs,
      progress = ifelse(progress, TRUE, FALSE)
    ))
    
  }
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  if(length(cohort) == 1) cat("\n")
  console(msg = "Wrapping up")
  
  # Transpose the output arrays to have data for each day as rows
  out.t <- transpose_array(input = out, cohortID = cohort, dates = date_seq) |> 
    purrr::set_names(nm = cohorts[id %in% cohort, abb])
  
  # Consolidate the data
  for (k in seq_along(cohort)) {
    if (cohort[k] %in% 4:5){
      out.t[[k]][["attrib"]] <- abind::abind(out.t[[k]][["attrib"]][[1]], out.t[[k]][["attrib"]][[2]], along = 2)
    }
    if (cohort[k] == 5){
      out.t[[k]][["E"]] <- abind::abind(out.t[[k]][["E"]][[1]], out.t[[k]][["E"]][[2]], along = 2)
      out.t[[k]][["kj"]] <- abind::abind(out.t[[k]][["kj"]][[1]], out.t[[k]][["kj"]][[2]], along = 2)
    } 
  }
  
  out.dt <- purrr::map(.x = out.t,
                       .f = ~abind::abind(
                         .x[["locs"]],
                         .x[["attrib"]],
                         .x[["stress"]],
                         .x[["E"]],
                         .x[["kj"]],
                         .x[["activ"]],
                         along = 2))
  
  outsim <- list()  
  outsim[["sim"]] <- consolidate(dtl = out.dt, 
                                 nsim = nsim,
                                 cnames = cohorts[id %in% cohort, name], 
                                 dates = date_seq, 
                                 months = monthIDs)
  
  outsim[["sim"]] <- purrr::map(.x = outsim[["sim"]], .f = ~{
    oncalvgrounds <- .x[day>0, list(entrySEUS = which.max(inseus),
                                    exitSEUS = 367-which.max(rev(inseus))), whale]
    oncalvgrounds[, tSEUS:= exitSEUS - entrySEUS]
    .x[oncalvgrounds, on = "whale", entrySEUS := i.entrySEUS]
    .x[oncalvgrounds, on = "whale", exitSEUS := i.exitSEUS]
    .x[oncalvgrounds, on = "whale", tSEUS := i.tSEUS]
  })
  
  # ............................................................
  # MORTALITY AND BODY CONDITION
  # ............................................................
  
  for(k in cohort){
    outsim[["sim"]][[cohorts[id == k, abb]]][, date_died:= date_died[which.max(date_died)], whale]
    if (k == 5) {
      outsim[["sim"]][[cohorts[id == k, abb]]][, date_died_calf:= date_died_calf[which.max(date_died_calf)], whale]
      outsim[["sim"]][[cohorts[id == k, abb]]][, dob:= dob[which.max(dob)], whale]
    }
  }
  
  gam.dt <- purrr::map(.x = seq_along(cohort), .f = ~ {
    
    dt_adjuv <-
      
      purrr::reduce(list(
        
        # Starting condition
        outsim[["sim"]][[.x]][day == 93, list(cohort = cohort, cohort_name = cohort_name, start_bc = bc), whale],
        
        # Finishing conditions
        outsim[["sim"]][[.x]][day == length(date_seq) - 1, list(end_bc = bc, alive, abort, strike, starve, died), whale],
        
        # Southward migration
        outsim[["sim"]][[.x]][, list(seus = max(inseus)), whale],
        
        # Date, day, region and coordinates of when death occurred if relevant
        # Add +1 to .SD call as need row number rather than day nunber
        outsim[["sim"]][[.x]][, .SD[unique(date_died) + 1], .SDcols = c("date", "day", "easting", "northing", "region"), whale]
      ), dplyr::left_join, by = "whale") |>
      dplyr::mutate(born = 1)
    
    dt_adjuv$event <- ifelse(dt_adjuv$alive == 1, "none", "death")
    
    if (cohort[.x] == 5) {
      
      calf_columns <- c(
        "cohort_calf",
        "cohort_name",
        "alive_calf",
        "born",
        "date",
        "day",
        "easting",
        "northing",
        "region"
      )
      
      dt_calves <- purrr::reduce(list(
        
        dplyr::bind_rows(
          
          # Births
          outsim[["sim"]][[.x]][, .SD[unique(dob) + 1], .SDcols = calf_columns, whale] |>
            dplyr::mutate(event = "birth") |>
            dplyr::rename(cohort = cohort_calf, alive = alive_calf),
          
          # Deaths
          outsim[["sim"]][[.x]][, .SD[unique(date_died_calf) + 1], .SDcols = calf_columns, whale] |>
            dplyr::mutate(event = "death") |>
            dplyr::rename(cohort = cohort_calf, alive = alive_calf) |>
            dplyr::mutate(
              date = replace(date, alive == 1, NA),
              day = replace(day, alive == 1, NA),
              easting = replace(easting, alive == 1, NA),
              northing = replace(northing, alive == 1, NA),
              region = replace(region, alive == 1, NA)
            )
        ) |> dplyr::select(-alive),
        
        # Body condition
        outsim[["sim"]][[.x]][, .SD[unique(dob) + 1], .SDcols = "bc_calf", whale] |> dplyr::rename(start_bc = bc_calf),
        
        # Finishing conditions
        outsim[["sim"]][[.x]][day == length(date_seq) - 1, 
                              list(bc_calf, alive_calf, abort, strike_calf, starve_calf, died_calf), whale] |>
          dplyr::rename(alive = alive_calf, end_bc = bc_calf, 
                        strike = strike_calf, starve = starve_calf, died = died_calf)
      ),
      dplyr::left_join,
      by = "whale"
      ) |> dplyr::filter(born > 0)
    } else {
      dt_calves <- data.table::data.table()
    }
    dt_calves$seus <- 0
    dplyr::bind_rows(dt_adjuv, dt_calves)
  }) |> data.table::rbindlist()
  
  
  
  # Add full cohort information
  gam.dt <- data.table::merge.data.table(x = gam.dt, y = cohorts, by.x = "cohort", by.y = "id", all.x = TRUE)
  
  # Reorder columns
  data.table::setcolorder(gam.dt, c("cohort", "name", "cohort_name", "class", "abb", "whale", "alive", "start_bc", "end_bc", "starve", "strike", "date", "day", "easting", "northing", "region"))
  
  # Update values
  gam.dt[cohort == 0, cohort_name:= "Calves (male, female)"]
  gam.dt[event == "death", cause_death:= ifelse(strike == 1, "strike",
                                                ifelse(starve == 1, "starve",
                                                       ifelse(died == 1, "natural", "other")))]
  
  # Update cause of death for calves that died as a result of their mothers dying
  mothers <- gam.dt[
    alive == 0 & cohort == 5,
    list(
      cause_mother = paste0(strike, starve),
      date_mother = date
    ), whale
  ]
  mothers[, cause_mother := switchv(cause_mother,
                                    "00" = "natural",
                                    "01" = "starve",
                                    "10" = "strike"
  )]
  
  gam.dt <- dplyr::left_join(gam.dt, mothers, by = "whale")
  
  gam.dt[alive == 0 & cohort == 0 & event == "death", 
         cause_death:=ifelse(!is.na(date_mother), 
                             ifelse(date_mother == date, paste0(cause_mother, " (female)"), cause_death), cause_death)]
  
  # if(5 %in% cohort){
  # 
  # gam.dt <- split(gam.dt, f = factor(gam.dt$cohort))
  # 
  # gam.dt[["0"]][event == "death"] <- gam.dt[["5"]][, c("whale", "cause_death")] |>
  # dplyr::rename(
  #   cause_mother = cause_death
  # ) |>
  # dplyr::left_join(x = gam.dt[["0"]][event == "death", ], by = "whale") |>
  # dplyr::mutate(cause_death = dplyr::case_when(
  #   !is.na(cause_mother) ~ cause_mother,
  #   .default = cause_death
  # )) |>
  # dplyr::select(-cause_mother)
  # 
  # gam.dt <- data.table::rbindlist(gam.dt)
  # 
  # }
  
  # gam.dt <- suppressWarnings(gam.dt |> dplyr::mutate(date = lubridate::as_date(date)))
  
  # ............................................................
  # Deaths
  # ............................................................
  
  locs.dead <- gam.dt[event == "death" & alive == 0 & born == 1]
  data.table::setorder(locs.dead, -cohort, whale)
  # if (nrow(locs.dead) == 0) locs.dead <- list(NULL)
  
  # ............................................................
  # Births
  # ............................................................
  
  if (5 %in% cohort) {
    locs.birth <- gam.dt[event == "birth" & cohort == 0 & born == 1]
  } else {
    locs.birth <- list(NULL)
  }
  
  # ............................................................
  # Abortions
  # ............................................................
  
  if(4 %in% cohort){
    abort.rate <- outsim[["sim"]][[cohorts[id == 4, abb]]][day > 0 & alive == 1, .(abort = max(abort)), whale]
  } else {
    abort.rate <- list(NULL)
  }
  
  console(msg = "Wrapping up", suffix = tickmark())
  
  #'...............................
  # Fit GAM models ----
  #'................................
  
  # https://stats.stackexchange.com/questions/403772/different-ways-of-modelling-interactions-between-continuous-and-categorical-pred
  
  console(msg = "Fitting terminal functions")
  
  surv.fit <- suppressWarnings(
    tryCatch(
      {
        mgcv::gam(alive ~ factor(cohort) + s(start_bc, by = factor(cohort), k = 4, bs = "tp"),
                  data = gam.dt[born == 1,],
                  method = "REML",
                  gamma = 1.4,
                  family = binomial(link = "logit")
        )
      },
      error = function(cond) {
        return(NA)
      }
    )
  )
  
  # Failsafes in the event that all animals in a cohort are dead
  # n.dead <- gam.dt[, .(alive = sum(alive)), cohort]
  # gam.alive <- gam.dt[gam.dt$alive == 1, ]
  # gam.dead <- gam.dt[cohort %in% n.dead[alive == 0, cohort]]
  # gam.dead[,start_bc:=0]
  # gam.dead[,end_bc:=0]
  # gam.bc <- rbindlist(list(gam.alive, gam.dead))
  
  bc.fit <- suppressWarnings(
    tryCatch(
      {
        mgcv::gam(end_bc ~ factor(cohort) + s(start_bc, by = factor(cohort), k = 4, bs = "tp"),
                  data = gam.dt[alive == 1,],
                  method = "REML",
                  gamma = 1.4,
                  family = betar(link = "logit"))
      },
      error = function(cond) {
        return(NA)
      }
    )
  )
  
  console(msg = "Fitting terminal functions", suffix = tickmark())
  
  # ............................................................
  # Add parameters to list output
  # ............................................................
  
  outsim$gam <- list(fit = list(surv = surv.fit, bc = bc.fit), 
                     dat = gam.dt)
  
  outsim[["dead"]] <- locs.dead
  outsim[["birth"]] <- locs.birth
  outsim[["abort"]] <- abort.rate
  # outsim[["max_bc"]] <- max_bc
  
  outsim$param <- list(nsim = nsim,
                       label = label,
                       cohort = cohort,
                       cohorts = cohorts,
                       date_seq = date_seq,
                       starvation = list(death = starvation.death,
                                         onset = starvation.onset,
                                         coefs = mortality_coefs),
                       nurse_cease = cease.nursing,
                       strike.scalar = strike.scalar,
                       prey.scalar = prey.scalar,
                       # step = step.size,
                       n.prop = n.prop,
                       ambient.dB = ambient.dB)
  
  outsim$stressors <- list(strike = vesselmaps.spdf,
                           noise = noisemaps.spdf)
  
  outsim$init <- list(month = init.month, 
                      xy = init.inds, 
                      dest = migration.df)
  
  outsim$run <- run_time
  if(is.null(scenario)){
    outsim$scenario <- get_scenarios(scenario = 0)
  } else {
    outsim$scenario <- scenario
  }
  
  # Add custom class and random seed attribute
  class(outsim) <- c("narwsim", class(outsim))
  class(outsim$init$xy) <- c("xyinits", class(outsim$init$xy))
  # attr(outsim, "seed") <- rdseed
  
  # Print run time
  cat("---------------------------------\n")
  cat(paste0("Done | Time elapsed: ", run_time), "\n")
  
  # Garbage collection
  gc()
  
  return(outsim)
}
