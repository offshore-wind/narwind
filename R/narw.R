#' Run the bioenergetic model
#'
#' Simulate right whale movements and behavior across a calendar year.
#' @export
#' @param nsim Number of simulated animals
#' @param n.prop Number of proposals used in the importance sampler (Michelot, 2019)
#' @import data.table
#' @importFrom foreach `%dopar%`
#' @importFrom doParallel registerDoParallel
#' @author Phil J. Bouchet
#' @seealso \code{\link{initialize}}
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' narw(10)
#' }

narw <- function(nsim = 1e3,
                 scenario = NULL,
                 n.cores = NULL,
                 progress = TRUE,
                 ...){
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values for optional arguments 
  init.month <- 10 # Start simulation in October
  cohortID <- 1:6
  stressors <- TRUE
  growth <- TRUE
  
  # Starvation threshold expressed as relative blubber mass
  # As per Pirotta et al. (2018) - to test extreme conditions of leanness
  starvation <- 0.05
  
  # Default values
  if(length(args) > 0) {
    
    if("cohortID" %in% names(args)) cohortID <- args[["cohortID"]] 
    if("stressors" %in% names(args)) stressors <- args[["stressors"]] 
    if("growth" %in% names(args)) growth <- args[["growth"]] 
    if("init.month" %in% names(args)) init.month <- args[["init.month"]]
    if("starvation" %in% names(args)) starvation <- args[["starvation"]]
  }
  
  # nsim = 10
  # cohortID = 5
  # scenario = NULL
  # init.month = 2
  # stressors = FALSE
  # growth = TRUE
  # n.cores = NULL
  # progress = TRUE

  cohortID <- cohortID[cohortID > 0]
  
  # Step size for movement model
  step.size = 120
  
  # Number of proposals when sampling new location
  # In line with sample sizes used in case studies from Michelot (2020)
  n.prop = 50
  
  # if(!init.month %in% 1:12) stop("<init.month> must be an integer between 1 and 12")
  # if(!"init.model" %in% ls(envir = .GlobalEnv)) stop("The model must first be initialized using <load_model>")
  if(any(!cohortID %in% 1:6)) stop("Unrecognized cohort")
  
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n")
  cat("\n")
  cat("     NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("                --- BIOENERGETIC MODEL ---\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n\n")
  
  cat("Starting up ...\n")
  
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

  # ............................................................
  # Dose-response
  # ............................................................
  
  spline_dr <- splinefun(x = doseresponse[,2], y = doseresponse[,1])
  
  # Use inverse transform sampling to generate 10,000 thresholds of response
  median_doseresponse <- spline_dr(runif(10000))
  
  # Checking that we get the dose-response back
  # y <- sapply(X = doseresponse[,1], FUN = function(x) sum(dose_n <= x) / 10000000)
  # plot(doseresponse[,1], doseresponse[,2], type = 'l')
  # lines(doseresponse[,1], y, col = "orange", lty = 2)
  
  # ............................................................
  # Load spatial support and regions
  # ............................................................
  
  # Convert from logical to numeric; also need to flip when converting from SpatialGridDataFrame to matrix
  geomap <- t(1*raster::as.matrix(density_support))
  
  # ............................................................
  # Load density surfaces
  # ............................................................
  
  maps <- lapply(density_narw, raster::as.matrix)
  if(!all(sapply(maps, class)[1,] %in% "matrix")) stop("Cannot find input of class <matrix>.")
  
  maps.weighted.seus <- lapply(density_weighted_seus, raster::as.matrix)
  maps.weighted.gsl <- lapply(density_weighted_gsl, raster::as.matrix)
  
  # ............................................................
  # Load prey, daylight, and stressor surfaces
  # ............................................................
  
  preymaps <- lapply(dummy_prey, raster::as.matrix)
  fishingmaps <- lapply(fishing_layer, raster::as.matrix)

  if(is.null(scenario)){
    vesselmaps <- lapply(dummy_vessels, raster::as.matrix)
    noisemaps <- lapply(dummy_noise, FUN = function(x) {
      n <- raster::as.matrix(x)
      n[!is.na(n)] <- 60
      n})
  } else {
    vesselmaps <- lapply(dummy_vessels, raster::as.matrix)
    noisemaps <- lapply(dummy_noise, raster::as.matrix)
  }
  
  daylightmaps <- lapply(daylight, raster::as.matrix)
  regionsmap <- raster::as.matrix(regions_m)
  
  
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
  # Reference information for maps
  # ............................................................
  
  coords <- sp::coordinates(density_narw[[1]])
  colnames(coords) <- c("x", "y")
  map_limits <- c(range(coords[,1]), range(coords[,2]))
  map_resolution <- density_narw[[1]]@grid@cellsize
  
  # ............................................................
  # Set up simulation parameters
  # ............................................................
  
  # Define time steps and sequence of density maps for each day
  # Use a non-leap year like 2021
  start.date <- paste0(2021, stringr::str_pad(init.month, width = 2, pad = "0"), "01")
  init.date <- paste0(2021, "-", stringr::str_pad(init.month, width = 2, pad = "0"), "-", "00")
  date_seq <- c(init.date, as.character(seq(from = lubridate::ymd(start.date), by = 'day', length.out = 365)))
  map_seq <- c(init.month, lubridate::month(date_seq[2:366]))
  
  # Migratory destinations
  migration.df <- purrr::map(.x = cohortID, 
                             .f = ~{data.frame(seus = prob_migration(nsim, "SEUS", .x),
                             gsl = prob_migration(nsim, "GSL", .x))}) |> purrr::set_names(nm = cohorts[id %in% cohortID, abb])
  
  init_maps <- purrr::map(.x = maps, .f = ~{
    p <- as.numeric(.x)
    p[!is.finite(p)] <- 0
    p})
  
  init_maps.seus <- purrr::map(.x = init_maps, .f = ~warp(.x, coords, -5000, -300, -1000, 3000))
  init_maps.gsl <- purrr::map(.x = init_maps, .f = ~warp(.x, coords, 1475, 5000, -1000, 1300))
  init_maps.main <- purrr::map(.x = init_maps, .f = ~warp(.x, coords))
  
  # Simulate initial locations for latent animals
  # Return cell ID for each animal (columns) x month (rows)
  init.inds <- purrr::map(
    .x = names(migration.df),
    .f = ~ initiate_xy(
      maplist = list(init_maps.seus, init_maps.gsl, init_maps.main, init_maps),
      coords = coords,
      migrate = migration.df[[.x]],
      init.month = init.month,
      nsim = nsim
    )
  ) |> purrr::set_names(nm = cohorts[id %in% cohortID, abb])
  
  
  # test_inits <- function(ind){
  #   xinit = matrix(data = coords[init.inds[[1]], 'x'], nrow = nrow(init.inds[[1]]))
  #   yinit = matrix(data = coords[init.inds[[1]], 'y'], nrow = nrow(init.inds[[1]]))
  #   plot(regions)
  #   plot(world, add = T, col = "grey")
  #   points(xinit[,ind], yinit[,ind], pch = 16, col = c(rep("#cb4154", 2), rep("#4199cb", 3), rep("orange",4), "#41cbb8", rep("#cb4154", 2)))
  #   legend("bottomright", legend=month.abb,
  #          fill=c(rep("#cb4154", 2), rep("#4199cb", 3), rep("orange",4), "#41cbb8", rep("#cb4154", 2))) 
  # }


  # init.inds <- purrr::map(.x = cohortID, 
  #   .f = ~init_xy(maps = maps, 
  #                 coords = coords, 
  #                 cohort.id = .x, 
  #                 nsim = nsim)) |> 
  #   purrr::set_names(nm = cohorts[id %in% cohortID, abb])
  
  # ............................................................
  # Run simulation
  # ............................................................

  if(length(cohortID) > 1){
    
    Ncores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
    
    # Set number of cores to use
    if(is.null(n.cores)){
      if(length(cohortID) <= (Ncores - 2)) n.cores <- length(cohortID) else n.cores <- Ncores - 2
    } else {
      if(n.cores > Ncores)
        stop("Insufficient cores.")
    }
    
    cat("Initializing parallel computing (", n.cores, " cores) ...\n", sep = "")
    
    isDarwin <- Sys.info()[['sysname']] == "Darwin"
    isWindows <- Sys.info()[['sysname']] == "Windows"
    
    # Trying snow
    cl <- snow::makeSOCKcluster(n.cores)
    doSNOW::registerDoSNOW(cl)
    on.exit(snow::stopCluster(cl))
    
    sink(tempfile())
    pb <- utils::txtProgressBar(max = length(cohortID), style = 3)
    sink()
    pbar <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = pbar)
    
    # # Create a parallel cluster and register the backend
    # cl <- suppressMessages(parallel::makeCluster(n.cores))
    # 
    # # After the function is run, shutdown the cluster.
    # on.exit(parallel::stopCluster(cl))
    # 
    # # Register parallel backend
    # doParallel::registerDoParallel(cl)   # Modify with any do*::registerDo*()
    
    cat("Running simulations ...")
    if(progress){cat("\n")}
    
    start.time <- Sys.time()
    
    # hashID <- hashID
    
    # Compute estimates
    out <- foreach::foreach(i = seq_along(cohortID), 
                            .options.snow = opts,
                            .packages = c("Rcpp"), 
                            .noexport = c("NARW_simulator")) %dopar% {
                              
                              # Set environments
                              .GlobalEnv$coords <- coords

                              Rcpp::sourceCpp("src/simtools.cpp") # Comment out for package use

                              # d <- maps
                              # if(cohortID[i] == 5) d <- maps.weighted else d <- maps
                              # if(cohortID[i] >0) d <- maps.weighted else d <- maps
                              
                              NARW_simulator(
                                cohortID = cohortID[i],
                                seus = migration.df[[i]][,1],
                                gsl = migration.df[[i]][,2],
                                densities = maps,
                                densities_seus = maps.weighted.seus,
                                densities_gsl = maps.weighted.gsl,
                                densitySeq = map_seq,
                                latentDensitySeq = 1:12,
                                prey = preymaps,
                                fishing = fishingmaps,
                                vessels = vesselmaps,
                                noise = noisemaps,
                                doseresp = median_doseresponse,
                                daylight = daylightmaps,
                                regions = regionsmap,
                                M = n.prop, # No.proposals in the importance sampler (Michelot, 2019)
                                stepsize = step.size / 2, # Movement model involves two half-steps (Michelot, 2019, 2020)
                                xinit = matrix(data = coords[init.inds[[i]], 'x'], nrow = nrow(init.inds[[i]])),
                                yinit = matrix(data = coords[init.inds[[i]], 'y'], nrow = nrow(init.inds[[i]])),
                                support = geomap,
                                limits = map_limits,
                                resolution = map_resolution,
                                stressors = stressors,
                                growth = growth,
                                starvation = starvation,
                                progress = FALSE)} # End foreach loop
    
    names(out) <- cohorts[id %in% cohortID, abb]
    
    close(pb) # Close progress bar
    
  } else {
    
    cat("Running simulations ...\n")

    # d <- maps
    # if(cohortID == 5) d <- maps.weighted else d <- maps
    # if(cohortID > 0) d <- maps.weighted else d <- maps
    
    # Simulated animals begin at the position of their corresponding latent animal for the starting month
    # This can be checked by running the below code
    # xpos <- apply(init.inds[[1]], 2, FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
    # ypos <- apply(init.inds[[1]], 2, FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
    
    start.time <- Sys.time()
    
    out <- list(NARW_simulator(
      cohortID = cohortID,
      seus = migration.df[[1]][,1],
      gsl = migration.df[[1]][,2],
      densities = maps,
      densities_seus = maps.weighted.seus,
      densities_gsl = maps.weighted.gsl,
      densitySeq = map_seq,
      latentDensitySeq = seq_along(density_narw),
      prey = preymaps,
      fishing = fishingmaps,
      vessels = vesselmaps,
      noise = noisemaps,
      doseresp = median_doseresponse,
      daylight = daylightmaps,
      regions = regionsmap,
      M = n.prop, # Number of proposals used in the importance sampler for movement (Michelot, 2019)
      stepsize = step.size / 2, # Movement framework involves two half-steps (Michelot 2019, 2020)
      xinit = matrix(data = coords[init.inds[[1]], 'x'], nrow = nrow(init.inds[[1]])),
      yinit = matrix(data = coords[init.inds[[1]], 'y'], nrow = nrow(init.inds[[1]])),
      support = geomap,
      limits = map_limits,
      resolution = map_resolution,
      stressors = stressors,
      growth = growth,
      starvation = starvation,
      progress = ifelse(progress, TRUE, FALSE)))
    
  }
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  # Transpose the output arrays to have data for each day as rows
  out.t <- transpose_array(input = out, cohortID = cohortID, dates = date_seq) |> 
    purrr::set_names(nm = cohorts[id %in% cohortID, abb])
  
  # Consolidate the data
  for (k in seq_along(cohortID)) {
    if (cohortID[k] %in% 4:5){
      out.t[[k]][["attrib"]] <- abind::abind(out.t[[k]][["attrib"]][[1]], out.t[[k]][["attrib"]][[2]], along = 2)
    }
    if (cohortID[k] == 5){
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
  outsim[["sim"]] <- consolidate(out.dt, nsim, cohorts[id %in% cohortID, name], date_seq, map_seq)
  
  # ............................................................
  # MORTALITY AND BODY CONDITION
  # ............................................................

  gam.dt <- purrr::map(.x = seq_along(cohortID), .f = ~ {

   outsim[["sim"]][[.x]][, death:= `if`(all(alive==0), 1, max(which(alive==1))), whale]
    
   dt_adjuv <- purrr::reduce(list(
     # Starting condition
     outsim[["sim"]][[.x]][day == 0, list(cohort = cohort, cohort_name = cohort_name, start_bc = bc), whale],
     # Finishing conditions
     outsim[["sim"]][[.x]][day == 365, list(end_bc = bc, alive, abort, strike, starve), whale],
     # Date, day, region and coordinates of when death occurred
     outsim[["sim"]][[.x]][, .SD[death], .SDcols = c("date", "day", "easting", "northing", "region"), whale] |> dplyr::distinct()
    ), dplyr::left_join, by = "whale") |> 
     dplyr::mutate(born = 1)

  if (cohortID[.x] == 5) {
    
    outsim[["sim"]][[.x]][, dob:= `if`(all(born==0), 1, min(which(born==1))), whale]
    outsim[["sim"]][[.x]][, death_calf:= `if`(all(alive_calf==0), dob, max(which(alive_calf==1))), whale]
    
    dt_calves <- purrr::reduce(list(
      # Starting condition
      outsim[["sim"]][[.x]][, .SD[dob], .SDcols = c("cohort_calf", "cohort_name", "bc_calf", "born"), whale] |>
        dplyr::distinct() |> 
        dplyr::rename(start_bc = bc_calf, cohort = cohort_calf),
      # Finishing conditions
      outsim[["sim"]][[.x]][day == 365, list(end_bc = bc, alive, abort, strike, starve), whale],
      # Date, day, region and coordinates of when death occurred
      outsim[["sim"]][[.x]][, .SD[death_calf], .SDcols = c("date", "day", "easting", "northing", "region"), whale] |> dplyr::distinct()
    ), dplyr::left_join, by = "whale") |> 
      dplyr::relocate(born, .after = dplyr::last_col())
    
  } else {
    dt_calves <- data.table::data.table()
  }
  dplyr::bind_rows(dt_adjuv, dt_calves)
}) |> data.table::rbindlist()
  
  # Add full cohort information
  gam.dt <- data.table::merge.data.table(x = gam.dt, y = cohorts, by.x = "cohort", by.y = "id", all.x = TRUE)
  data.table::setcolorder(gam.dt,
  c("cohort", "cohort_name", "name", "class", "abb", "whale", "alive", "start_bc", "end_bc", "starve", "strike", "date", "day", "easting", "northing", "region"))
  gam.dt[cohort == 0, cohort_name:= "Calves (male, female)"]
  gam.dt[,cause_death:= ifelse(strike == 1, "strike", ifelse(starve == 1, "starve", "NA"))]
  
  # ............................................................
  # Deaths
  # ............................................................
  
  locs.dead <- gam.dt[alive == 0 & born == 1]
  data.table::setorder(locs.dead, -cohort, whale)

  # ............................................................
  # Births
  # ............................................................
  
  if (5 %in% cohortID) {
   locs.birth <- gam.dt[cohort == 0 & born == 1]
  #   locs.birth <-
  #     outsim[["sim"]][[cohorts[id == 5, abb]]][day > 0, .SD[ifelse(which.max(born) == 1, 0, which.max(born))],
  #       .SDcols = c("date", "day", "cohort_calf", "easting", "northing", "region"), whale
  #     ] |>
  #     dplyr::mutate(event = "birth", date = lubridate::as_date(date))
    if (nrow(locs.birth) == 0) locs.birth <- list(NULL)
  } else {
    locs.birth <- list(NULL)
  }

  # ............................................................
  # Abortions
  # ............................................................
  
  if(4 %in% cohortID){
    abort.rate <- outsim[["sim"]][[cohorts[id == 4, abb]]][day > 0 & alive == 1, .(abort = max(abort)), whale]
  } else {
    abort.rate <- list(NULL)
  }
  
  # ............................................................
  # Fit GAM models
  # ............................................................
  
  # https://stats.stackexchange.com/questions/403772/different-ways-of-modelling-interactions-between-continuous-and-categorical-pred
  
  # print(gam.dt)
  
  # Standard GAM formulation
  # surv.fit <- suppressWarnings(
  #   tryCatch(
  #     {mgcv::gam(alive ~ factor(cohort) + s(start_bc, by = factor(cohort), bs = "tp"),
  #         data = gam.dt,
  #         method = "REML",
  #         family = binomial("logit"))},
  #     error = function(cond) {
  #       return(NA)
  #     }))
  
  # Shape-constrained GAM formulation
  surv.fit <- 
    suppressWarnings(
    tryCatch(
      {
        scam::scam(alive ~ factor(cohort) + s(start_bc, by = factor(cohort), bs = "mpi"),
                 data = gam.dt,
                 family = binomial("logit"))
        },
      error = function(cond) {
        return(NA)
      }))
  
  # In case all animals in a cohort are dead
  n.dead <- gam.dt[, .(alive = sum(alive)), cohort]

  gam.alive <- gam.dt[gam.dt$alive == 1, ]
  gam.dead <- gam.dt[cohort %in% n.dead[alive == 0, cohort]]
  gam.dead[,start_bc:=0.001]
  gam.dead[,end_bc:=0.001]
  
  gam.alive <- rbindlist(list(gam.alive, gam.dead))
                              
  # Standard GAM formulation
  # bc.fit <- suppressWarnings(
  #   tryCatch(
  #     {mgcv::gam(end_bc ~ factor(cohort) + s(start_bc, by = factor(cohort), bs = "tp", k = 3),
  #                data = gam.alive,
  #                method = "REML",
  #                family = "Gamma")},
  #     error = function(cond) {
  #       return(NA)
  #     }))
  
  # Shape-constrained GAM formulation
  bc.fit <- 
    suppressWarnings(
    tryCatch(
      {
        scam::scam(end_bc ~ factor(cohort) + s(start_bc, by = factor(cohort), bs = "mpi"),
                 data = gam.alive,
                 family = "Gamma")
        },
      error = function(cond) {
        return(NA)
      }))
  
  # Spline interpolation of the fitted GAM relationships
  mbc_preds <- splinefun(x = gam_gest$mod[,2], y = gam_gest$fitted.values)
  
  surv_preds <- suppressWarnings(
    tryCatch(
      {purrr::map(.x = unique(cohorts$id), 
                  .f = ~{
                    x <- seq(0,1,by = 0.01)
                    p <- tryCatch(
                      {predict(object = surv.fit, 
                                 newdata = data.frame(start_bc = x, cohort = .x), 
                                 type = "response", 
                                 newdata.guaranteed = TRUE)
                      }, error = function(cond){return(rep(0, length(x)))})
                    splinefun(x = x, y = p)}) |> 
          purrr::set_names(nm = unique(cohorts$id))},
      error = function(cond) {
        return(NA)
      }))
  
  bc_preds <- suppressWarnings(
    tryCatch(
      {purrr::map(.x = unique(cohorts$id), 
                  .f = ~{
                    x <- seq(0,1,by = 0.01)
                    p <- tryCatch(
                      {predict(object = bc.fit, 
                                 newdata = data.frame(start_bc = x, cohort = .x), 
                                 type = "response", 
                                 newdata.guaranteed = TRUE)
                      }, error = function(cond){return(rep(0, length(x)))})
                    splinefun(x = x, y = p)}) |> 
          purrr::set_names(nm = unique(cohorts$id))},
      error = function(cond) {
        return(NA)
      }))
  
  # Separate GAMS for each cohort
  
  # gam.dtl <- split(gam.dt, f = factor(gam.dt$cohort))
  # 
  # gamfit <- suppressWarnings(purrr::map(.x = unique(cohorts$id), .f = ~ {
  #   gam.dat <- gam.dtl[[as.character(.x)]]
  #   gam.out <- list(surv = NA, bc = NA)
  #   surv.fit <- tryCatch({
  #       surv = mgcv::gam(alive ~ s(start_bc),
  #                        data = gam.dat,
  #                        method = "REML",
  #                        family = binomial("logit"))},
  #       error = function(cond) {return(NA)})
  #   bc.fit <- tryCatch({  bc = mgcv::gam(end_bc ~ s(start_bc, k = 3),
  #                      data = gam.dat[gam.dat$alive == 1, ],
  #                      family = "Gamma",
  #                      method = "REML")},
  #              error = function(cond) {return(NA)})
  #   gam.out[["surv"]] <- surv.fit
  #   gam.out[["bc"]] <- bc.fit
  #   gam.out
  #   }) |> purrr::set_names(unique(cohorts$id)) |> 
  #     tibble::enframe() |> 
  #     dplyr::rename(cohort = name, gam = value)
  # )
  # gamfit$surv <- purrr::map(.x = gamfit$gam, .f = ~.x[["surv"]])
  # gamfit$bc <- purrr::map(.x = gamfit$gam, .f = ~.x[["bc"]])
  # gamfit$gam <- NULL
  # 
  # all.fitted <- sum(purrr::map_lgl(.x = gamfit$surv, .f = ~sum(class(.x)=="logical"))) + 
  #   sum(purrr::map_lgl(.x = gamfit$bc, .f = ~sum(class(.x)=="logical"))) == 0
  # 
  # if(!all.fitted) simwarn <- "Warning: Some terminal functions could not be fitted. Insufficient data available." else simwarn <- NULL
  
  # ............................................................
  # Add parameters to list output
  # ............................................................
  
  outsim$gam <- list(fit = list(surv = surv.fit, bc = bc.fit), 
                     dat = gam.dt, 
                     pred = list(bc_gest = mbc_preds, bc = bc_preds, surv = surv_preds))
  
  outsim[["dead"]] <- locs.dead
  outsim[["birth"]] <- locs.birth
  outsim[["abort"]] <- abort.rate

  outsim$param <- list(nsim = nsim, 
                       cohortID = cohortID,
                       cohorts = cohorts)

  outsim$init <- list(month = init.month, xy = init.inds, dest = migration.df)
  outsim$run <- run_time
  if(is.null(scenario)){
  outsim$scenario <- list(scenario)
  } else {
  outsim$scenario <- scenario
  }
  
  class(outsim) <- c("narwsim", class(outsim))
  class(outsim$init$xy) <- c("xyinits", class(outsim$init$xy))
  cat("\nDone!\n")
  cat(paste0("Time elapsed: ", run_time))
  gc()
  return(outsim)
}