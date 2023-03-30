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
                 cohortID = 1:6,
                 scenario = NULL,
                 init.month = 2,
                 step.size = 120,
                 stressors = TRUE,
                 growth = TRUE,
                 n.cores = NULL,
                 n.prop = 50, # In line with sample sizes used in case studies from Michelot (2020)
                 progress = TRUE){
  
  # nsim = 10
  # cohortID = 5
  # scenario = NULL
  # init.month = 2
  # step.size = 120
  # stressors = FALSE
  # growth = TRUE
  # n.cores = NULL
  # n.prop = 50
  # progress = TRUE

  if(!init.month %in% 1:12) stop("<init.month> must be an integer between 1 and 12")
  if(!"init.model" %in% ls(envir = .GlobalEnv)) stop("The model must first be initialized using <load_model>")
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

  # # For simulator
  # cohort.names <- c("Juveniles (male)", # 1
  #                   "Juveniles (female)", # 2
  #                   "Adults (male)", # 3
  #                   "Adults (female, pregnant)", # 4
  #                   "Adults (female, lactating)", # 5
  #                   "Adults (female, resting)") # 6
  # 
  # # For population model
  # if(5 %in% cohortID) all.cohorts <- c("Calves (male + female)", cohort.names) else all.cohorts <- cohort.names
  # age.cat <- c("Calves", "Juveniles", "Juveniles", "Adults", "Adults", "Adults", "Adults")
  # 
  # cohort.ab <- abbreviate(tolower(cohort.names), minlength = 6)
  # cohort.ab <- gsub(pattern = "j\\(", replacement = "jv\\(", x = cohort.ab)
  # cohort.ab <- unname(gsub(pattern = "a\\(", replacement = "ad\\(", x = cohort.ab))
  
  # age.minmax <- matrix(data = c(rep(c(1,9), 2), rep(c(9,69), 4)), nrow = 6, ncol = 2, byrow = TRUE)
  # row.names(age.minmax) <- cohort.ab
  
  # cohort.names <- cohort.names[cohortID]
  # cohort.ab <- cohort.ab[cohortID] 
  
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
  maps.weighted <- lapply(density_weighted, raster::as.matrix)
  if(!all(sapply(maps, class)[1,] %in% "matrix")) stop("Cannot find input of class <matrix>.")
  
  # ............................................................
  # Load prey, daylight, and stressor surfaces
  # ............................................................
  
  preymaps <- lapply(dummy_prey, raster::as.matrix)
  fishingmaps <- lapply(dummy_fishing, raster::as.matrix)
  vesselmaps <- lapply(dummy_vessels, raster::as.matrix)
  noisemaps <- lapply(dummy_noise, raster::as.matrix)
  daylightmaps <- lapply(daylight, raster::as.matrix)
  regionsmap <- raster::as.matrix(regions_m)
  
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
  start.date <- paste0(lubridate::year(Sys.Date()), stringr::str_pad(init.month, width = 2, pad = "0"), "01")
  init.date <- paste0(lubridate::year(Sys.Date()), "-", stringr::str_pad(init.month, width = 2, pad = "0"), "-", "00")
  date_seq <- c(init.date, as.character(seq(from = lubridate::ymd(start.date), by = 'day', length.out = 365)))
  map_seq <- c(init.month, lubridate::month(date_seq[2:366]))
  
  # Simulate initial locations for latent animals
  # Return cell ID for each animal (columns) x month (rows)
  init.inds <- purrr::map(.x = cohortID, 
    .f = ~init_xy(maps = maps, maps.weighted = maps.weighted, coords = coords, cohort.id = .x, nsim = nsim)) |> 
    purrr::set_names(nm = cohorts[id %in% cohortID, abb])
  
  # ............................................................
  # Run simulation
  # ............................................................

  if(length(cohortID) > 1){
    
    Ncores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)-1
    
    # Set number of cores to use
    if(is.null(n.cores)){
      n.cores <- min(4, length(cohortID), Ncores)
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
    
    # # Compute estimates
    out <- foreach::foreach(i = seq_along(cohortID), 
                            .options.snow = opts,
                            .packages = c("Rcpp"), 
                            .noexport = c("NARW_simulator")) %dopar% {
                              
                              # Set environments
                              .GlobalEnv$coords <- coords

                              Rcpp::sourceCpp("src/simtools.cpp")

                              if(cohortID[i] == 5) d <- maps.weighted else d <- maps
                              if(cohortID[i] >0) d <- maps.weighted else d <- maps
                              
                              NARW_simulator(
                                cohortID = cohortID[i],
                                densities = d,
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
                                progress = FALSE)} # End foreach loop
    
    names(out) <- cohorts[id %in% cohortID, abb]
    
    close(pb) # Close progress bar
    
  } else {
    
    cat("Running simulations ...\n")

    if(cohortID == 5) d <- maps.weighted else d <- maps
    if(cohortID > 0) d <- maps.weighted else d <- maps
    
    # Simulated animals begin at the position of their corresponding latent animal for the starting month
    # This can be checked by running the below code
    # xpos <- apply(init.inds[[1]], 2, FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
    # ypos <- apply(init.inds[[1]], 2, FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
    
    start.time <- Sys.time()
    
    out <- list(NARW_simulator(
      cohortID = cohortID,
      densities = d,
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
  outsim[["sim"]] <- consolidate(out.dt, nsim, cohorts[id %in% cohortID, name], date_seq)
  
  # ............................................................
  # Mortality and start/end body condition
  # ............................................................

  gam.dt <- purrr::map(.x = seq_along(cohortID), .f = ~ {
    
    dt_adjuv <- data.table::merge.data.table(x = dplyr::left_join(
      x = outsim[["sim"]][[.x]][day == 0, list(cohort = unique(cohortID), cohort_name = unique(cohort_name), start_bc = bc), whale],
      y = outsim[["sim"]][[.x]][day == 365, list(end_bc = bc, alive), whale], by = "whale"
    ), y = outsim[["sim"]][[.x]][day > 0, .SD[ifelse(which.min(alive) == 1, 0, which.min(alive))],
                                 .SDcols = c("date", "day", "cohort", "easting", "northing", "region", "bc", "strike"), whale
    ], by = c("cohort", "whale"), all.x = TRUE)
    
    if (cohortID[.x] == 5) {
      dt_calves <- data.table::merge.data.table(
        x = dplyr::left_join(
          x = outsim[["sim"]][[.x]][day == 0, list(cohort = unique(cohort_calf), cohort_name = unique(cohort_name), start_bc = bc_calf), whale],
          y = outsim[["sim"]][[.x]][day == 365, list(end_bc = bc_calf, alive_calf), whale], by = "whale"
        ) |> dplyr::rename_with(~ tolower(gsub("_calf", "", .x, fixed = TRUE))),
        y = outsim[["sim"]][[.x]][day > 0, .SD[ifelse(which.min(alive_calf) == 1, 0, which.min(alive_calf))],
                                  .SDcols = c("date", "day", "cohort_calf", "easting", "northing", "region", "bc_calf", "strike"), whale
        ] |>
          dplyr::rename_with(~ tolower(gsub("_calf", "", .x, fixed = TRUE))), by = c("cohort", "whale"), all.x = TRUE
      )
    } else {
      dt_calves <- data.table::data.table()
    }
    
    dplyr::bind_rows(dt_adjuv, dt_calves) 
  }) |> data.table::rbindlist()
  
  # Add full cohort information
  gam.dt <- data.table::merge.data.table(x = gam.dt, y = cohorts, by.x = "cohort", by.y = "id", all.x = TRUE)
  data.table::setcolorder(gam.dt, c("cohort", "cohort_name", "name", "class", "abb", "whale", "alive", "start_bc", "end_bc", "bc", "strike",
                                    "date", "day", "easting", "northing", "region")) 
  
  # locs.dead <- purrr::map(.x = seq_along(cohortID), .f = ~ {
  # 
  #   dt1 <- outsim[["sim"]][[.x]][day > 0, .SD[ifelse(which.min(alive) == 1, 0, which.min(alive))],
  #             .SDcols = c("date", "day", "cohort", "easting", "northing", "region", "bc", "strike"), whale]
  # 
  #   if(cohortID[.x] == 5){
  #   dt2 <- outsim[["sim"]][[.x]][day > 0, .SD[ifelse(which.min(alive_calf) == 1, 0, which.min(alive_calf))],
  #        .SDcols = c("date", "day", "cohort_calf", "easting", "northing", "region", "bc_calf", "strike"), whale] |>
  #       dplyr::rename_with(~ tolower(gsub("_calf", "", .x, fixed = TRUE)))
  #   } else {
  #     dt2 <- data.table::data.table()
  #   }
  #     dplyr::bind_rows(dt1, dt2)
  #   }) |> data.table::rbindlist()

  # Only retain dead animals
  locs.dead <- gam.dt[alive == 0]
  locs.dead[, cause_death:=ifelse(strike == 1, "strike", "starve")]
  data.table::setorder(locs.dead, whale, -cohort)
  outsim[["dead"]] <- locs.dead

  # ............................................................
  # Fit GAM models
  # ............................................................
  
  # https://stats.stackexchange.com/questions/403772/different-ways-of-modelling-interactions-between-continuous-and-categorical-pred
  
  # print(gam.dt)
  
  surv.fit <- suppressWarnings(
    tryCatch(
      {mgcv::gam(alive ~ factor(cohort) + s(start_bc, by = factor(cohort)),
          data = gam.dt,
          method = "REML",
          family = binomial("logit"))},
      error = function(cond) {
        return(NA)
      }))
  
  bc.fit <- suppressWarnings(
    tryCatch(
      {mgcv::gam(end_bc ~ factor(cohort) + s(start_bc, by = factor(cohort)),
                 data = gam.dt[gam.dt$alive == 1, ],
                 method = "REML",
                 family = "Gamma")},
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
  
  # outsim$gam <- list(fit = gamfit, dt = gam.dt)
  outsim$gam <- list(fit = list(surv = surv.fit, bc = bc.fit), dat = gam.dt)
  outsim$param <- list(nsim = nsim, 
                       cohortID = cohortID,
                       cohorts = cohorts)
  outsim$init <- list(month = init.month, xy = init.inds)
  outsim$run <- run_time
  
  class(outsim) <- c("narwsim", class(outsim))
  class(outsim$init$xy) <- c("xyinits", class(outsim$init$xy))
  cat("\nDone!\n")
  cat(paste0("Time elapsed: ", run_time))
  gc()
  return(outsim)
}