#' Run the bioenergetic model
#'
#' Prepare data required for model runs. 
#' @export
#' @param nsim Number of simulated animals
#' @param n.prop Number of proposals used in the importance sampler (Michelot, 2019)
#'
#' @importFrom foreach `%dopar%`
#' @importFrom doParallel registerDoParallel
#' @author Phil J. Bouchet
#' @seealso \code{\link{initialize}}
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' run_model()
#' }

run_model <- function(nsim = 1e3,
                      init.month = 2,
                      step.size = 120,
                      cohort.id = 1:6,
                      stressors = TRUE,
                      growth = TRUE,
                      n.cores = NULL,
                      n.prop = 50, # In line with sample sizes used in case studies from Michelot (2020)
                      progress = TRUE){

  # nsim = 1
  # init.month = 2
  # step.size = 120
  # cohort.id = 5
  # n.cores = NULL
  # n.prop = 50
  # stressors = FALSE
  # growth = FALSE
  # progress = TRUE

  if(!init.month %in% 1:12) stop("<init.month> must be an integer between 1 and 12")
  if(!"init.model" %in% ls(envir = .GlobalEnv)) stop("The model must first be initialized using <load_model>")
  if(any(!cohort.id %in% 1:6)) stop("Unrecognized cohort")
  
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
  
  cohort.names <- c("Juveniles (male)", # 1
                    "Juveniles (female)", # 2
                    "Adults (male)", # 3
                    "Adults (female, pregnant)", # 4
                    "Adults (female, lactating)", # 5
                    "Adults (female, resting)") # 6
  
  cohort.ab <- abbreviate(tolower(cohort.names), minlength = 6)
  cohort.ab <- gsub(pattern = "j\\(", replacement = "jv\\(", x = cohort.ab)
  cohort.ab <- unname(gsub(pattern = "a\\(", replacement = "ad\\(", x = cohort.ab))
  
  # age.minmax <- matrix(data = c(rep(c(1,9), 2), rep(c(9,69), 4)), nrow = 6, ncol = 2, byrow = TRUE)
  # row.names(age.minmax) <- cohort.ab
  
  cohort.names <- cohort.names[cohort.id]
  cohort.ab <- cohort.ab[cohort.id] 
  
  # ............................................................
  # Dose-response
  # ............................................................
  
  spline_dr <- splinefun(x = doseresponse[,2], y = doseresponse[,1])
  median_doseresponse <- spline_dr(runif(10000)) # Using inverse transform sampling
  
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
  # Check inputs
  # ............................................................
  
  if(!all(sapply(maps, class)[1,] %in% "matrix")) stop("Cannot find input of class <matrix>.")
  
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
  
  start.date <- paste0('2022', stringr::str_pad(init.month, width = 2, pad = "0"), '01')
  date_seq <- seq(from = lubridate::ymd(start.date), by = 'day', length.out = 365)
  map_seq <- lubridate::month(date_seq)
  
  # Simulate initial locations for latent animals
  # Return cell ID for each animal (columns) x month (rows)
  init.inds <- purrr::map(.x = cohort.id, 
    .f = ~init_xy(maps = maps, maps.weighted = maps.weighted, coords = coords, cohort.id = .x, nsim = nsim)) |> 
    purrr::set_names(nm = cohort.ab)
  
  # ............................................................
  # Run simulation
  # ............................................................
  
  # TODO: estimate spatial step size
  
  if(length(cohort.id) > 1){
    
    Ncores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)-1
    
    # Set number of cores to use
    if(is.null(n.cores)){
      n.cores <- min(4, length(cohort.names), Ncores)
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
    pb <- utils::txtProgressBar(max = length(cohort.names), style = 3)
    sink()
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    # # Create a parallel cluster and register the backend
    # cl <- suppressMessages(parallel::makeCluster(n.cores))
    # 
    # # After the function is run, shutdown the cluster.
    # on.exit(parallel::stopCluster(cl))
    # 
    # # Register parallel backend
    # doParallel::registerDoParallel(cl)   # Modify with any do*::registerDo*()
    
    cat("Running simulations ...\n")
    
    start.time <- Sys.time()
    
    # # Compute estimates
    out <- foreach::foreach(i = seq_along(cohort.id), 
                            .options.snow = opts,
                            .packages = c("Rcpp"), 
                            .noexport = c("MovementSimulator")) %dopar% {
                              
                              # Set environments
                              .GlobalEnv$coords <- coords

                              Rcpp::sourceCpp("src/simtools.cpp")

                              if(cohort.id[i] == 5) d <- maps.weighted else d <- maps
                              if(cohort.id[i] >0) d <- maps.weighted else d <- maps
                              
                              MovementSimulator(
                                cohortID = cohort.id[i],
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
    
    names(out) <- cohort.ab
    
    close(pb) # Close progress bar
    
  } else {
    
    cat("Running simulations ...\n")

    if(cohort.id == 5) d <- maps.weighted else d <- maps
    if(cohort.id >0) d <- maps.weighted else d <- maps
    
    # Simulated animals begin at the position of their corresponding latent animal for the starting month
    # This can be checked by running the below code
    # xpos <- apply(init_inds, 2, FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
    # ypos <- apply(init_inds, 2, FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
    
    start.time <- Sys.time()
    
    out <- list(MovementSimulator(
      cohortID = cohort.id,
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
  
  # Transpose the output array to have data for each day as rows
  
  out.t <- transpose_array(input = out, cohortID = cohort.id, dates = date_seq) |> 
    purrr::set_names(nm = cohort.ab)

  for (k in seq_along(cohort.id)) {
    if (cohort.id[k] %in% 4:5) out.t[[k]][["attrib"]] <- abind::abind(out.t[[k]][["attrib"]][[1]], out.t[[k]][["attrib"]][[2]], along = 2)
    if (cohort.id[k] == 5){
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
 
 # for (k in seq_along(cohort.id)) {
 #   outsim[["sim"]][[cohort.ab[k]]] <-
 #     abind::abind(
 #       out.t[[k]][["locs"]],
 #       out.t[[k]][["attrib"]],
 #       out.t[[k]][["stress"]],
 #       out.t[[k]][["E"]],
 #       out.t[[k]][["kj"]],
 #       out.t[[k]][["activ"]],
 #       along = 2)
 # }

  outsim <- list()  
  outsim[["sim"]] <- purrr::map2(.x = out.dt,
                                 .y = cohort.names, 
                                 .f = ~{
      dt <- data.table::data.table(day = rep(1:365, times = nsim), 
                                   date = rep(date_seq, times = nsim),
                                   whale = rep(1:nsim, each = 365))
        a <- cbind(array2dt(.x), dt)
        a$region <- sort(regions$region)[a$region]
        a$cohort_name <- .y
        data.table::setcolorder(a, c((ncol(a)-3):(ncol(a)-1), 1:(ncol(a)-4))) 
      
    }
  )
  
  # ............................................................
  # Fit GAM models
  # ............................................................
  
  pop.dt <- purrr::map(.x = outsim[["sim"]], .f = ~{
    dplyr::left_join(x = .x[day == 1, list(cohort = unique(cohort), cohort_name = unique(cohort_name), start = bc), whale],
                     y = .x[day == 365, list(end = bc, alive), whale], by = "whale")
    }) |> do.call(what = rbind)
  
  
  
  # ............................................................
  # Add parameters to list output
  # ............................................................
  
  outsim$param <- list(nsim = nsim, 
                       cohort.id = cohort.id,
                       cohort.names = cohort.names,
                       cohort.ab = cohort.ab)
  
  outsim$init <- list(month = init.month, xy = init.inds)
  outsim$run <- run_time
  
  class(outsim) <- c("narwsim", class(outsim))
  class(outsim$init$xy) <- c("xyinits", class(outsim$init$xy))
  cat("\nDone!\n")
  cat(paste0("Time elapsed: ", run_time))
  gc()
  return(outsim)
}