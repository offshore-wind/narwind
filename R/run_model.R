#' Run the bioenergetic model
#'
#' Prepare data required for model runs. 
#'
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
                      init.month = 9,
                      step.size = 120,
                      cohorts = FALSE,
                      n.cores = NULL,
                      n.prop = 50,
                      show.progress = TRUE){
  
  # nsim = 10
  # init.month = 9
  # step.size = 79.8
  # cohorts = FALSE
  # n.cores = NULL
  # prey = NULL
  # show.progress = TRUE
  
  if(!init.month %in% 1:12) stop("<init.month> must be an integer between 1 and 12")
  if(!"init.model" %in% ls(envir = .GlobalEnv)) stop("The model must first be initialized")
  
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
  
  if(cohorts){
    
    cohort.names <- c("Calves",
                      "Juveniles (male)",
                      "Juveniles (female)",
                      "Adults (male)",
                      "Adults (female, pregnant)",
                      "Adults (female, lactating)",
                      "Adults (female, resting)")
    
    cohort.ab <- abbreviate(tolower(cohort.names), minlength = 6)
    cohort.ab <- gsub(pattern = "j\\(", replacement = "jv\\(", x = cohort.ab)
    cohort.ab <- unname(gsub(pattern = "a\\(", replacement = "ad\\(", x = cohort.ab))
    
  } else {
    
    cohort.names <- "North Atlantic right whales"
    cohort.ab <- abbreviate(tolower(cohort.names), minlength = 4)
    
  }
  
  # Calculate maximum Manhattan distance for given step size
  # D <- max(sin(seq(0,90)*pi/180) * step.size + cos(seq(0,90)*pi/180) * step.size)
  
  # ............................................................
  # Load spatial support and regions
  # ............................................................
  
  # Convert from logical to numeric; also need to flip when converting from SpatialGridDataFrame to matrix
  geomap <- t(1*raster::as.matrix(density_support))
  
  # ............................................................
  # Load density surfaces
  # ............................................................
  
  maps <- lapply(density_narw, raster::as.matrix)
  
  # ............................................................
  # Load prey surfaces
  # ............................................................
  
  preymaps <- lapply(prey_maps, raster::as.matrix)
  
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
  
  # ............................................................
  # Run simulation
  # ............................................................
  
  # TODO: estimate spatial step size
  
  start.time <- Sys.time()
  
  # if(method == "geodesic") geo <- TRUE else geo <- FALSE
  
  if(cohorts){
    
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
    
    # # Compute estimates
    locs <- foreach::foreach(i = seq_along(cohort.names), 
                             # .options.snow = opts,
                             .packages = c("Rcpp"), 
                             .noexport = c("simAnnualCoupled")) %dopar% {
                               
                               Rcpp::sourceCpp("src/simtools.cpp")
                               
                               # Simulate initial locations for latent animals
                               # Return cell ID for each animal (columns) x month (rows)
                               init_inds <- do.call(rbind, lapply(maps, function(m) {
                                 p <- as.numeric(m)
                                 p[!is.finite(p)] <- 0
                                 sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
                               }))
                               
                               simAnnualCoupled(
                                 geodesic = TRUE,
                                 densities = maps,
                                 densitySeq = map_seq,
                                 latentDensitySeq = 1:12,
                                 prey = preymaps,
                                 M = n.prop, # Number of proposals used in the importance sampler (Michelot, 2019)
                                 stepsize = step.size / 2, # Movement framework involves two steps (Michelot 2019)
                                 xinit = matrix(data = coords[init_inds, 'x'], nrow = nrow(init_inds)),
                                 yinit = matrix(data = coords[init_inds, 'y'], nrow = nrow(init_inds)),
                                 support = geomap,
                                 limits = map_limits,
                                 resolution = map_resolution,
                                 progress = TRUE)} # End foreach loop
    
    names(locs) <- cohort.ab
    
    close(pb)
    
  } else {
    
    cat("Running simulations ...\n")
    
    # Simulate initial locations for latent animals
    # Return cell ID for each animal (columns) x month (rows)
    init_inds <- do.call(rbind, lapply(maps, function(m) {
      p <- as.numeric(m)
      p[!is.finite(p)] <- 0
      sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
    }))
    
    locs <- list(simAnnualCoupled(
      geodesic = TRUE,
      densities = maps,
      densitySeq = map_seq,
      latentDensitySeq = seq_along(density_narw),
      prey = preymaps,
      M = n.prop, # Number of proposals used in the importance sampler for movement (Michelot, 2019)
      stepsize = step.size / 2, # Movement framework involves two steps (Michelot 2019)
      xinit = matrix(data = coords[init_inds, 'x'], nrow = nrow(init_inds)),
      yinit = matrix(data = coords[init_inds, 'y'], nrow = nrow(init_inds)),
      support = geomap,
      limits = map_limits,
      resolution = map_resolution,
      progress = ifelse(show.progress, TRUE, FALSE)))

   names(locs) <- cohort.ab

  }
  
  # Transpose the output array to have x,y coordinates for each day as rows
  outsim <- purrr::map(.x = locs, .f = ~{
    .x <- aperm(.x, c(3,1,2))
    dimnames(.x)[[1]] <- format(date_seq)
    dimnames(.x)[[3]] <- paste0("sim", seq_len(dim(.x)[3]))
    # dimnames(.x)[[3]] <- stringi::stri_rand_strings(n = dim(.x)[3], pattern = "[[A-Z][0-9]]", length = 8)
    .x
  })
  
  outsim <- list(locs = outsim)
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  # ............................................................
  # Validate sampling scheme
  # ............................................................
  
  # tidy version of map rasters
  # raster_id <- map_seq[1]
  # df <- as.data.frame(raster::raster(map.spgrd[[raster_id]]), xy = TRUE)
  
  # TODO: do data format checks on df here?
  # scroll across longitude first, so we're in row-major order
  # head(coords)
  # we start with (xmin,ymax) and work our way over.  the cells are all center-offset,
  # so the bounding box is not the correct xmin value
  # map.spgrd[[raster_id]]@bbox
  
  # ............................................................
  # Add parameters to list output
  # ............................................................
  
  outsim$param <- list(cohorts = cohorts,
                       cohort.names = cohort.names,
                       cohort.ab = cohort.ab)
  
  class(outsim) <- c("narwsim", class(outsim))
  cat("\nDone!\n")
  cat(paste0("Time elapsed: ", run_time))
  
  # suppressWarnings(rm("init.model", envir = .GlobalEnv))
  
  return(outsim)
}