#' Predictions of right whale abundance
#' 
#' Takes an object of class \code{narwsim} (as returned by \code{\link[narwind]{narw}}) and produces forecasts of right whale population size (by cohort) over a time horizon of interest. Prediction uncertainty is estimated from both process variance (i.e., replicate projections) and parameter uncertainty (i.e., statistical uncertainty in the relationships between individual body condition and health/survival, respectively). The latter is considered only if replicate coefficients have been sampled from the posteriors of the fitted survival and body condition GAMs using \code{\link[narwind]{augment}}.
#'
#' @param ... One or more objects of class \code{narwsim}.
#' @param n Integer. Number of replicate projections. Defaults to \code{100}.
#' @param yrs Integer. Time horizon, specified either as the desired number of years (from current) or the desired target end year. Defaults to \code{35}, which is commensurate with the average lifespan of a typical wind farm.
#' @param lpmatrix Logical. If \code{TRUE}, prediction uncertainty includes parameter uncertainty.  
#' @param piling Integer. Year of construction. Defaults to \code{1}, such that piling occurs on the first year of the projection, followed by O&M for the remainder.
#' @param progress Logical. If \code{TRUE}, a progress bar is shown during execution. Defaults to \code{FALSE}
#' @return A list of class \code{narwproj}.
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' p <- predict(m)
#' }
predict.narwsim <- function(...,
                            n = 100,
                            yrs = 35,
                            lpmatrix = TRUE,
                            piling = 1,
                            progress = TRUE) {

  cat("--------------------------------------------------------------------------------\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("\n")
  cat("                 NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("                            *** POPULATION MODEL ***\n")
  cat("\n")
  cat("--------------------------------------------------------------------------------\n")
  cat("--------------------------------------------------------------------------------\n\n")
  
  # Adapted from original code by Scott Creel
  # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
  
  #'------------------------------------------------------------
  # Function checks ----
  #'------------------------------------------------------------
  
  args <- list(...)
  
  # Identify and extract <narwsim> objects 
  which.obj <- which(purrr::map_lgl(.x = args, .f = ~inherits(.x, "narwsim")))
  if(length(which.obj)==0) stop("No object of class <narwsim> found.")
  obj <- args[which.obj]
  
  # Check whether posterior samples are available
  if(any(purrr::map_lgl(.x = obj, .f = ~is.null(.x$post)))) warning("Posterior samples for terminal functions not available!\nRun the augment() function to estimate variance through posterior simulation.", sep = "")

  # Check whether terminal functions are available
  if(sum(sapply(X = obj, FUN = function(o){
    sum(purrr::map_lgl(.x = o$gam$fit, .f = ~any(is.null(unlist(.x)))))}) > 0)){
    stop("Sample sizes insufficient to reliably estimate survival and health functions. Re-run narw() with a higher <nsim>.") 
  }
  
  # Check that all 6 population cohorts have been simulated
  if(any(!purrr::map_dbl(.x = obj, .f = ~length(.x$param$cohort)) == 6)) stop("Missing cohorts! Cannot generate projections")
  
  if(length(obj) > 3) stop("Too many objects - maximum of 3 scenarios (baseline, construction, O&M) allowed.")
  
  if(yrs <= 0) stop("<yrs> must be a positive integer")
  if(piling <= 0) stop("<piling> must be a positive integer")
  
  # If only one object is supplied, must either baseline conditions or O&M
  if(length(obj) == 1){
    # Change <piling> if only baseline has been provided
    if(obj[[1]]$scenario$phase == 0 & piling < yrs) piling <- yrs + 1 
    if(obj[[1]]$scenario$phase == 1) stop("Missing O&M scenario")
    # Change <piling> if only O&M has been provided
    if(obj[[1]]$scenario$phase == 2 & piling < yrs) piling <- yrs + 1 
  }
  
  # If multiple objects are supplied
  if(length(obj) > 1 & sum(duplicated(purrr::map_dbl(.x = obj, .f = ~.x$scenario$phase))) > 0) stop("Duplicate scenarios detected.")
  if(length(obj) > 1 & piling < yrs){
    if(piling > 1 & all(purrr::map_dbl(.x = obj, .f = ~.x$scenario$phase) %in% 1:2)) stop("Missing baseline scenario")
    if(all(purrr::map_dbl(.x = obj, .f = ~.x$scenario$phase) %in% c(0,2))) stop("Missing construction scenario")
    if(all(purrr::map_dbl(.x = obj, .f = ~.x$scenario$phase) %in% 0:1)) stop("Missing O&M scenario")
  }
  
  # Current year
  current.yr <- lubridate::year(lubridate::now())
  
  # User can specify the <yrs> argument either as a number of years from now or as a target end year
  if(yrs > 2000) yrs <- yrs - current.yr
  if(piling > 2000) piling <- piling - current.yr
  
  # Cohort information (augmented to separate out male/female calves)
  # cohortID <- obj[[1]]$param$cohort
  cohorts.proj <- obj[[1]]$param$cohorts |> dplyr::slice(1) |> 
    dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |> 
    dplyr::bind_rows(obj[[1]]$param$cohorts) |> 
    dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
  cohorts.proj <- cohorts.proj[c(1,3,5,2,4,6,7,8)]
  
  # Attributes to monitor during projection
  mat.attribs <- c("alive", "cohort", "female", "age", "length", "tot_mass", "lean_mass", "bc",
                   "p_surv", "min_bc", "time_rest", "birth", "reprod")
  
  #'------------------------------------------------------
  # PHASES ----
  #'------------------------------------------------------
  
  # Define schedule of wind farm development
  phases <- sapply(X = obj, FUN = function(o) o$scenario$phase)
  phase.names <- c("base", "const", "ops")[phases+1]
  schedule <- rep(0, yrs+1)
  if(piling < yrs){
    schedule[piling + 1] <- 1
    schedule[(piling + 2): (yrs+1)] <- 2
  }
  names(schedule) <- paste0("yr ", 0:yrs)
  
  #'------------------------------------------------------
  # TERMINAL FUNCTIONS ----
  #'------------------------------------------------------
  
  # Fitted models
  mbc.model <- gam_gest
  survbc.model <- purrr::map(.x = obj, .f = ~ list("surv" = .x$gam$fit$surv, "bc" = .x$gam$fit$bc)) |>
    purrr::set_names(nm = phase.names)
  
  # Inverse link functions
  ilink.mbc <- family(mbc.model)$linkinv
  ilink.survbc <- purrr::map(.x = obj[1], 
                             .f = ~list("surv" = family(.x$gam$fit$surv)$linkinv,
                                        "bc" = family(.x$gam$fit$bc)$linkinv))[[1]]

  # Posterior samples
  mcmc <- purrr::map2(.x = obj, 
                      .y = survbc.model,
                      .f = ~{
    
    # If posterior samples are not available
    if(is.null(.x$post)){
      
      list(
        surv = matrix(data = coef(.y$surv), nrow = 1, ncol = length(coef(.y$surv))),
        bc = matrix(data = coef(.y$bc), nrow = 1, ncol = length(coef(.y$bc))),
        mbc = matrix(data = coef(mbc.model), nrow = 1, ncol = length(coef(mbc.model))),
        n = 1
      )
      
    } else { # If posterior samples are available
      
      list(
        surv = .x$post$samples$surv,
        bc = .x$post$samples$bc,
        mbc = .x$post$samples$mbc,
        n = .x$post$nsamples
      )
      
    }
    
  }) |> purrr::set_names(nm = phase.names)
  
  #'------------------------------------------------------
  # Abortion rate ----
  #'------------------------------------------------------
  
  abort.rate <- purrr::map_dbl(.x = obj, .f = ~sum(.x$abort[, abort])/.x$param$nsim)
  names(abort.rate) <- phase.names
  
  #'------------------------------------------------------
  # INITIALIZATION ----
  #'------------------------------------------------------
  
  cat("Date:", as.character(Sys.Date()), "\n")
  now.time <- Sys.time() |> 
    stringr::str_split(" ") |> 
    purrr::map_chr(2)
  cat("Time:", now.time, "\n\n")
  
  cat("––– Projections:\n\n")
  
  cat("+ Replicates: N =", formatC(n, big.mark = ","), "\n")
  cat("+ Horizon:", yrs, "years\n")
  cat("+ Parameter uncertainty:", ifelse(lpmatrix, "Yes", "No"), "\n")
  cat("\n––– Timeline:\n\n")
  proj_timeline(schedule)
  cat("––– Execution:\n\n")
  
  # Define initial population vector using predictions from NOAA model
  N_0 <- N0[["yr2019"]]
  names(N_0) <- cohorts.proj[, name]
  
  console(msg = "Initializing")
  
  # Initial cohort vector
  cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts.proj), FUN = function(x) rep(cohorts.proj$id[x], each = N_0[x])))
  
  # Sterile females - 4% never reproduce (based on PCoMS estimate)
  nonrep <- 0.04
  reprod.fem <- rep(1, sum(N_0))
  reprod.fem[sample(
    x = which(cohort.vec == 6), # Assigning 0 to resting females only
    size = round(nonrep * sum(N_0[6:8]), 0) # Based on total number of mature females in the population
  )] <- 0
  
  # Number of individuals in each cohort
  narw.pop <- array(
    data = NA, c(n, yrs + 1, nrow(cohorts.proj)),
    dimnames = list(
      paste0("prj ", 1:n),
      paste0("yr ", 0:yrs),
      cohorts.proj$name
    )
  )
  narw.pop[, 1, ] <- rep(N_0, each = n)
  
  # Total population size
  tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj", 1:n), paste0("yr ", 0:yrs)))
  tot.pop[, 1] <- sum(N_0)
  
  # Initialise array to store numbers of births
  births <- array(
    data = NA, c(yrs + 1, n),
    dimnames = list(
      paste0("yr ", 0:yrs),
      paste0("prj ", 1:n)
    )
  )
  
  maxbc <- find_maxBC()
  
  #'------------------------------------------------------
  # PROJECTIONS ----
  #'------------------------------------------------------
  
  # This uses nested loops. 
  # The prj loop (outermost loop) replicates the projection <n> times.
  # The i loop is next, and steps across all years of projection from an initial population vector.
  
  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "Running projections [:bar] :percent eta: :eta",
    total = n, 
    clear = TRUE, 
    width = 80
  )
  
  console(msg = "Initializing", suffix = tickmark())
  
  start.time <- Sys.time()
  
  narw.individuals <- vector(mode = "list", length = n)
  
  for(prj in 1:n){
    
    # Update progress bar
    if(progress) pb$tick() 
    
    # Initial development phase
    initial.phase <- phase.names[schedule[1] + 1]
    
    # Create matrices and initialize them
    # rows: years <yrs>
    # columns: <attributes>
    # layers: individuals <n>
    # 4th dimension: replicate projections -> later converted to list
    narw.indiv <- array(
      data = NA, c(yrs + 1, length(mat.attribs), sum(N_0)),
      dimnames = list(
        paste0("yr ", 0:yrs), mat.attribs,
        paste0("whale ", 1:(sum(N_0)))
      )
    )
    
    # Alive and population cohort
    narw.indiv[1, "alive", ] <- 1
    narw.indiv[1, "cohort", ] <- cohort.vec
    
    # Sex
    narw.indiv[, "female", 1:N_0[1]] <- 0                                      #  -- Calves (male)
    fem <- which(cohort.vec == 0)                                              #  -- Calves (female)
    narw.indiv[, "female", fem[fem > N_0[1]]] <- 1                             #  -- Juveniles and adults
    narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6))] <- 1
    narw.indiv[, "female", which(cohort.vec %in% c(1, 3))] <- 0
    
    # Age
    ages <- start_age_vec(cohort.vec)
    narw.indiv[1, "age", ] <- ages
    
    # Total body length
    lengths <- age2length_vec(ages)
    narw.indiv[1, "length", ] <- lengths
    
    # Lean mass
    mass <- length2mass_vec(lengths)
    narw.indiv[1, "lean_mass", ] <- mass
    
    # Initial body condition - sampled from values obtained after
    # a 3-month burn-in within the spatial simulation
    bc_init <- purrr::set_names(unique(obj[[schedule[1] + 1]]$gam$dat$cohort)) |> 
      purrr::map(.f = ~obj[[schedule[1] + 1]]$gam$dat[cohort == .x, start_bc])
    bc <- sapply(X = cohort.vec, FUN = function(x){ sample(x = bc_init[[as.character(x)]], size = 1, replace = TRUE)})
    narw.indiv[1, "bc", ] <- bc
    
    # Total mass
    narw.indiv[1, "tot_mass", ] <- mass / (1-bc)
    
    #' NOTE: ---------------------------------------------------------------------------------------
    # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001).
    # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
    # This corresponds to a Normal (7.7, 1.45)
    # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
    # 
    # Stewart et al. 2022 -- 
    # The degree to which the energetic reserves of females are depleted during lactation may govern
    # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).
    #' --------------------------------------------------------------------------------------------
    
    # Non-reproductive females
    narw.indiv[1, "reprod", ] <- reprod.fem
    
    # Reserves needed to leave resting state and initiate pregnancy
    
    if(lpmatrix){

      Xp <- predict(mbc.model, data.frame(mass = narw.indiv[1, "tot_mass", ]), type = "lpmatrix")
      
      mbc.sample <- sample(
        x = seq_len(mcmc[[initial.phase]]$n),
        size = nrow(Xp),
        # size = ifelse(mcmc[[initial.phase]]$n == 1, 1, nrow(Xp)),
        replace = TRUE)
      
      mbc.curves <- mcmc[[initial.phase]]$mbc[mbc.sample,]
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * ilink.mbc(rowSums(Xp*mbc.curves))  
    } else {
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * predict(mbc.model, data.frame(mass = narw.indiv[1, "tot_mass", ]), type = "response")
    }

    
    # Time spent in resting state
    narw.indiv[1, "time_rest", ] <- (cohort.vec == 6)

    # Calving events
    narw.indiv[1, "birth", ] <- as.numeric(cohort.vec == 4)
    
    if(lpmatrix){
    
    # Survival probability
    Xp <- predict(survbc.model[[initial.phase]]$surv, data.frame(start_bc = bc, cohort = cohort.vec), type = "lpmatrix")
    surv.sample <- sample(
      x = seq_len(mcmc[[initial.phase]]$n),
      size = nrow(Xp),
      replace = TRUE)
    
    surv.curves <- mcmc[[initial.phase]]$surv[surv.sample,]
    narw.indiv[1, "p_surv", ] <- ilink.survbc[["surv"]](rowSums(Xp*surv.curves))
    
    } else {
      
      narw.indiv[1, "p_surv", ] <- predict(survbc.model[[initial.phase]]$surv, data.frame(start_bc = bc, cohort = cohort.vec), type = "response")
      
    }

    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      if(tot.pop[prj, i-1] > 0) {
        
        # Current phase of development
        current.phase <- phase.names[schedule[i]+1]
        
        # New samples from the posteriors of model coefficients
        post.sample <- sample(
          x = seq_len(mcmc[[current.phase]]$n),
          size = nrow(Xp),
          replace = TRUE)
        
        bc.curves <- mcmc[[current.phase]]$bc[post.sample,]
        mbc.curves <- mcmc[[current.phase]]$mbc[post.sample,]
        
        #' ----------------------------
        # ALIVE
        #' ----------------------------
        # Determine whether the animal survived
        # alive <- rbinom(n = length(ps), size = 1, prob = ps) * (narw.indiv[i-1, "age", ] <=69)
        alive <- rbinom(
          n = dim(narw.indiv)[3], 
          size = 1,
          prob = (narw.indiv[i - 1, "alive", ] * narw.indiv[i - 1, "p_surv", ])
        ) * (narw.indiv[i - 1, "age", ] <= 69)
        
        narw.indiv[i, "alive", ] <- alive
        
        #' ----------------------------
        # AGE & SEX
        #' ----------------------------
        
        # Increment age
        narw.indiv[i, "age", ] <- alive * (narw.indiv[i-1, "age", ] + 1)
        
        # Sex remains the same
        narw.indiv[i, "female", ] <- narw.indiv[i-1, "female", ]
        
        #' ----------------------------
        # COHORTS
        #' ----------------------------
        
        # Maturity - transitions between cohorts
        narw.indiv[i, "cohort", ] <-
          increment_cohort(
            alive = alive,
            cohort = narw.indiv[i-1, "cohort", ],
            age = narw.indiv[i, "age", ],
            female = narw.indiv[i, "female", ],
            bc = narw.indiv[i-1, "bc", ],
            min_bc = narw.indiv[i-1, "min_bc", ],
            reprod = narw.indiv[i-1, "reprod", ],
            abort = abort.rate[current.phase])
        
        #' ----------------------------
        # GROWTH
        #' ----------------------------
        
        # Increment length
        narw.indiv[i, "length", ] <- alive * age2length_vec(narw.indiv[i, "age", ])
        
        # Increment lean mass
        narw.indiv[i, "lean_mass", ] <- alive * length2mass_vec(narw.indiv[i, "length", ])
        
        # Predict new body condition from current body condition
        if (lpmatrix) {
          Xp <- predict(survbc.model[[current.phase]]$bc,
            data.frame(
              start_bc = narw.indiv[i - 1, "bc", ],
              cohort = narw.indiv[i - 1, "cohort", ]
            ),
            type = "lpmatrix"
          )

          narw.indiv[i, "bc", ] <- alive * clamp(ilink.survbc$bc(rowSums(Xp * bc.curves)), maxbc)
        } else {
          narw.indiv[i, "bc", ] <- alive * clamp(predict(survbc.model[[current.phase]]$bc,
            data.frame(
              start_bc = narw.indiv[i - 1, "bc", ],
              cohort = narw.indiv[i - 1, "cohort", ]
            ),
            type = "response"
          ), maxbc)
        }
         
        # Increment total mass
        narw.indiv[i, "tot_mass", ] <- alive * narw.indiv[i, "lean_mass", ] / (1 - narw.indiv[i, "bc", ])
        
        #' ----------------------------
        # REPRODUCTION
        #' ----------------------------
        
        # Reproductive females
        narw.indiv[i, "reprod", ] <- alive * narw.indiv[i-1, "reprod", ]
     
        # Time spent in resting phase
        narw.indiv[i, "time_rest", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i-1, "time_rest", ] + 1)
        
        # Minimum body condition needed to successfully bring fetus to term without starving
        # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
        if (lpmatrix) {
          Xp <- predict(mbc.model, data.frame(mass = narw.indiv[i, "tot_mass", ]), type = "lpmatrix")
          narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * ilink.mbc(rowSums(Xp * mbc.curves))
        } else {
          narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * predict(mbc.model, data.frame(mass = narw.indiv[i, "tot_mass", ]), type = "response")
        }
          
        # Birth of new calf, conditional on the mother being alive and in pregnant state
        narw.indiv[i, "birth", ] <- alive * (narw.indiv[i, "cohort", ] == 4)
        new.births <- sum(narw.indiv[i, "birth", ])
        births[i, prj] <- new.births
        
        # Number of non-reproductive females in the population
        Nr <- nrow(t(narw.indiv[i, ,])[t(narw.indiv[i, ,])[,"alive"] == 1 & t(narw.indiv[i, ,])[,"female"] == 1 & t(narw.indiv[i, ,])[,"reprod"] == 0, ])
        
        # Total abundance of females
        Nt <- nrow(t(narw.indiv[i, ,])[t(narw.indiv[i, ,])[,"alive"] == 1 & t(narw.indiv[i, ,])[,"female"] == 1 & t(narw.indiv[i, ,])[,"cohort"] %in% c(2,4,5,6), ])
        
        if(new.births > 0){
          
          # Determine sex of newborn calves
          # Note: This needs to be calculated here as maintaining a % of non-reproductive females
          # in the population requires knowing the number of female calves being born
          # Sex ratio based on mean estimated "probability of being female" reported in Linden (2023)
          # https://www.fisheries.noaa.gov/s3/2023-10/TM314-508-0.pdf
          new.females <- rbinom(new.births, 1, 0.46)
          
          # Total abundance of females given new births
          Nt <- Nt + sum(new.females)
          
          # Newborn females that are non-reproductive
          Nb <- min(sum(new.females), max(0, round(nonrep * Nt, 0) - Nr))
          
          newborns.nonreprod <- rep(1, new.births)
          if(Nb > 0 & sum(new.females) > 0) newborns.nonreprod[sample(which(new.females == 1), size = Nb)] <- 0
          
          # Create array for newborns
          new.calves <- array(
            data = NA, c(yrs + 1, length(mat.attribs), new.births),
            dimnames = list(
              paste0("yr ", 0:yrs), mat.attribs,
              paste0("whale ", ((dim(narw.indiv)[3]+1):(dim(narw.indiv)[3] + new.births)))
            )
          )
          
          new.calves[i,,] <- add_calf(new.births, mat.attribs, new.females, newborns.nonreprod)
          narw.indiv <- abind::abind(narw.indiv, new.calves, along = 3)
          
        }
        
        #' ----------------------------
        # SURVIVAL
        #' ----------------------------
        
        if(lpmatrix){
        # Predict survival probability based on body condition
        Xp <- predict(survbc.model[[current.phase]]$surv,
         data.frame(
           start_bc = narw.indiv[i, "bc", ],
           cohort = narw.indiv[i, "cohort", ]
         ), type = "lpmatrix")

        post.sample <- sample(
          x = seq_len(mcmc[[current.phase]]$n),
          size = nrow(Xp),
          replace = TRUE)
        
        surv.curves <- mcmc[[current.phase]]$surv[post.sample,]
        narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * ilink.survbc$surv(rowSums(Xp*surv.curves))
        } else {
          narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * predict(survbc.model[[current.phase]]$surv,
                                                                          data.frame(
                                                                            start_bc = narw.indiv[i, "bc", ],
                                                                            cohort = narw.indiv[i, "cohort", ]
                                                                          ), type = "response")
        }
        
        # Formerly using sapply which is much slower
        # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.surv(Xp[r,] %*% surv.curves[r,])}) 

        #' ----------------------------
        # TOTALS
        #' ----------------------------
        
        # Number of  in each cohort
        # Calves (male)
        narw.pop[prj, i, 1] <- sum((narw.indiv[i, "cohort", ] == 0) *
                                     (narw.indiv[i, "female", ] == 0) *
                                     (narw.indiv[i, "alive", ] == 1))
        
        # Juveniles and adults (male)
        narw.pop[prj, i, 2] <- sum((narw.indiv[i, "cohort", ] == 1) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 3] <- sum((narw.indiv[i, "cohort", ] == 3) * (narw.indiv[i, "alive", ] == 1))
        
        # Calves (female)
        narw.pop[prj, i, 4] <- sum((narw.indiv[i, "cohort", ] == 0) *
                                     (narw.indiv[i, "female", ] == 1) *
                                     (narw.indiv[i, "alive", ] == 1))
        
        # Juvenile and reproductive adults (female)
        narw.pop[prj, i, 5] <- sum((narw.indiv[i, "cohort", ] == 2) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 6] <- sum((narw.indiv[i, "cohort", ] == 4) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 7] <- sum((narw.indiv[i, "cohort", ] == 5) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 8] <- sum((narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i, "alive", ] == 1))
        
        # Total population size
        tot.pop[prj, i] <- sum(narw.indiv[i, "alive", ], na.rm = TRUE)
        
        } else {
          
        narw.indiv[i, "alive", ] <- narw.indiv[i, "cohort", ] <- narw.indiv[i, "female", ] <- 
          narw.indiv[i, "age", ] <- narw.indiv[i, "length", ] <- narw.indiv[i, "tot_mass", ] <- 
          narw.indiv[i, "lean_mass", ] <- narw.indiv[i, "bc", ] <- narw.indiv[i, "p_surv", ] <- 
          narw.indiv[i, "min_bc", ] <- narw.indiv[i, "time_rest", ] <- narw.indiv[i, "birth", ] <- 
          narw.indiv[i, "reprod", ] <- 0
        
        narw.pop[prj, i, 1:8] <- 0
        
        tot.pop[prj, i] <- 0
          
        } # End totpop >0
    } # End years
    
    narw.individuals[[prj]] <- narw.indiv
    
  } # End projections
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  #'------------------------------------------------------
  # SUMMARY ----
  #'------------------------------------------------------

  console(msg = "Running projections", suffix = tickmark())
  console(msg = "Summarizing outputs")
  
  # **** Survival probability ****
  
  prjpreds <- purrr::map(.x = narw.individuals, .f = ~ {
    out <- apply(.x, 3, function(x) x[x[, "alive"] == 1 & !is.na(x[, "alive"]), c("cohort", "p_surv", "bc", "min_bc"), drop = FALSE])
    if(inherits(out, "list")) out <- do.call(rbind, out)
  }) |> do.call(what = rbind) |> 
    data.table::as.data.table()
  
  # **** Total births ****
  
  births.df <- purrr::map(.x = 1:n, .f = ~{
    m <- matrix(rowSums(narw.individuals[[.x]][2:(yrs+1),"birth",], na.rm = TRUE), ncol = 1)
    colnames(m) <- .x
    m
  }) |> do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |>
    dplyr::select(prj, year, birth) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |> 
    dplyr::mutate(prj = as.numeric(prj))
  
  # **** Total deaths ****
  
  deaths.df <- purrr::map(.x = 1:n, .f = ~{
    m <- matrix(apply(X = narw.individuals[[.x]][2:(yrs+1),"alive",],
                      MARGIN = 1,
                      FUN = function(x) {
                        r <- x[!is.na(x)]
                        r <- sum(r == 0)
                        r
                      }), ncol = 1)
    colnames(m) <- .x
    m
  }) |> do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
    dplyr::select(prj, year, death) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |> 
    dplyr::mutate(prj = as.numeric(prj))
  
  # **** Number of births per female ****
  
  births.per.female <- purrr::map(.x = narw.individuals, .f = ~ {
    # For each projection
    out <- apply(.x, 3, function(x) {
      # Filter data to retain instances of pregnant females being alive
      fem <- x[x[, "alive"] == 1 & x[, "female"] == 1 & x[, "cohort"] == 4 & x[, "reprod"] == 1, , drop = FALSE]
      # Remove NAs (from new individuals born during the simulation)
      fem <- fem[complete.cases(fem), , drop = FALSE]
      # Calculate the number of births 
      ifelse(nrow(fem) == 0, 0, sum(fem[, "birth"], na.rm = TRUE))
    })
    if(inherits(out, "list")) out <- do.call(c, out)
    out
  }) |> do.call(what = c)
  
  # **** Time spent in resting phase ****
  
  time.resting <- purrr::map(.x = narw.individuals, .f = ~ {
      lapply(X = seq_len(dim(.x)[3]), FUN = function(d) {
        tr <- .x[,,d]
        tr <- tr[tr[,"alive"] == 1 & tr[,"female"]==1 & tr[, "cohort"] == 6 & tr[, "reprod"] == 1,,drop =FALSE]
        tr <- tr[, "time_rest"]
        tr[tr == 0] <- NA
        unname(purrr::map_dbl(.x = split_byNA(tr), .f = ~max(.x)))
      }) |> do.call(what = c)
    }) |> do.call(what = c)

  # **** Inter-birth interval ****

  inter.birth <- purrr::map(.x = narw.individuals, .f = ~ {
    lapply(1:dim(.x)[3], FUN = function(i) {
      x <- .x[,,i]
      # Filter out NA records for individuals born during the simulation
      x <- x[complete.cases(x),,drop = FALSE]
      # Exclude juveniles
      birth.record <- x[x[, "alive"] == 1 & x[, "female"] == 1 & x[, "cohort"] > 2 & x[, "reprod"] == 1, , drop = FALSE]
      birth.record <- birth.record[, "birth"]
      if (length(birth.record) == 0 | all(birth.record == 0) | all(is.na(birth.record))) {
        NA
      } else {
        # Extract timeline of births
        birth.record <- birth.record[tapply(seq_along(birth.record), birth.record, min)[2]:tapply(seq_along(birth.record), birth.record, max)[2]]
        birth.record[birth.record == 1] <- NA
        unname(purrr::map_dbl(.x = split_byNA(birth.record), .f = ~ length(.x)))
      }
    }) |> do.call(what = c)
  }) |> do.call(what = c)
  
  # **** Non-reproductive females ****
  
  nonreprod.females <- purrr::map(.x = seq_len(n), .f = ~ {

    # Non-reproductive females (alive) by year
    nonrep <- lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      rep.status <- rep(0, yrs+1)
      rep.status[narw.individuals[[.x]][,,d][, "female"] == 1 & narw.individuals[[.x]][,,d][, "alive"] == 1 & narw.individuals[[.x]][,,d][, "reprod"] == 0] <- 1
      rep.status
    }) |> do.call(what = rbind) |> 
      colSums()
    
    # Total females (by year)
    totfem <- rowSums(narw.pop[.x,,5:8])
    totfem[1] <- sum(N_0[6:8])
    out <- nonrep/totfem
    # Correct cases where population has gone extinct - which leads to NAs here
    out[is.na(out)] <- 0
    out
    
  }) |> do.call(what = rbind)
  
  # narw.out <- purrr::map(.x = 1:n, .f = ~{
  #   reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
  #     dplyr::mutate(attr = mat.attribs[attr]) |> 
  #     tidyr::pivot_wider(names_from = attr, values_from = value) |> 
  #     dplyr::mutate(prj = .x) |> 
  #     dplyr::relocate(prj, .before = yr)
  # }) |> do.call(what = rbind) |> 
  #   data.table::data.table()
  
  # narw.out <- narw.out[is.finite(rowSums(narw.out)),]
  
  
# microbenchmark::microbenchmark(
# e1 = { 
#   
#   m1 <- apply(narw.pop, c(2,3), FUN = function(x) mean(x, na.rm = TRUE))
#   l1 <- apply(narw.pop, c(2,3), FUN = function(x) quantile(x, 0.025, na.rm = TRUE))
#   u1 <- apply(narw.pop, c(2,3), FUN = function(x) quantile(x, 0.975, na.rm = TRUE))
#   
#   test <- data.frame(year = current.yr + as.numeric(gsub("yr ", "", 0:yrs)))
#   
#   out <- purrr::map(.x = cohorts.proj$name,
#              .f = ~{
#               cbind(test, .x, m1[, .x], l1[, .x], u1[, .x])
#              }) |> do.call(what = "rbind")
#  
  
  # 
  # test <- data.frame(year = current.yr + as.numeric(gsub("yr ", "", 0:yrs)),
  #                         cohort = "Adults (female, resting)",
  #                         mean = apply(narw.pop, c(2), FUN = function(x) mean(x, na.rm = TRUE)),
  #                         lwr = apply(narw.pop, c(2), FUN = function(x) quantile(x, 0.025, na.rm = TRUE)),
  #                         uppr = apply(narw.pop, c(2), FUN = function(x) quantile(x, 0.975, na.rm = TRUE)))


# },
# e2 = { },
# times = 10
# )
  
  narw.df <- purrr::map(.x = cohorts.proj$name, .f = ~{
    tmp <- narw.pop[,,.x] |>
      tibble::as_tibble() |> 
      tibble::rownames_to_column(var = "prj")
    if(n == 1){
      tmp |>
        dplyr::rename(year = "prj") |> 
        dplyr::mutate(prj = 1, year = as.numeric(year) - 1) |>
        dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::rename(N = "value") |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x)) |> 
        dplyr::relocate(prj, .before = year)
    } else {
      tmp |>
        tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
        dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x))
    }
  }) |> do.call(what = rbind) |> 
    data.table::data.table()
 
  tot.df <- tibble::as_tibble(tot.pop) |> 
    tibble::rownames_to_column(var = "prj") |> 
    tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
    dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
    dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
  
  narw.conf <- narw.df[
    , list(
      mean = mean(N, na.rm = TRUE),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>
    dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)"
    )))
  
  tot.conf <- tot.df[
    , list(
      mean = mean(N, na.rm = TRUE),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)",
      "North Atlantic right whales"
    )))

  console(msg = "Summarizing outputs", suffix = tickmark())
  cat("---------------------------------\n")
  cat(paste0("Done | Time elapsed: ", run_time), "\n")
  
  outproj <- list(param = list(label = obj[[1]]$param$label,
                               yrs = yrs,
                               n = n,
                               current.yr = current.yr,
                               abort = abort.rate,
                               nonrep = nonrep,
                               phases = phase.names,
                               schedule = schedule),
                  init = list(N_0 = N_0,
                              reprod.fem = reprod.fem),
                  dat = list(ind = narw.individuals,
                             birth = list(tot = births.df, 
                                          perfemale = births.per.female, 
                                          inter = inter.birth),
                             death = deaths.df,
                             health = prjpreds,
                             nonrepfem = nonreprod.females,
                             rest = time.resting,
                             pop = narw.pop,
                             tot = tot.pop),
                  prj = list(proj = rbind(narw.df, tot.df) |> dplyr::mutate(prj = as.numeric(prj)), 
                             mean = rbind(narw.conf, tot.conf)))
  class(outproj) <- c("narwproj", class(outproj))
  return(outproj)
  
}
