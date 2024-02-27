#' Predictions of right whale abundance
#' 
#' Takes an object of class \code{narwsim} (as returned by \code{\link[narwind]{narw}}) and produces forecasts of right whale population size (by cohort) over a time horizon of interest. Prediction uncertainty is estimated from both process variance (i.e., replicate projections) and parameter uncertainty (i.e., statistical uncertainty in the relationships between individual body condition and health/survival, respectively). The latter is considered only if replicate coefficients have been sampled from the posteriors of the fitted survival and body condition GAMs using \code{\link[narwind]{augment}}.
#'
#' @param ... One or more objects of class \code{narwsim}.
#' @param n Integer. Number of replicate projections. Defaults to \code{100}.
#' @param yrs Integer. Time horizon, specified either as the desired number of years (from current) or the desired target end year. Defaults to \code{35}, which is commensurate with the average lifespan of a typical wind farm.
#' @param param Logical. If \code{TRUE}, prediction variance includes parameter uncertainty.  
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
                            param = TRUE,
                            piling = 1,
                            survival_prob = NULL,
                            min_gest = NULL,
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
  if(any(purrr::map_lgl(.x = obj, .f = ~is.null(.x$post)))){
    param <- FALSE
    warning("Posterior samples for terminal functions not available!\nRun the augment() function to estimate variance through posterior simulation.", sep = "")
  }

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
  
  # The initial population vector from NOAA is for 2019, so need to adjust the time horizon
  burn.in <- current.yr - 2019
  yrs <- yrs + burn.in
  piling <- piling + burn.in
  
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

  survbc.model <- purrr::map(.x = obj, .f = ~ list("surv" = .x$gam$fit$surv, "bc" = .x$gam$fit$bc)) |>
    purrr::set_names(nm = phase.names)
  
  # Inverse link functions
  ilink.mbc <- family(gam_gest)$linkinv
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
        # mbc = matrix(data = coef(gam_gest), nrow = 1, ncol = length(coef(gam_gest))),
        n = 1
      )
      
    } else { # If posterior samples are available
      
      list(
        surv = .x$post$samples$surv,
        bc = .x$post$samples$bc,
        # mbc = .x$post$samples$mbc,
        # Take only the median curve here
        # mbc = matrix(data = coef(gam_gest), nrow = 1, ncol = length(coef(gam_gest))),
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
  
  date_time()
  
  cat("––– Projections:\n\n")
  
  cat("+ Replicates: N =", formatC(n, big.mark = ","), "\n")
  cat("+ Burn-in:", burn.in, "years\n")
  cat("+ Horizon:", yrs - burn.in, "years\n")
  cat("+ Parameter uncertainty:", ifelse(param, "Yes", "No"), "\n")
  cat("\n––– Timeline:\n\n")
  proj_timeline(schedule, burn.in)
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
  
  # # Initialise array to store numbers of births
  # births <- array(
  #   data = NA, c(yrs + 1, n),
  #   dimnames = list(
  #     paste0("yr ", 0:yrs),
  #     paste0("prj ", 1:n)
  #   )
  # )
  
  # Initialise array to store numbers of abortion
  abortions <- array(
    data = NA, c(yrs + 1, n),
    dimnames = list(
      paste0("yr ", 0:yrs),
      paste0("prj ", 1:n)
    )
  )
  
  # Maximum allowable body condition
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
    # if(param){
    # 
    if(!is.null(min_gest)){
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * min_gest
    } else {
      Xp <- predict(gam_gest, data.frame(mass = narw.indiv[1, "tot_mass", ]), type = "lpmatrix")
      #       
      #       # mbc.sample <- sample(
      #       #   x = seq_len(mcmc[[initial.phase]]$n),
      #       #   size = nrow(Xp),
      #       #   replace = TRUE)
      #       
      #       mbc.curves <- mcmc[[initial.phase]]$mbc[1,]
      #       narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * ilink.mbc(rowSums(Xp*mbc.curves))  
      #       
      #     } else {
      
      
      
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * ilink.mbc(Xp %*% coef(gam_gest))
    }
      
    # }
    
    # Time spent in resting state
    narw.indiv[1, "time_rest", ] <- (cohort.vec == 6)

    # Calving events
    narw.indiv[1, "birth", ] <- as.numeric(cohort.vec == 4)
    
    if(param){
    
    # Survival probability
    Xp <- predict(survbc.model[[initial.phase]]$surv, data.frame(start_bc = bc, cohort = cohort.vec), type = "lpmatrix")
    surv.sample <- sample(
      x = seq_len(mcmc[[initial.phase]]$n),
      size = nrow(Xp),
      replace = TRUE)
    
    surv.curves <- mcmc[[initial.phase]]$surv[surv.sample,]
    if(!is.null(survival_prob)){
      narw.indiv[1, "p_surv", ] <- survival_prob
    } else {
      narw.indiv[1, "p_surv", ] <- ilink.survbc[["surv"]](rowSums(Xp*surv.curves))
    }

    
    } else {
      
      if(!is.null(survival_prob)){
        narw.indiv[1, "p_surv", ] <- survival_prob
      } else {
        narw.indiv[1, "p_surv", ] <- predict(survbc.model[[initial.phase]]$surv, 
                                             data.frame(start_bc = bc, cohort = cohort.vec), type = "response")
      }

      
    }

    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      if(tot.pop[prj, i-1] > 0) {
        
        # Current phase of development
        current.phase <- phase.names[schedule[i]+1]
        
        if(param){
          
          # New samples from the posteriors of model coefficients
          post.sample <- sample(
            x = seq_len(mcmc[[current.phase]]$n),
            size = nrow(Xp),
            replace = TRUE)
          
          bc.curves <- mcmc[[current.phase]]$bc[post.sample,]
          # mbc.curves <- mcmc[[current.phase]]$mbc[post.sample,]
        }
        
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
        
        # Compare females that are in resting phase at time t with females that were in pregnant phase at time t-1
        n.abortions <- sum((narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i-1, "cohort", ] == 4))
        
        # Record abortions
        abortions[i, prj] <- n.abortions
        
        #' ----------------------------
        # GROWTH
        #' ----------------------------
        
        # Increment length
        narw.indiv[i, "length", ] <- alive * age2length_vec(narw.indiv[i, "age", ])
        
        # Increment lean mass
        narw.indiv[i, "lean_mass", ] <- alive * length2mass_vec(narw.indiv[i, "length", ])
        
        # Predict new body condition from current body condition
        if (param) {
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
        
        if(!is.null(min_gest)){
          narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * min_gest
        } else {
        
        # Minimum body condition needed to successfully bring fetus to term without starving
        # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
        # if (param) {
          Xp <- predict(gam_gest, data.frame(mass = narw.indiv[i, "tot_mass", ]), type = "lpmatrix")
        #   narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * ilink.mbc(rowSums(Xp * mbc.curves))
        # } else {
          narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * 
            ilink.mbc(Xp %*% coef(gam_gest))
        # }
        }
          
        # Birth of new calf, conditional on the mother being alive and in pregnant state
        narw.indiv[i, "birth", ] <- alive * (narw.indiv[i, "cohort", ] == 4)
        new.births <- sum(narw.indiv[i, "birth", ])
        # births[i, prj] <- new.births
        
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
        
        if(param){
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
        
        if(!is.null(survival_prob)){
          narw.indiv[i, "p_surv", ] <- survival_prob
        } else {
          narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * ilink.survbc$surv(rowSums(Xp*surv.curves))
        }
        
 
        } else {
          
          if(!is.null(survival_prob)){
            narw.indiv[i, "p_surv", ] <- survival_prob
          } else {
          narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * predict(survbc.model[[current.phase]]$surv,
                                                                          data.frame(
                                                                            start_bc = narw.indiv[i, "bc", ],
                                                                            cohort = narw.indiv[i, "cohort", ]
                                                                          ), type = "response")
          }
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
  
  # Survival probability ----
  
  prjpreds <- purrr::map(.x = seq_len(n), .f = ~ {
    out <- apply(narw.individuals[[.x]], 3, function(x) x[x[, "alive"] == 1 & !is.na(x[, "alive"]), c("cohort", "p_surv", "bc", "min_bc"), drop = FALSE])
    if(inherits(out, "list")) out <- do.call(rbind, out)
     out <- cbind(out, .x, rep(0:yrs, length.out = nrow(out)))
     colnames(out)[5:6] <- c("prj", "year")
     out
  }) |> do.call(what = rbind) |> 
    data.table::as.data.table()
  prjpreds[, year:= 2019 + year]
  
  # Total births ----
  
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
  births.df[, year:= 2019 + year]
  
  # Total deaths ----
  
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
  
  deaths.df[, year:= 2019 + year]
  n.deaths <- deaths.df[, list(n = c(death[1], diff(death))), prj]
  deaths.df$death <- n.deaths$n
  
  # Number of births per female ----
  
  births.per.female <- purrr::map(.x = seq_len(n), .f = ~ {
    lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      data.table::data.table(
        prj = .x, whale = d,
        nbirths = sum(narw.individuals[[.x]][, "birth", d], na.rm = TRUE)
      )
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
    

  # Time spent in resting phase ----
  
  time.resting <- purrr::map(.x = seq_len(n), .f = ~ {
    lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      x <- narw.individuals[[.x]][,,d]
      x <- x[x[,"reprod"] == 1,, drop = FALSE]
      if(nrow(x)<=1){
        NULL
      } else {
        tr <- split_NAzero(x[, "time_rest"])
        if (length(tr) == 0) {
          NULL
        } else {
          tr <- tr |>
            unname() |>
            sapply(FUN = function(x) x[length(x)])
          
          data.table::data.table(
            prj = .x, whale = d,
            year = 2019 + as.numeric(gsub("yr ", "", names(tr))),
            t_rest = tr)
        }
      }
    }) |> data.table::rbindlist()
 }) |> data.table::rbindlist()

  # Inter-birth interval ----

  inter.birth <- purrr::map(.x = seq_len(n), .f = ~ {
    lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      x <- narw.individuals[[.x]][,,d]
      birth.record <- x[x[, "alive"] == 1, "birth"]
      if (length(birth.record) == 0 | sum(birth.record, na.rm = TRUE) <=1) {
        NULL
      } else {
        birth.record <- birth.record[tapply(seq_along(birth.record), birth.record, min)[2]:tapply(seq_along(birth.record), birth.record, max)[2]]
        birth.record[birth.record == 1] <- NA
        data.table::data.table(
          prj = .x, whale = d,
          ibi = sapply(split_NA(birth.record), length) + 1
        )
      }
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
  
  # Non-reproductive females ----
  
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
  
  # Abortions ----
 abort.df <- as.data.frame.table(abortions[2:nrow(abortions), ], responseName = "value") |>
   dplyr::rename(year = Var1, prj = Var2, n_abort = value) |>
   dplyr::mutate(
     year = gsub(pattern = "yr ", replacement = "", year),
     prj = gsub(pattern = "prj ", replacement = "", prj)
   ) |>
   data.table::as.data.table() |>
   dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year)))
  
 # Pregnant females ----
 preg.fem <- as.data.frame.table(narw.pop[, , cohorts.proj[id == 4, name]]) |>
   dplyr::rename(year = Var2, prj = Var1, n_pregfem = Freq) |>
   dplyr::mutate(
     year = gsub(pattern = "yr ", replacement = "", year),
     prj = gsub(pattern = "prj ", replacement = "", prj)
   ) |>
   data.table::as.data.table() |>
   dplyr::filter(year > 0) |>
   dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year)))
 
  narw.df <- purrr::map(.x = cohorts.proj$name, .f = ~{
    tmp <- narw.pop[,,.x] |>
      tibble::as_tibble() |> 
      tibble::rownames_to_column(var = "prj")
    if(n == 1){
      tmp |>
        dplyr::rename(year = "prj") |> 
        dplyr::mutate(prj = 1, year = as.numeric(year) - 1) |>
        dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::rename(N = "value") |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x)) |> 
        dplyr::relocate(prj, .before = year)
    } else {
      tmp |>
        tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
        dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x))
    }
  }) |> do.call(what = rbind) |> 
    data.table::data.table()
 
  # Population totals ----
  
  tot.df <- tibble::as_tibble(tot.pop) |> 
    tibble::rownames_to_column(var = "prj") |> 
    tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
    dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year))) |> 
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
  
  # Output list ----
  
  outproj <- list(param = list(label = obj[[1]]$param$label,
                               yrs = yrs - burn.in,
                               burn = burn.in,
                               n = n,
                               cohorts = cohorts.proj,
                               current.yr = current.yr,
                               abort = abort.rate,
                               nonrep = nonrep,
                               phases = phase.names,
                               schedule = schedule),
                  init = list(N_0 = N_0,
                              cohort = cohort.vec,
                              reprod.fem = reprod.fem),
                  dat = list(ind = narw.individuals,
                             birth = list(tot = births.df, 
                                          perfemale = births.per.female, 
                                          inter = inter.birth),
                             death = deaths.df,
                             health = prjpreds,
                             abort = abort.df,
                             nonrepfem = nonreprod.females,
                             pregfem = preg.fem,
                             rest = time.resting,
                             pop = narw.pop,
                             tot = tot.pop),
                  proj = list(tbl = rbind(narw.df, tot.df) |> dplyr::mutate(prj = as.numeric(prj)), 
                              mean = rbind(narw.conf, tot.conf)))
  class(outproj) <- c("narwproj", class(outproj))
  return(outproj)
  
}
