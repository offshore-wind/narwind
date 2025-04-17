#' Predictions of right whale abundance
#'
#' Takes an object of class \code{narwsim} (as returned by \code{\link[narwind]{narw}}) and produces forecasts of right whale population size (by cohort) over a time horizon of interest. Prediction uncertainty is estimated from both process variance (i.e., replicate projections) and parameter uncertainty (i.e., statistical uncertainty in the relationships between individual body condition and health/survival, respectively). The latter is considered only if replicate coefficients have been sampled from the posteriors of the fitted survival and body condition GAMs using \code{\link[narwind]{augment}}.
#'
#' @param ... One or more objects of class \code{narwsim}.
#' @param n Integer. Number of replicate projections. Defaults to \code{100}.
#' @param yrs Integer. Time horizon, specified either as the desired number of years (from current) or the desired target end year. Defaults to \code{35}, which is commensurate with the average lifespan of a typical wind farm.
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

  #' ------------------------------------------------------------
  # Function checks ----
  #' ------------------------------------------------------------

  # For package development only
  param <- TRUE
  noaa <- FALSE
  survival_prob <- NULL
  min_gest <- NULL

  args <- list(...)

  # Identify and extract <narwsim> objects
  which.obj <- which(purrr::map_lgl(.x = args, .f = ~ inherits(.x, "narwsim")))
  if (length(which.obj) == 0) stop("No object of class <narwsim> found.")
  obj <- args[which.obj]

  # Check whether posterior samples are available
  if (any(purrr::map_lgl(.x = obj, .f = ~ is.null(.x$post)))) {
    param <- FALSE
    warning("Posterior samples for terminal functions not available!\nRun the augment() function to estimate variance through posterior simulation.", sep = "")
  }

  # Check whether terminal functions are available
  if (sum(sapply(X = obj, FUN = function(o) {
    sum(purrr::map_lgl(.x = o$gam$fit, .f = ~ any(is.null(unlist(.x)))))
  }) > 0)) {
    stop("Sample sizes insufficient to reliably estimate survival and health functions. Re-run narw() with a higher <nsim>.")
  }

  # Check that all 6 population cohorts have been simulated
  if (any(!purrr::map_dbl(.x = obj, .f = ~ length(.x$param$cohort)) == 6)) stop("Missing cohorts! Cannot generate projections")

  if (length(obj) > 3) stop("Too many objects - maximum of 3 scenarios (baseline, construction, O&M) allowed.")

  if (yrs <= 0) stop("<yrs> must be a positive integer")
  if (piling <= 0) stop("<piling> must be a positive integer")

  # If only one object is supplied, must be baseline conditions
  if (length(obj) == 1) {
    # Change <piling> if only baseline has been provided
    if (obj[[1]]$scenario$phase == 0 & piling < yrs) piling <- yrs + 1
    if (obj[[1]]$scenario$phase == 1) stop("Missing O&M scenario")
    if (obj[[1]]$scenario$phase == 2) stop("Missing baseline scenario")
    # Change <piling> if only O&M has been provided
    # if(obj[[1]]$scenario$phase == 2 & piling < yrs) piling <- yrs + 1
  }

  phases <- sapply(X = obj, FUN = function(o) o$scenario$phase)

  # If multiple objects are supplied
  if (length(obj) > 1 & sum(duplicated(phases)) > 0) stop("Duplicate scenarios detected.")

  if (length(obj) > 1 & piling < yrs) {
    if (piling > 1 & all(phases %in% 1:2)) stop("Missing baseline scenario")
    if (all(phases %in% c(0, 2))) piling <- yrs + 1
    if (all(phases %in% 0:1)) stop("Missing O&M scenario")
  }

  # Current year
  current.yr <- lubridate::year(lubridate::now())

  # User can specify the <yrs> argument either as a number of years from now or as a target end year
  if (yrs > 2000) yrs <- yrs - current.yr
  if (piling > 2000) piling <- piling - current.yr

  # The initial population vector from NOAA is for 2019, so need to adjust the time horizon
  burn.in <- current.yr - 2019
  yrs <- yrs + burn.in
  piling <- piling + burn.in

  # Cohort information (augmented to separate out male/female calves)
  # cohortID <- obj[[1]]$param$cohort
  cohorts.proj <- obj[[1]]$param$cohorts |>
    dplyr::slice(1) |>
    dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |>
    dplyr::bind_rows(obj[[1]]$param$cohorts) |>
    dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
  cohorts.proj <- cohorts.proj[c(1, 3, 5, 2, 4, 6, 7, 8)]

  # Attributes to monitor during projection
  mat.attribs <- c(
    "alive", "cohort", "female", "age", "length", "tot_mass", "lean_mass", "bc",
    "p_surv", "min_bc", "time_rest", "birth", "reprod"
  )

  #' ------------------------------------------------------
  # PHASES ----
  #' ------------------------------------------------------

  # Define schedule of wind farm development
  phase.names <- c("base", "const", "ops")[phases + 1]
  schedule <- rep(0, yrs + 1)
  if (piling < yrs) {
    schedule[piling + 1] <- 1
    schedule[(piling + 2):(yrs + 1)] <- 2
  } else if (piling > yrs) {
    if (sum(phases) > 0) {
      schedule[(burn.in + 1):(yrs + 1)] <- 2
    }
  }
  names(schedule) <- paste0("yr ", 0:yrs)

  #' ------------------------------------------------------
  # TERMINAL FUNCTIONS ----
  #' ------------------------------------------------------

  # Fitted models
  survbc.model <- purrr::map(.x = obj, .f = ~ list("surv" = .x$gam$fit$surv, "bc" = .x$gam$fit$bc)) |>
    purrr::set_names(nm = phase.names)

  survbc.model_calf <- purrr::map(.x = obj, .f = ~ list("surv" = .x$gam$fit_calf$surv, "bc" = .x$gam$fit_calf$bc)) |>
    purrr::set_names(nm = phase.names)

  # Inverse link functions
  ilink.mbc <- family(gam_gest)$linkinv

  ilink.survbc <- purrr::map(
    .x = obj[1],
    .f = ~ list(
      "surv" = family(.x$gam$fit$surv)$linkinv,
      "bc" = family(.x$gam$fit$bc)$linkinv
    )
  )[[1]]

  ilink.survbc_calf <- purrr::map(
    .x = obj[1],
    .f = ~ list(
      "surv" = family(.x$gam$fit_calf$surv)$linkinv,
      "bc" = family(.x$gam$fit_calf$bc)$linkinv
    )
  )[[1]]

  # Posterior samples
  mcmc <- purrr::map2(
    .x = obj,
    .y = survbc.model,
    .f = ~ {

      # If posterior samples are not available
      if (is.null(.x$post)) {
        list(
          surv = matrix(data = coef(.y$surv), nrow = 1, ncol = length(coef(.y$surv))),
          bc = matrix(data = coef(.y$bc), nrow = 1, ncol = length(coef(.y$bc))),
          n = 1
        )
      } else { # If posterior samples are available

        list(
          surv = .x$post$samples$all$surv,
          bc = .x$post$samples$all$bc,
          n = .x$post$nsamples
        )
      }
    }
  ) |> purrr::set_names(nm = phase.names)

  mcmc_calf <- purrr::map2(
    .x = obj,
    .y = survbc.model_calf,
    .f = ~ {

      # If posterior samples are not available
      if (is.null(.x$post)) {
        list(
          surv = matrix(data = coef(.y$surv), nrow = 1, ncol = length(coef(.y$surv))),
          bc = matrix(data = coef(.y$bc), nrow = 1, ncol = length(coef(.y$bc))),
          n = 1
        )
      } else { # If posterior samples are available

        list(
          surv = .x$post$samples$calf$surv,
          bc = .x$post$samples$calf$bc,
          n = .x$post$nsamples
        )
      }
    }
  ) |> purrr::set_names(nm = phase.names)


  # Theta parameters from the Beta GAMs of body condition
  theta <- purrr::map2(
    .x = survbc.model,
    .y = survbc.model_calf,
    .f = ~ c(
      .x$bc$family$getTheta(),
      .y$bc$family$getTheta()
    )
  ) |>
    purrr::set_names(nm = phase.names)

  # Theta parameter from the Beta GAM of min body condition needed for gestation
  theta_gest <- gam_gest$family$getTheta()

  #' ------------------------------------------------------
  # Abortion rate ----
  #' ------------------------------------------------------

  abort.rate <- purrr::map_dbl(.x = obj, .f = ~ sum(.x$abort[, abort]) / .x$param$nsim)
  names(abort.rate) <- phase.names

  #' ------------------------------------------------------
  # INITIALIZATION ----
  #' ------------------------------------------------------

  date_time()

  cat("––– Projections:\n\n")

  cat("+ Replicates: N =", formatC(n, big.mark = ","), "\n")
  cat("+ Burn-in:", burn.in, "years\n")
  cat("+ Horizon:", yrs - burn.in, "years [from current]\n")
  # cat("+ Parameter uncertainty:", ifelse(param, "Yes", "No"), "\n")
  cat("\n––– Timeline:\n\n")
  proj_timeline(schedule, burn.in)
  cat("––– Execution:\n\n")

  # Define initial population vector
  # Adjust predictions from NOAA model to reflect a more stable age/stage structure
  # This is achieved by running adjust_popvec() on a test projection (noaa = TRUE)
  if (noaa) {
    # Note that the sum of cohort abundances in N_0 differs slightly
    # from the median population size under NOAA's most recent PVA
    # We therefore correct N_0 here to ensure that abundance estimates match
    N_0 <- N0[["yr2019"]]
    N_0[N_0 > 0] <- N_0[N_0 > 0] - 1
    N_0[N_0 > 200] <- N_0[N_0 > 200] - 1
  } else {
    N_0 <- c(9, 45, 168, 6, 39, 17, 15, 57)
  }
  names(N_0) <- cohorts.proj[, name]

  # Check if there are the same numbers of calves and lactating females
  if (!sum(N_0[1], N_0[4]) == N_0[7]) stop("Initial conditions do not match!")

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

  # Initialize array to store numbers of abortion
  abortions <- array(
    data = NA, c(yrs + 1, n),
    dimnames = list(
      paste0("yr ", 0:yrs),
      paste0("prj ", 1:n)
    )
  )

  # Maximum allowable body condition
  maxbc <- find_maxBC()

  #' ------------------------------------------------------
  # PROJECTIONS ----
  #' ------------------------------------------------------

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

  schedule.phases <- schedule + 1
  if (sum(phases) == 2) schedule.phases[schedule.phases == 3] <- 2

  for (prj in 1:n) {

    # Update progress bar
    if (progress) pb$tick()

    # Initial development phase
    initial.phase <- phase.names[schedule.phases[1]]

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
    narw.indiv[, "female", 1:N_0[1]] <- 0 #  -- Calves (male)
    fem <- which(cohort.vec == 0) #  -- Calves (female)
    narw.indiv[, "female", fem[fem > N_0[1]]] <- 1 #  -- Juveniles and adults
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

    # Initial body condition - from Christiansen et al. (2020) and (2022)
    # Values obtained using the init_bc() function
    bc_dist <- data.table::data.table(
      cohort = unique(obj[[schedule.phases[1]]]$gam$dat$cohort),
      mean = c(0.35884621, 0.2904040, 0.2904040, 0.27617455, 0.27617455, 0.30989147, 0.27617455),
      sd = c(0.07788818, 0.0681619, 0.0681619, 0.06969136, 0.06969136, 0.06917943, 0.06969136)
    )

    bc <- sapply(X = cohort.vec, FUN = function(x) {
      clamp(rnorm(1, mean = bc_dist[cohort == x, mean], sd = bc_dist[cohort == x, sd]), maxbc)
    })

    narw.indiv[1, "bc", ] <- bc

    lac.f <- sample(which(cohort.vec == 5))

    # Overwrite starting BC values with predictions from the calf model
    if (param) {
      Xp_bc0 <- predict(survbc.model_calf[[initial.phase]]$bc,
        newdata = data.frame(start_bc = bc[lac.f]),
        type = "lpmatrix"
      )

      bc.sample0 <- sample(
        x = seq_len(mcmc_calf[[initial.phase]]$n),
        size = nrow(Xp_bc0),
        replace = TRUE
      )

      bc.curves0 <- mcmc_calf[[initial.phase]]$bc[bc.sample0, ]

      bc_mu0 <- as.numeric(clamp(ilink.survbc_calf$bc(rowSums(Xp_bc0 * bc.curves0)), maxbc))
      bc_var0 <- bc_mu0 * (1 - bc_mu0) / (1 + theta[[initial.phase]][2])
      beta_pars0 <- estBetaParams(bc_mu0, bc_var0)

      narw.indiv[1, "bc", cohort.vec == 0] <- clamp(rbeta(nrow(beta_pars0), beta_pars0[, "shape1"], beta_pars0[, "shape2"]), maxbc)
    } else {
      narw.indiv[1, "bc", cohort.vec == 0] <- predict(survbc.model_calf[[initial.phase]]$bc,
        data.frame(start_bc = bc[lac.f]),
        type = "response"
      )
    }

    # narw.indiv[1, "bc_mother", cohort.vec==0] <- bc[lac.f]

    # Total mass
    narw.indiv[1, "tot_mass", ] <- mass / (1 - bc)

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
    if (!is.null(min_gest)) {
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * min_gest
    } else {
      Xp_gest0 <- predict(gam_gest, data.frame(mass = narw.indiv[1, "tot_mass", ]), type = "lpmatrix")
      
      gest.sample0 <- sample(
        x = seq_len(nrow(gam_gestation$post)),
        size = nrow(Xp_gest0),
        replace = TRUE
      )
      
      gest.curves0 <- gam_gestation$post[gest.sample0, ]
      minbc_mu0 <- as.numeric(ilink.mbc(rowSums(Xp_gest0 * gest.curves0)))
      
      # minbc_mu <- as.numeric(ilink.mbc(Xp_gest0 %*% coef(gam_gest))) # Predicted mean from Beta GAM
      minbc_var0 <- minbc_mu0 * (1 - minbc_mu0) / (1 + theta_gest) # Estimated variance from Beta GAM
      beta_pars0 <- estBetaParams(minbc_mu0, minbc_var0) # Get shape parameters from mean and variance
      narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * rbeta(nrow(beta_pars0), beta_pars0[, "shape1"], beta_pars0[, "shape2"])
    }

    # Time spent in resting state
    narw.indiv[1, "time_rest", ] <- (cohort.vec == 6)

    # Calving events
    narw.indiv[1, "birth", ] <- as.numeric(cohort.vec == 4)

    if (!is.null(survival_prob)) {
      narw.indiv[1, "p_surv", ] <- survival_prob
    } else {
      if (param) {

        # Survival probability
        Xp_surv0 <- predict(survbc.model[[initial.phase]]$surv,
          data.frame(start_bc = bc, cohort = cohort.vec),
          type = "lpmatrix"
        )
        surv.sample <- sample(
          x = seq_len(mcmc[[initial.phase]]$n),
          size = nrow(Xp_surv0),
          replace = TRUE
        )

        surv.curves <- mcmc[[initial.phase]]$surv[surv.sample, ]
        narw.indiv[1, "p_surv", ] <- ilink.survbc[["surv"]](rowSums(Xp_surv0 * surv.curves))
      } else {
        narw.indiv[1, "p_surv", ] <- predict(survbc.model[[initial.phase]]$surv,
          data.frame(start_bc = bc, cohort = cohort.vec),
          type = "response"
        )
      }
    }



    # Overwrite values for calves based on calf model
    if (!is.null(survival_prob)) {
      narw.indiv[1, "p_surv", cohort.vec == 0] <- survival_prob
    } else {
      if (param) {
        Xp_surv0_calf <- predict(survbc.model_calf[[initial.phase]]$surv, data.frame(start_bc = bc[lac.f]), type = "lpmatrix")
        surv.sample <- sample(
          x = seq_len(mcmc_calf[[initial.phase]]$n),
          size = nrow(Xp_surv0_calf),
          replace = TRUE
        )

        surv.curves <- mcmc_calf[[initial.phase]]$surv[surv.sample, ]
        narw.indiv[1, "p_surv", cohort.vec == 0] <- ilink.survbc_calf[["surv"]](rowSums(Xp_surv0_calf * surv.curves))
      } else {
        narw.indiv[1, "p_surv", cohort.vec == 0] <- predict(survbc.model_calf[[initial.phase]]$surv,
          data.frame(start_bc = bc[lac.f]),
          type = "response"
        )
      }
    }


    #' ------------------------------------------------------
    # Loop over years
    #' ------------------------------------------------------

    for (i in 2:(yrs + 1)) {
      if (tot.pop[prj, i - 1] > 0) {

        # Current phase of development
        current.phase <- phase.names[schedule.phases[i]]

        if (param) {

          # New samples from the posteriors of model coefficients
          post.sample <- sample(
            x = seq_len(mcmc[[current.phase]]$n),
            size = dim(narw.indiv)[3],
            replace = TRUE
          )

          bc.curves <- mcmc[[current.phase]]$bc[post.sample, ]
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
        narw.indiv[i, "age", ] <- alive * (narw.indiv[i - 1, "age", ] + 1)

        # Sex remains the same
        narw.indiv[i, "female", ] <- narw.indiv[i - 1, "female", ]

        #' ----------------------------
        # COHORTS
        #' ----------------------------

        # Maturity - transitions between cohorts
        narw.indiv[i, "cohort", ] <-
          increment_cohort(
            alive = alive,
            cohort = narw.indiv[i - 1, "cohort", ],
            age = narw.indiv[i, "age", ],
            female = narw.indiv[i, "female", ],
            bc = narw.indiv[i - 1, "bc", ],
            min_bc = narw.indiv[i - 1, "min_bc", ],
            reprod = narw.indiv[i - 1, "reprod", ],
            abort = abort.rate[current.phase]
          )

        # Compare females that are in resting phase at time t with females that were in pregnant phase at time t-1
        n.abortions <- sum((narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i - 1, "cohort", ] == 4))

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
          Xp_BC <- predict(survbc.model[[current.phase]]$bc,
            data.frame(
              start_bc = narw.indiv[i - 1, "bc", ],
              cohort = narw.indiv[i - 1, "cohort", ]
            ),
            type = "lpmatrix"
          )

          bc_mu <- as.numeric(clamp(ilink.survbc$bc(rowSums(Xp_BC * bc.curves)), maxbc))
          bc_var <- bc_mu * (1 - bc_mu) / (1 + theta[[current.phase]][1])
          beta_pars <- estBetaParams(bc_mu, bc_var)
          narw.indiv[i, "bc", ] <- alive * clamp(rbeta(nrow(beta_pars), beta_pars[, "shape1"], beta_pars[, "shape2"]), maxbc)
        } else {
          narw.indiv[i, "bc", ] <- alive * clamp(
            predict(survbc.model[[current.phase]]$bc,
              data.frame(
                start_bc = narw.indiv[i - 1, "bc", ],
                cohort = narw.indiv[i - 1, "cohort", ]
              ),
              type = "response"
            ), maxbc
          )
        }

        # Increment total mass
        narw.indiv[i, "tot_mass", ] <- alive * narw.indiv[i, "lean_mass", ] / (1 - narw.indiv[i, "bc", ])

        #' ----------------------------
        # REPRODUCTION
        #' ----------------------------

        # Reproductive females
        narw.indiv[i, "reprod", ] <- alive * narw.indiv[i - 1, "reprod", ]

        # Time spent in resting phase
        narw.indiv[i, "time_rest", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i - 1, "time_rest", ] + 1)

        if (!is.null(min_gest)) {
          narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * min_gest
        } else {

          # Minimum body condition needed to successfully bring fetus to term without starving
          # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
          Xp_gest <- predict(gam_gest, data.frame(mass = narw.indiv[i, "tot_mass", ]), type = "lpmatrix")

          gest.sample <- sample(
            x = seq_len(nrow(gam_gestation$post)),
            size = nrow(Xp_gest),
            replace = TRUE
          )
          
          gest.curves <- gam_gestation$post[gest.sample, ]
          minbc_mu <- as.numeric(ilink.mbc(rowSums(Xp_gest * gest.curves)))
          
          # minbc_mu <- as.numeric(ilink.mbc(Xp_gest %*% coef(gam_gest))) # Predicted mean from Beta GAM
          
          minbc_var <- minbc_mu * (1 - minbc_mu) / (1 + theta_gest) # Estimated variance from Beta GAM
          beta_pars <- estBetaParams(minbc_mu, minbc_var) # Get shape parameters from mean and variance
          narw.indiv[i, "min_bc", ] <- clamp(alive * (narw.indiv[i, "cohort", ] == 6) *
            rbeta(nrow(beta_pars), beta_pars[, "shape1"], beta_pars[, "shape2"]), maxbc)
        }

        # Birth of new calf, conditional on the mother being alive and in pregnant state
        narw.indiv[i, "birth", ] <- alive * (narw.indiv[i, "cohort", ] == 4)
        new.births <- sum(narw.indiv[i, "birth", ])

        # Number of non-reproductive females in the population
        Nr <- nrow(t(narw.indiv[i, , ])[t(narw.indiv[i, , ])[, "alive"] == 1 & t(narw.indiv[i, , ])[, "female"] == 1 & t(narw.indiv[i, , ])[, "reprod"] == 0 & t(narw.indiv[i, , ])[, "cohort"] %in% c(2, 4, 5, 6), ])

        # Total abundance of females
        Nt <- nrow(t(narw.indiv[i, , ])[t(narw.indiv[i, , ])[, "alive"] == 1 & t(narw.indiv[i, , ])[, "female"] == 1 & t(narw.indiv[i, , ])[, "cohort"] %in% c(2, 4, 5, 6), ])

        
        if (new.births > 0) {

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
          if (Nb > 0 & sum(new.females) > 0) newborns.nonreprod[sample(which(new.females == 1), size = Nb)] <- 0
          
          # Create array for newborns
          new.calves <- array(
            data = NA, c(yrs + 1, length(mat.attribs), new.births),
            dimnames = list(
              paste0("yr ", 0:yrs), mat.attribs,
              paste0("whale ", ((dim(narw.indiv)[3] + 1):(dim(narw.indiv)[3] + new.births)))
            )
          )

          # BC of new calves
          if (param) {
            Xp_bc_calf <- predict(survbc.model_calf[[current.phase]]$bc,
              newdata = data.frame(start_bc = narw.indiv[i, "bc", alive * narw.indiv[i, "cohort", ] == 4]),
              type = "lpmatrix"
            )

            bc.sample_calf <- sample(
              x = seq_len(mcmc_calf[[current.phase]]$n),
              size = new.births,
              replace = TRUE
            )

            bc.curves_calf <- mcmc_calf[[current.phase]]$bc[bc.sample_calf, ]

            bc_mu_calf <- as.numeric(clamp(ilink.survbc_calf$bc(rowSums(Xp_bc_calf * bc.curves_calf)), maxbc))
            bc_var_calf <- bc_mu_calf * (1 - bc_mu_calf) / (1 + theta[[current.phase]][2])
            beta_pars_calf <- estBetaParams(bc_mu_calf, bc_var_calf)

            bcpreds_calf <- clamp(rbeta(nrow(beta_pars_calf), beta_pars_calf[, "shape1"], beta_pars_calf[, "shape2"]), maxbc)
          } else {
            bcpreds_calf <- predict(survbc.model_calf[[current.phase]]$bc,
              data.frame(start_bc = narw.indiv[i, "bc", alive * narw.indiv[i, "cohort", ] == 4]),
              type = "response"
            )
          }


          # Survival probability of new calves

          if (!is.null(survival_prob)) {
            survpreds_calf <- survival_prob
          } else {
            if (param) {
              Xp_surv_calf <- predict(survbc.model_calf[[current.phase]]$surv,
                data.frame(start_bc = narw.indiv[i, "bc", alive * narw.indiv[i, "cohort", ] == 4]),
                type = "lpmatrix"
              )

              surv.sample_calf <- sample(
                x = seq_len(mcmc_calf[[current.phase]]$n),
                size = new.births,
                replace = TRUE
              )

              surv.curves_calf <- mcmc_calf[[current.phase]]$surv[surv.sample_calf, ]

              survpreds_calf <- ilink.survbc_calf[["surv"]](rowSums(Xp_surv_calf * surv.curves_calf))
            } else {
              survpreds_calf <- predict(survbc.model_calf[[current.phase]]$surv,
                data.frame(start_bc = narw.indiv[i, "bc", narw.indiv[i, "cohort", ] == 4]),
                type = "response"
              )
            }
          }

          new.calves[i, , ] <- add_calf(
            new.births,
            mat.attribs,
            new.females,
            newborns.nonreprod,
            bcpreds_calf,
            survpreds_calf
          )
          
          narw.indiv <- abind::abind(narw.indiv, new.calves, along = 3)
        }

        
        #' ----------------------------
        # SURVIVAL
        #' ----------------------------

        if (!is.null(survival_prob)) {
          narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * survival_prob
        } else {
          if (param) {
            # Predict survival probability based on body condition
            XpSU <- predict(survbc.model[[current.phase]]$surv,
              data.frame(
                start_bc = narw.indiv[i, "bc", ],
                cohort = narw.indiv[i, "cohort", ]
              ),
              type = "lpmatrix"
            )

            post.sample <- sample(
              x = seq_len(mcmc[[current.phase]]$n),
              size = nrow(XpSU),
              replace = TRUE
            )

            surv.curves <- mcmc[[current.phase]]$surv[post.sample, ]
            narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * ilink.survbc$surv(rowSums(XpSU * surv.curves))
          } else {
            narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * predict(survbc.model[[current.phase]]$surv,
              data.frame(
                start_bc = narw.indiv[i, "bc", ],
                cohort = narw.indiv[i, "cohort", ]
              ),
              type = "response"
            )
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

  #' ------------------------------------------------------
  # SUMMARY ----
  #' ------------------------------------------------------

  console(msg = "Running projections", suffix = tickmark())
  console(msg = "Summarizing outputs")

  survival <- purrr::map(.x = seq_len(n), .f = ~ {
    narw.tbl <- as.data.frame.table(narw.individuals[[.x]]) |>
      dplyr::rename(year = Var1, param = Var2, whale = Var3, value = Freq) |>
      dplyr::mutate(
        year = as.numeric(gsub("yr ", "", year)),
        whale = as.numeric(gsub("whale ", "", whale))
      ) |>
      tidyr::pivot_wider(names_from = param, values_from = value) |>
      dplyr::mutate(prj = .x) |>
      dplyr::select(prj, year, whale, cohort, alive, female) |>
      dplyr::mutate(year = 2019 + year) |>
      data.table::as.data.table()

    # Calculate survival by shifting the alive column up one row
    narw.tbl[, survival := dplyr::lead(alive), list(prj, whale)]

    # Set survival to NA when alive is NA
    narw.tbl[is.na(alive), ]$survival <- NA
    narw.tbl <- narw.tbl[alive == 1 & !is.na(alive) & !is.na(survival), ]

    narw.tbl[
      year < max(narw.tbl$year) & !is.na(cohort),
      list(
        alive = sum(survival, na.rm = TRUE),
        total = .N,
        surv = sum(survival, na.rm = TRUE) / .N
      ),
      list(prj, year, cohort)
    ]
  }) |> data.table::rbindlist()


  # Health ----

  prjpreds <- purrr::map(.x = seq_len(n), .f = ~ {

    # Extract data
    out <- apply(narw.individuals[[.x]], 3, function(x) {
      x <- cbind(x, dplyr::lead(x[, "bc"]))
      colnames(x)[ncol(x)] <- "end_bc"
      x[x[, "alive"] == 1 & !is.na(x[, "alive"]),
        c("cohort", "end_bc", "min_bc", "female"),
        drop = FALSE
      ]
    })

    # Add whale IDs
    whaleID <- purrr::map(.x = seq_len(dim(narw.individuals[[.x]])[3]), .f = ~
      rep(.x, each = nrow(out[[.x]])))

    # And corresponding years
    yy <- purrr::map(.x = seq_len(dim(narw.individuals[[.x]])[3]), .f = ~
      as.numeric(gsub(pattern = "yr ", replacement = "", x = rownames(out[[.x]]))))

    out <- purrr::pmap(.l = list(.x, yy, whaleID, out), function(first, second, third, fourth) {
      cbind(prj = first, year = second, whale = third, fourth)
    })

    if (inherits(out, "list")) out <- do.call(rbind, out)
    out
  }) |>
    do.call(what = rbind) |>
    data.table::as.data.table() |>
    dplyr::rename(bc = end_bc)
  prjpreds[, year := 2019 + year]
  prjpreds <- prjpreds[order(prj, whale, year)]

  # Total births ----

  births.df <- purrr::map(.x = 1:n, .f = ~ {
    m <- matrix(rowSums(narw.individuals[[.x]][2:(yrs + 1), "birth", ], na.rm = TRUE), ncol = 1)
    colnames(m) <- .x
    m
  }) |>
    do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |>
    dplyr::select(prj, year, birth) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |>
    dplyr::mutate(prj = as.numeric(prj))
  births.df[, year := 2019 + year]

  # Total deaths ----
  deaths.df <- purrr::map(.x = 1:n, .f = ~ {
    m <- matrix(apply(
      X = narw.individuals[[.x]][2:(yrs + 1), "alive", ],
      MARGIN = 1,
      FUN = function(x) {
        r <- x[!is.na(x)]
        r <- sum(r == 0)
        r
      }
    ), ncol = 1)
    colnames(m) <- .x
    m
  }) |>
    do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
    dplyr::select(prj, year, death) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |>
    dplyr::mutate(prj = as.numeric(prj))

  deaths.df[, year := 2019 + year]
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
      x <- narw.individuals[[.x]][, , d]
      x <- x[x[, "reprod"] == 1, , drop = FALSE]
      if (nrow(x) <= 1) {
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
            t_rest = tr
          )
        }
      }
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()

  # Inter-birth interval ----

  inter.birth <- purrr::map(.x = seq_len(n), .f = ~ {
    lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      x <- narw.individuals[[.x]][, , d]
      birth.record <- x[x[, "alive"] == 1, "birth"]
      if (length(birth.record) == 0 | sum(birth.record, na.rm = TRUE) <= 1) {
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
      rep.status <- rep(0, yrs + 1)
      rep.status[narw.individuals[[.x]][, , d][, "female"] == 1 & narw.individuals[[.x]][, , d][, "alive"] == 1 & narw.individuals[[.x]][, , d][, "reprod"] == 0] <- 1
      rep.status
    }) |>
      do.call(what = rbind) |>
      colSums()

    # Total females (by year)
    totfem <- rowSums(narw.pop[.x, , 5:8])
    totfem[1] <- sum(N_0[6:8])
    out <- nonrep / totfem
    # Correct cases where population has gone extinct - which leads to NAs here
    out[is.na(out)] <- 0
    out
  }) |> do.call(what = rbind)

  # Abortions ----
  abort.df <- as.data.frame.table(abortions[2:nrow(abortions), , drop = FALSE], responseName = "value") |>
    dplyr::rename(year = Var1, prj = Var2, n_abort = value) |>
    dplyr::mutate(
      year = gsub(pattern = "yr ", replacement = "", year),
      prj = gsub(pattern = "prj ", replacement = "", prj)
    ) |>
    data.table::as.data.table() |>
    dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year)))

  # Pregnant females ----
  preg.fem <- as.data.frame.table(narw.pop[, , cohorts.proj[id == 4, name], drop = FALSE]) |>
    dplyr::rename(year = Var2, prj = Var1, n_pregfem = Freq) |>
    dplyr::mutate(
      year = gsub(pattern = "yr ", replacement = "", year),
      prj = gsub(pattern = "prj ", replacement = "", prj)
    ) |>
    data.table::as.data.table() |>
    dplyr::filter(year > 0) |>
    dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year)))

  narw.df <- purrr::map(.x = cohorts.proj$name, .f = ~ {
    tmp <- narw.pop[, , .x] |>
      tibble::as_tibble() |>
      tibble::rownames_to_column(var = "prj")
    if (n == 1) {
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
  }) |>
    do.call(what = rbind) |>
    data.table::data.table()

  # Population totals ----

  tot.df <- tibble::as_tibble(tot.pop) |>
    tibble::rownames_to_column(var = "prj") |>
    tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |>
    dplyr::mutate(year = 2019 + as.numeric(gsub("yr ", "", year))) |>
    dplyr::mutate(cohort = "North Atlantic right whales") |>
    data.table::data.table()

  narw.conf <- narw.df[
    , list(
      mean = mean(N, na.rm = TRUE),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>
    dplyr::rename(cohortID = cohort) |>
    dplyr::mutate(cohort = factor(cohortID, levels = c(
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
  ] |>
    dplyr::rename(cohortID = cohort) |>
    dplyr::mutate(cohort = factor(cohortID, levels = c(
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

  # Minimum expected population size ----

  min_popsize <- function(array) {
    out <- apply(X = array, MARGIN = 1, FUN = cummin)
    out <- apply(X = out, MARGIN = 1, FUN = mean)
    return(out)
  }

  minpop <- min_popsize(tot.pop)
  minpop <- data.table::data.table(year = 0:yrs, minpop = minpop) |>
    dplyr::mutate(year = 2019 + year)

  # Quasi-extinction risk ----

  # This is the  probability that the number of reproductive females would fall below thresholds of 10,
  # 50, or 100 animals by any point in time. However, a quasi-extinction threshold of 100
  # is not informative because the population is already below this value at the start.
  # Similarly, a quasi-extinction threshold of 10 mature females is also not very informative
  # because the probability is uniformly low for all scenarios and time points.

  qr <- 50

  quasi_ext <- function(array, threshold) {
    out <- apply(X = array, MARGIN = 1, FUN = cummin)
    out <- apply(X = out, MARGIN = 1, FUN = function(x) sum(x <= threshold) / length(x))
    return(out)
  }

  quasi.females <- quasi_ext(rowSums(narw.pop[, , cohorts.proj[id >= 4, name], drop = FALSE], dims = 2), qr)
  quasi.females <- data.table::data.table(year = 0:yrs, quasi = quasi.females) |>
    dplyr::mutate(year = 2019 + year)

  # IUCN ----
  # Calculated over 100 years
  if (yrs - burn.in >= 100) {
    iucn.threshold <- c(0.3, 0.5, 0.8)
    iucn.p <- purrr::map_dbl(.x = iucn.threshold, .f = ~ {
      out <- as.numeric(((sum(N_0) - apply(X = tot.pop[, 1:100], MARGIN = 1, min)) / sum(N_0)) > .x)
      sum(out) / length(out)
    })
    names(iucn.p) <- paste0(100 * iucn.threshold, "%")
  } else {
    iucn.p <- list(NULL)
  }

  # Label projections
  if (length(obj) == 1) {
    out.label <- obj[[1]]$param$label
  } else if (length(obj) > 1) {
    out.label <- obj[[2]]$param$label
  }

  # Output list ----

  outproj <- list(
    param = list(
      label = out.label,
      start.year = current.yr,
      yrs = yrs - burn.in,
      burn = burn.in,
      n = n,
      cohorts = cohorts.proj,
      abort = abort.rate,
      nonrep = nonrep,
      phases = phase.names,
      schedule = schedule,
      run = run_time
    ),
    init = list(
      N_0 = N_0,
      cohort = cohort.vec,
      reprod.fem = reprod.fem
    ),
    dat = list(
      ind = narw.individuals,
      birth = list(
        tot = births.df,
        perfemale = births.per.female,
        inter = inter.birth
      ),
      death = deaths.df,
      health = prjpreds,
      survival = survival,
      abort = abort.df,
      nonrepfem = nonreprod.females,
      pregfem = preg.fem,
      rest = time.resting,
      pop = narw.pop,
      tot = tot.pop
    ),
    proj = list(
      tbl = rbind(narw.df, tot.df) |> dplyr::mutate(prj = as.numeric(prj)),
      mean = rbind(narw.conf, tot.conf),
      quasi = quasi.females,
      minpop = minpop,
      iucn.p = iucn.p
    )
  )
  class(outproj) <- c("narwproj", class(outproj))
  return(outproj)
}
