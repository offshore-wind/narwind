#' Forecast NARW abundance
#'
#' @param obj Object returned by narw()
#' @param n Number of replicate projections
#' @param yrs Number of years for the projection
#' 
#' @export
#'
#' @author Phil J. Bouchet
#'
predict.narwsim <- function(obj,
                            n = 100,
                            # scenario,
                            yrs = 35,
                            ...) {

  if(sum(suppressWarnings(purrr::map_lgl(.x = obj$gam$pred, .f = ~any(is.na(.x)))) > 0)) 
    stop("Insufficient sample size. Cannot make predictions.") 
  
  if(length(obj$param$cohortID) < 6) stop("Missing cohorts. Cannot make predictions.")
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values
  if(length(args) == 0){
    progress <- TRUE
  } else {
    if("progress" %in% names(args)) progress <- args[["progress"]] else progress <- TRUE
  }
  
  if(yrs > 2000){
    yrs <- yrs - lubridate::year(lubridate::now())
  }
  
  cat("Initializing ...\n")
  
  # if(is.null(obj$gam)) stop("Insufficient data available. Cannot proceed with population projections.")
  # if(!identical(cohortID, 1:6)) stop("Missing cohorts in input <obj>. Cannot proceed with population projections.")
  # if(length(obj$gam$fit$surv$xlevels[[1]]) < 6 | length(obj$gam$fit$bc$xlevels[[1]]) < 6) stop("Missing factor levels in input <obj>. Cannot proceed with population projections.")
  
  # plogis("link" predictions + error)
  
  # Adapted from original code by Scott Creel
  # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
  # 
  # Prediction intervals
  # https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html
  
  # Population estimate as per 2022 NARW report card is 340 (+/- 7).
  
  # test <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = seq(0,0.6, 0.01)), type = "link", se.fit = TRUE)
  # plot(seq(0,0.6, 0.01), plogis(test$fit), type = "l")
  # lines(seq(0,0.6, 0.01), plogis(Reduce("+", test)), lty = 2)
  # lines(seq(0,0.6, 0.01), plogis(Reduce("-", test)), lty = 2)
  # test2 <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = 0.2), type = "response")
  # abline(v = 0.2)
  # abline(h = tesct2)
  
  cohortID <- obj$param$cohortID
  cohorts <- obj$param$cohorts |> dplyr::slice(1) |> 
    dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |> 
    dplyr::bind_rows(obj$param$cohorts) |> 
    dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
  cohorts <- cohorts[c(1,3,5,2,4,6,7,8)]
  
  # Attributes to monitor during projection
  mat.attribs <- c("alive", "cohort", "female", "age", "length", "tot_mass", "lean_mass", "bc",
                   "p_surv", "min_bc", "rest", "time_rest", "birth", "reprod")
  
  # Current year
  current.yr <- lubridate::year(lubridate::now())
  
  # Extract terminal functions
  mod <- obj$gam$fit
  mod[["gest"]] <- gam_gest
  
  #'------------------------------------------------------
  # GAM PARAMETERS
  #'------------------------------------------------------
  
  mbc_preds <- obj$gam$pred$bc_gest
  bc_preds <- obj$gam$pred$bc
  surv_preds <- obj$gam$pred$surv
  
  #'------------------------------------------------------
  # Abortion rate
  #'------------------------------------------------------
  
  abort.rate <- sum(obj$abort[, abort])/obj$param$nsim
  
  #'------------------------------------------------------
  # INITIALIZATION
  #'------------------------------------------------------
  
  # Define initial population vector
  N0 <- c(0, 7, 212, 0, 71, 1, 7, 61)
  names(N0) <- cohorts[, name]
  
  cat("Setting up ...\n")
  
  # Initial cohort vector
  cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts), FUN = function(x) rep(cohorts$id[x], each = N0[x])))
  
  # Reproductive females - 4% never reproduce
  reprod.fem <- rbinom(n = sum(N0), size = 1, prob = ifelse(cohort.vec == 6, (1 - 0.04), 1))
  
  # Number of individuals in each cohort
  narw.pop <- array(
    data = NA, c(n, yrs + 1, nrow(cohorts)),
    dimnames = list(
      paste0("prj ", 1:n),
      paste0("yr ", 0:yrs),
      cohorts$name
    )
  )
  narw.pop[, 1, ] <- rep(N0, each = n)
  
  # Total population size
  tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj", 1:n), paste0("yr ", 0:yrs)))
  tot.pop[, 1] <- sum(N0)
  
  births <- array(
    data = NA, c(yrs + 1, n),
    dimnames = list(
      paste0("yr ", 0:yrs),
      paste0("prj ", 1:n)
    )
  )
  
  #'------------------------------------------------------
  # RUN PROJECTIONS
  #'------------------------------------------------------
  
  # This uses nested loops. 
  # The prj loop (outermost loop) replicates the projection <n> times.
  # The i loop is next, and steps across all years of projection from an initial population vector.
  
  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n, clear = FALSE, width = 80
  )
  
  cat("Running projections ...\n")
  start.time <- Sys.time()
  
  narw.individuals <- vector(mode = "list", length = n)
  
  for(prj in 1:n){
    
    if(progress) pb$tick() # Update progress bar
    
    # Create matrices and initialize them
    # rows: years <yrs>
    # columns: <attributes>
    # layers: individuals <n>
    # 4th dimension: replicate projection -> later converted to list
    narw.indiv <- array(
      data = NA, c(yrs + 1, length(mat.attribs), sum(N0)),
      dimnames = list(
        paste0("yr ", 0:yrs), mat.attribs,
        paste0("whale ", 1:(sum(N0)))
      )
    )
    
    # Alive and population cohort
    narw.indiv[1, "alive", ] <- 1
    narw.indiv[1, "cohort", ] <- cohort.vec
    
    # Sex
    #  -- Calves (male)
    narw.indiv[, "female", 1:N0[1]] <- 0
    #  -- Calves (female)
    fem <- which(cohort.vec == 0)
    narw.indiv[, "female", fem[fem > N0[1]]] <- 1
    #  -- Juveniles and adults
    narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6))] <- 1
    narw.indiv[, "female", which(cohort.vec %in% c(1, 3))] <- 0
    
    # Age
    ages <- start_age_vec(cohort.vec)
    narw.indiv[1, "age", ] <- ages
    
    # Total body length
    lengths <- age2length_vec(ages)
    narw.indiv[1, "length", ] <- lengths
    
    mass <- length2mass_vec(lengths)
    narw.indiv[1, "lean_mass", ] <- mass
    
    # Body conditon
    bc <- start_bcondition_vec(cohort.vec)
    narw.indiv[1, "bc", ] <- bc
    
    # Total mass
    narw.indiv[1, "tot_mass", ] <- mass / (1-bc)
    
    # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
    # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
    # This corresponds to a Normal (7.7, 1.45)
    # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
    # 
    # Stewart et al. 2022 -- 
    # The degree to which the energetic reserves of females are depleted during lactation may govern
    # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).

    # Non-reproductive females
    narw.indiv[1, "reprod", ] <- reprod.fem
    
    # Reserves needed to leave resting state and initiate pregnancy
    narw.indiv[1, "min_bc", ] <- mbc_preds(narw.indiv[1, "tot_mass", ]) * (cohort.vec == 6)
    
    narw.indiv[1, "time_rest", ] <- (cohort.vec == 6)
      
    narw.indiv[1, "rest", ] <- 
      ifelse(reprod.fem == 0, 1, (cohort.vec == 6) * (1 - (bc >= narw.indiv[1, "min_bc", ])))

    # Calving events
    narw.indiv[1, "birth", ] <- as.numeric(cohort.vec == 4)
    
    # Survival
    narw.indiv[1, "p_surv", ] <- (surv_preds[["0"]](bc) * (cohort.vec == 0) +
                                    surv_preds[["1"]](bc) * (cohort.vec == 1) +
                                    surv_preds[["2"]](bc) * (cohort.vec == 2) +
                                    surv_preds[["3"]](bc) * (cohort.vec == 3) +
                                    surv_preds[["4"]](bc) * (cohort.vec == 4) +
                                    surv_preds[["5"]](bc) * (cohort.vec == 5) +
                                    surv_preds[["6"]](bc) * (cohort.vec == 6))

    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      if(tot.pop[prj, i-1] > 0) {
        
        # alive <- narw.indiv[i-1, "alive", ] * (narw.indiv[i-1, "age", ] <=69)
        
        #' ----------------------------
        # ALIVE
        #' ----------------------------
        # Determine whether the animal survived
        # alive <- rbinom(n = length(ps), size = 1, prob = ps) * (narw.indiv[i-1, "age", ] <=69)
        alive <- rbinom(n = dim(narw.indiv)[3], size = 1, 
                        prob = (narw.indiv[i-1, "alive", ] * narw.indiv[i-1, "p_surv", ])) * (narw.indiv[i-1, "age", ] <=69)
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
          alive * increment_cohort(
            cohort = narw.indiv[i-1, "cohort", ],
            age = narw.indiv[i, "age", ],
            female = narw.indiv[i, "female", ],
            rest = narw.indiv[i-1, "rest", ],
            abort = abort.rate
          )
        
        #' ----------------------------
        # GROWTH
        #' ----------------------------
        
        # Increment length
        narw.indiv[i, "length", ] <- alive * age2length_vec(narw.indiv[i, "age", ])
        
        # Increment lean mass
        narw.indiv[i, "lean_mass", ] <- alive * length2mass_vec(narw.indiv[i, "length", ])
        
        # Predict new body condition from current body condition
        narw.indiv[i, "bc", ] <-
          alive * (bc_preds[["0"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 0) +
                     bc_preds[["1"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 1) +
                     bc_preds[["2"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 2) +
                     bc_preds[["3"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 3) +
                     bc_preds[["4"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 4) +
                     bc_preds[["5"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 5) +
                     bc_preds[["6"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i, "cohort", ] == 6))
        
        # Increment total mass
        narw.indiv[i, "tot_mass", ] <- alive * narw.indiv[i, "lean_mass", ] / (1 - narw.indiv[i, "bc", ])
        
        #' ----------------------------
        # REPRODUCTION
        #' ----------------------------
        
        narw.indiv[i, "reprod", ] <- alive * narw.indiv[i-1, "reprod", ]
        
        # Resting state
        narw.indiv[i, "rest", ] <-
          (narw.indiv[i, "reprod", ] == 1) * alive * (narw.indiv[i, "cohort", ] == 6) *
          (1 - (narw.indiv[i, "bc", ] >= narw.indiv[i-1, "min_bc", ]))
                 
        # Which animals are juvenile females that are ready to start reproducing
        # juvenile.females.ofage <- 
        #   (narw.indiv[i-1,"cohort", ] == 2) * (narw.indiv[i-1, "age", ] >= 9) * (narw.indiv[i-1, "bc", ] >= narw.indiv[i-1, "min_bc", ])
        
        # Time spent in resting phase
        narw.indiv[i, "time_rest", ] <-
          alive * (narw.indiv[i, "cohort", ] == 6) * 
          narw.indiv[i, "rest", ] * (narw.indiv[i-1, "time_rest", ] + 1)
        
        # Minimum body condition needed to successfully bring fetus to term without starving
        # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
        narw.indiv[i, "min_bc", ] <-
          alive * mbc_preds(narw.indiv[i, "tot_mass", ]) * (narw.indiv[i, "cohort", ] == 6) 
        # * narw.indiv[i, "rest", ]
        
        # Birth of new calf, conditional on the mother being alive and in pregnant state
        narw.indiv[i, "birth", ] <- alive * (narw.indiv[i, "cohort", ] == 4)
        new.births <- sum(narw.indiv[i, "birth", ])
        births[i, prj] <- new.births
        
        if(new.births > 0){
          
          new.calves <- array(
            data = NA, c(yrs + 1, length(mat.attribs), new.births),
            dimnames = list(
              paste0("yr ", 0:yrs), mat.attribs,
              paste0("whale ", ((dim(narw.indiv)[3]+1):(dim(narw.indiv)[3] + new.births)))
            )
          )
          
          new.calves[i,,] <- add_calf(new.births, mat.attribs)
          narw.indiv <- abind::abind(narw.indiv, new.calves, along = 3)
        }
        
        #' ----------------------------
        # SURVIVAL
        #' ----------------------------
        
        # Predict survival probability based on body condition
        ps <- narw.indiv[i, "alive", ] *
          (surv_preds[["0"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 0) +
             surv_preds[["1"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 1) +
             surv_preds[["2"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 2) +
             surv_preds[["3"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 3) +
             surv_preds[["4"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 4) +
             surv_preds[["5"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 5) +
             surv_preds[["6"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 6))
        
          # 0.95 * narw.indiv[i, "alive", ]

        narw.indiv[i, "p_surv", ] <- ps
        
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
        
      } # End totpop >0
    } # End years
    
    narw.individuals[[prj]] <- narw.indiv
    
  } # End projections
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  # cat("Processing outputs ...\n")
  
  # narw.out <- purrr::map(.x = 1:n, .f = ~{
  #   reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
  #     dplyr::mutate(attr = mat.attribs[attr]) |> 
  #     tidyr::pivot_wider(names_from = attr, values_from = value) |> 
  #     dplyr::mutate(prj = .x) |> 
  #     dplyr::relocate(prj, .before = yr)
  # }) |> do.call(what = rbind) |> 
  #   data.table::data.table()
  
  # narw.out <- narw.out[is.finite(rowSums(narw.out)),]
  
  narw.df <- purrr::map(.x = cohorts$name, .f = ~{
    narw.pop[,,.x] |> 
      tibble::as_tibble() |> 
      tibble::rownames_to_column(var = "prj") |> 
      tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
      dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
      dplyr::mutate(cohort = stringr::str_to_sentence(.x))
  }) |> do.call(what = rbind) |> data.table::data.table()
  
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
  
  # Find 95% confidence intervals on final population size
  cat("Final population size:\n")
  final.pop <- unname(tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)])
  cat("N = ", round(final.pop[1],0), " (95% CI: ", round(final.pop[2],0), "–", round(final.pop[3],0), ")\n", sep = "")
  
  cat(paste0("Time elapsed: ", run_time))
  cat("\n")
  
  outproj <- list(param = list(yrs = yrs,
                               n = n),
                  dat = list(ind = narw.individuals,
                             birth = births,
                             pop = narw.pop,
                             tot = tot.pop),
                  prj = list(final = c(round(final.pop[1],0), round(final.pop[2],0), round(final.pop[3],0)),
                             proj = rbind(narw.df, tot.df) |> dplyr::mutate(prj = as.numeric(prj)), 
                             mean = rbind(narw.conf, tot.conf)))
  
  class(outproj) <- c("narwproj", class(outproj))
  return(outproj)
  
}
