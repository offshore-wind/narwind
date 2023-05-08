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
                            yrs = 35,
                            n = 100,
                            popr = 1,
                            do.plot = FALSE,
                            seed = 125897,
                            ...) {
  
  set.seed(seed)
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values
  if(length(args) == 0){
    spline <- TRUE
    progress <- TRUE
  } else {
    if("spline" %in% names(args)) spline <- args[["spline"]] else spline <- TRUE
    if("progress" %in% names(args)) progress <- args[["progress"]] else progress <- TRUE
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
  mat.attribs <- c("alive", "cohort", "female", "age", "length", "length_a", "length_b", "length_c",
                   "tot_mass", "lean_mass", "bc", "mass_a", "mass_b", "p_surv", "min_bc", "trest", "t2calf", "birth")
  
  # Current year
  current.yr <- lubridate::year(lubridate::now())
  
  # Extract terminal functions
  mod <- obj$gam$fit
  mod[["gest"]] <- gam_gest
  
  #'------------------------------------------------------
  # GAM PARAMETERS
  #'------------------------------------------------------
  
  mbc_preds <- obj$gam$pred$min_bc
  bc_preds <- obj$gam$pred$bc
  surv_preds <- obj$gam$pred$surv

  # plot(seq(0,1,by = 0.01), bc_preds[["5"]](seq(0,1,by = 0.01)), ylim = c(0,1))
  # plot(seq(0,1,by = 0.01), surv_preds[["0"]](seq(0,1,by = 0.01)), ylim = c(0,1))
  # plot(seq(0,1,by = 0.01), surv_preds[["4"]](seq(0,1,by = 0.01)), ylim = c(0,1))
  # plot(seq(0,1,by = 0.01), surv_preds[["2"]](seq(0,1,by = 0.01)), ylim = c(0,1))
  
  #'------------------------------------------------------
  # INITIALIZATION
  #'------------------------------------------------------
  
  # Define initial population vector
  N0 <- c(2, 5, 212, 2, 69, 1, 7, 61)
  names(N0) <- cohorts[, name]
  totn <- sum(N0)* (1 + popr)
  
  # Create matrices and initialize them
  # rows: years <yrs>
  # columns: <attributes>
  # layers: individuals <n>
  # 4th dimension: replicate projection -> later converted to list
  narw.indiv <- array(data = NA, c(yrs + 1, length(mat.attribs), totn, n), 
                      dimnames = list(paste0("yr ", 0:yrs), 
                                      mat.attribs,
                                      paste0("whale ", 1:totn),
                                      paste0("prj ", 1:n)))
  
  cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts), FUN = function(x) rep(cohorts$id[x], each = N0[x])))
  
  animals <- 1:sum(N0)
  
  # Alive and population cohort
  narw.indiv[1, "alive", animals, ] <- 1
  narw.indiv[1, "cohort", animals, ] <- rep(cohort.vec, n)
  
  # Sex
  #  -- Calves (male)
  narw.indiv[, "female", 1:N0[1], ] <- 0
  #  -- Calves (female)
  fem <- which(cohort.vec == 0)
  narw.indiv[, "female", fem[fem > N0[1]], ] <- 1
  #  -- Juveniles and adults
  narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6)), ] <- 1
  narw.indiv[, "female", which(cohort.vec %in% c(1, 3)), ] <- 0
  
  # Age
  ages <- start_age_vec(rep(cohort.vec, n))
  narw.indiv[1, "age", animals, ] <- ages
  
  # Total body length
  l.params <- agL_vec(ages)
  lengths <- age2length_vec(ages, l.params)
  narw.indiv[1, "length", animals, ] <- lengths
  
  narw.indiv[, "length_a", animals, ] <- rep(l.params[, 1], each = yrs + 1)
  narw.indiv[, "length_b", animals, ] <- rep(l.params[, 2], each = yrs + 1)
  narw.indiv[, "length_c", animals, ] <- rep(l.params[, 3], each = yrs + 1)
  
  # Total mass
  m.params <- mL(n * sum(N0))
  mass <- length2mass_vec(lengths, m.params, FALSE)
  narw.indiv[, "mass_a", animals, ] <- rep(m.params[, 1], each = yrs + 1)
  narw.indiv[, "mass_b", animals, ] <- rep(m.params[, 2], each = yrs + 1)
  narw.indiv[1, "tot_mass", animals, ] <- mass
  
  # Body conditon
  narw.indiv[1, "bc", animals, ] <- start_bcondition_vec(rep(cohort.vec, n))
  
  # lean mass
  narw.indiv[1, "lean_mass", animals, ] <- mass - (narw.indiv[1, "bc", animals, ] * mass)
  
  # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
  # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
  # This corresponds to a Normal (7.7, 1.45)
  # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
  # 
  # Stewart et al. 2022 -- 
  # The degree to which the energetic reserves of females are depleted during lactation may govern
  # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).
  t2calf <- (rep(cohort.vec, n) == 6) * random_int(sum(N0) * n)
  narw.indiv[1, "t2calf", animals, ] <- t2calf
  narw.indiv[1, "trest", animals, ] <- as.numeric(ifelse(t2calf == 0, 13, 1)) * (as.numeric(narw.indiv[1, "cohort", animals, ]) == 6)
  
  if (!spline) {
    narw.indiv[1, "min_bc", animals, ] <- predict_m(
      model = mod,
      values = as.vector(narw.indiv[1, "tot_mass", animals, ]),
      prediction = "gest"
    ) * as.vector(narw.indiv[1, "cohort", animals, ] == 6)
  } else {
    narw.indiv[1, "min_bc", animals, ] <- mbc_preds(narw.indiv[1, "tot_mass", animals, ]) * (narw.indiv[1, "cohort", animals, ] == 6)
  }

  narw.indiv[1, "birth", animals, ] <- ifelse(narw.indiv[1, "trest", animals, ] == 13 & narw.indiv[1, "t2calf", animals, ] == 0, 1, 0)
  narw.indiv[1, "p_surv", animals, ] <- 1
  
  #' ---------------------------
  # IMPORTANT 
  #' ---------------------------
  # Turn array into a list
  narw.indiv <- purrr::array_branch(narw.indiv, 4)

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
  
  for(prj in 1:n){
    
    if(progress) pb$tick() # Update progress bar
    
    animals <- 1:sum(N0)
    
    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      current.dat <- as.matrix(narw.indiv[[prj]][i-1, , animals])
      alive <- current.dat["alive", animals] * (current.dat["age", animals] <=69)
      
      #' ----------------------------
      # SURVIVAL
      #' ----------------------------
      # Predict survival probability based on body condition
      
      if(!spline){

        ps <- alive * predict_m(model = mod, cohort = current.dat["cohort",animals], 
                                values = current.dat["bc",animals], prediction = "surv")

      } else {

        ps <- alive * (surv_preds[["0"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 0) +
                         surv_preds[["1"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 1) +
                         surv_preds[["2"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 2) +
                         surv_preds[["3"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 3) +
                         surv_preds[["4"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 4) +
                         surv_preds[["5"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 5) +
                         surv_preds[["6"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 6))
      }

      # ps <- 1
      narw.indiv[[prj]][i, "p_surv", animals] <- ps
      
      # Determine whether the animal survived
      alive <- rbinom(n = animals, size = 1, prob = ps) * (current.dat["age", animals] <=69)
      narw.indiv[[prj]][i, "alive", animals] <- alive
      
      # Sex remains the same
      narw.indiv[[prj]][i, "female", animals] <- current.dat["female", animals]
      
      #' ----------------------------
      # GROWTH
      #' ----------------------------
      
      # Increment age
      narw.indiv[[prj]][i, "age", animals] <- alive * (current.dat["age", animals] + 1)
      
      # Increment length
      newLp <- agL_vec(animals)
      
      narw.indiv[[prj]][i,"length_a", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,1], current.dat["length_a", animals])
      
      narw.indiv[[prj]][i,"length_b", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,2], current.dat["length_b", animals])
      
      narw.indiv[[prj]][i,"length_c", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,3], current.dat["length_c", animals])

      narw.indiv[[prj]][i, "length", animals] <-
        alive * age2length_vec(
          narw.indiv[[prj]][i, "age", animals],
          t(narw.indiv[[prj]][i, c("length_a", "length_b", "length_c"), animals])
        )
      
      # Increment lean mass
     narw.indiv[[prj]][i, "lean_mass", animals] <- alive * length2mass_vec(narw.indiv[[prj]][i, "length", animals],
      t(narw.indiv[[prj]][i, c("mass_a", "mass_b"), animals]), lean = TRUE)
     
      # Predict new body condition from current body condition
      if (!spline) {
        narw.indiv[[prj]][i, "bc", animals] <- alive * predict_m(
          model = mod, cohort = current.dat["cohort", animals],
          values = current.dat["bc", animals], prediction = "bc")
      } else {
        narw.indiv[[prj]][i, "bc", animals] <-
          alive * (bc_preds[["0"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 0) +
                     bc_preds[["1"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 1) +
                     bc_preds[["2"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 2) +
                     bc_preds[["3"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 3) +
                     bc_preds[["4"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 4) +
                     bc_preds[["5"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 5) +
                     bc_preds[["6"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 6))
      }
      
      # Increment total mass
      narw.indiv[[prj]][i, "tot_mass", animals] <-
        alive * narw.indiv[[prj]][i, "lean_mass", animals] / (1 - narw.indiv[[prj]][i, "bc", animals])
      
      #' ----------------------------
      # REPRODUCTION
      #' ----------------------------
      
      # Which animals are resting females?
      rest.females <- (current.dat["cohort", animals] == 6)
      
      # Which animals are juvenile females that are ready to start reproducing
      juvenile.females.ofage <- (current.dat["cohort", animals] == 2) * (current.dat["age", animals] >= 8)
      
      # Which animals calved in previous step?
      prev.births <- current.dat["birth", animals]
      
      newt2calf <- ifelse(prev.births == 1, random_int(sum(prev.births), lwr = 1), 
                          ifelse(current.dat["t2calf", animals] == 0, 0, 
                                 current.dat["t2calf", animals] - 1))
      
      # Time spent in resting state - only incremented if calving event hasn't occurred, otherwise reset
      narw.indiv[[prj]][i, "t2calf", animals] <- alive * rest.females * newt2calf
      
      # Years until next calving event
      narw.indiv[[prj]][i, "trest", animals] <- 
        alive * (rest.females | juvenile.females.ofage) * 
        ifelse(narw.indiv[[prj]][i - 1, "trest", animals] == 13, 1, current.dat["trest", animals] + 1)
      
      # Minimum body condition needed to successfully bring fetus to term without starving
      # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
      
      if (!spline) {
        narw.indiv[[prj]][i, "min_bc", animals] <-
          alive * predict_m(model = mod, values = narw.indiv[[prj]][i, "tot_mass", animals], prediction = "gest") * rest.females
      } else {
        narw.indiv[[prj]][i, "min_bc", animals] <-
          alive * mbc_preds(narw.indiv[[prj]][i, "tot_mass", animals]) * rest.females
      }
      
      # Birth of new calf, conditional on the mother being alive, in pregnant state
      narw.indiv[[prj]][i, "birth", animals] <- alive * (current.dat["cohort", animals] == 4)
      
      # Maturity - transitions between cohorts
      narw.indiv[[prj]][i, "cohort", animals] <-
        alive * increment_cohort(
          cohort = narw.indiv[[prj]][i - 1, "cohort", animals],
          age = narw.indiv[[prj]][i, "age", animals],
          female = narw.indiv[[prj]][i, "female", animals],
          bc = narw.indiv[[prj]][i, "bc", animals],
          min_bc = narw.indiv[[prj]][i, "min_bc", animals],
          trest = narw.indiv[[prj]][i, "trest", animals],
          t2calf = narw.indiv[[prj]][i, "t2calf", animals])
      
      new.births <- sum(narw.indiv[[prj]][i, "birth", animals])
      
      if(new.births > 0){
        narw.indiv[[prj]][i, , (max(animals)+1):(max(animals)+new.births)] <- add_calf(n = new.births, attr = mat.attribs)
        animals <- 1:(length(animals) + new.births)
      }
      
      # Number of animals in each cohort
      # Calves (male)
      narw.pop[prj, i, 1] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
                                   (narw.indiv[[prj]][i, "female", animals] == 0) *
                                   (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Juveniles and adults (male)
      narw.pop[prj, i, 2] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 1) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 3] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 3) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Calves (female)
      narw.pop[prj, i, 4] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
                                   (narw.indiv[[prj]][i, "female", animals] == 1) *
                                   (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Juvenile and reproductive adults (female)
      narw.pop[prj, i, 5] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 2) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 6] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 4) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 7] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 5) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 8] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 6) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Total population size
      tot.pop[prj, i] <- sum(narw.indiv[[prj]][i, "alive", animals], na.rm = TRUE)
      
    } # End years
  } # End projections
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  cat("Processing outputs ...\n")
  
  narw.out <- purrr::map(.x = 1:n, .f = ~{
    reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
      dplyr::mutate(attr = mat.attribs[attr]) |> 
      tidyr::pivot_wider(names_from = attr, values_from = value) |> 
      dplyr::mutate(prj = .x) |> 
      dplyr::relocate(prj, .before = yr)
  }) |> do.call(what = rbind) |> 
    data.table::data.table()
  
  narw.out <- narw.out[is.finite(rowSums(narw.out)),]
  
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
  
  # births.df <- purrr::map(.x = 1:n, .f = ~{
  #   m <- matrix(rowSums(narw.indiv[[.x]][2:(yrs+1),"birth",], na.rm = TRUE), ncol = 1)
  #   colnames(m) <- .x
  #   m
  # }) |> do.call(what = cbind) |> 
  #   tibble::as_tibble() |> 
  #   tibble::rownames_to_column(var = "year") |> 
  #   dplyr::mutate(year = as.numeric(year)) |> 
  #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |> 
  #   dplyr::select(prj, year, birth) |> 
  #   dplyr::arrange(prj, year) |> 
  #   data.table::data.table()
  # 
  # deaths.df <- purrr::map(.x = 1:n, .f = ~{
  #   m <- matrix(apply(X = narw.indiv[[.x]][2:(yrs+1),"alive",],
  #                     MARGIN = 1,
  #                     FUN = function(x) {
  #     r <- x[!is.na(x)]
  #     r <- sum(r == 0)
  #     r
  #     }), ncol = 1)
  #   colnames(m) <- .x
  #   m
  # }) |> do.call(what = cbind) |>
  #   tibble::as_tibble() |>
  #   tibble::rownames_to_column(var = "year") |>
  #   dplyr::mutate(year = as.numeric(year)) |>
  #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
  #   dplyr::select(prj, year, death) |>
  #   dplyr::arrange(prj, year) |>
  #   data.table::data.table()
  
  narw.conf <- narw.df[
  , list(
    mean = mean(N),
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
      mean = mean(N),
      lwr = quantile(N, 0.025),
      uppr = quantile(N, 0.975)
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
      "Adults (female, resting)",
      "North Atlantic right whales"
    )))
  
  p1 <- plot_projection(narw.df, narw.conf)
  p2 <- plot_projection(tot.df, tot.conf)
  
  if(do.plot){
    print(p1)
    print(p2)
  }
  
  # Find 95% confidence intervals on final population size
  cat("Final population size:\n")
  final.pop <- unname(tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.5, 0.025, 0.975))])
  cat("N = ", round(final.pop[1],0), " (95% CI: ", round(final.pop[2],0), "–", round(final.pop[3],0), ")\n", sep = "")
  
  cat(paste0("Time elapsed: ", run_time))
  cat("\n")

  return(list(dat = narw.out,
              out = list(df = rbind(narw.df, tot.df), 
                         ci = rbind(narw.conf, tot.conf),
                         plot = list(p1, p2))))
  
}
