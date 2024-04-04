# Define methods
#' @export
augment <- function (x, ...) UseMethod("augment", x)
#' @export
export <- function(x, ...) UseMethod("export")
#' @export
animate <- function (x, ...) UseMethod("animate", x)

# PACKAGE DEVELOPMENT ------------------------------------------------------

# ++ [FUNCTION] Load all non-essential internal functions [used only for package development]
load_functions <- function(){
  functions.files <- list.files("inst/functions", full.names = TRUE)
  for(j in functions.files) source(j)
  cat("Done!\n")
}

# ++ [FUNCTION] Source C++ code and required functions
source_pkg <- function(){
  Rcpp::compileAttributes()
  Rcpp::sourceCpp("src/simtools.cpp")
  suppressWarnings(source("R/narw.R"))
}

# MODEL ------------------------------------------------------

# ++ [FUNCTION] Generate scenario objects for use in the simulator
# ++ [PARAM] scenario –– Integer. Unique scenario ID, from 1 to 3. A value of 0 indicates baseline conditions.
# ++ [RETURN] <narwscenario> object
get_scenarios <- function(scenario = 1){
  
  # BOEM offshore wind scenarios
  if(scenario > 0){
  
  # vparams <- vessel_transits[, list(speed_knt = unique(speed_knt), Nvessels = sum(Nvessels)), list(category, windfarm, routeID, phase)]
  vessel_transits <- targets::tar_read(vessel_transits)
  turbines <- targets::tar_read(turbines)
  vessel_routes <- targets::tar_read(vessel_routes)
  
  vparams <- vessel_transits[phase == ifelse(scenario < 3, "Construction", "Operation and maintenance")]
  
  if(scenario == 1){
    start.month <- c(1,2,1)
    start.day = c(15,1,1)
  } else if(scenario == 2){
    start.month <- c(7,9,5)
    start.day = c(15,15,1)
  } else {
    start.month <- c(1,1,1)
    start.day = c(1,1,1)
  }
  
  phase <- ifelse(scenario <=2, 1, 2)
  
  out <-
    list(phase = phase,
      locs = data.table::as.data.table(turbines[[paste0("scenario_0", scenario)]]),
      routes = vessel_routes,
      vessels = vparams,
      start.month = start.month,
      start.day = start.day,
      piles.per.day = 1,
      ambient = 80,
      sourceLvL = 220,
      lowerdB = c(0,10,0)[scenario],
      logrange = 15,
      absorb = 1.175
    )
  
  } else {
    
    # Baseline conditions
    
    out <-
      list(phase = 0,
           locs = NULL,
           routes = NULL,
           vessels = NULL,
           start.month = NULL,
           start.day = NULL,
           piles.per.day = NULL,
           ambient = NULL,
           sourceLvL = NULL,
           lowerdB = NULL,
           logrange = NULL,
           absorb = NULL
      )
    
  }
  
  class(out) <- c("narwscenario", class(out))
  return(out)
}

# ++ [FUNCTION] Re-scale density surfaces within pre-defined bounds.
# ++ [PARAM] p –– Density surface, converted to matrix format. 
# ++ [PARAM] coords –– Two-column matrix of raster cell coordinates.
# ++ [PARAM] lower –– Numeric. Northing of the southernmost boundary.
# ++ [PARAM] upper –– Numeric. Northing of the northermost boundary.
# ++ [PARAM] left –– Numeric. Easting of the westernmost boundary.
# ++ [PARAM] right –– Numeric. Easting of the westernmost boundary.
# ++ [NOTE] Used in allocation of initial coordinates for simulated animals 
# ++ [RETURN] Modified matrix.
warp <- function(p, coords, lower = 700, upper = 850, left = 500, right = 750){
  
  inside <- which(coords[,2] > lower & coords[,2] < upper & coords[,1] > left & coords[,1] < right)
  outside <- which(!1:nrow(coords) %in% inside)
  N.outside <- sum(p[outside])
  p[outside] <- 0
  p[inside][p[inside] > 0] <- p[inside][p[inside] > 0] + (N.outside / length(p[inside][p[inside] > 0]))
  return(p)
}

# ++ [FUNCTION] Generate initial coordinates for simulated animals
# ++ [PARAM] maplist –– List of density surfaces, weighted and unweighted.
# ++ [PARAM] coords –– Two-column matrix of raster cell coordinates.
# ++ [PARAM] migrate –– Two-column matrix indicating whether individuals
# migrate to the SEUS/GSL as returned by prob_migration.
# ++ [PARAM] init.month –– Integer. Start month for the simulation.
# ++ [PARAM] nsim –– Integer. Number of simulated animals.
# ++ [RETURN] Initial coordinates for each animal
initiate_xy <- function(maplist,
                         coords,
                         migrate,
                         init.month,
                         nsim){
  
  out <- do.call(rbind, lapply(seq_along(maplist[[1]]), function(mo) { # For each month of the simulation year
    
    ncells <- seq_along(maplist[[1]][[mo]])
    
    sampled.xy <- purrr::map(.x = maplist, .f = ~{
      sample(x = ncells, size = nsim, replace = TRUE, prob = .x[[mo]])
    })
    
    purrr::map_int(.x = 1:nsim, .f = ~ { # For each individual
      
      # Start of simulation off Southern New England
      if (mo == init.month) {
        sampled.xy[["sne"]][.x]
        
        # Migration to the calving grounds between November and February
      } else if (mo %in% c(1:2, 11:12) & migrate[.x, 1] == 1) {
        sampled.xy[["seus"]][.x]
        
        # Migration to the Gulf of St Lawrence between June and September
      } else if (mo %in% c(6:9) & migrate[.x, 2] == 1) {
        sampled.xy[["gsl"]][.x]
        
        # Last trimester of simulation for lactating females
        # and spring season for all agents
      } else if (mo %in% c(3:5, 13:15)) {
        sampled.xy[["gom"]][.x]
        
      } else {
        sampled.xy[["gom"]][.x]
      }
      
    }) # End purrr
  }))
  
  row.names(out) <- names(maplist[[1]])
  return(out)
}

# ++ [FUNCTION] Generate initial coordinates for simulated animals
# ++ [PARAM] maplist –– List of density surfaces, weighted and unweighted.
# ++ [PARAM] coords –– Two-column matrix of raster cell coordinates.
# ++ [PARAM] migrate –– Two-column matrix indicating whether individuals
# migrate to the SEUS/GSL as returned by prob_migration.
# ++ [PARAM] init.month –– Integer. Start month for the simulation.
# ++ [PARAM] nsim –– Integer. Number of simulated animals.
# ++ [RETURN] Initial coordinates for each animal
# initiate_xy <- function(maplist,
#                         coords,
#                         migrate,
#                         init.month,
#                         nsim){
#   
#   out <- do.call(rbind, lapply(seq_along(maplist[[1]]), function(mo) {
#     
#     sampled.xy <- purrr::map(.x = maplist, .f = ~{
#       sample(x = 1:283077, size = nsim, replace = TRUE, prob = .x[[mo]])
#     })
#     
#     purrr::map_int(.x = 1:nsim, .f = ~ {
#       if (mo %in% c(1:2, 11:12) & migrate[.x, 1] == 1) {
#         sampled.xy[[1]][.x]
#       } else if (mo %in% c(6:9) & migrate[.x, 2] == 1) {
#         sampled.xy[[2]][.x]
#       } else if(mo == init.month){
#         sampled.xy[[3]][.x]
#       } else {
#         sampled.xy[[4]][.x]
#       }
#     })
#     
#   }))
#   row.names(out) <- names(maplist[[1]])
#   return(out)
# }

# ++ [FUNCTION] Extract and format outputs returned by the simulator
# ++ [PARAM] input –– List object as returned by <NARW_simulator>.
# ++ [PARAM] cohortID –– Integer. Cohort ID number.
# ++ [PARAM] dates –– Character vector. Simulation dates, as returned by <get_dates>.
# ++ [RETURN] data.table in which each row represents data for each day of the simulation.
transpose_array <- function(input, cohortID, dates) {
  
  out.array <- purrr::map2(.x = input, .y = cohortID, .f = ~ {
    
    tp <- lapply(X = .x, FUN = function(l) {
      tmp <- aperm(l, c(3, 1, 2))
      if(dim(tmp)[3] > 1){
      dimnames(tmp)[[1]] <- format(dates)
      dimnames(tmp)[[3]] <- paste0("whale.", seq_len(dim(tmp)[3]))
      }
      tmp
    })
    
    outl <- list()
    outl[["locs"]] <- tp[["locs"]]
    outl[["stress"]] <- tp[["stressors"]]
    outl[["activ"]] <- tp[["activity"]]
    
    if(.y == 5){ # Lactating females / calves
      
      outl[["E"]] <- list(adjuv = tp[["E"]], calves = tp[["E_calves"]])
      
      tp[grepl("fetus", names(tp))] <- NULL  
      outl[["inits"]] <- list(adjuv = tp[["inits"]], calves = tp[["inits_calves"]])
      outl[["attrib"]] <- list(adjuv = tp[["attrib"]], calves = tp[["attrib_calves"]])
      outl[["kj"]] <- list(adjuv = do.call(abind::abind, list(tp[["in.kj"]], tp[["out.kj"]], along = 2)),
                           calves = do.call(abind::abind, list(tp[["in.kj_calves"]], tp[["out.kj_calves"]], along = 2)))
      
    } else if(.y == 4){ # Pregnant females
      
      outl[["E"]] <- tp[["E"]]
      outl[["inits"]] <- tp[["inits"]]
      tp[grepl("calves", names(tp))] <- NULL  
      outl[["attrib"]] <- list(adjuv = tp[["attrib"]], fetus = tp[["attrib_fetus"]])
      outl[["kj"]] <- do.call(abind::abind, list(tp[grepl("kj", names(tp))], along = 2))
      

    } else {
      
      outl[["E"]] <- tp[["E"]]
      outl[["inits"]] <- tp[["inits"]]
      tp[grepl("fetus", names(tp))] <- NULL  
      tp[grepl("calves", names(tp))] <- NULL  
      outl[["attrib"]] <- tp[["attrib"]]
      outl[["kj"]] <- do.call(abind::abind, list(tp[grepl("kj", names(tp))], along = 2))
    }
    
    outl
    
  })
  
  return(out.array)
}

# ++ [FUNCTION] Consolidate simulator outputs into a single data.table
# ++ [PARAM] dtl –– List of data.table containing simulator outputs.
# ++ [PARAM] nsim –– Integer. Number of simulated animals.
# ++ [PARAM] cnames –– Character vector. Cohort names.
# ++ [PARAM] dates –– Character vector. Simulation dates, as returned by <get_dates>.
# ++ [PARAM] months -- Integer vector, indicating the month of the year for each day of the simulation.
# ++ [RETURN] list of data.table objects.
consolidate <- function(dtl, nsim, cnames, dates, months){
  purrr::map2(.x = dtl,
              .y = cnames, 
              .f = ~{
                dt <- data.table::data.table(day = rep(0:(length(dates)-1), times = nsim), 
                                             date = rep(dates, times = nsim),
                                             month = rep(months, times = nsim),
                                             whale = rep(1:nsim, each = length(dates)))
                a <- cbind(array2dt(.x), dt)
                a$region <- sort(regions$region)[a$region]
                a$cohort_name <- .y
                data.table::setcolorder(a, c((ncol(a)-4):(ncol(a)-1), 1:(ncol(a)-5))) 
              })
}

# ++ [FUNCTION] Reshape a 3D array into a 2D table
# ++ [NOTE] Used in <export.narwproj>
# ++ [RETURN] 2D matrix.
reshape_array <- function(a, value, ..., .id = NULL) 
{
  qs <- rlang::quos(...)
  if (missing(value)) {
    evalue <- rlang::sym("var")
  }
  else {
    evalue <- rlang::enquo(value)
  }
  len <- length(qs)
  d <- dim(a)
  if (len > 0) {
    dimnames <- purrr::map(qs, rlang::quo_name) |> purrr::as_vector()
  }
  else {
    dimnames <- paste0("dim_", 1:length(d))
  }
  l <- list()
  for (i in 1:length(d)) {
    l[[i]] <- 1:d[i]
  }
  names(l) <- dimnames
  tidy <- expand.grid(l)
  tidy[[rlang::quo_name(evalue)]] <- a[as.matrix(tidy)]
  if (!is.null(.id)) 
    tidy[[.id]] <- rlang::expr_name(a)
  return(tidy)
}

# ++ [FUNCTION] Calculate daily probability from yearly probability
# ++ [PARAM] yearly_p –– Numeric. Annual mortality rate 
# Estimate daily p(dead) to match a target annual mortality rate
# ++ [NOTE] If the daily probability of a mortality event occurring is p, then the daily probability of
# survival is (1-p). The probability of survival over a year is then given by (1-p)^365.
# Therefore the probability of a mortality event occurring at least once during the year
# (annual mortality rate) is p_year = 1-((1-p)^365). It follows that p_day = 1-((1-p_year)^(1/365))
# ++ [RETURN] Numeric value
yearly2daily <- function(yearly_p = 0.02){
  daily_p = 1-((1-yearly_p)^(1/365))
  return(daily_p)
}

# ++ [FUNCTION] Custom starvation mortality curve
# ++ [NOTE] Starvation probability is taken to increase exponentially with decreasing body condition
# # ++ [PARAM] starvation_death –– Numeric. Body condition at which certain death occurs.
# ++ [PARAM] starvation_onset –– Numeric. Body condition at which the starvation process begins.
# ++ [PARAM] exp.factor –– Integer. Scaling of exponential decay.
# ++ [PARAM] do.plot –– Logical. if TRUE, the resulting curves are also plotted.
# ++ [RETURN] data.frame including annual and daily probabilities
mortality_curve <- function(starvation_death = 0.005,
                            starvation_onset = 0.05,
                            exp.factor = 5,
                            do.plot = FALSE) {
  
  x <- seq(starvation_death, starvation_onset, 0.001)
  df <- data.frame(bc = x, p_annual = exp(-exp.factor * seq(0, 1, length.out = length(x))))
  mod <- lm(log(df$p_annual)~df$bc)
  names(mod$coefficients)[2] <- "bc"
  df$p <- exp(mod$coefficients[1] + mod$coefficients[2] * x)
  
  # Convert from annual to daily probability
  # df$p_daily <- sapply(X = df$p_annual, FUN = yearly2daily)

  if (do.plot) {
   
    # tag.xy <- c(0.93, .93)
    
    # Daily
    # p1 <- ggplot2::ggplot(data = rbind(data.frame(bc = 0, p_annual = 1, p_daily = 1),df), 
    #                       aes(x = bc, y = p_daily)) +
    #   ggplot2::geom_line() +
    #   theme_narw() +
    #   labs(
    #     y = "Daily mortality rate"
    #   ) +
    #   ggplot2::coord_cartesian(ylim = c(0, 0.01)) +
    #   ggplot2::scale_x_continuous(breaks = pretty(c(0, starvation_onset), n = 5))
    # p1 <- add_label(p1, tag.xy[1], tag.xy[2])
    
    # Annual
    p2 <- ggplot2::ggplot(data = rbind(data.frame(bc = 0, p_annual = 1, p = 1), df), 
                          aes(x = bc, y = p_annual)) +
      ggplot2::geom_line(aes(x = bc, y = p), col = "black") +
      theme_narw() +
      labs(
        y = "Annual mortality rate",
        x = "Body condition (relative reserve mass)"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::scale_x_continuous(breaks = pretty(c(0, starvation_onset), n = 5))
    # p2 <- add_label(p2, tag.xy[1], tag.xy[2])

    # p <- (p1 + p2) & xlab(NULL) & ggplot2::theme(plot.margin = margin(5.5, 5.5, 10, 5.5))
    # p <- patchwork::wrap_elements(panel = p) +
    #   ggplot2::labs(tag = "Body condition (relative reserve mass, %)", ) +
    #   theme(
    #     plot.tag = ggplot2::element_text(size = rel(1)),
    #     plot.tag.position = "bottom"
    #   )
    print(p2)
  }

  if(do.plot) return(p2) else return(mod$coefficients)
}

# ++ [FUNCTION] Determine maximum theoretical body condition based on feeding effort curves
# ++ [PARAM] cohort –– Integer. Unique cohort ID.
# ++ [PARAM] threshold –– Numeric. Feeding effort threshold above which feeding is negligible.
# ++ [RETURN] Body condition.
find_maxBC <- function(cohort = 4, threshold = 0.05) {
  if (cohort == 4) {
    myfunc <- function(par) {
      (threshold - feeding_effort(30, 0.6, par))^2
    }
  } else {
    myfunc <- function(par) {
      (threshold - feeding_effort_vec(10, 0.3674076, par))^2
    }
  }
  optimize(myfunc, c(0, 1), tol = 0.0001)$minimum
}

# ++ [FUNCTION] Adjust population vector with a more stable age structure
# ++ [PARAM] obj –– Population projection object as returned by <predict.narwsim>
adjust_popvec <- function(obj){
  
  if(!inherits(obj, "narwproj")) stop("Object must be of class <narwproj>")
  
  pop <- obj$proj$tbl[!cohort == "North Atlantic right whales", ]
  totpop <- obj$proj$tbl[cohort == "North Atlantic right whales", ] |> 
    dplyr::rename(tot = N)
  pop <- dplyr::left_join(pop, totpop[, c("prj", "year", "tot")], by = c("prj", "year"))
  pop$perc <- pop$N/pop$tot
  
  ggplot2::ggplot(data = pop[cohort == "Adults (female, pregnant)"], aes(year, perc, group = prj)) +
    ggplot2::geom_line() +
    theme_narw() +
    ylab("Percentage of pregnant females in the population") +
    xlab("")
  
  propsum <- pop[year == 2059, list(mean_perc = round(mean(perc), 3)), cohort] |>
    dplyr::mutate(N_adj = round(sum(obj$init$N_0) * mean_perc)) |>
    dplyr::left_join(y = data.frame(N_0 = obj$init$N_0) |> 
                       tibble::rownames_to_column() |> 
                       dplyr::rename(cohort = rowname), by = "cohort")
  
  cat("N_0:", sum(propsum$N_0), "\n")
  cat("N_adj:", sum(propsum$N_adj), "\n\n")
  print(propsum$N_adj)
  return(propsum)
  
}


# ++ [FUNCTION] Helper function to determine the mean +/- SD of initial
#  body condition distributions for different NARW cohorts, based on Fredrik
#  Christiansen's published results.
# ++ [PARAM] pop -- Logical. If FALSE (the default), returns values for the agent-based model. If TRUE, returns
# values for the stochastic population model.
# ++ [PARAM] calves –– Mean +/- SE body condition index of calves at birth (in %) from Christiansen et al.
# ++ [PARAM] juveniles –– Mean +/- SE body condition index of juveniles (in %) from Christiansen et al.
# ++ [PARAM] adults –– Mean +/- SE body condition index of adults (in %) from Christiansen et al.
# ++ [PARAM] lactating –– Mean +/- SE body condition index of lactating females (in %) from Christiansen et al.
# ++ [RETURN] List of mean +/- SD values
init_bc <- function(pop = FALSE,
                    obj = NULL,
                    calves = c(-0.61, 14.64272),
                    juveniles = c(-13.1, 11.6),
                    adults = c(-16.7, 12.49),
                    lactating = c(-9.4, 11.75755)){
  
  if(pop){
    
    if(is.null(obj)) stop("Missing projection object")
    if(!inherits(obj, "narwproj")) stop("Input object must be of class <narwproj>")
    
    max.year <- max(obj$dat$health$year)
    tmp <- obj$dat$health[year == max.year - 1,]

    out <- tmp[, list(mean = mean(bc, na.rm = TRUE), sd = sd(bc, na.rm = TRUE)), list(cohort)] |> 
      dplyr::left_join(x = obj$param$cohorts, by = c("id" = "cohort")) |> 
      dplyr::arrange(id)
    
  } else {
    
    set.seed(1358961)
    
    # From Christiansen et al. (2022) J Physiol
    # The birth condition value is -0.0061 (or -0.61%) with a SE of 0.0049 (or 0.49%)
    # The sample size is: n = 893
    # This gives an SD = sqrt(n)*0.49 = 14.64272%
    # 
    # From Christiansen et al. (2020) MEPS
    # immature NARWs (mean = −13.1%, SE = 2.9, n = 16, SD = 11.6)
    # adults (mean = −16.7%, SE = 2.0, n = 39, SD = 12.49)
    # lactating NARW females (mean = −9.4%, SE = 4.8, n = 6, SD = 11.75755)
    
    # Calculate relative blubber mass from body condition indices (mean +/- 2SE)
    # using a spline function <bc_index> fitted to data from Christiansen et al. (2022) -- Figure 7D
    get_bc <- function(v){
      out <- sapply(
        X = 1:100000,
        FUN = function(x) {
          rtnorm(
            mean = v[1] / 100,
            sd = v[2] / 100, 
            low = -0.3,
            high = 0.6
          )}) |> bc_index()
      return(c(mean(out), sd(out)))
    }
    
    out <- purrr::map(.x = list(calves, juveniles, adults, lactating),
                      .f = ~get_bc(.x)) |> 
      purrr::set_names(nm = "Calves", "Juveniles", "Adults", "Females (lactating)")
    
  }
  return(out)
  
}

# From <broman> package
switchv <- function (EXPR, ...) 
{
  result <- EXPR
  for (i in seq(along = result)) result[i] <- switch(EXPR[i], ...)
  result
}

# NOISE ------------------------------------------------------

# ++ [FUNCTION] Sum overlapping noise surfaces
# ++ [PARAM] spl –– List of input surfaces
# ++ [PARAM] raster –– Logical. Set to TRUE if inputs are in raster form.
# ++ [NOTE] To calculate the resultant sound pressure level (SPL) .
# of multiple sources, the SPLs must be added logarithmically
# ++ [RETURN] Output surface.
add_SPL <- function(spl, raster = FALSE){
  if(raster){
    10*log(Reduce("+", lapply(X = spl, FUN = function(x) 10^(x/10))), base = 10)
  } else {
    10*log(sum(10^(spl/10)), base = 10)
  }
}

# ++ [FUNCTION] Create a uniform dummy raster
# ++ [PARAM] value –– Value to assign to the raster cells. If NULL, defaults to 0.
# ++ [NOTE] The dummy raster has the same extent and resolution as the NARW density surfaces.
# ++ [RETURN] Output raster
dummy_raster <- function(value = NULL){
  dummy.r <- density_support |> raster::raster()
  dummy.r <- dummy.r - 1
  dummy.r[dummy.r < 0] <- NA
  if(!is.null(value)) dummy.r[dummy.r == 0] <- value
  return(dummy.r)
}

# ++ [FUNCTION] Calculate the amount of acoustic transmission for a given range from a noise source
# ++ [PARAM] r –– Numeric. Range in km.
# ++ [PARAM] logfac –– Numeric. Log-range coefficient.
# ++ [PARAM] a –– Numeric. Absorption coefficient.
# ++ [NOTE] Assumes a simple propagation model in which transmission loss (TL) 
# depends on log-range (logfac) and frequency-specific absorption (a).
# ++ [RETURN] Numeric value.
TL <- function(r, logfac = 15, a = 1.175){
  # r in km
  # alpha is in dB/km
  loss <- logfac * log10(r*1000)
  loss[loss < 0] <- 0
  loss <- loss + a * r
  return(loss)
}

# ++ [FUNCTION] Calculate the range at which a desired received sound level is attained.
# ++ [PARAM] target.dB –– Numeric. Target received sound level in dB.
# ++ [PARAM] SL –– Numeric. Source level in dB.
# ++ [PARAM] logfac –– Numeric. Log-range coefficient.
# ++ [PARAM] a –– Numeric. Absorption coefficient.
# ++ [PARAM] mitigation –– Numeric. Magnitude of noise attenuation from noise abatement systems (in dB).
# ++ [NOTE] Assumes a simple propagation model in which transmission loss (TL) 
# depends on log-range (logfac) and frequency-specific absorption (a).
# ++ [RETURN] Range in km.
dB2km <- function(target.dB, SL = 220, logfac = 15, a = 1.175, mitigation = 0){
  opt.fun <- function(r, SL, target.L, logfac){
    return(((SL-mitigation)-TL(r, logfac = logfac, a = a)-target.L)^2)}
  out <- optimize(f = opt.fun, interval = c(0,30000), SL = SL, target.L = target.dB, logfac = logfac)
  return(out$minimum)
}

# ++ [FUNCTION] Calculate the sound level received at a specified range from a point source.
# ++ [PARAM] r –– Numeric. Range in km.
# ++ [PARAM] SL –– Numeric. Source level in dB.
# ++ [PARAM] logfac –– Numeric. Log-range coefficient.
# ++ [PARAM] a –– Numeric. Absorption coefficient.
# ++ [PARAM] mitigation –– Numeric. Magnitude of noise attenuation from noise abatement systems (in dB).
# ++ [NOTE] Assumes a simple propagation model in which transmission loss (TL) 
# depends on log-range (logfac) and frequency-specific absorption (a).
# ++ [RETURN] Range in km.
km2dB <- function(r, SL = 220, logfac = 15, a = 1.175, mitigation = 0){
  return((SL-mitigation)-TL(r, logfac = logfac, a = a))
}

# POSTERIOR SAMPLING -----------------------------------------------------

# ++ [FUNCTION] Posterior simulation with gam fits
# ++ [PARAM] b –– A fitted GAM model object.
# ++ [PARAM] ns –– Integer. Number of samples to generate.
# ++ [PARAM] burn –– Integer. Length of the initial burn-in period to discard.
# ++ [PARAM] t.df –– Integer. Degrees of freedom for static multivariate t proposal. Lower for heavier tailed proposals.
# ++ [PARAM] rw.scale –– Numeric. actor by which to scale posterior covariance matrix when generating random walk proposals. Negative or non finite to skip the random walk step.
# ++ [PARAM] thin –– Integer. Retain only every <thin> samples.
# ++ [NOTE] From mgcv::gam.mh –– Includes a fix by Len Thomas
# ++ [RETURN] Matrix of sampled coefficients.
# mgcv::gam.mh with fix via Len Thomas
gam.mh <- function(b, ns = 10000, burn = 1000, t.df = 40, rw.scale = .25, thin = 1) {
  ## generate posterior samples for fitted gam using Metroplois Hastings sampler
  ## alternating fixed proposal and random walk proposal, both based on Gaussian
  ## approximation to posterior...
  if (inherits(b, "bam")) stop("not usable with bam fits")
  beta <- coef(b)
  Vb <- vcov(b)
  X <- model.matrix(b)
  burn <- max(0, burn)
  prog <- interactive()
  iprog <- 0
  di <- floor((ns + burn) / 100)
  if (prog) {
    prg <- txtProgressBar(
      min = 0, max = ns + burn, initial = 0,
      char = "=", width = NA, title = "Progress", style = 3
    )
  }
  bp <- rmvt(ns + burn, beta, Vb, df = t.df) ## beta proposals
  bp[1, ] <- beta ## Don't change this after density step!!
  lfp <- dmvt(t(bp), beta, Vb, df = t.df) ## log proposal density

  rw <- is.finite(rw.scale) && rw.scale > 0
  if (rw) {
    R <- chol(Vb)
    step <- rmvn(ns + burn, beta * 0, Vb * rw.scale) ## random walk steps (mgcv::rmvn)
  }
  u <- runif(ns + burn)
  us <- runif(ns + burn) ## for acceptance check
  bs <- bp
  j <- 1
  accept <- rw.accept <- 0
  lpl0 <- lpl(b, bs[1, ], X)
  for (i in 2:(ns + burn)) { ## MH loop
    ## first a static proposal...
    lpl1 <- lpl(b, bs[i, ], X)
    if (u[i] < exp(lfp[j] - lfp[i] + lpl1 - lpl0)) {
      lpl0 <- lpl1
      accept <- accept + 1
      j <- i ## row of bs containing last accepted beta
    } else {
      bs[i, ] <- bs[i - 1, ]
    }
    ## now a random walk proposal...
    if (rw) {
      lpl1 <- lpl(b, bs[i, ] + step[i, ], X)
      if (us[i] < exp(lpl1 - lpl0)) { ## accept random walk step
        lpl0 <- lpl1
        j <- i
        bs[i, ] <- bs[i, ] + step[i, ]
        rw.accept <- rw.accept + 1
        # lfp[i] <- dmvt(bs[i,],beta,Vb,df=4,R=R) ## have to update static proposal density
        # FIX via LJT 5/10/20
        lfp[i] <- dmvt(bs[i, ], beta, Vb, df = t.df, R = R)
      }
    }
    if (i == burn) accept <- rw.accept <- 0
    if (prog && i %% di == 0) setTxtProgressBar(prg, i)
  } ## MH loop
  if (burn > 0) bs <- bs[-(1:burn), ]
  if (thin > 1) bs <- bs[seq(1, ns, by = thin), ]
  if (prog) {
    setTxtProgressBar(prg, i)
    cat("\n")
    cat("fixed acceptance = ", accept / ns, "  RW acceptance = ", rw.accept / ns, "\n")
  }
  list(bs = bs, rw.accept = rw.accept / ns, accept = accept / ns)
} ## gam.mh

# other parts of mcmc.r in mgcv that are not exported
## t-distribution stuff for mgcv.
## (c) Simon N. Wood (2020)

## some useful densities (require mgcv::rmvn)...

## simulate multivariate t variates
rmvt <- function(n, mu, V, df) {
  y <- rmvn(n, mu * 0, V)
  v <- rchisq(n, df = df)
  t(mu + t(sqrt(df / v) * y))
}

r.mvt <- function(n, mu, V, df) rmvt(n, mu, V, df)

## multivariate t log density...
dmvt <- function(x, mu, V, df, R = NULL) {
  p <- length(mu)
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R), x - mu)
  k <- -sum(log(diag(R))) - p * log(df * pi) / 2 + lgamma((df + p) / 2) - lgamma(df / 2)
  k - if (is.matrix(z)) (df + p) * log1p(colSums(z^2) / df) / 2 else (df + p) * log1p(sum(z^2) / df) / 2
}

d.mvt <- function(x, mu, V, df, R = NULL) dmvt(x, mu, V, df, R)

# Taken from mgcv source code
# some functions to extract important components of joint density from fitted gam
## evaluate penalty for fitted gam, possibly with new beta
# patched to include parapen
bSb <- function(b, beta = coef(b)) {
  bSb <- k <- 0
  sp <- if (is.null(b$full.sp)) b$sp else b$full.sp ## handling linked sp's

  # the parapen bits
  # need to do something clever with L at some point
  if (!is.null(b$paraPen)) {
    for (i in 1:length(b$paraPen$S)) {
      k <- k + 1
      # get indices
      ii <- b$paraPen$off[i]
      ii <- ii:(ii + ncol(b$paraPen$S[[i]]) - 1)
      # add to penalty
      bSb <- bSb + sp[b$paraPen$full.sp.names[i]] *
        (t(beta[ii]) %*% b$paraPen$S[[i]] %*% beta[ii])
    }
  }


  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + sp[k] * (t(beta[ii]) %*% b$smooth[[i]]$S[[j]] %*% beta[ii])
      }
    }
  }


  bSb
} ## bSb

devg <- function(b, beta = coef(b), X = model.matrix(b)) {
  ## evaluate the deviance of a fitted gam given possibly new coefs, beta
  ## for general families this is simply -2*log.lik
  if (inherits(b$family, "general.family")) {
    -2 * b$family$ll(b$y, X, beta, b$prior.weights, b$family, offset = b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y, b$family$linkinv(X %*% beta + b$offset), b$prior.weights))
  }
} ## devg

lpl <- function(b, beta = coef(b), X = model.matrix(b)) {
  ## log joint density for beta, to within uninteresting constants
  -(devg(b, beta, X) / b$sig2 + bSb(b, beta) / b$sig2) / 2
}

# SPATIAL DATA ------------------------------------------------------

# ++ [FUNCTION] Spatial coordinate systems (projected + geographic)
# ++ [PARAM] latlon –– Logical. If TRUE, returns a lat/lon CRS, otherwise uses
# an Albers Equal Area projection
# ++ [RETURN] Object of class CRS
narw_crs <- function(latlon = FALSE){
  if(latlon) sp::CRS("+init=epsg:4326") else sp::CRS("+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs") 
}

# ++ [FUNCTION] Returns the extent coordinates of an input SpatialGridDataFrame
# ++ [PARAM] spgdf –– Input SpatialGridDataFrame.
# ++ [RETURN] Numeric vector
get_limits <- function(spgdf){
  out <- sp::coordinates(spgdf)
  return(c(range(out[,1]), range(out[,2])))
}

# ++ [FUNCTION] Add projected coordinates to a data.frame containing lat/lon
# ++ [PARAM] dat –– Input data.frame.
# ++ [PARAM] CRS.obj –– Projected coordinate system to use.
# ++ [RETURN] Same data.frame with two new columns "x" and "y"
add_xy <- function(dat, CRS.obj = narw_crs()){
  
  tmp <- sp::SpatialPoints(
    coords = dat[, c("longitude", "latitude")],
    proj4string = sp::CRS("+proj=longlat")) |>
    sp::spTransform(CRSobj = CRS.obj)
  
    xy <- sp::coordinates(tmp) |> as.data.frame()
    names(xy) <- c("x", "y")
    cbind(dat, xy)
 
}

# ++ [FUNCTION] Add geographic coordinates (lat/lon) to a data.frame containing projected coordinates
# ++ [PARAM] dat –– Input data.frame.
# ++ [PARAM] CRS.obj –– Projected coordinate system to use.
# ++ [RETURN] Same data.frame with two new columns "long" and "lat"
add_latlon <- function(dat, CRS.obj = narw_crs()){
  
  tmp <- sp::SpatialPoints(
    coords = dat[, c("easting", "northing")],
    proj4string = narw_crs()) |>
    sp::spTransform(CRSobj = sp::CRS("+proj=longlat"))
  
  xy <- sp::coordinates(tmp) |> as.data.frame()
  names(xy) <- c("long", "lat")
  cbind(dat, xy)
  
}

# PLOTTING ------------------------------------------------------

# ++ [FUNCTION] Custom ggplot2 theme
# ++ [PARAM] vertical –– Logical. If TRUE, axis labels are rotated 90 degrees.
theme_narw <- function(vertical = FALSE, grey = TRUE, bbox = FALSE){
  
  # font <- "Georgia"   # Assign font family up front
  
  if(grey) p <- ggplot2::theme_grey() else p <- ggplot2::theme_bw()
  
  p %+replace%    # Replace elements we want to change
    
    ggplot2::theme(
      
      # Axis elements
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.text.y = ggplot2::element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), 
                                          angle = ifelse(vertical, 90, 0),
                                          vjust = ifelse(vertical, 0.5, 0),
                                          hjust = ifelse(vertical, 0.5, 0)),
      axis.text.x = ggplot2::element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), angle = 90),
      axis.title.x = ggplot2::element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
      
      # Panel elements
      panel.border = ggplot2::element_rect(colour = ifelse(bbox, "black", "transparent"), 
                                           fill = NA, linewidth = 0.75),
      
      # Facet elements
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text = ggplot2::element_text(colour = 'white', size = 12),
      strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "cm")),
      strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm"), angle = -90)
      
    )
}

# ++ [FUNCTION] Plot rasters and SpatialGridDataFrames using ggplot
# ++ [PARAM] r –– Input raster or SpatialGridDataFrame.
# ++ [PARAM] duke –– Logical. Whether to match the formatting of Jason Roberts' density maps.
# ++ [PARAM] positive –– Logical. If TRUE, only retain positive values.
# ++ [PARAM] quantile –– Logical. If TRUE, colour breaks are based on quantiles of the data.
# ++ [PARAM] breaks –– Numeric vector defining custom colours breaks, if necessary.
# ++ [PARAM] title –– Character. Plot title.
# ++ [PARAM] xmin –– Numeric. Lower limit for the x-axis.
# ++ [PARAM] xmax –– Numeric. Upper limit for the x-axis.
# ++ [PARAM] ymin –– Numeric. Lower limit for the y-axis.
# ++ [PARAM] ymax –– Numeric. Upper limit for the y-axis.
plot_raster <- function(r, 
                        duke = FALSE, 
                        positive = FALSE, 
                        quantile = FALSE, 
                        breaks = NULL, 
                        title = NULL,
                        xmin = NULL,
                        xmax = NULL,
                        ymin = NULL,
                        ymax = NULL){
  
  if(class(r) == "SpatialGridDataFrame") r <- raster::raster(r)
  
  dat <- raster::as.data.frame(r, xy = TRUE)
  names(dat)[3] <- "Nhat"
  dat <- dat[complete.cases(dat),]
  
  if(positive) dat <- dat |> dplyr::filter(Nhat > 0)
  
  if (is.null(breaks)) {
    if (quantile) {
      colour.breaks <- c(0, quantile(dat$Nhat, seq(0, 1, 0.1)))
    } else {
      if (duke) {
        colour.breaks <- colour_breaks(dat)
      } else {
        colour.breaks <- c(0, rgeoda::natural_breaks(10, dat[, "Nhat", drop = FALSE]), max(dat$Nhat))
      }
    }
  } else {
    colour.breaks <- breaks
  }
  
  colour.breaks <- unique(colour.breaks)
  
  all.breaks <- data.frame(Ncol = cut(colour.breaks, breaks = colour.breaks, include.lowest = TRUE))
  all.breaks$mapcol <- pals::viridis(nrow(all.breaks))
  
  dat <- suppressWarnings(dat |> 
    dplyr::mutate(Ncol = cut(Nhat, breaks = colour.breaks, include.lowest = TRUE)) |> 
    dplyr::mutate(Ncol = factor(Ncol)) |> 
    dplyr::left_join(y = all.breaks, by = "Ncol"))
  
  levels(dat$Ncol) <-  gsub(pattern = ",", replacement = " – ", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\[|\\]", replacement = "", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\(|\\)", replacement = "", levels(dat$Ncol))
  
  r.dat <- raster::rasterFromXYZ(dat[, c("x", "y", "Nhat")], res = raster::res(r), crs = narw_crs())
  
  if(all(is.null(c(xmin, xmax, ymin, ymax)))){
    
    x.lim <- raster::extent(r.dat)[1:2]
    y.lim <- raster::extent(r.dat)[3:4]
    
  } else {
    
    x.lim <- c(xmin, xmax)
    y.lim <- c(ymin, ymax)
    
  }
  
  world.sf <- sf::st_as_sf(world)
  
  ocean.colour <- "#A9C9D2"
  land.colour <- "#E1D4B2"
  land.outline <- "#75684A"
  country.colour <- "#E9E2D0"
  state.colour <- "#43391F"
  grid.line.colour <- "#D4EAEF"
  country.text <- "#43391F"
  scale.colour <- "#35484E"
  ocean.feature <- "#6899A7"
  
  p <- ggplot2::ggplot(data = dat) +  
    ggplot2::geom_tile(aes(x,y, fill = Ncol)) +
    ggplot2::geom_sf(data = world.sf, fill = country.colour, color = "black", size = 0.25) +
    ggplot2::geom_sf(data = gis$map$states, fill = land.colour, col = land.outline, size = 0.15) +
    # ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    ylab("") + xlab("") +
    theme_narw(bbox = TRUE) + 
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::coord_sf(xlim = x.lim, ylim = y.lim, expand = FALSE) +
    ggplot2::scale_fill_manual(values = pals::viridis(length(levels(dat$Ncol))),
                               guide = guide_legend(reverse = TRUE),
                               drop = FALSE,
                               na.value = "transparent") +
    ggplot2::annotate("text", x = 1200, y = 0,
                      label = "ATLANTIC \nOCEAN",
                      fontface = "bold.italic", 
                      color = "white", size = 4) +
    ggplot2::annotate("text", x = -100, y = 300,
                      label = "UNITED \nSTATES",
                      fontface = "bold",
                      color = country.text, size = 4) +
    ggplot2::annotate("text", x = 0, y = 1400,
                      label = "CANADA",
                      fontface = "bold",
                      color = country.text, size = 4)
  
  if(!is.null(title)) p <- p + ggplot2::ggtitle(title)
  
  return(p)
  
}

# ++ [FUNCTION] Replicate colour breaks used in Jason Robert's density maps
# ++ [PARAM] dat –– Input data.frame containing density values
# ++ [RETURN] Vector of breaks
colour_breaks <- function(dat){
  colour.breaks <- 25*c(0,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3,10)/100
  # colour.breaks <- c(colour.breaks, seq(max(colour.breaks), ceiling(max(dat$Nhat, na.rm = TRUE)), length.out = 5))
  colour.breaks <- round(colour.breaks[!duplicated(colour.breaks)],3)
  return(colour.breaks)
}

# UTILITIES ------------------------------------------------------

compareSim <- function(...){

   args <- list(...)
  
   allcombs <- combn(seq_along(args), 2)

   match.test <- lapply(X = 1:ncol(allcombs), FUN = function(x){
     
     purrr::map2(.x = args[[allcombs[1,x]]]$sim,
                 .y = args[[allcombs[2,x]]]$sim,
                 .f = ~{
                  out <- data.table:::all.equal.data.table(target = .x[, list(easting, northing)], 
                                                     current = .y[, list(easting, northing)], tolerance = 1e-3)
                  if(is.character(out)) FALSE else out
                     
                 }) |> do.call(what = c)
     })
   
   do.match <- do.call(match.test, what = c)
   return(do.match)
   
    # match.test <- purrr::set_names(x = cohorts[id %in% cohort, abb]) |>
    #   purrr::map(.f = ~{
    #     exclude.cols <- c("alive", "bc", "mass", "fatmass", "leanmass",
    #                       "abort", "delta_fat", "delta_m", "noise_resp", "noise_lvl", "strike_risk",
    #                       "E_tot", "E_in", "E_out", "feed_effort", "LC", "t_feed", "t_rest_nurse")
    #     dt1 <- outsim$sim[[.x]][, -exclude.cols, with = FALSE]
    #     dt2 <- pair$sim[[.x]][, -exclude.cols, with = FALSE]
    #     data.table:::all.equal.data.table(target = dt1, current = dt2, tolerance = 1e-3)
    #   }) |> tibble::enframe() |>
    #   dplyr::rename(cohort = name, paired = value) |>
    #   tidyr::unnest(cols = c(paired)) |>
    #   dplyr::mutate(paired = !is.character(paired)) |>
    #   data.table::as.data.table()
    # 
    # if(all(match.test$paired)){
    #   console(msg = "Match testing", suffix = tickmark())
    # } else {
    #   console(msg = "Match testing", suffix = crossmark())
    # }

}



# ++ [FUNCTION] Round any number to desired accuracy
# ++ [PARAM] x –– Numeric vector
# ++ [PARAM] accuracy –– Number to round to
# ++ [PARAM] f –– Rounding function
roundany <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

# ++ [FUNCTION] Print current date and time
date_time <- function(){
  cat("Date:", as.character(Sys.Date()), "\n")
  now.time <- Sys.time() |> 
    format('%Y-%m-%d %H:%M:%S') |> 
    stringr::str_split(" ") |> 
    purrr::map_chr(2)
  timezone <- lubridate::tz(lubridate::ymd_hms(Sys.time()))
  cat("Time: ", now.time, " (", timezone, ")\n\n", sep = "")
}

# ++ [FUNCTION] Run basic checks on simulator outputs
# ++ [PARAM] obj –– An object of class <narwsim>, as returned by <narw>
# ++ [RETURN] Prints warnings if errors are detected
check <- function(obj){
  
  max_bc <- purrr::map(.x = obj[["sim"]], .f = ~max(.x$bc)) |> do.call(what = c)
  
  if(any(is.na(max_bc))) warning("NA values found in body condition column")
  if(any(max_bc[!is.na(max_bc)] > find_maxBC())) warning("Body condition value(s) > max")
  if(any(max_bc[!is.na(max_bc)] < 0)) warning("Body condition value(s) < 0")
  
  time_feed <- purrr::map(.x = obj[["sim"]], .f = ~max(.x$t_feed)) |> do.call(what = c)
  time_rest <- purrr::map(.x = obj[["sim"]], .f = ~max(.x$t_rest_nurse)) |> do.call(what = c)
  time_travel <- purrr::map(.x = obj[["sim"]], .f = ~max(.x$t_travel)) |> do.call(what = c)
  time_tot <- purrr::map(.x = obj[["sim"]], .f = ~{
    time_total <- .x[, list(tot = sum(t_feed,t_travel, t_rest_nurse)), list(day, whale)]
    max(time_total$tot)
  }) |> do.call(what = c)
  
  if(any(is.na(time_tot))) warning("NA values found in activity budgets")
  if(any(time_tot[!is.na(time_tot)] > 24)) warning("Budget > 24 hrs")
  if(any(time_feed[!is.na(time_feed)] < 0) | 
     any(time_rest[!is.na(time_rest)] < 0) |
     any(time_travel[!is.na(time_travel)] < 0)
     ) warning("Budget < 0")
  
  swim_speed <- purrr::map(.x = obj[["sim"]], .f = ~max(.x$swimspeed)) |> do.call(what = c)
  if(any(is.na(swim_speed))) warning("NA values found in swimming speeds")
  if(any(is.infinite(swim_speed))) warning("Inf found in swimming speeds")
  
}

# ++ [FUNCTION] Add a label to a simulation
# ++ [PARAM] obj –– An object of class <narwsim>, as returned by <narw>
# ++ [PARAM] lb –– Character. Label to assign to the object
label <- function(obj, lb){
  if(!inherits(obj, "narwsim") & !inherits(obj, "narwproj")) 
    stop("Input object must be of class <narwsim> or <narwproj>")
  obj$param$label <- lb
  return(obj)
}

# ++ [FUNCTION] Text-based timeline of wind farm development activities
# ++ [PARAM] schedule –– Vector of integers indicating which phase of development takes place in which year
proj_timeline <- function(schedule, burn){
  
  operator <- " -----> "
  
  schedule <- schedule[(burn+1):length(schedule)]
  
  # Baseline
  if(any(schedule == 0)){
    if(all(schedule == 0)) base.msg <- paste0("Baseline (yr 1 to ", length(schedule) - 1, ")") else base.msg <- NULL
  } else {
    base.msg <- NULL
  }
  
  # Construction phase
  if(any(schedule == 1)){
    const.msg <- paste0(if(!is.null(base.msg)) operator, "Construction (yr ", min(which(schedule == 1) - 1), ")")
  } else {
    const.msg <- NULL
  }
  
  # O&M phase
  if(any(schedule == 2)){
    if(all(schedule == 2)) om.msg <- paste0("Operation & maintenance (yr 1 to ", length(schedule) - 1, ")") else om.msg <- paste0(operator, "Operation & maintenance (yr ", min(which(schedule == 2) - 1), " to ", max(which(schedule == 2) - 1), ")") 
  } else {
    om.msg <- NULL
  }
  
  cat(base.msg, const.msg, om.msg, sep = "")
  cat("\n\n")
}

# ++ [FUNCTION] Split vector at NA values and zeroes
# ++ [PARAM] x –– Input vector
split_NAzero <- function(x){
  idx <- 1 + cumsum( is.na(x) | x == 0)
  not.na <- !is.na(x) & !x == 0
  split(x[not.na], idx[not.na])
}

# ++ [FUNCTION] Split vector at NA values
# ++ [PARAM] x –– Input vector
split_NA <- function(x){
  idx <- 1 + cumsum(is.na(x))
  not.na <- !is.na(x)
  split(x[not.na], idx[not.na])
}

# ++ [FUNCTION] Simple tick mark
tickmark <- function(){
  return("[ok]")
  # return("\U2714")
}

# ++ [FUNCTION] Simple cross mark
crossmark <- function(){
  return("[x]")
  # return("\U2717")
}

# ++ [FUNCTION] Print current status in console 
console <- function(msg, suffix = NULL){
  if(is.null(suffix)){
    cat("\r", paste0(msg, " ..."), sep = "")
  } else {
    if(nchar(suffix) > 1){
      blank <- paste0(rep(" ", max(nchar(suffix)-3,2)), collapse = "")
    } else {
      blank <- "  "
    }
    cat("\r", paste0(msg, " ", suffix, blank), sep = "")
    cat("\n")
  }
}

# ++ [FUNCTION] Generate date sequence for the simulation
# ++ [PARAM] start.month –– Integer indicating the month in which the simulation begins.
# ++ [PARAM] ndays –– Integer. Duration of the simulation in days.
# ++ [PARAM] strip –– Logical. If TRUE, remove year from output date.
# ++ [PARAM] gantt –– Logical. If TRUE, update year for use in Shiny Gantt chart.
# ++ [RETURN] A character vector of dates
get_dates <- function(start.month = 10, ndays = 457, strip = TRUE, gantt = FALSE){
  
  thisyear <- lubridate::year(lubridate::now())
  
  if(gantt){
    
    date_seq <- seq(lubridate::ymd(paste0(thisyear, "-01-01")), lubridate::ymd(paste0(thisyear, "-12-31")), by = '1 day')
    if(sum(grepl(pattern = "02-29", x = date_seq))>0) date_seq <- date_seq[-which(grepl(pattern = "02-29", x = date_seq))]
    
  } else {
    
    # Use 2021 as non-leap year, as AIS data were obtained for 2019, which is a non-leap year as well
    start.date <- lubridate::as_date(paste0(2021, "-", stringr::str_pad(start.month, width = 2, pad = "0"), "-01"))
    date_seq <- seq(start.date, by = 'day', length.out = ndays)
    if(strip){
      date_seq <- format(as.Date(date_seq), "%d-%m") # Strip year
      date_seq <- c(paste0("00-", start.month), date_seq)
    } else {
      date_seq <- c(paste0("2021-", start.month, "-00"), as.character(date_seq))
    }
  }
  return(date_seq)
}

# ++ [FUNCTION] Fix paths to vignette assets (images)
fix_paths_vignettes <- function(vignette.name = "narwind"){
  vignette.path <- file.path("./docs/articles", paste0(vignette.name, ".html"))
  tx  <- readLines(vignette.path)
  tx_mod  <- gsub(pattern = "../../../../../My%20Drive/Documents/git/narwind/articles/", replace = "", x = tx)
  writeLines(tx_mod, con = vignette.path)
}

# ++ [FUNCTION] Simple plots of simulator outputs
# ++ [PARAM] obj –– An object of class <narwsim>, as returned by <narw>.
# ++ [PARAM] cohort –– Cohort label
# ++ [PARAM] whaleID –– Individual ID number.
# ++ [PARAM] calf –– Logical. If TRUE, includes calf parameters.
inspect <- function(obj, cohort = "ad(f,l)", whaleID = 1, calf = FALSE){
  pdf(file = "allplots.pdf")
  par(mfrow = c(3,3))
  if(calf){
    nn <- names(m$sim[[1]])[7:112][grepl(pattern = "calf", names(m$sim[[1]])[7:112])]
  } else {
    nn <- names(m$sim[[1]])[7:112]
  }
  for(i in nn){
    dat <- obj[["sim"]][[cohort]][whale == whaleID, c("day", i), with = FALSE]
    plot(dat, xlab = "Day", ylab = "", main = i, type = "l")
  }
  par(mfrow = c(1,1))
  dev.off()
}

# ++ [FUNCTION] Format data.tables for use in <summary> by adding %
# ++ [PARAM] dt –– Input data.table
# ++ [PARAM] relative –– Logical. If TRUE, percentages are calculated relative to class totals.
# ++ [PARAM] N –– Total to use.
# ++ [PARAM] direction –– One of "row" or "col". Indicates how totals are calculated.
format_dt <- function(dt, relative = FALSE, N, direction = "col") {
  if (nrow(dt) > 0) {
    if (relative) {
      dt |>
        janitor::adorn_percentages(denominator = direction) |>
        janitor::adorn_pct_formatting() |>
        janitor::adorn_ns()
    } else {
      dt |>
        janitor::adorn_totals("col") |> 
        dplyr::mutate(Total = N - Total) |> 
        janitor::untabyl() |> 
        janitor::adorn_percentages(denominator = "row") |>
        janitor::adorn_pct_formatting() |>
        janitor::adorn_ns() |> 
        dplyr::select(-Total)

    }
  } else {
    dt
  }
}

# ++ [FUNCTION] Convert a 3D array to a data.table
# ++ [PARAM] a –– Input array
array2dt <- function(a){
  y <- aperm(a, c(1, 3, 2))
  dim(y) <- c(prod(dim(a)[-2]), dim(a)[2])
  y <- data.table::data.table(y)
  names(y) <- colnames(a)
  return(y)
}