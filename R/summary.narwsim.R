#' Summary
#'
#' Summary information
#' @import data.table
#' @export
#' 
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' animals <- run_model(10)
#' plot(animals)
#' }

summary.narwsim <- function(obj, 
                            what = "all", 
                            plot = FALSE,
                            ...){
  
  options(pillar.sigfig = 7)
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default optional arguments
  cohortID <- obj$param$cohortID
  whaleID <- 1:obj$param$nsim
  
  # Default values
  if(length(args) > 0){
    if("cohortID" %in% names(args)) cohortID <- args[["cohortID"]]
    if("whaleID" %in% names(args)) whaleID <- args[["whaleID"]]
  }
  
 # Preamble ---------------------------------------------------------------
  
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")
  
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n")
  cat("\n")
  cat("     NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("                --- MODEL SUMMARY ---\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n\n")
  
  cat("BOF: Bay of Fundy (lower, upper)\n")
  cat("CABOT: Cabot Strait\n")
  cat("CCB: Cape Cod Bay\n")
  cat("GOM: Gulf of Maine and Georges Bank\n")
  cat("GSL: Gulf of St Lawrence\n")
  cat("MIDA: Mid-Atlantic\n")
  cat("SCOS: Scotian Shelf\n")
  cat("SEUS: South-east United States\n")
  cat("SNE: Southern New England\n\n")

  n.ind <- obj$param$nsim
  cohorts <- obj$param$cohorts

  cohort.ab <- cohorts[id %in% cohortID, abb]
  cohort.names <- cohorts[id %in% cohortID, name]
  init.month <- obj$init$month
  
  sim <- purrr::map(.x = obj$sim[cohort.ab],
                    .f = ~.x[whale %in% whaleID & cohort %in% cohortID, ]) |> 
    purrr::set_names(nm = cohort.ab)

  gamdat <- obj$gam$dat[whale %in% whaleID,]
  dead.df <- obj$dead[whale %in% whaleID,]
  
  if(5 %in% cohortID){
    gamdat <- gamdat[cohort %in% c(0,cohortID),]
    dead.df <- dead.df[cohort %in% c(0,cohortID),]
  } else {
    gamdat <- gamdat[cohort %in% cohortID,]
    dead.df <- dead.df[cohort %in% cohortID,]
  }
  
  cat("=============================================================\n")
  cat("SIMULATIONS\n")
  cat("=============================================================\n\n")
  
  cat("No. animals:", format(n.ind, big.mark = ","), "\n\n")
  cat("Cohort(s)\n")
  cat("----------\n")
  for (h in cohortID) cat(cohorts[id==h, abb], ": ", cohorts[id==h, name], "\n", sep = "")
  cat("\n")
  cat("Simulation start:", month.name[init.month], "\n")
  cat("\n")
  
  # cat("% SEUS:", sum(sim[[1]][day > 0, list(sum(unique(seus))), whale][, 2]) / n.ind, "\n\n")

  # Health ---------------------------------------------------------------
    
  if("health" %in% what | what == "all"){
    
    cat("=============================================================\n")
    cat("HEALTH\n")
    cat("=============================================================\n")
    
    cat("\n+++++++++++ Mortality +++++++++++")
    
    n.dead <- gamdat[whale %in% whaleID & cohort %in% cohortID,
                          list(born = sum(born == 1),
                          alive = sum(alive == 1)), name] |>
      dplyr::mutate(dead = born - alive) |> 
      dplyr::select(-born) |> 
      format_dt(direction = "row") |> 
      dplyr::rename(cohort = name)
    
    if(5 %in% cohortID){
    n.dead.calves <- gamdat[whale %in% whaleID & cohort ==0,
                     list(alive = sum(alive == 1)), name] |>
      dplyr::mutate(dead = nrow(gamdat[whale %in% whaleID & cohort ==0 & event == "birth"]) - alive) |> 
      format_dt(direction = "row") |> 
      dplyr::rename(cohort = name)
    
    print(knitr::kable(data.table::rbindlist(list(n.dead, n.dead.calves)), format = "simple"))
    
    } else {
      
      print(knitr::kable(n.dead, format = "simple"))
    
    }
    
    # Mortality by region
    
    if(5 %in% cohortID){
      
      dead.df[cohort == 0, ] <- dead.df |> 
        dplyr::filter(cohort == 0) |> 
        dplyr::rowwise() |> 
        dplyr::mutate(cause_death = ifelse(sum(starve, strike, died) == 0, 
                      paste0(cause_death, " (female)"), cause_death)) |>
        dplyr::ungroup()
      
    }
    
    n.dead.region <- dead.df[, .N, list(abb, region, cause_death)] |> 
      # dplyr::mutate(strike = dplyr::coalesce(strike, 0),
      #               starve = dplyr::coalesce(starve, 0),
      #               other = dplyr::coalesce(other, 0),
      #               region = dplyr::coalesce(region, "region")) |> 
      # tidyr::pivot_longer(!c(abb, region), names_to = "cause_death", values_to = "count") |> 
      tidyr::pivot_wider(names_from = abb, values_from = N) |> 
      {\(.) {replace(.,is.na(.),0)}}() |> 
      dplyr::arrange(cause_death, region)
    
    if(nrow(n.dead.region)> 0){
      
      n.dead.region <- split(n.dead.region, f = factor(n.dead.region$cause_death)) |> 
        purrr::map(.f = ~format_dt(.x, direction = "col"))
      
      purrr::walk(.x = n.dead.region,
                  .f = ~ print(knitr::kable(.x, format = "simple")))
    
    }
    cat("\n")
    
    if(4 %in% cohortID){
      
      cat("\n+++++++++++ Pregnancies +++++++++++\n") 
      
      cat("\n")
      cat("Abortion rate: ", 100 * as.numeric(obj$abort[whale %in% whaleID, .(abort = sum(abort))]) / nrow(obj$abort), 
          "% (", as.numeric(obj$abort[whale %in% whaleID, .(abort = sum(abort))]), ")\n", sep = "")
    }
    
    if (5 %in% cohortID) {
      
      cat("\n+++++++++++ Births +++++++++++\n\n")

      if(!is.null(obj$birth[[1]])){
        
        birth.df <- obj$birth[whale %in% whaleID, ]
        n.births <- nrow(birth.df)
        cat("No. births:", paste0(n.births, " (", 100 * nrow(birth.df) / n.ind, "%)\n"))
        
        dob <- birth.df[, .(mean_DOB = mean(date), min_DOB = min(date), max_DOB = max(date))]
        
        cat("Mean:", lubridate::day(dob$mean_DOB), 
            as.character(lubridate::month(dob$mean_DOB, label = TRUE, abbr = TRUE)), "\n")
        cat("Range:", lubridate::day(dob$min_DOB), 
            as.character(lubridate::month(dob$min_DOB, label = TRUE, abbr = TRUE)), "–",
            lubridate::day(dob$max_DOB), 
            as.character(lubridate::month(dob$max_DOB, label = TRUE, abbr = TRUE)), "\n")
        cat("\n")
      } else {
        
        cat("No. births: 0% (0)\n\n")
        
      }
    }
    
    ## Body condition ---------------------------------------------------------------
    
    # .....................................................
    # BODY CONDITION
    # .....................................................
    
    if(plot){
      
      bodycondition.df <- purrr::map(.x = sim, .f = ~{
        cc <- unique(.x$cohort)
        if(cc[cc>0] == 5){
          dplyr::select(.x, whale, day, bc, bc_calf, cohort, cohort_name) |> 
            tidyr::pivot_longer(!c("whale", "day", "cohort", "cohort_name"), names_to = "animal", values_to = "bc")
        } else {
          dplyr::select(.x, whale, day, bc, cohort, cohort_name) |> 
            tidyr::pivot_longer(!c("whale", "day", "cohort", "cohort_name"), names_to = "animal", values_to = "bc")
        }}) |> 
        do.call(what = rbind) |> 
        dplyr::mutate(animal = gsub("bc_calf", "Calves", animal)) |> 
        dplyr::mutate(animal = gsub("bc", "Adults", animal))
      
      bodycondition.plot <- 
        purrr::map(.x = cohortID, .f = ~{
          ggplot2::ggplot(data = bodycondition.df |> dplyr::filter(cohort == .x),
                          aes(x = day, y = bc, group = whale)) +
            ggplot2::geom_path(col = "grey70") +
            {if(!.x == 5) ggplot2::facet_wrap(vars(cohort_name), scales = 'free') } +
            {if(.x == 5) ggplot2::facet_grid(vars(animal), vars(cohort_name), scales = 'free') } +
            theme_narw() +
            ylab("Body condition") +
            xlab("Day") +
            ggplot2::scale_x_continuous(breaks = pretty(0:365, n = 10)) +
            ggplot2::scale_y_continuous(
              limits = ~ c(0, ceiling(max(.x))),
              breaks = ~ pretty(.x, 5),
              expand = c(0, 0)) +
            ggplot2::geom_hline(yintercept = 0.05, col = "#cb4154")
        })
      
      if(length(cohortID) == 1){
        plot(bodycondition.plot[[1]])
      } else {
        bcp <- patchwork::wrap_plots(bodycondition.plot, 
                              nrow = ifelse(length(cohortID) == 2, 1, 2), 
                              ncol = ifelse(length(cohortID) == 2, 1, 3))
        print(bcp)
      }

      # Growth ---------------------------------------------------------------
      
      # # .....................................................
      # # GROWTH
      # # .....................................................
      # 
      # body_growth_l <-
      #  purrr::map(
      #    .x = cohortID,
      #    .f = ~ growth_curve(param = "length", 
      #                        obj = obj, 
      #                        cohortID = .x, 
      #                        whaleID = whaleID,
      #                        ylabel = "Body length (m)")
      #  )
      # 
      # if (4 %in% cohortID) {
      #   body_growth_l <-
      #     append(body_growth_l, 
      #            list(growth_curve(param = "fetus_l", 
      #                              obj = obj, 
      #                              cohortID = 4, 
      #                              whaleID = whaleID,
      #                              ylabel = "Body length (m)")))
      # }
      # 
      # if (5 %in% cohortID) {
      #   body_growth_l <-
      #     append(body_growth_l, 
      #            list(growth_curve(param = "length_calf", 
      #                              obj = obj, 
      #                              cohortID = 5, 
      #                              whaleID = whaleID,
      #                              ylabel = "Body length (m)")))
      # }
      # 
      # growth.plot <- patchwork::wrap_plots(purrr::discard(.x = body_growth_l, .p = ~is.list(.x) & length(.x) == 1))
      # print(growth.plot)
      
      } # End plot = TRUE
  } # End what == health
  
  cat("\n")
  
  move.out <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~{
      all.d <- sim[[.x]][day > 0, .SD[1:(ifelse(which.min(alive)==1, 365, which.min(alive)-1))], 
                         .SDcols = c("d_travel", "region", "day", "seus", "gsl"), whale]
      all.d[, d_migr := sum(d_travel), whale]
      all.d$cohort <- .x
      all.d[]
    })
  
  # Movements ---------------------------------------------------------------
  
  if("movements" %in% what | what == "all"){
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================")
  
  locs.by.area <- purrr::map(.x = cohort.ab, .f = ~
    dplyr::right_join(move.out[[.x]][, .N, region],
      data.table::data.table(region = regions$region, country = regions$country),
      by = "region"
    ) |> dplyr::mutate(N = dplyr::coalesce(N, 0), cohort = .x)) |> 
    do.call(what = rbind)

  locs.by.country <- locs.by.area |> 
    dplyr::group_by(cohort, country) |> 
    dplyr::summarise(n = sum(N), .groups = "keep") |> 
    dplyr::ungroup() |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    format_dt()

  locs.by.region <- locs.by.area |> 
    dplyr::group_by(cohort, region) |> 
    dplyr::summarise(n = sum(N), .groups = "keep") |> 
    dplyr::ungroup() |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    format_dt()
    
  print(knitr::kable(locs.by.region, format = "simple"))
  print(knitr::kable(locs.by.country, format = "simple"))
  cat("\n")
  
  cat("=============================================================\n")
  cat("MOVEMENTS (km)\n")
  cat("=============================================================\n")
  
  if(plot){
  par(mfrow = c(1,2))
  purrr::walk(.x = seq_along(move.out), .f = ~{
    
    # Step lengths
    hist(move.out[[.x]]$d_travel, freq = FALSE, main = cohort.names[.x], xlab = "Step length (km)", breaks = 20)
    lines(density(move.out[[.x]]$d_travel, adjust = 2), col = "orange", lwd = 1.5)
    boxplot(d_travel~region, data = move.out[[.x]], ylab = "Step lengths (km)", xlab = "")
    
  })
  par(mfrow = c(1,1))
  }
  
  #'---------------------------------------
  # Summary of SEUS/GSL visitation by cohort
  #'---------------------------------------
  journey.seus <- purrr::map(.x = cohort.ab, .f =~{

    visits <- dplyr::left_join(move.out[[.x]][, list(target_SEUS = unique(seus), reach_SEUS = as.numeric(sum(grepl(pattern = "SEUS", x = region)) > 0)), whale], move.out[[.x]][, list(target_GSL = unique(gsl), reach_GSL = as.numeric(sum(grepl(pattern = "GSL", x = region)) > 0)), whale], by = "whale") |> dplyr::select(-whale) |> dplyr::mutate_if(is.numeric, list(~factor(., levels = c("0", "1"))))
    
    visits |> janitor::tabyl(target_SEUS, reach_SEUS) |> 
      dplyr::rename(SEUS = target_SEUS, reach_not = '0', reach = '1') |> 
      dplyr::mutate(cohort = .x)
      
  }) |> do.call(what = rbind) |> 
    dplyr::relocate(cohort, .before = "SEUS")
  
  journey.seus.all <- data.table::data.table(journey.seus)
  journey.seus.all <- journey.seus.all[,.(reach_not = sum(reach_not), reach = sum(reach)), SEUS]
  journey.seus.all <- format_dt(journey.seus.all, direction = "row")
  journey.seus[, 2:4] <- format_dt(data.frame(journey.seus[, 2:4]), direction = "row")
  
  journey.gsl <- purrr::map(.x = cohort.ab, .f =~{
    
    visits <- dplyr::left_join(move.out[[.x]][, list(target_SEUS = unique(seus), reach_SEUS = as.numeric(sum(grepl(pattern = "SEUS", x = region)) > 0)), whale], move.out[[.x]][, list(target_GSL = unique(gsl), reach_GSL = as.numeric(sum(grepl(pattern = "GSL", x = region)) > 0)), whale], by = "whale") |> dplyr::select(-whale) |> dplyr::mutate_if(is.numeric, list(~factor(., levels = c("0", "1"))))
   
    gsl <- visits |> janitor::tabyl(target_GSL, reach_GSL) |> 
      dplyr::rename(GSL = target_GSL, reach_not = '0', reach = '1') |> 
      dplyr::mutate(cohort = .x)
    
  }) |> do.call(what = rbind) |> 
    dplyr::relocate(cohort, .before = "GSL")
  
  journey.gsl.all <- data.table::data.table(journey.gsl)
  journey.gsl.all <- journey.gsl.all[,.(reach_not = sum(reach_not), reach = sum(reach)), GSL]
  journey.gsl.all <- format_dt(journey.gsl.all, direction = "row")
  journey.gsl[, 2:4] <- format_dt(data.frame(journey.gsl[, 2:4]), direction = "row")
  
  cat("\n+++++++++++ Migratory destinations (by cohort) +++++++++++")
  print(knitr::kable(journey.seus, format = "simple"))
  print(knitr::kable(journey.gsl, format = "simple"))
  
  cat("\n+++++++++++ Migratory destinations (all individuals) +++++++++++")
  print(knitr::kable(journey.seus.all, format = "simple"))
  print(knitr::kable(journey.gsl.all, format = "simple"))
  
  #'---------------------------------------
  # Summary statistics for daily step lengths
  #'---------------------------------------
  step.df <- purrr::map(.x = cohort.ab, .f = ~
    move.out[[.x]][, list(
      min = round(min(d_travel), 1),
      mean = round(mean(d_travel), 1),
      max = round(max(d_travel), 1),
      sd = round(sd(d_travel), 1)
    ), ] |>
      dplyr::mutate(cohort = .x)) |>
    do.call(what = rbind) |>
    dplyr::relocate(cohort, .before = min) |>
    dplyr::mutate(step = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |>
    dplyr::select(-min, -max, -mean, -sd)

  #'---------------------------------------
  # Summary statistics for overall migration distance
  #'---------------------------------------
  migr.df <- purrr::map(.x = cohort.ab, .f = ~
    move.out[[.x]][, list(
      min = format(round(min(d_migr), 0), big.mark = ","), 
      mean = format(round(mean(d_migr), 0), big.mark = ","), 
      max = format(round(max(d_migr), 0), big.mark = ","), 
      sd = format(round(sd(d_migr), 0), big.mark = ",")), ] |>
      dplyr::mutate(cohort = .x)) |>
    do.call(what = rbind) |>
    dplyr::relocate(cohort, .before = min) |>
    dplyr::mutate(migration = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |>
    dplyr::select(-min, -max, -mean, -sd)
 
  cat("\n+++++++++++ Step lengths and migration distances +++++++++++")
  print(knitr::kable(dplyr::left_join(step.df, migr.df, by = "cohort"), format = "simple"))

  cat("\n")
  
  }
  
  # Habitat use ---------------------------------------------------------------
  
  if("habitatuse" %in% what | what == "all"){
  
  cat("=============================================================\n")
  cat("HABITAT USE\n")
  cat("=============================================================\n\n")
  
  # Number of individuals visiting each region 
  
  indiv_reg.df <- purrr::map(.x = cohort.ab, .f = ~{
    move.out[[.x]] |> 
    dplyr::group_by(region) |> 
    dplyr::summarise(n = dplyr::n_distinct(whale)) |> 
    dplyr::left_join(x = tibble::tibble(region = sort(regions$region)), by = "region") |> 
    dplyr::mutate(n = dplyr::coalesce(n, 0), cohort = .x)
    }) |> do.call(what = rbind) |> 
    dplyr::mutate(n = glue::glue("{format(round(100*n/n.ind, 1), nsmall = 1)}% ({n})")) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n)
    
  days.df <- purrr::map(.x = cohort.ab, .f = ~{
    move.out[[.x]][, .N, list(region, whale)] |> 
      (\(.) .[, list(
                min = min(N),
                mean = round(mean(N), 1),
                max = max(N),
                sd = round(sd(N), 1)
              ), region])() |>  
      dplyr::left_join(x = tibble::tibble(region = sort(regions$region)), by = "region") |> 
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~dplyr::coalesce(.x, 0))) |> 
      dplyr::mutate(cohort = .x)
    }) |> do.call(what = rbind) |>  
      dplyr::mutate(days = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |>
      dplyr::select(-min, -max, -mean, -sd) |> 
      tidyr::pivot_wider(names_from = cohort, values_from = days)

  # Number of regions per individual
  
  reg.per.ind.df <- purrr::map(.x = cohort.ab, .f = ~{
    dt1 <- move.out[[.x]][, .(nreg = uniqueN(region)), whale]
    dt1[,.(n = .N, cohort = .x), nreg]  
  }) |> data.table::rbindlist() |> 
    dplyr::rename(No.regions = nreg) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    dplyr::arrange(No.regions) |> 
    format_dt()
  
  # reg.per.ind.df <- purrr::map(.x = cohort.ab, .f = ~{
  #   move.out[[.x]] |> 
  #     dplyr::group_by(whale) |> 
  #     dplyr::summarise(nreg = dplyr::n_distinct(region)) |> 
  #     janitor::tabyl(nreg) |> 
  #     dplyr::mutate(cohort = .x)
  #   }) |> 
  #   do.call(what = rbind) |> 
  #   dplyr::select(-percent) |>
  #   dplyr::rename(No.regions = nreg) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
  #   format_dt()
  
  cat("+++++++++++ Number of animals visiting each region (N = ", n.ind, ") +++++++++++", sep = "")
  print(knitr::kable(indiv_reg.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Days spent in each region +++++++++++")
  print(knitr::kable(days.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Total number of regions visited +++++++++++")
  print(knitr::kable(reg.per.ind.df, format = "simple"))
 
  }
  
  # Behavior ---------------------------------------------------------------
  
  if("behavior" %in% what | what == "all"){
  
  cat("\n")
  cat("=============================================================\n")
  cat("ACTIVITY BUDGETS\n")
  cat("=============================================================\n")
 
  tot.hrs <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ mean(sim[[.x]][
      alive == 1 & day > 0,
      list(total = sum(t_travel, t_feed, t_rest, t_nurse)), list(whale, day)]$total)) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name, total_hrs = value)

  if(sum(do.call(c, purrr::map(.x = sim, .f = ~c(unlist(.x[, .(t_travel, t_rest, t_nurse, t_feed)]))))>24) > 0){
    warning("Inconsistent time allocation in activity budgets")
  } 
  
   activity.df <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][alive == 1, list(cohort = .x, region, alive, t_travel, t_feed, t_rest, t_nurse, feed)]) |> 
     do.call(what = rbind)
      
   activ.summary.other <- activity.df[,list(
     `travel (hrs)` = paste0(round(mean(t_travel), 2), " (±", round(sd(t_travel), 2), ")"),
     `rest (hrs)` = paste0(round(mean(t_rest), 2), " (±", round(sd(t_rest), 2), ")"),
     `nurse (hrs)` = paste0(round(mean(t_nurse), 2), " (±", round(sd(t_nurse), 2), ")")
   ), list(cohort, region)]
   
   activ.summary.feed <- activity.df[feed == 1,
     list(`feed (hrs)` = paste0(round(mean(t_feed), 2), " (±", round(sd(t_feed), 2), ")")), list(cohort, region)]
   
   dt <- expand.grid(cohort.ab, sort(regions$region)) |> 
     dplyr::rename(cohort = Var1, region = Var2) |> 
     data.table::as.data.table()
   
   activ.summary <- purrr::reduce(list(
     activ.summary.other,
     activ.summary.feed, dt
   ), dplyr::left_join, by = c("cohort", "region")) |>
     dplyr::mutate_at(vars(matches("hrs")), ~ dplyr::coalesce(.x, "0 (± 0)"))
   
     # dplyr::left_join(activ.summary.other, activ.summary.feed, by = c("cohort", "region")) |> 
     # dplyr::left_join(x = dt, by = c("cohort", "region")) |>

  purrr::walk2(.x = cohort.ab, .y = cohort.names, .f = ~{
   cat("\n+++++++++++ ", .y, " +++++++++++")
   print(knitr::kable(activ.summary[cohort == .x, .SD, .SDcols = !c("cohort")], format = "simple"))
 })
 
  activ.feed <- activity.df |>
    dplyr::filter(feed == 1) |> 
    dplyr::select(-feed) |> 
    tidyr::pivot_longer(!c("region", "cohort"), names_to = "behav", values_to = "hrs") |> 
    dplyr::mutate(behav = gsub("t_", "", behav)) |> 
    dplyr::right_join(y = tibble::as_tibble(expand.grid(cohort = cohort.ab, region = regions$region, behav = c("travel", "feed", "rest", "nurse"))), by = c("region","cohort", "behav")) |> 
    dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |> 
    dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region)) |> 
    dplyr::mutate(behav = stringr::str_to_upper(behav)) |> 
    dplyr::filter(behav == "FEED")
  
  activ.other <- activity.df |> 
    dplyr::select(-feed) |> 
    tidyr::pivot_longer(!c("region", "cohort"), names_to = "behav", values_to = "hrs") |> 
    dplyr::mutate(behav = gsub("t_", "", behav)) |> 
    dplyr::right_join(y = tibble::as_tibble(expand.grid(cohort = cohort.ab, region = regions$region, behav = c("travel", "feed", "rest", "nurse"))), by = c("region","cohort", "behav")) |> 
    dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |> 
    dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region)) |> 
    dplyr::mutate(behav = stringr::str_to_upper(behav)) |> 
    dplyr::filter(!behav == "FEED")
  
  activ.out <- rbind(activ.feed, activ.other)
  
  if(plot){
  for(k in cohort.ab){
  tmp <- activ.out |> dplyr::mutate(cohort = k)
  p <- ggplot2::ggplot(tmp, aes(x = factor(region), y = hrs)) + 
    ggplot2::geom_boxplot(fill = "darkgrey") +
    ggplot2::facet_wrap(. ~ factor(behav), scales = 'free', ncol = 2) +
    ylab("Time spent (hrs)") +
    xlab("") + theme_narw() +
    labs(title = cohort.names[which(cohort.ab==k)]) +
    scale_y_continuous(breaks = seq(0,24,4), limits = c(0,24))
  suppressWarnings(print(p))
  } 
  }
  
  }
 

  # Stressors ---------------------------------------------------------------
  
  if(any(c("stressors", "strike", "gear", "noise", "other") %in% what) | what == "all"){
    
    cat("\n")
    cat("=============================================================\n")
    cat("STRESSORS\n")
    cat("=============================================================\n")

    ## Entanglements ---------------------------------------------------------------
    
    if("gear" %in% what){
    
    #'---------------------------------------
    # Entanglement rates
    #'---------------------------------------
    
    entgl.rate <- purrr::set_names(cohort.ab) |>
      purrr::map(.f = ~ sim[[.x]][
        day > 0,
        list(entangled = max(is_entgl),
               n_events = length(unique(entgl_d[entgl_d>0]))), 
        whale]) |>
      tibble::enframe() |>
      tidyr::unnest(cols = c(value)) |>
      dplyr::rename(cohort = name) |> 
      data.table::as.data.table()
    
    if(5 %in% cohortID){
    entgl.rate.calf <- 
      sim[[cohorts[id==5,abb]]][day > 0,
        list(entangled = max(is_entgl_calf),
             n_events = length(unique(entgl_d_calf[entgl_d_calf>0])),
             born = max(born)), 
        whale] |> 
      dplyr::mutate(cohort = "c(m,f)") |> 
      dplyr::relocate(cohort, .before = "whale") |> 
      dplyr::filter(born == 1) |> 
      dplyr::select(-born)
    
    entgl.rate <- data.table::rbindlist(list(entgl.rate, entgl.rate.calf))
    
    }
  
    entgl.overall <- suppressWarnings(entgl.rate |> 
      janitor::tabyl(var1 = entangled)) |> 
      dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
      dplyr::select(entangled, rate) |> 
      dplyr::mutate(entangled = ifelse(entangled == 0, "no", "yes"))
    
    entgl.overall.calf <- suppressWarnings(entgl.rate.calf |> 
      janitor::tabyl(var1 = entangled)) |> 
      dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
      dplyr::select(entangled, rate) |> 
      dplyr::mutate(entangled = ifelse(entangled == 0, "no", "yes"))
    
    #'---------------------------------------
    # Entanglement rates by cohort
    #'---------------------------------------
    
    entgl.bycohort <- suppressWarnings(entgl.rate |> 
      janitor::tabyl(var1 = cohort, var2 = entangled)) |> 
      format_dt(direction = "row") |> 
      dplyr::rename_at(dplyr::vars(starts_with("0")), ~"not entangled") |> 
      dplyr::rename_at(dplyr::vars(starts_with("1")), ~"entangled") 

    #'---------------------------------------
    # Rope position
    #'---------------------------------------
    
    entgl.head <- suppressWarnings(purrr::set_names(cohort.ab) |>
      purrr::map(.f = ~ sim[[.x]][
        day > 0 & is_entgl == 1,
        list(anterior = max(entgl_head)), whale]) |>
      tibble::enframe() |>
      tidyr::unnest(cols = c(value)) |>
      dplyr::rename(cohort = name) |> 
      data.table::as.data.table())
    
    if(5 %in% cohortID){
      
      entgl.head.calf <- suppressWarnings(sim[[cohorts[id==5,abb]]][
          day > 0 & is_entgl_calf == 1,
          list(anterior = max(entgl_head_calf), born = max(born)), whale] |> 
        dplyr::mutate(cohort = "c(m,f)") |> 
        dplyr::relocate(cohort, .before = "whale") |> 
        dplyr::filter(born == 1) |> 
        dplyr::select(-born))
      
      entgl.head <- data.table::rbindlist(list(entgl.head, entgl.head.calf))
    }
    
    head.overall <- suppressWarnings(entgl.head |> 
      janitor::tabyl(var1 = anterior)) |> 
      dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
      dplyr::rename(position = anterior) |> 
      dplyr::select(position, rate) |> 
      dplyr::mutate(position = ifelse(position == 0, "body", "head"))
    
    head.overall.calf <- suppressWarnings(entgl.head.calf |>
      janitor::tabyl(var1 = anterior)) |> 
      dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
      dplyr::rename(position = anterior) |> 
      dplyr::select(position, rate) |> 
      dplyr::mutate(position = ifelse(position == 0, "body", "head"))

    head.bycohort <- suppressMessages(entgl.head |>
                                        janitor::tabyl(var1 = cohort, var2 = anterior)) |>
                                        format_dt(direction = "row") |>
                                        dplyr::rename_at(dplyr::vars(starts_with("0")), ~"body") |>
                                        dplyr::rename_at(dplyr::vars(starts_with("1")), ~"head")
    
    #'---------------------------------------
    # Entanglement events
    #'---------------------------------------
    
    entgl.events <- entgl.rate[n_events > 0, list(
      `No. entanglements` =
        paste0(
          round(mean(n_events), 2),
          " (±",
          round(sd(n_events), 2),
          ") [", min(n_events), "–", max(n_events), "]"
        )
    ), cohort]
    
    entgl.df <- 
      purrr::reduce(list(entgl.bycohort, 
                         data.table::data.table(cohort = entgl.bycohort$cohort, ` ` = "|"),
                         head.bycohort, 
                         data.table::data.table(cohort = entgl.bycohort$cohort, `  ` = "|"),
                         entgl.events), dplyr::left_join, by = 'cohort')
    
    #'---------------------------------------
    # Entanglement severity
    #'---------------------------------------
    
    entgl.sev <- purrr::set_names(cohort.ab) |>
      purrr::map(.f = ~ sim[[.x]][
        day > 0 & is_entgl == 1,
        list(entgl = sum(is_entgl),
             minor = sum(entgl_sev == 0),
             moderate = sum(entgl_sev == 1),
             severe = sum(entgl_sev == 2)), whale]) |>
      tibble::enframe() |>
      tidyr::unnest(cols = c(value)) |>
      dplyr::rename(cohort = name) |> data.table::as.data.table()
    
    if(5 %in% cohortID){
      
      entgl.sev.calf <- suppressWarnings(sim[[cohorts[id==5,abb]]][
        day > 0 & is_entgl_calf == 1,
        list(entgl = sum(is_entgl_calf),
             minor = sum(entgl_sev_calf == 0),
             moderate = sum(entgl_sev_calf == 1),
             severe = sum(entgl_sev_calf == 2),
             born = max(born)), whale] |>
        dplyr::mutate(cohort = "c(m,f)") |> 
        dplyr::relocate(cohort, .before = "whale") |> 
        dplyr::filter(born == 1) |> 
        dplyr::select(-born))
      
      entgl.sev <- data.table::rbindlist(list(entgl.sev, entgl.sev.calf))
    }
   
    #'---------------------------------------
    # Entanglement duration
    #'---------------------------------------

    entgl.durations <-
      suppressWarnings(entgl.sev[
        , list(
          `minor (days)` = paste0(
            ifelse(all(minor == 0), 0, round(mean(minor[minor>0]), 0)),
            " (±", 
            ifelse(all(minor == 0) | is.na(sd(minor[minor>0])), 0, round(sd(minor[minor>0]), 0)),
            ") [", 
            ifelse(all(minor == 0), 0, min(minor[minor>0])),
            "–", 
            ifelse(all(minor==0), 0, max(minor[minor>0])), "]"
          ),
          `moderate (days)` = paste0(
            ifelse(all(moderate == 0), 0, round(mean(moderate[moderate>0]), 0)),
            " (±", 
            ifelse(all(moderate == 0) | is.na(sd(moderate[moderate>0])), 0, round(sd(moderate[moderate>0]), 0)),
            ") [", 
            ifelse(all(moderate == 0), 0, min(moderate[moderate>0])),
            "–", 
            ifelse(all(moderate==0), 0, max(moderate[moderate>0])), "]"
          ),
          `severe (days)` = paste0(
            ifelse(all(severe == 0), 0, round(mean(severe[severe>0]), 0)),
            " (±", 
            ifelse(all(severe == 0) | is.na(sd(severe[severe>0])), 0, round(sd(severe[severe>0]), 0)),
            ") [", 
            ifelse(all(severe == 0), 0, min(severe[severe>0])),
            "–", 
            ifelse(all(severe==0), 0, max(severe[severe>0])), "]"
          )
        ),
        cohort
      ])
    
    if(plot){
      
      if(nrow(entgl.sev) == 0){
        warning("No entanglement events recorded.")
      } else {
      
      events.df <- 
        tidyr::pivot_longer(data = entgl.sev[, c("cohort", "minor", "moderate", "severe")],
                            !cohort, names_to = "severity", values_to = "days") |> 
        dplyr::mutate(title = "Entanglement severity") |> 
        dplyr::mutate(severity = stringr::str_to_sentence(severity)) |> 
        dplyr::left_join(y = cohorts[, c("name", "abb")], by = c("cohort" = "abb"))
      
      print(ggplot2::ggplot(events.df, aes(x = factor(severity), y = days)) + 
        ggplot2::geom_boxplot(fill = "darkgrey") +
        ggplot2::facet_grid(name ~ factor(title)) +
        ylab("Time spent (hrs)") +
        xlab("") + 
        ylab("Duration (days)") + 
        theme_narw() + 
        ggplot2::scale_y_continuous(breaks = pretty(range(events.df$days), n = 10)))
      
      }
      
    }
    
  gear_risk <- purrr::set_names(cohort.ab) |>
      purrr::map(.f = ~ sim[[.x]][
        alive == 1 & day > 0,
        list(mean = round(mean(gear_risk), 3),
             sd = round(sd(gear_risk), 3),
             min = round(min(gear_risk), 3),
             max = round(max(gear_risk), 3))] |> 
          dplyr::mutate(`p(entangled)` = paste0(format(mean, big.mark = ","),
                                " (±", format(sd, big.mark = ","), ") [",
                                format(min, big.mark = ","), "–",
                                format(max, big.mark = ","), "]"))) |>
      tibble::enframe() |>
      tidyr::unnest(cols = c(value)) |>
      dplyr::rename(cohort = name) |> 
    dplyr::select(cohort, "p(entangled)")
  
  
  cat("\n+++++++++++ Entanglements +++++++++++")
  
  # if(length(cohortID) > 1 | (length(cohortID) == 1 & all(cohortID == 5))){
  #   print(knitr::kable(entgl.bycohort, format = "simple"))
  #   if(nrow(head.bycohort) > 0) print(knitr::kable(head.bycohort, format = "simple"))
  # } else {
  #   print(knitr::kable(entgl.overall, format = "simple"))
  #   print(knitr::kable(head.overall, format = "simple"))
  # }
  
  print(knitr::kable(entgl.overall, format = "simple"))
  print(knitr::kable(head.overall, format = "simple"))
  print(knitr::kable(entgl.bycohort, format = "simple"))
  print(knitr::kable(purrr::reduce(list(entgl.events, gear_risk), dplyr::left_join, by = 'cohort'), format = "simple"))
  print(knitr::kable(entgl.durations, format = "simple"))
  
    }
  
  ## Vessel strikes ---------------------------------------------------------------
  
  if("strike" %in% what){
    
  strike.rate <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      day > 0 & alive == 0,
      list(strike = max(strike)), whale]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    data.table::as.data.table()
  
  if(5 %in% cohortID){
    
  strike.rate.calf <- 
    suppressWarnings(sim[["ad(f,l)"]][day > 0, list(strike = max(strike_calf), born = max(born)), whale] |>
    dplyr::mutate(cohort = "c(m,f)") |> 
    dplyr::relocate(cohort, .before = "whale") |> 
    dplyr::filter(born == 1) |> 
    dplyr::select(-born))
  
  strike.rate <- data.table::rbindlist(list(strike.rate, strike.rate.calf))
  
  }

  strike.overall <- suppressWarnings(
    strike.rate |> 
    janitor::tabyl(var1 = strike)) |> 
    dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
    dplyr::select(strike, rate) |> 
    dplyr::mutate(strike = ifelse(strike == 0, "no", "yes"))

  
  # strike.overall.calf <- suppressWarnings(strike.rate.calf |> 
  #   janitor::tabyl(var1 = strike)) |> 
  #   dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
  #   dplyr::select(strike, rate) |> 
  #   dplyr::mutate(strike = ifelse(strike == 0, "no", "yes"))
  
  # By cohort
  strike.bycohort <- suppressWarnings(suppressMessages(strike.rate |> 
    janitor::tabyl(var1 = cohort, var2 = strike)) |> 
    format_dt(direction = "row") |> 
    dplyr::rename_at(dplyr::vars(starts_with("0")), ~"not struck") |> 
    dplyr::rename_at(dplyr::vars(starts_with("1")), ~"struck"))
  
  strike_risk <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      alive == 1 & day > 0,
      list(mean = round(mean(strike_risk), 3),
           sd = round(sd(strike_risk), 3),
           min = round(min(strike_risk), 3),
           max = round(max(strike_risk), 3))]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    dplyr::mutate(`p(strike)` = paste0(format(mean, big.mark = ","),
                                " (±", format(sd, big.mark = ","), ") [",
                                format(min, big.mark = ","), "–",
                                format(max, big.mark = ","), "]")) |> 
    dplyr::select(cohort, "p(strike)")
  
  cat("\n\n+++++++++++ Vessel strikes +++++++++++")
  
  print(knitr::kable(strike.overall, format = "simple"))
  print(knitr::kable(strike.bycohort, format = "simple"))
  print(knitr::kable(strike_risk, format = "simple"))
 
  }
  
  ## Mortality from other sources --------------------------------------------

 if("other" %in% what){
    
  mortality <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      day > 0,
      list(dead = max(died)), whale]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    data.table::as.data.table()
  
  if(5 %in% cohortID){
    
    mortality.calf <- 
      suppressWarnings(sim[["ad(f,l)"]][day > 0, list(dead = max(died_calf), born = max(born)), whale] |>
                         dplyr::mutate(cohort = "c(m,f)") |> 
                         dplyr::relocate(cohort, .before = "whale") |> 
                         dplyr::filter(born == 1) |> 
                         dplyr::select(-born))
    
    mortality.rate <- data.table::rbindlist(list(mortality, mortality.calf))
    
  }
  
  mortality.overall <- suppressWarnings(mortality.rate |> janitor::tabyl(var1 = dead)) |> 
    dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
    dplyr::select(dead, rate) |>
    dplyr::rename(mortality = dead) |> 
    dplyr::mutate(mortality = ifelse(mortality == 0, "alive", "dead"))
  
  mortality.overall.calf <- suppressWarnings(mortality.calf |> janitor::tabyl(var1 = dead)) |> 
    dplyr::mutate(rate = paste0(round(100 * percent, 1), "% (", n, ")")) |> 
    dplyr::select(dead, rate) |> 
    dplyr::rename(mortality = dead) |> 
    dplyr::mutate(mortality = ifelse(mortality == 0, "alive", "dead"))
  
  # By cohort
  mortality.bycohort <- suppressWarnings(suppressMessages(mortality.rate |> 
                                                         janitor::tabyl(var1 = cohort, var2 = dead)) |> 
                                        format_dt(direction = "row") |> 
                                        dplyr::rename_at(dplyr::vars(starts_with("0")), ~"alive") |> 
                                        dplyr::rename_at(dplyr::vars(starts_with("1")), ~"dead"))
  
  cat("\n\n+++++++++++ Other sources of mortality +++++++++++")
  
  print(knitr::kable(mortality.overall, format = "simple"))
  print(knitr::kable(mortality.bycohort, format = "simple"))
 
 }
  
  ## Noise ---------------------------------------------------------------
  
    if("noise" %in% what){
    
  noise_lvl <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      alive == 1 & day > 0,
      list(mean = round(mean(noise_lvl), 3),
           sd = round(sd(noise_lvl), 3),
           min = round(min(noise_lvl), 3),
           max = round(max(noise_lvl), 3))]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    dplyr::mutate(`noise level` = paste0(format(mean, big.mark = ","),
                                " (±", format(sd, big.mark = ","), ") [",
                                format(min, big.mark = ","), "–",
                                format(max, big.mark = ","), "]")) |> 
    dplyr::select(cohort, `noise level`)
  
  dB_thresh <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      alive == 1 & day > 0,
      list(mean = round(mean(dB_thresh), 1),
           sd = round(sd(dB_thresh), 1),
           min = round(min(dB_thresh), 1),
           max = round(max(dB_thresh), 1))]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    dplyr::mutate(`response threshold` = paste0(format(mean, big.mark = ","),
                                         " (±", format(sd, big.mark = ","), ") [",
                                         format(min, big.mark = ","), "–",
                                         format(max, big.mark = ","), "]")) |> 
    dplyr::select(cohort, `response threshold`)
     
  noise.rate <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      day > 0 & alive == 1, list(noise_resp = sum(noise_resp)), whale]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    data.table::as.data.table()
  
  beh.resp <- noise.rate[, list(mean = round(mean(noise_resp), 1),
                    sd = round(sd(noise_resp), 1),
                    min = round(min(noise_resp), 1),
                    max = round(max(noise_resp), 1)), cohort] |> 
    dplyr::mutate(`Responses (No. days)` = paste0(format(mean, big.mark = ","),
                                                " (±", format(sd, big.mark = ","), ") [",
                                                format(min, big.mark = ","), "–",
                                                format(max, big.mark = ","), "]")) |> 
    dplyr::select(cohort, `Responses (No. days)`)
  
  cat("\n\n+++++++++++ Pile-driving noise +++++++++++")
  print(knitr::kable(purrr::reduce(list(noise_lvl, dB_thresh, beh.resp), dplyr::left_join, by = 'cohort'), format = "simple"))
    }
  }
  
  # Energy budgets ---------------------------------------------------------------
  
  if("energy" %in% what | what == "all"){
    
    cat("\n")
    cat("=============================================================\n")
    cat("ENERGY BUDGETS (MJ per day) \n")
    cat("=============================================================")
    
    ## Intake ---------------------------------------------------------------
    
  Ein <- purrr::set_names(cohort.ab) |>
      purrr::map(.f = ~ sim[[.x]][
        alive == 1 & day > 0,
        list(mean = round(mean(E_in), 1),
             sd = round(sd(E_in), 1),
             min = round(min(E_in), 1),
             max = round(max(E_in), 1))]) |>
      tibble::enframe() |>
      tidyr::unnest(cols = c(value)) |>
      dplyr::rename(cohort = name) |> 
      dplyr::mutate(E_in = paste0(format(mean, big.mark = ","),
                                  " (±", format(sd, big.mark = ","), ") [",
                                  format(min, big.mark = ","), "–",
                                  format(max, big.mark = ","), "]"))
   
  if(5 %in% cohortID){
    
    Ein.calf <- suppressWarnings(sim[["ad(f,l)"]][
        alive_calf == 1 & day > 0,
        list(mean = round(mean(E_in_calf), 1),
             sd = round(sd(E_in_calf), 1),
             min = round(min(E_in_calf), 1),
             max = round(max(E_in_calf), 1))] |>
      dplyr::mutate(cohort = "c(m,f)") |> 
      dplyr::relocate(cohort, .before = "mean") |> 
      dplyr::mutate(E_in = paste0(format(mean, big.mark = ","),
                                  " (±", format(sd, big.mark = ","), ") [",
                                  format(min, big.mark = ","), "–",
                                  format(max, big.mark = ","), "]")))
    
    Ein <- data.table::rbindlist(list(Ein, Ein.calf))
    
  }
  
    ## Expenditure ---------------------------------------------------------------
    
  Eout <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      alive == 1 & day > 0,
      list(mean = round(mean(E_out), 1),
           sd = round(sd(E_out), 1),
           min = round(min(E_out), 1),
           max = round(max(E_out), 1))]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
    dplyr::mutate(E_out = paste0(format(mean, big.mark = ","),
                                " (±", format(sd, big.mark = ","), ") [",
                                format(min, big.mark = ","), "–",
                                format(max, big.mark = ","), "]"))
  
  if(5 %in% cohortID){
    
    Eout.calf <- suppressWarnings(sim[["ad(f,l)"]][
      alive_calf == 1 & day > 0,
      list(mean = round(mean(E_out_calf), 1),
           sd = round(sd(E_out_calf), 1),
           min = round(min(E_out_calf), 1),
           max = round(max(E_out_calf), 1))] |>
      dplyr::mutate(cohort = "c(m,f)") |> 
      dplyr::relocate(cohort, .before = "mean") |> 
      dplyr::mutate(E_out = paste0(format(mean, big.mark = ","),
                                  " (±", format(sd, big.mark = ","), ") [",
                                  format(min, big.mark = ","), "–",
                                  format(max, big.mark = ","), "]")))
    
    Eout <- data.table::rbindlist(list(Eout, Eout.calf))
    
  }
  
    ## Balance ---------------------------------------------------------------
    
  Ebalance <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][
      alive == 1 & day > 0,
      list(negative = sum(E_tot < 0),
           positive = sum(E_tot > 0)), 
           whale]) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name) |> 
      data.table::as.data.table()
    
    if(5 %in% cohortID){
      
    Ebalance.calf <- suppressWarnings(sim[["ad(f,l)"]][
        alive_calf == 1 & day > 0,
        list(negative = sum(E_tot_calf < 0),
             positive = sum(E_tot_calf > 0)), whale] |>
      dplyr::mutate(cohort = "c(m,f)") |> 
      dplyr::relocate(cohort, .before = "whale"))
    
    Ebalance <- data.table::rbindlist(list(Ebalance, Ebalance.calf))
    
    }
  
  Ebalance <- janitor::adorn_percentages(Ebalance[, c(1,3,4)], "row") |> 
    dplyr::rename(balance_neg = negative, balance_pos = positive) |> 
    dplyr::select(balance_neg, balance_pos) |> 
    cbind(Ebalance)

  Ebalance <- Ebalance[, list(
    negative = paste0(
      100 * round(mean(balance_neg), 3), "% (±",
      100 * round(sd(balance_neg), 3),
      ") [",
      100 * round(min(balance_neg), 3),
      "–",
      100 * round(max(balance_neg), 3),
      "]"
    ),
    positive = paste0(
      100 * round(mean(balance_pos), 3), "% (±",
      100 * round(sd(balance_pos), 3),
      ") [",
      100 * round(min(balance_pos), 3),
      "–",
      100 * round(max(balance_pos), 3),
      "]"
    )
  ), cohort]
  
  print(knitr::kable(purrr::reduce(list(Ein[, c("cohort", "E_in")],
                                        Eout[, c("cohort", "E_out")]), dplyr::left_join, by = 'cohort') |> 
                       dplyr::rename(Energy_intake = E_in,
                                     Energy_expenditure = E_out), format = "simple"))
  
  print(knitr::kable(Ebalance |> dplyr::rename(Deficit = negative, 
                                     Surplus = positive), format = "simple"))
    
  }
   
}
