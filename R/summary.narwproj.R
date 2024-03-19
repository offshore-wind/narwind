#' Population projection summary
#'
#' Prints a text-based summary of outputs from the stochastic population model to the R console. 
#'
#' @param obj An object of class \code{narwproj}.
#'
#' @import data.table
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' m <- augment(m)
#' prj <- predict(m)
#' summary(prj)
#' }

summary.narwproj <- function(obj){
  
  options(pillar.sigfig = 7)
  
 # Preamble ---------------------------------------------------------------
  
  if(!inherits(obj, "narwproj")) stop("Object must be of class <narwproj>")
  
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n")
  cat("\n")
  cat("     NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("          *** POPULATION MODEL SUMMARY ***\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n\n")

  # Extract projection parameters
  current.yr <- obj$param$start.year
  yrs <- obj$param$yrs
  burn <- obj$param$burn
  N <- obj$param$n
  N_0 <- obj$init$N_0
  
  # Extract results
  tot.df <- obj$proj$tbl
  births.per.female <- obj$dat$birth$perfemale
  tot.births <- obj$dat$birth$tot
  tot.deaths <- obj$dat$death
  inter.birth <- obj$dat$birth$inter
  time.resting <- obj$dat$rest
  nonreprod.females <- obj$dat$nonrepfem

  schedule <- data.table::data.table(phase = obj$param$phases[obj$param$schedule + 1]) |> 
    dplyr::mutate(year = 0:(obj$param$yrs+burn)) |> 
    dplyr::mutate(year = 2019 + year)
  
  cat("Replicates: N =", N, "\n")
  cat("Projection horizon:", yrs, "years\n\n")
  cat("Timeline:\n")
  proj_timeline(obj$param$schedule, burn)
  
  cat("=============================================================\n")
  cat("ABUNDANCE\n")
  cat("=============================================================\n\n")
  
  cat("Initial population size:\n")
  cat("N =", sum(N_0))
  init.pop <- tibble::tibble(cohort = names(N_0), N = N_0) |> 
    janitor::adorn_totals()
 
  print(knitr::kable(init.pop, format = "simple"))
  cat("\n")
   
  # Find 95% confidence intervals on final population size
  final.pop <- tot.df[year == current.yr + yrs, list(N = quantile(N, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)), cohort]
  final.pop[, value := rep(c("median", "lower", "upper"), length.out = nrow(final.pop))]
  final.pop <- final.pop |> tidyr::pivot_wider(names_from = value, values_from = N) |> data.table::data.table()
  final.pop[, N:= paste0(format(round(median,0), big.mark = ","), " (95% CI: ", format(round(lower,0), big.mark = ","), " – ", format(round(upper,0), big.mark = ","), ")")]
  
  cat("Final population size:\n")
  cat("N =", final.pop[cohort == "North Atlantic right whales", N])
  
  final.pop <- tibble::tibble(final.pop[, list(cohort, N)]) |> 
    dplyr::mutate(N = trimws(stringr::str_replace_all(N, "\\h+", " ")))
  
  print(knitr::kable(final.pop[!final.pop$cohort == "North Atlantic right whales", ], format = "simple"))
  
  cat("\n=============================================================\n")
  cat("HEALTH & MORTALITY\n")
  cat("=============================================================\n\n")
  
  # sampled.ind <- sample(seq_len(max(obj$dat$health$whale)), size = 500)
  
  # BC plots based on cohort each animal starts in
  init.cohort <- obj$dat$health[, list(init.cohort = cohort[1]),list(prj, whale)]

  bc.df <- dplyr::left_join(obj$dat$health[!is.na(bc) & bc > 0 & prj %in% sample(1:obj$param$n, size = 10)], init.cohort, by = c("whale", "prj")) |> 
    dplyr::left_join(y = obj$param$cohorts[id>0, list(id,name)], by = c("init.cohort" = "id")) |> 
    dplyr::mutate(name = ifelse(init.cohort > 0, name, ifelse(female == 0, "Calves (male)", "Calves (female)"))) |> 
    dplyr::mutate(name = factor(name, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)"
    )))
  
  bc.series <- ggplot2::ggplot(data = bc.df, aes(x = year, y = bc, group = interaction(whale, prj))) +
    ggplot2::geom_path(alpha = 0.15) +
    ggplot2::facet_wrap(vars(name), scales = "fixed") +
    scale_y_continuous(limits = c(0,0.6)) +
    theme_narw() + xlab("") + ylab("Body condition (%)")
  
  suppressWarnings(print(bc.series))
  
  # Extract health and survival data
  health.df <- obj$dat$health |> 
    dplyr::filter(!is.na(bc)) |> 
    dplyr::select(-prj, -year, -whale, -female) |> 
    tidyr::pivot_longer(!cohort, names_to = "param", values_to = "value") |> 
    dplyr::mutate(row = dplyr::row_number()) |> 
    data.table::as.data.table()
   
  health.df <- health.df[!row %in% health.df[param == "min_bc" & value == 0, row]]
  
  survival.df <- obj$dat$survival |> 
    dplyr::select(surv, cohort)
  
  # Calculate summary statistics
  health.tbl <- health.df[, list(min = min(value),
                               max = max(value),
                               mean = mean(value),
                               sd = sd(value)), list(cohort, param)]
  
  health.tbl[param == "bc", p := paste0(round(mean, 3), " ± ",
                                                   round(sd, 3), " [",
                                                   round(min, 3), "–",
                                                   round(max, 3), "]")]
  
  health.tbl[param == "min_bc", p := paste0(round(mean, 3), " ± ",
                                        round(sd, 3), " [",
                                        round(min, 3), "–",
                                        round(max, 3), "]")]
  
  health.tbl <- health.tbl |> 
    dplyr::select(cohort, param, p) |> 
    tidyr::pivot_wider(names_from = "param", values_from = "p") |> 
    tidyr::replace_na(list(min_bc = "-")) |> 
    dplyr::rename(`body condition` = bc, `min condition (gestation)` = min_bc) |> 
    dplyr::arrange(cohort) |> 
    dplyr::left_join(x = cohorts[, c("id", "abb")], by = c("id" = "cohort")) |> 
    dplyr::select(-id) |> 
    dplyr::rename(cohort = abb)

  survival.tbl <- survival.df[, list(min = min(surv),
                                     max = max(surv),
                                     mean = mean(surv),
                                     sd = sd(surv)), list(cohort)]
  
  survival.tbl[, p := paste0(round(mean, 3), " ± ",
                             round(sd, 3), " [",
                             round(min, 3), "–",
                             round(max, 3), "]")]
  
  survival.tbl <- survival.tbl |> 
    dplyr::select(cohort, p) |> 
    dplyr::rename(`survival` = p) |> 
    dplyr::arrange(cohort) |> 
    dplyr::left_join(x = cohorts[, c("id", "abb")], by = c("id" = "cohort")) |> 
    dplyr::select(-id) |> 
    dplyr::rename(cohort = abb)
  
  health.print <- trimws(knitr::kable(dplyr::left_join(survival.tbl, health.tbl, by = "cohort"),
                                      format = "simple"))     
  cat(health.print, sep = "\n")
  
  cat("\n=============================================================\n")
  cat("FECUNDITY\n")
  cat("=============================================================\n\n")
  
  # Combine data on births and deaths
  plot.df <- rbind(
    tot.deaths |> dplyr::mutate(param = "No. deaths") |> dplyr::rename(N = death),
    tot.births |> dplyr::mutate(param = "No. births") |> dplyr::rename(N = birth)
  )
  
  p <- ggplot2::ggplot(data = plot.df) +
    ggplot2::geom_line(aes(x = year, y = N, group = prj), col = "grey50", alpha = 0.5) +
    ggplot2::geom_smooth(aes(x = year, y = N), 
                         fill = "#0098A0", colour = "#0098A0", method = "loess", formula = 'y ~ x') +
    ggplot2::facet_wrap(vars(param), scales = "free_y") +
    theme_narw() +
    labs(x = "", y = "")
  
  calving.events <- tibble::tibble(`Calving events` = c(
    paste0(
      "Per year: N = ",
      round(mean(tot.births$birth, na.rm = TRUE), 2), " ± ",
      round(sd(tot.births$birth, na.rm = TRUE), 2), " [",
      min(tot.births$birth, na.rm = TRUE), "–",
      max(tot.births$birth, na.rm = TRUE), "]"
    ),
    paste0("Per female: N = ", paste0(
      round(mean(births.per.female$nbirths, na.rm = TRUE), 2), " ± ",
      round(sd(births.per.female$nbirths, na.rm = TRUE), 2), " [",
      min(births.per.female$nbirths, na.rm = TRUE), "–",
      max(births.per.female$nbirths, na.rm = TRUE), "]"
    ))
  ))
  
  resting.period <- tibble::tibble(`Resting phase` = c(
    paste0("t(rest): ", paste0(
      median(time.resting$t_rest, na.rm = TRUE), " ± ",
      round(sd(time.resting$t_rest, na.rm = TRUE), 0), " [",
      min(time.resting$t_rest, na.rm = TRUE), "–",
      max(time.resting$t_rest, na.rm = TRUE), "]"
    )),
    paste0("t(inter-birth): ", paste0(
      median(inter.birth$ibi, na.rm = TRUE), " ± ",
      round(sd(inter.birth$ibi, na.rm = TRUE), 0), " [",
      min(inter.birth$ibi, na.rm = TRUE), "–",
      max(inter.birth$ibi, na.rm = TRUE), "]"
    ))
  ))
      
  calve.print <- trimws(knitr::kable(calving.events, format = "simple"))   
  rest.print <- trimws(knitr::kable(resting.period, format = "simple"))
  cat(calve.print, sep = "\n")
  cat("\n")
  cat(rest.print, sep = "\n")
  cat("\n")
  
  nrf <- nonreprod.females |> 
    tibble::as_tibble() |> 
    dplyr::mutate(proj = dplyr::row_number()) |> 
    tidyr::pivot_longer(!proj, names_to = "year", values_to = "nonrep") |> 
    dplyr::mutate(year = as.numeric(gsub("yr ", "", year))) |> 
    data.table::as.data.table()
  
  nrf.tbl <- nrf[, list(mean = round(100*mean(nonrep), 2),
                        sd = round(100*sd(nonrep), 2),
                        min = round(100*min(nonrep), 2),
                        max = round(100*max(nonrep), 2)), ] |> 
    dplyr::mutate(`Non-reproductive females` = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |> 
    dplyr::select(-mean, -sd, -min, -max) |> 
    dplyr::pull(`Non-reproductive females`)
  
  cat("Non-reproductive females:\n", as.character(nrf.tbl), " %", sep ="")
  cat("\n\n")
  
  # Abortion rates
  abrt.proj <- dplyr::left_join(obj$dat$pregfem, obj$dat$abort, by = c("year", "prj")) |> 
    dplyr::filter(n_pregfem > 0) |> 
    dplyr::mutate(year = as.integer(year)) |> 
    dplyr::mutate(rate = n_abort / n_pregfem) |> 
    dplyr::filter(year >= current.yr) |> 
    dplyr::left_join(y = schedule, by = "year")
  
  abrt.summary <- abrt.proj[,list(mean = round(100*mean(rate),2),
                                  sd = round(100*sd(rate),2),
                                  min = round(100*min(rate),2),
                                  max = round(100*max(rate),2)), phase] |> 
    dplyr::mutate(`Projection (%)` = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |> 
    dplyr::select(-mean, -sd, -min, -max)
  
  abrt.sim <- tibble::enframe(obj$param$abort) |> 
    dplyr::mutate(value = 100*value) |> 
    dplyr::rename(phase = name, `Simulation (%)` = value)
  
  abrt.out <- dplyr::left_join(abrt.sim, abrt.summary, by = "phase")
  
  abrt.out$phase <- obj$param$phases
  abrt.out$phase <- sapply(X = abrt.out$phase, FUN = function(x) switch(x,
                       "base" = "Baseline",
                       "const" = "Construction",
                       "ops" = "Operation & maintenance"))
  abrt.out <- abrt.out |> dplyr::rename(Abortion = phase)
  abrt.print <- trimws(knitr::kable(abrt.out, format = "simple"))

  cat(abrt.print, sep = "\n")
  cat("\n")
  
  cat("\n=============================================================\n")
  cat("POPULATION VIABILITY \n")
  cat("=============================================================")
  
  pm <- ggplot2::ggplot(data = obj$proj$minpop) +
    ggplot2::geom_line(aes(x = year, y = minpop)) +
    theme_narw() +
    labs(x = "", y = "Expected minimum population size") +
    scale_x_continuous(breaks = pretty(obj$proj$minpop$year, n = 5))
  
  pq <- ggplot2::ggplot(data = obj$proj$quasi) +
    ggplot2::geom_line(aes(x = year, y = quasi)) +
    theme_narw() +
    labs(x = "", y = "Quasi-extinction risk") +
    scale_x_continuous(breaks = pretty(obj$proj$quasi$year, n = 5))
  

  
  pq.tbl <- obj$proj$quasi[year %in% c(obj$param$start.year + pretty(1:obj$param$yrs, n = 10)),] |> 
    dplyr::rename(`p(quasi-extinct)` = quasi)
  pq.tbl <- trimws(knitr::kable(pq.tbl, format = "simple"))
  print(pq.tbl)
  
  cat("\n")
  if(obj$param$yrs >=100){
  p.decline <- obj$proj$iucn.p |> 
    tibble::enframe() |> 
    dplyr::rename(`pop decline` = name, `prob` = value)
  p.decline <- trimws(knitr::kable(p.decline, format = "simple"))
  cat(p.decline, sep = "\n")
  }
  
  # Plots
  print(p / pm + pq)
  
  }
