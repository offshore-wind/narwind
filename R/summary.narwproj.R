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
  current.yr <- obj$param$current.yr
  yrs <- obj$param$yrs
  N <- obj$param$n
  N_0 <- obj$init$N_0
  
  # Extract results
  tot.df <- obj$proj$tbl
  births.per.female <- obj$dat$birth$perfemale
  tot.births <- obj$dat$birth$tot
  inter.birth <- obj$dat$birth$inter
  time.resting <- obj$dat$rest
  nonreprod.females <- obj$dat$nonrepfem

  cat("Replicates: N =", N, "\n")
  cat("Projection horizon:", yrs, "years\n\n")
  cat("Timeline:\n")
  proj_timeline(obj$param$schedule)
  
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
  
  health.df <- obj$dat$health |> 
    tidyr::pivot_longer(!cohort, names_to = "param", values_to = "value") |> 
    dplyr::mutate(row = dplyr::row_number()) |> 
    data.table::as.data.table()
  
  health.df <- health.df[!row %in% health.df[param == "min_bc" & value == 0, row]]
  
  health.tbl <- health.df[, list(min = min(value),
                               max = max(value),
                               mean = mean(value),
                               sd = sd(value)), list(cohort, param)]
  
  health.tbl[param == "p_surv", p := paste0(round(mean, 3), " ± ",
                                                 round(sd, 3), " [",
                                                 round(min, 2), "–",
                                                 round(max, 2), "]")]
  
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
    dplyr::rename(`p(survival)` = p_surv, `body condition` = bc, `min condition (gestation)` = min_bc) |> 
    dplyr::arrange(cohort) |> 
    dplyr::left_join(x = cohorts[, c("id", "abb")], by = c("id" = "cohort")) |> 
    dplyr::select(-id) |> 
    dplyr::rename(cohort = abb)

  health.print <- trimws(knitr::kable(health.tbl, format = "simple"))     
  cat(health.print, sep = "\n")
  
  cat("\n=============================================================\n")
  cat("FECUNDITY\n")
  cat("=============================================================\n\n")
  
  plot(tot.births[prj == 1, (list(year, birth))], type = "l", col = "grey", ylim = range(tot.births$birth),
       xlab = "", ylab = "Total No. births")
  for(i in 2:prj$param$n){
    lines(tot.births[prj == i, (list(year, birth))], col = "grey")
  }
  
  # cat("Calving events [per year]:\n")
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
  
  abrt <- tibble::enframe(obj$param$abort) |> 
    dplyr::mutate(value = 100*value) |> 
    dplyr::rename(phase = name, `Abortion (%)` = value)
  abrt$phase <- obj$param$phases
  abrt$phase <- sapply(X = abrt$phase, FUN = function(x) switch(x,
                       "base" = "Baseline",
                       "const" = "Construction",
                       "ops" = "Operation & maintenance"))
  abrt.print <- trimws(knitr::kable(abrt, format = "simple"))
  cat(abrt.print, sep = "\n")
  cat("\n")
  
  }
