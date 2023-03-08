#' Plot
#'
#' Make plots
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

summary.narwsim <- function(obj){
  
  format_dt <- function(dt, direction = "col"){
    dt |>  janitor::adorn_percentages(denominator = direction) |>
      janitor::adorn_pct_formatting() |> 
      janitor::adorn_ns()
  }
  
  options(pillar.sigfig = 7)
  
  gg.opts <- ggplot2::theme(axis.text = element_text(size = 10, color = "black"),
                 axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
                 axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                 axis.title = element_text(size = 12),
                 axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                 axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                 strip.background = element_rect(fill = "grey20"),
                 strip.text = element_text(colour = 'white', size = 12))
  
  if(!"narwsim" %in% class(obj)) stop("Input must be of class <narwsim>")
  
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
  
  n.ind <- dim(obj$locs[[1]])[3]
  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  cohort.id <- obj$param$cohort.id
  init.month <- obj$init$month
  
  # Compile all data
  sim <- purrr::pmap(
    list(x = obj$sim, y = obj$locs, z = cohort.names),
    function(x,y,z) {
      a <- cbind(array2dt(y), array2dt(x))
      a$region <- sort(regions$region)[a$region]
      a$cohort_name <- z
      add_whale(a, n.ind = n.ind)
    }
  )
  
  cat("=============================================================\n")
  cat("SIMULATIONS\n")
  cat("=============================================================\n\n")
  
  cat("No. animals:", format(n.ind, big.mark = ","), "\n")
  cat("Cohort(s):", cohort.names, "\n")
  cat("Simulation start:", month.name[init.month], "\n")
  cat("\n")
  
  cat("=============================================================\n")
  cat("MORTALITY & HEALTH\n")
  cat("=============================================================")

  dead.dt <- purrr::map(.x = sim, .f = ~{
        .x[, list(row = which.min(alive)), whale] |>
      dplyr::mutate(dead = as.numeric(row > 1))
    }) |> tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name)

  n.dead <- dead.dt |> 
    dplyr::group_by(cohort) |> 
    dplyr::summarise(dead = sum(dead), alive = n.ind - sum(dead)) |> 
    dplyr::ungroup() |> 
    format_dt()
  
  # Breakdown of deaths by region
  locs.dead <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][, .SD[ifelse(which.min(alive)==1,0,which.min(alive))], .SDcols = c("easting", "northing", "region", "bc", "strike"), whale])
  
  n.dead.region <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~ {
      if(nrow(locs.dead[[.x]]) > 0){
      locs.dead[[.x]] |> 
        dplyr::mutate(starve = ifelse(bc < 0.05, 1, 0)) |>
        dplyr::count(region, strike, starve) |> 
        dplyr::mutate(strike = strike * n, starve = starve * n) |> 
        dplyr::group_by(region) |>
        dplyr::summarise(strike = sum(strike), starve = sum(starve)) |> 
        dplyr::ungroup() |> 
        format_dt()
      } else {
        NULL
      }
    }) |> tibble::enframe() |> 
    dplyr::rename(cohort = name) |> 
    tidyr::unnest(cols = c(value))
 
  print(knitr::kable(n.dead, format = "simple"))
  print(knitr::kable(n.dead.region, format = "simple"))
  cat("\n")
 
  
  bodyc.plot <- ggplot2::ggplot(data = do.call(sim, what = rbind), aes(x = day, y = bc, group = whale)) +
    ggplot2::geom_line(col = "grey70") +
    ggplot2::facet_wrap(. ~ factor(cohort_name), scales = 'free', ncol = 2) +
    gg.opts +
    ylab("Body condition") +
    xlab("Day of the year") +
    ggplot2::scale_x_continuous(breaks = c(0, seq(5, 365, by = 15))) +
    ggplot2::scale_y_continuous(breaks = seq(0,1, 0.025)) +
    ggplot2::geom_hline(yintercept = 0.05, col = "#cb4154")
  
  print(bodyc.plot)
  
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================")
  
  move.out <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~{
      all.d <- sim[[.x]][, .SD[1:(ifelse(which.min(alive)==1, 365, which.min(alive)-1))], .SDcols = c("d_travel", "region"), whale]
      all.d[, d_migr := sum(d_travel), whale]
      all.d$cohort <- .x
      all.d[]
    })
  
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
  
  par(mfrow = c(1,2))
  purrr::walk(.x = seq_along(move.out), .f = ~{
    
    # Step lengths
    hist(move.out[[.x]]$d_travel, freq = FALSE, main = cohort.names[.x], xlab = "Step length (km)", breaks = 20)
    lines(density(move.out[[.x]]$d_travel, adjust = 2), col = "orange", lwd = 1.5)
    boxplot(d_travel~region, data = move.out[[.x]], ylab = "Step lengths (km)", xlab = "")
    
  })
  par(mfrow = c(1,1))
  
  
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
    dplyr::mutate(step = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |>
    dplyr::select(-min, -max, -mean, -sd)

  migr.df <- purrr::map(.x = cohort.ab, .f = ~
    move.out[[.x]][, list(
      min = format(round(min(d_migr), 0), big.mark = ","), 
      mean = format(round(mean(d_migr), 0), big.mark = ","), 
      max = format(round(max(d_migr), 0), big.mark = ","), 
      sd = format(round(sd(d_migr), 0), big.mark = ",")), ] |>
      dplyr::mutate(cohort = .x)) |>
    do.call(what = rbind) |>
    dplyr::relocate(cohort, .before = min) |>
    dplyr::mutate(migration = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |>
    dplyr::select(-min, -max, -mean, -sd)
 
  print(knitr::kable(dplyr::left_join(step.df, migr.df, by = "cohort"), format = "simple"))
  cat("\n")
  
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
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    format_dt()
    
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
      dplyr::mutate(days = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |>
      dplyr::select(-min, -max, -mean, -sd) |> 
      tidyr::pivot_wider(names_from = cohort, values_from = days)

  # Number of regions per individual
  
  reg.per.ind.df <- purrr::map(.x = cohort.ab, .f = ~{
    move.out[[.x]] |> 
      dplyr::group_by(whale) |> 
      dplyr::summarise(nreg = dplyr::n_distinct(region)) |> 
      janitor::tabyl(nreg) |> 
      dplyr::mutate(cohort = .x)}) |> 
    do.call(what = rbind) |> 
    dplyr::select(-percent) |>
    dplyr::rename(No.regions = nreg) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    format_dt()
  
  cat("+++++++++++ Number of animals visiting each region (N = ", n.ind, ") +++++++++++", sep = "")
  print(knitr::kable(indiv_reg.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Days spent in each region +++++++++++")
  print(knitr::kable(days.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Total number of regions visited +++++++++++")
  print(knitr::kable(reg.per.ind.df, format = "simple"))
 
  cat("\n\n")
  cat("=============================================================\n")
  cat("ACTIVITY BUDGETS\n")
  cat("=============================================================\n")
 
  tot.hrs <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ mean(sim[[.x]][
      alive == 1,
      list(total = sum(t_travel, t_feed, t_rest, t_nurse)), row_id]$total)) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name, total_hrs = value)
  cat("Total:", mean(tot.hrs$total_hrs), "hrs\n")
  if (mean(tot.hrs$total_hrs) < 24) warning("Inconsistent time allocation in activity budgets")
  
   activity.df <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][alive == 1, list(cohort = .x, region, alive, t_travel, t_feed, t_rest, t_nurse, feed)]) |> 
     do.call(what = rbind)
      
   activ.summary.other <- activity.df[,list(
     `travel (hrs)` = paste0(round(mean(t_travel), 2), " (± ", round(sd(t_travel), 2), ")"),
     `rest (hrs)` = paste0(round(mean(t_rest), 2), " (± ", round(sd(t_rest), 2), ")"),
     `nurse (hrs)` = paste0(round(mean(t_nurse), 2), " (± ", round(sd(t_nurse), 2), ")")
   ), list(cohort, region)]
   
   activ.summary.feed <- activity.df[feed == 1,
     list(`feed (hrs)` = paste0(round(mean(t_feed), 2), " (± ", round(sd(t_feed), 2), ")")), list(cohort, region)]
   
   dt <- expand.grid(cohort.ab, sort(regions$region)) |> 
     dplyr::rename(cohort = Var1, region = Var2) |> 
     data.table::as.data.table()
   
   activ.summary <- dplyr::left_join(activ.summary.other, activ.summary.feed, by = c("cohort", "region")) |> 
     dplyr::left_join(x = dt, by = c("cohort", "region")) |>
     dplyr::mutate_at(vars(matches("hrs")), ~ dplyr::coalesce(.x, "0 (± 0)"))

  purrr::walk2(.x = cohort.ab, .y = cohort.names, .f = ~{
   cat("\n+++++++++++ ", .y, " +++++++++++")
   print(knitr::kable(activ.summary[cohort == .x], format = "simple"))
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
  
  for(k in cohort.ab){
  tmp <- activ.out |> dplyr::mutate(cohort = k)
  p <- ggplot2::ggplot(tmp, aes(x = factor(region), y = hrs)) + 
    ggplot2::geom_boxplot(fill = "darkgrey") +
    ggplot2::facet_wrap(. ~ factor(behav), scales = 'free', ncol = 2) +
    ylab("Time spent (hrs)") +
    xlab("") + gg.opts +
    labs(title = cohort.names[which(cohort.ab==k)]) +
    scale_y_continuous(breaks = seq(0,24,4), limits = c(0,24))
  suppressWarnings(print(p))
  } 
  
}
