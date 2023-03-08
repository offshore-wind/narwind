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
  
  
  # locations <- obj$locs
  
  n.ind <- dim(obj$locs[[1]])[3]
  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  cohort.id <- obj$param$cohort.id
  init.month <- obj$init$month
  
  # Compile all data
  sim <- purrr::map2(
    .x = obj$sim,
    .y = obj$locs,
    .f = ~ {
      a <- cbind(array2dt(.y), array2dt(.x))
      a$region <- sort(regions$region)[a$region]
      add_whale(a, n.ind = n.ind)
    }
  )
  
  # locations <- purrr::map(.x = sim, .f = ~.x[, .(easting, northing, region)])
  
  cat("=============================================================\n")
  cat("SIMULATIONS\n")
  cat("=============================================================\n\n")
  
  cat("No. animals:", format(n.ind, big.mark = ","), "\n")
  cat("Cohort(s):", cohort.names, "\n")
  cat("Simulation start:", month.name[init.month], "\n")
  cat("\n")
  
  # Dynamic scoping to make sure geodesic calculations are fast
  environment(summary_geo) <- .GlobalEnv
  
  cat("=============================================================\n")
  cat("HEALTH\n")
  cat("=============================================================")
  
  # row.alive <- purrr::map(.x = cohort.ab, .f = ~sim[[.x]][, list(cohort = .x, row = which.min(alive) - 1), whale]) |>
  #   do.call(what = rbind)

  dead.dt <- purrr::map(.x = sim, .f = ~{
        .x[, list(row = which.min(alive)), whale] |>
      dplyr::mutate(dead = as.numeric(row < 365))
    }) |> tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name)

  n.dead <- dead.dt |> 
    dplyr::group_by(cohort) |> 
    dplyr::summarise(dead = sum(dead), alive = n.ind - sum(dead)) |> 
    dplyr::ungroup() |> 
    format_dt()

  # row.alive <- purrr::map(.x = seq_along(cohort.ab), .f = ~{
  #   dat <- obj[["sim"]][[cohort.ab[.x]]]
  #   if(cohort.id[.x] %in% c(4,5)) dat <- dat[[1]]
  #   out <- sapply(seq_len(n.ind), function(u) which.min(unname(dat[ , "alive", , drop = FALSE][,,u])) - 1)
  #   out[out==0] <- 365
  #   out}) |> purrr::set_names(nm = cohort.ab)
  
  
  # n.dead <- purrr::map(.x = row.alive, .f = ~ sum(.x < 365)) |>
  #   tibble::enframe() |>
  #   tidyr::unnest(cols = c(value)) |>
  #   dplyr::rename(cohort = name, dead = value) |>
  #   dplyr::mutate(alive = n.ind - dead)
  
  # if (percent) {
  #   n.dead <- n.dead |>
  #     dplyr::mutate(
  #       dead = paste(round(100 * dead / sum(n.ind), 1), "%"),
  #       alive = paste(round(100 * alive / sum(n.ind), 1), "%")
  #     )}
  
  print(knitr::kable(n.dead, format = "simple"))
  
  # Breakdown of deaths by region
  locs.dead <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ sim[[.x]][, .SD[which.min(alive)], .SDcols = c("easting", "northing", "region", "bc", "strike"), whale])
  
  n.dead.region <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~ {
      locs.dead[[.x]] |> 
        dplyr::mutate(starve = ifelse(bc < 0.05, 1, 0)) |>
        dplyr::count(region, strike, starve) |> 
        dplyr::mutate(strike = strike * n, starve = starve * n) |> 
        dplyr::group_by(region) |>
        dplyr::summarise(strike = sum(strike), starve = sum(starve)) |> 
        dplyr::ungroup() |> 
        format_dt()
    }) |> tibble::enframe() |> 
    dplyr::rename(cohort = name) |> 
    tidyr::unnest(cols = c(value))
   
  
  # locs.dead <- purrr::map(.x = cohort.ab, .f = ~{
  #   lapply(X = seq_along(row.alive[[.x]]), FUN = function(a) obj[["locs"]][[.x]][row.alive[[.x]][a],,a]) |>  do.call(what = rbind) |>
  #     tibble::as_tibble() |> dplyr::mutate(trackID = paste0("whale.", dplyr::row_number())) |>
  #     dplyr::slice(which(row.alive[[.x]] < 365))
  # }) |> purrr::set_names(nm = cohort.ab)
  
  # n.dead.region <- purrr::set_names(cohort.ab) |> purrr::map(.f = ~ {
  #   region <- locs.dead[[.x]]$region
  #   out <- janitor::tabyl(region) |>
  #     dplyr::mutate(region = sort(regions$region)[region]) |>
  #     # janitor::adorn_totals() |>
  #     dplyr::select(-percent)
  #   if (percent) {
  #     out <- out |>
  #       dplyr::mutate(percent = paste(round(100 * n / sum(n.ind), 1), "%")) |>
  #       dplyr::select(-n)
  #   }
  #   # out <- format_table(out)
  #   row.names(out) <- NULL
  #   out
  # }) |> tibble::enframe() |> 
  #   dplyr::rename(cohort = name) |> 
  #   tidyr::unnest(cols = c(value))
  
  print(knitr::kable(n.dead.region, format = "simple"))
  
  cat("\n")
  
  # # Set data after death to NA
  # obj[["locs"]] <- purrr::map2(.x = obj[["locs"]], .y = row.alive, .f = ~{
  #   for(i in seq_len(n)){
  #     if(.y[i] < 365){
  #       r <- seq(.y[i]+1,365)
  #       .x[r,,i] <- NA
  #     }}
  #   .x
  # })
  # 
  # obj[["sim"]] <- purrr::map2(.x = obj[["sim"]], .y = row.alive, .f = ~{
  #   for(i in seq_len(n)){
  #     if(.y[i] < 365){
  #       r <- seq(.y[i]+1,365)
  #       .x[r,,i] <- NA
  #     }}
  #   .x
  # })
  # 
 
  
  # ............................................................
  # Load spatial support
  # ............................................................
  
  # if(geodesic){
  # coords <- sp::coordinates(density_support)
  # colnames(coords) <- c("x", "y")
  # map_limits <- c(range(coords[,1]), range(coords[,2]))
  # map_resolution <- raster::res(raster::raster(density_support))
  # geomap <- t(1*raster::as.matrix(density_support))
  # }
  
  # ............................................................
  # Calculate distances traveled and tally regions visited
  # ............................................................
  
  # dist_euclid <- function(df){
  #   obj.ind <- df[, c("easting", "northing")]
  #   obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
  #   dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
  #   dd <- dd[!is.na(dd)]
  #   unname(dd)
  # }
  # 
  # dist_euclid <- function(df){
  #   obj.ind <- df[, c("easting", "northing")]
  #   obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
  #   dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
  #   dd <- dd[!is.na(dd)]
  #   unname(dd)
  # }
  # 
  # out <- purrr::set_names(cohort.ab) |>
  #   purrr::map2(.y = cohort.names, .f = ~{
  # 
  #     # ............................................................
  #     # Calculate Euclidean distances between pairs of sequential points
  #     # ............................................................
  # 
  #     dist.euclid <- sapply(seq(n.ind), function(x) dist_euclid(locations[[.x]][,,x]), simplify = FALSE)
  # 
  #     dist.euclid.region <- lapply(seq(n.ind), function(x) cbind(locations[[.x]][,,x], c(dist.euclid[[x]], NA))) |>
  #       do.call(what = rbind) |>
  #       data.frame() |>
  #       tibble::as_tibble() |>
  #       dplyr::rename_at(vars(starts_with("V")), function(x) "step") |>
  #       dplyr::mutate(region = sort(regions$region)[region]) |>
  #       dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |>
  #       dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region))
  # 
  # 
  #     # dist.euclid.region <- sort(regions$region)[dist.euclid.region$region])
  # 
  #     # dist.euclid.region <- dist.euclid.region |> dplyr::mutate(region = ifelse(sort(regions$region)[dist.euclid.region$region])
  # 
  #     # dist.euclid <- lapply(seq(n.ind), function(x){
  #     #   obj.ind <- locations[[.x]][,c("easting", "northing"),x]
  #     #   obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
  #     #   dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
  #     #   dd <- dd[!is.na(dd)]
  #     #   unname(dd)
  #     # })
  # 
  #     # ............................................................
  #     # Calculate geodesic distances between pairs of sequential points
  #     # ............................................................
  # 
  #     if(geodesic){
  #       dist.geo <- summary_geo(obj, map_limits, map_resolution, geomap)
  #     } else {
  #       dist.geo <- NA
  #     }
  # 
  #     # ............................................................
  #     # Total migration distances
  #     # ............................................................
  # 
  #     tot.dist.euclid <- purrr::map_dbl(.x = dist.euclid, .f = ~sum(.x))
  #     if(geodesic) tot.dist.geo <- purrr::map_dbl(.x = dist.geo, .f = ~sum(.x)) else tot.dist.geo <- NULL
  # 
  #     # ............................................................
  #     # Tally points by region
  #     # ............................................................
  # 
  #     obj.sp <- lapply(seq(dim(locations[[.x]])[3]), function(x){
  #       sp::SpatialPoints(coords = locations[[.x]][ , c("easting", "northing"), x], proj4string = narw_crs())
  #     })
  # 
  #     locations.per.region <- sp::over(regions, do.call(rbind, obj.sp), returnList = TRUE) |>
  #       purrr::set_names(nm = regions$region)
  # 
  #     pr.ind <- purrr::map(.x = obj.sp, .f = ~sp::over(regions, .x, returnList = TRUE)) |>
  #       purrr::map(.f = ~purrr::set_names(x = .x, nm = regions$region))
  # 
  #     regions.per.ind <- purrr::map_dbl(.x = pr.ind, .f = ~length(purrr::discard(.x = .x, .p = function(x) length(x)==0)))
  # 
  #     days.per.regions <- purrr::map(sort(regions$region), .f = ~{
  #       sapply(X = pr.ind, FUN = function(x) length(x[[.x]]))
  #     }) |> purrr::set_names(sort(regions$region))
  # 
  #     ind.per.regions <- purrr::map(.x = days.per.regions, .f = ~as.numeric(sum(.x>0)))
  # 
  #     locations.per.country <-
  #       list("Canada" = length(unlist(locations.per.region[c("GSL", "SCOS", "CABOT", "BOF_lower", "BOF_upper")])),
  #            "U.S" = length(unlist(locations.per.region[c("CCB", "MIDA", "SNE", "SEUS", "GOM")])))
  # 
  #     locations.per.region <- purrr::map(.x = locations.per.region, .f = ~length(.x))
  # 
  #     # Return outputs
  # 
  #     list(dist = list(euclid = dist.euclid,
  #                      euclid.region = dist.euclid.region,
  #                      geo = dist.geo),
  #          totdist = list(euclid = tot.dist.euclid,
  #                         geo = tot.dist.geo),
  #          n = list(locreg = locations.per.region,
  #                   locctry = locations.per.country,
  #                   dayreg = days.per.regions,
  #                   indreg = ind.per.regions,
  #                   regind = regions.per.ind))
  #   }) # End purrr loop
  
  # ............................................................
  # OUTPUTS
  # ............................................................
  
  move.out <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~{
      all.d <- sim[[.x]][, .SD[1:(which.min(alive)-1)], .SDcols = c("d_travel", "region"), whale]
      all.d[, d_migr := sum(d_travel), whale]
      all.d$cohort <- .x
      all.d[]
    })
  
  par(mfrow = c(1,2))

  purrr::walk(.x = seq_along(move.out), .f = ~{

    # Step lengths -- overall
    hist(move.out[[.x]]$d_travel, freq = FALSE, main = cohort.names[.x], xlab = "Step length (km)", breaks = 20)
    lines(density(move.out[[.x]]$d_travel, adjust = 2), col = "orange", lwd = 1.5)

    # Step lengths -- by region
    boxplot(d_travel~region, data = move.out[[.x]], ylab = "Step lengths (km)", xlab = "")

  })
  
  # par(mfrow = c(1,2))
  # 
  # purrr::walk(.x = seq_along(out), .f = ~{
  #   
  #   deucl <- unlist(out[[cohort.ab[.x]]][["dist"]][["euclid"]])
  #   dreg <- out[[cohort.ab[.x]]][["dist"]][["euclid.region"]]
  #   
  #   # Step lengths -- overall
  #   hist(deucl, freq = FALSE, main = cohort.names[.x], xlab = "Step length (km)", breaks = 20)
  #   lines(density(deucl, adjust = 2), col = "orange", lwd = 1.5)
  #   
  #   # Step lengths -- by region
  #   boxplot(step~region, data = dreg, ylab = "Step lengths (km)", xlab = "")
  #   
  # })
  
  par(mfrow = c(1,1))

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
    
    
    
    # locs.by.area |> 
    # dplyr::group_by(cohort, region) |> 
    # dplyr::summarise(n = sum(N), .groups = "keep") |> 
    # tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    # janitor::adorn_totals( where = "row")
    # (\(.) (if(percent) janitor::adorn_totals(., where = "row") |> 
    #          dplyr::mutate(across(where(is.numeric), ~ paste(round(100 * .x / Total, 1), "%"))) |> 
    #          dplyr::select(-Total)
    # ))() 
  
  # locs.by.region <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   tibble::enframe(.x$n$locreg) |>
  #     dplyr::mutate(value = unlist(value)) |>
  #     dplyr::rename(region = name, n = value) |>
  #     dplyr::mutate(cohort = .y) |>
  #     dplyr::relocate(cohort, .after = region) |>
  #     dplyr::mutate(percent = round(100 * n / sum(n), 1))
  # }) |> do.call(what = rbind)
  
  # if (percent) {
  #   locs.by.region <- locs.by.region |>
  #     dplyr::select(-n) |>
  #     tidyr::pivot_wider(names_from = cohort, values_from = percent) |>
  #     janitor::adorn_totals() |>
  #     dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  # } else {
  #   locs.by.region <- locs.by.region |>
  #     dplyr::select(-percent) |>
  #     tidyr::pivot_wider(names_from = cohort, values_from = n) |>
  #     janitor::adorn_totals()
  # }
  # 
  # locs.by.region <- format_table(locs.by.region) |> dplyr::select(-country)
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================")
  
  print(knitr::kable(locs.by.region, format = "simple"))
  print(knitr::kable(locs.by.country, format = "simple"))
  cat("\n")
  
  # if(percent){
  #   locs.by.country <- locs.by.country |> 
  #     dplyr::select(-n) |> 
  #     dplyr::arrange(country) |>
  #     tidyr::pivot_wider(names_from = cohort, values_from = percent) |> 
  #     janitor::adorn_totals() |>
  #     dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  # } else {
  #   locs.by.country <- locs.by.country |> 
  #     dplyr::select(-percent) |> 
  #     tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
  #     dplyr::arrange(country) |>
  #     janitor::adorn_totals()
  # }
  # 
  # locs.by.country <- format_table(locs.by.country)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  cat("=============================================================\n")
  cat("MOVEMENTS (km)\n")
  cat("=============================================================\n")
  
  # cat("\n+++++++++++ Daily steps +++++++++++", sep = "")
  
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
    
  
  # dist.df <- purrr::map2(.x = move.out, .y = cohort.ab, .f = ~{
  #   fivenum(.x$d_travel)[c(1,3,5)] |> 
  #     tibble::as_tibble(.name_repair = "minimal") |> 
  #     dplyr::mutate(param = c("min", "mean", "max")) |> 
  #     dplyr::bind_rows(tibble::tibble(value = sd(.x$d_travel), param = "sd", cohort = .y)) |> 
  #     dplyr::rename(distance = value) |> 
  #     dplyr::mutate(distance = round(distance, 1)) |> 
  #     dplyr::mutate(cohort = .y)
  # }) |> do.call(what = rbind) |> 
  #   tidyr::pivot_wider(names_from = param, values_from = distance)
  # 
  # dist.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   fivenum(unlist(.x$dist$euclid))[c(1,3,5)] |> 
  #     tibble::as_tibble(.name_repair = "minimal") |> 
  #     dplyr::mutate(param = c("min", "mean", "max")) |> 
  #     dplyr::bind_rows(tibble::tibble(value = sd(unlist(.x$dist$euclid)), param = "sd", cohort = .y)) |> 
  #     dplyr::rename(distance = value) |> 
  #     dplyr::mutate(distance = round(distance, 1)) |> 
  #     dplyr::mutate(cohort = .y)
  # }) |> do.call(what = rbind) |> 
  #   tidyr::pivot_wider(names_from = param, values_from = distance)
  
  # geodesic.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   fivenum(unlist(.x$dist$geo))[c(1,3,5)] |> 
  #     tibble::as_tibble(.name_repair = "minimal") |> 
  #     dplyr::mutate(param = c("min", "mean", "max")) |> 
  #     dplyr::bind_rows(tibble::tibble(value = sd(unlist(.x$dist$geo)), param = "sd", cohort = .y)) |> 
  #     dplyr::rename(distance = value) |> 
  #     dplyr::mutate(distance = round(distance, 1)) |> 
  #     dplyr::mutate(cohort = .y)
  # }) |> do.call(what = rbind) |> 
  #   tidyr::pivot_wider(names_from = param, values_from = distance)
  # 
  # dist.type <- "euclidean"
  
  # if(geodesic){
  #   dist.type <-  c(dist.type, "geodesic")
  #   dist.df <- rbind(dist.df, geodesic.df)
  # }
  # 
  # dist.df <- dist.df |> 
    # dplyr::mutate(type = dist.type) |> 
    # dplyr::relocate(type, .before = cohort) |> 
    # dplyr::mutate(step = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
    # dplyr::select(-min, -max, -mean, -sd)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  # print(knitr::kable(dist.df, format = "simple"))
  # cat("\n")
  # 
  # cat("\n+++++++++++ Migration +++++++++++", sep = "")
  
  # tot.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   tibble::tibble(d = c(.x$totdist$euclid, .x$totdist$geo), 
  #                  type = rep(dist.type, each = length(.x$totdist$euclid))) |> 
  #     dplyr::mutate(cohort = .y) |> 
  #     dplyr::group_by(type, cohort) |> 
  #     dplyr::summarise(min = format(round(min(d),0), big.mark = ","), 
  #                      max = format(round(max(d),0), big.mark = ","), 
  #                      mean = format(round(mean(d),0), big.mark = ","),
  #                      sd = format(round(sd(d),0), big.mark = ","), .groups = 'drop')
  # }) |> do.call(what = rbind) |> 
  #   dplyr::mutate(migration = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
  #   dplyr::select(-min, -max, -mean, -sd)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
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

      
  #     
  #     dplyr::group_by(region) |> 
  #     dplyr::summarise(n = dplyr::n_distinct(whale)) |> 
  #     dplyr::left_join(x = tibble::tibble(region = sort(regions$region)), by = "region") |> 
  #     dplyr::mutate(n = dplyr::coalesce(n, 0), cohort = .x)
  # }) |> do.call(what = rbind) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
  #   format_dt()
  
  # indreg.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   tibble::enframe(.x$n$indreg) |> 
  #     dplyr::mutate(value = unlist(value)) |> 
  #     dplyr::rename(n = value, region = name) |> 
  #     dplyr::mutate(cohort = .y)
  # }) |> do.call(what = rbind)
  # 
  # if(percent){
  #   indreg.df <- indreg.df |> 
  #   dplyr::mutate(perc = round(100 * n / n.ind)) |> 
  #   dplyr::select(-n) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = perc) |> 
  #   dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  # } else {
  #   indreg.df <- indreg.df |> 
  #     tidyr::pivot_wider(names_from = cohort, values_from = n) 
  # }
  # 
  # Days spent in each region
  
  # days.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #   
  #   purrr::map(.x = .x$n$dayreg, .f = ~{
  #     if(all(.x == 0)) 0 else .x[.x >0]
  #   }) |> tibble::enframe() |> 
  #     dplyr::rename(region = name, days = value) |> 
  #     dplyr::group_by(region) |> 
  #     dplyr::mutate(cohort = .y) |> 
  #     dplyr::group_by(region, cohort) |> 
  #     dplyr::summarise(min = min(unlist(days)),
  #                      max = max(unlist(days)),
  #                      mean = round(mean(unlist(days)),1),
  #                      sd = round(sd(unlist(days)), 1), .groups = 'drop')
  # }) |> do.call(what = rbind) |> 
  #   dplyr::mutate(days = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
  #   dplyr::select(-min, -max, -mean, -sd) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = days) 

  
  # days.min <- days.df |> dplyr::select(region, cohort, min) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = min)
  # 
  # days.max <- days.df |> dplyr::select(region, cohort, max) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = max) 
  # 
  # days.mean <- days.df |> dplyr::select(region, cohort, mean) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = mean) 
  # 
  # days.sd <- days.df |> dplyr::select(region, cohort, sd) |> 
  #   tidyr::pivot_wider(names_from = cohort, values_from = sd) 

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

    
  
  # regind.df <- suppressWarnings(purrr::map2(.x = out, .y = cohort.ab, .f = ~{
  #     tmp <- .x$n$regind |> 
  #     janitor::tabyl() |> dplyr::mutate(percent = 100 * percent)
  #     names(tmp)[1] <- "No.regions"
  #     tmp |>  dplyr::mutate(cohort = .y)
  #     })
  #     ) |> do.call(what = rbind)
  # 
  # if(percent){
  #   regind.df <- regind.df |> 
  #     dplyr::select(-n) |> 
  #     tidyr::pivot_wider(names_from = cohort, values_from = percent) |> 
  #     janitor::adorn_totals() |>
  #     dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  # } else {
  #   regind.df <- regind.df |> 
  #     dplyr::select(-percent) |> 
  #     tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
  #     janitor::adorn_totals()
  # }
  # 
  # regind.df <- format_table(regind.df)
  

  #   fivenum(.x$n$regind)[c(1,3,5)] |> 
  #     tibble::as_tibble(.name_repair = "minimal") |> 
  #     dplyr::mutate(param = c("min", "median", "max")) |> 
  #     dplyr::rename(n = value) |> 
  #     dplyr::mutate(cohort = .y)
  # }) |> do.call(what = rbind) |> 
  #   tidyr::pivot_wider(names_from = param, values_from = n) |> 
  #   dplyr::mutate(N = paste0(median, " [", min, " – ", max, "]")) |> 
  #   dplyr::select(-min, -max, -median)
  
  
  # ............................................................
  # Print outputs
  # ............................................................
  
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
  
  # activ <- purrr::set_names(cohort.ab) |>
  #   purrr::map(.f = ~{
  #     lapply(X = seq_len(n.ind), FUN = function(i){
  #       if(cohort.id %in% c(4,5)){
  #         tmp <- cbind(obj[["locs"]][[.x]][, "region", i, drop = FALSE],
  #                      obj[["sim"]][[.x]][[1]][, c("alive", "t_travel", "t_feed", "t_rest", "t_nurse", "feed"), i])
  #       } else {
  #         tmp <- cbind(obj[["locs"]][[.x]][, "region", i, drop = FALSE],
  #                      obj[["sim"]][[.x]][, c("alive", "t_travel", "t_feed", "t_rest", "t_nurse", "feed"), i])
  #       }
  #       colnames(tmp)[1] <- "region"
  #       tmp
  #     }) |> do.call(what = rbind)
  #   })
  # 
  # activ.df <- tibble::enframe(activ) |>
  #   dplyr::rename(cohort = name, data = value) |>
  #   dplyr::mutate(data = purrr::map(.x = data, .f = ~{
  #     tibble::as_tibble(.x)
  #   })) |> tidyr::unnest(cols = c(data)) |>
  #   dplyr::mutate(region = sort(regions$region)[region]) |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(total = sum(t_travel, t_feed, t_nurse, t_rest)) |>
  #     dplyr::ungroup() |>
  #   dplyr::filter(alive == 1)
  # 
  # 
  # 
  # cat("Total:", mean(activ.df$total), "hrs\n")

  tot.hrs <- purrr::set_names(cohort.ab) |>
    purrr::map(.f = ~ mean(sim[[.x]][
      alive == 1,
      list(total = sum(t_travel, t_feed, t_rest, t_nurse)), row_id
    ]$total)) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name, total_hrs = value)
  cat("Total:", mean(tot.hrs$total_hrs), "hrs\n")
  if (mean(tot.hrs$total_hrs) < 24) warning("Inconsistent time allocation in activity budgets")
  
  # activ.df |> dplyr::mutate(ind = 1:nrow(activ.df)) |> dplyr::filter(total < 24)
  
  # which(activ.df$total < 24)
 
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
     


   
   
   
   
    #   # Traveling, resting, and nursing
    #   dt1 <- sim[[.x]][
    #     alive == 1,
    #     list(
    #       `travel (hrs)` = paste0(round(mean(t_travel), 2), " (± ", round(sd(t_travel), 2), ")"),
    #       `rest (hrs)` = paste0(round(mean(t_rest), 2), " (± ", round(sd(t_rest), 2), ")"),
    #       `nurse (hrs)` = paste0(round(mean(t_nurse), 2), " (± ", round(sd(t_nurse), 2), ")")
    #     ), region
    #   ]
    #   
    #   # Feeding
    #   dt2 <- sim[[.x]][
    #     alive == 1 & feed == 1,
    #     list(`feed (hrs)` = paste0(round(mean(t_feed), 2), " (± ", round(sd(t_feed), 2), ")")), region
    #   ]
    #   dplyr::left_join(dt1, dt2, by = "region") |>
    #     dplyr::relocate(`feed (hrs)`, .before = "rest (hrs)") |>
    #     dplyr::left_join(x = tibble::tibble(region = sort(regions$region)), by = "region") |>
    #     dplyr::mutate_at(vars(matches("hrs")), ~ dplyr::coalesce(.x, "0 (± 0)"))
    # })
    # 
    # 
  
 # activ.df <- purrr::set_names(cohort.ab) |>
 #  purrr::map(.f = ~ {
 #    
 #    # Traveling, resting, and nursing
 #    dt1 <- sim[[.x]][
 #      alive == 1,
 #      list(
 #        `travel (hrs)` = paste0(round(mean(t_travel), 2), " (± ", round(sd(t_travel), 2), ")"),
 #        `rest (hrs)` = paste0(round(mean(t_rest), 2), " (± ", round(sd(t_rest), 2), ")"),
 #        `nurse (hrs)` = paste0(round(mean(t_nurse), 2), " (± ", round(sd(t_nurse), 2), ")")
 #      ), region
 #    ]
 # 
 #    # Feeding
 #    dt2 <- sim[[.x]][
 #      alive == 1 & feed == 1,
 #      list(`feed (hrs)` = paste0(round(mean(t_feed), 2), " (± ", round(sd(t_feed), 2), ")")), region
 #    ]
 #    dplyr::left_join(dt1, dt2, by = "region") |>
 #      dplyr::relocate(`feed (hrs)`, .before = "rest (hrs)") |>
 #      dplyr::left_join(x = tibble::tibble(region = sort(regions$region)), by = "region") |>
 #      dplyr::mutate_at(vars(matches("hrs")), ~ dplyr::coalesce(.x, "0 (± 0)"))
 #  })

  purrr::walk2(.x = cohort.ab, .y = cohort.names, .f = ~{
   cat("\n+++++++++++ ", .y, " +++++++++++")
   print(knitr::kable(activ.summary[cohort == .x], format = "simple"))
 })
 
  #  
  # 
  # travel.df <- activ.df |> dplyr::group_by(cohort, region) |>
  #   dplyr::summarise(mean = mean(t_travel), sd = sd(t_travel), .groups = "drop") |>
  #   dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |>
  #   dplyr::select(-mean, -sd) |>
  #   tidyr::pivot_wider(names_from = cohort, values_from = value) |>
  #   dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |>
  #   dplyr::arrange(region)
  # 
  # feed.df <- activ.df |> dplyr::filter(feed == 1) |> dplyr::select(-feed) |> dplyr::group_by(cohort, region) |>
  #   dplyr::summarise(mean = mean(t_feed), sd = sd(t_feed), .groups = "drop") |>
  #   dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |>
  #   dplyr::select(-mean, -sd) |>
  #   tidyr::pivot_wider(names_from = cohort, values_from = value) |>
  #   dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |>
  #   dplyr::arrange(region)
  # 
  # if(ncol(feed.df) == 1) feed.df[, cohort.ab] <- 0
  # 
  # rest.df <- activ.df |> dplyr::group_by(cohort, region) |>
  #   dplyr::summarise(mean = mean(t_rest), sd = sd(t_rest), .groups = "drop") |>
  #   dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |>
  #   dplyr::select(-mean, -sd) |>
  #   tidyr::pivot_wider(names_from = cohort, values_from = value) |>
  #   dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |>
  #   dplyr::arrange(region)
  # 
  # nurse.df <- activ.df |> dplyr::group_by(cohort, region) |>
  #   dplyr::summarise(mean = mean(t_nurse), sd = sd(t_nurse), .groups = "drop") |>
  #   dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |>
  #   dplyr::select(-mean, -sd) |>
  #   tidyr::pivot_wider(names_from = cohort, values_from = value) |>
  #   dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |>
  #   dplyr::arrange(region)
  # 
  # activ.list <- purrr::set_names(cohort.ab) |>  purrr::map(.f = ~{
  #   tmp <- dplyr::bind_cols(travel.df[, c("region", .x)], feed.df[, .x], rest.df[,.x], nurse.df[,.x], .name_repair = "minimal")
  #   names(tmp)[2:5] <- c("travel (hrs)", "feed (hrs)", "rest (hrs)", "nurse (hrs)")
  #   tmp
  # })
  # 
  # purrr::walk2(.x = cohort.ab, .y = cohort.names, .f = ~{
  #   cat("\n+++++++++++ ", .y, " +++++++++++")
  #   print(knitr::kable(activ.list[[.x]], format = "simple"))
  # })
  # 
  

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
    xlab("") +
    ggplot2::theme(axis.text = element_text(size = 10, color = "black"),
                   axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
                   axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                   axis.title = element_text(size = 14),
                   axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                   axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                   strip.background = element_rect(fill = "grey20"),
                   strip.text = element_text(colour = 'white')) +
    labs(title = cohort.names[which(cohort.ab==k)]) +
    scale_y_continuous(breaks = seq(0,24,4), limits = c(0,24))
  suppressWarnings(print(p))
  } 
  
}
