#' Plot
#'
#' Make plots
#'
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
                            percent = TRUE,
                            show.legend = FALSE,
                            geodesic = FALSE){
  
  # obj = m
  # percent = TRUE
  # show.legend = FALSE
  # geodesic = FALSE

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
  
  locations <- obj$locs
  n.ind <- dim(locations[[1]])[3]
  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  cohort.id <- obj$param$cohort.id
  init.month <- obj$param$init$month
  
  if(show.legend){
    cat("BOF: Bay of Fundy (lower, upper)\n")
    cat("CABOT: Cabot Strait\n")
    cat("CCB: Cape Cod Bay\n")
    cat("GOM: Gulf of Maine and Georges Bank\n")
    cat("GSL: Gulf of St Lawrence\n")
    cat("MIDA: Mid-Atlantic\n")
    cat("SCOS: Scotian Shelf\n")
    cat("SEUS: South-east United States\n")
    cat("SNE: Southern New England\n\n")
  }
  
  cat("=============================================================\n")
  cat("SIMULATIONS\n")
  cat("=============================================================\n\n")
  
  cat("N animals:", format(n.ind, big.mark = ","), "\n")
  cat("Cohort(s):", cohort.names, "\n")
  cat("Initialization:", month.name[init.month], "\n")
  cat("\n")
  
  # Dynamic scoping to make sure geodesic calculations are fast
  environment(summary_geo) <- .GlobalEnv
  
  cat("=============================================================\n")
  cat("HEALTH\n")
  cat("=============================================================")

  row.alive <- purrr::map(.x = seq_along(cohort.ab), .f = ~{
    dat <- obj[["sim"]][[cohort.ab[.x]]]
    if(cohort.id[.x] %in% c(4,5)) dat <- dat[[1]]
    out <- sapply(seq_len(n.ind), function(u) which.min(unname(dat[ , "alive", , drop = FALSE][,,u])) - 1)
    out[out==0] <- 365
    out}) |> purrr::set_names(nm = cohort.ab)
  
  n.dead <- purrr::map(.x = row.alive, .f = ~ sum(.x < 365)) |>
    tibble::enframe() |>
    tidyr::unnest(cols = c(value)) |>
    dplyr::rename(cohort = name, dead = value) |>
    dplyr::mutate(alive = n.ind - dead)
  if (percent) {
    n.dead <- n.dead |>
      dplyr::mutate(
        dead = paste(round(100 * dead / sum(n.ind), 1), "%"),
        alive = paste(round(100 * alive / sum(n.ind), 1), "%")
      )}
  
  print(knitr::kable(n.dead, format = "simple"))
  
  locs.dead <- purrr::map(.x = cohort.ab, .f = ~{
    lapply(X = seq_along(row.alive[[.x]]), FUN = function(a) obj[["locs"]][[.x]][row.ind[a],,a]) |>  do.call(what = rbind) |> 
      tibble::as_tibble() |> dplyr::mutate(trackID = paste0("whale.", dplyr::row_number())) |> 
      dplyr::slice(which(row.alive[[.x]] < 365))
  }) |> purrr::set_names(nm = cohort.ab)
  
  n.dead.region <- purrr::set_names(x = cohort.ab) |> purrr::map(.f = ~ {
    region <- locs.dead[[.x]]$region
    out <- janitor::tabyl(region) |>
      dplyr::mutate(region = sort(regions$region)[region]) |>
      # janitor::adorn_totals() |>
      dplyr::select(-percent)
    if (percent) {
      out <- out |>
        dplyr::mutate(percent = paste(round(100 * n / sum(n.ind), 1), "%")) |>
        dplyr::select(-n)
    }
    # out <- format_table(out)
    row.names(out) <- NULL
    out
  }) |> tibble::enframe() |> 
    dplyr::rename(cohort = name) |> 
    tidyr::unnest(cols = c(value))
  
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
  
  if(geodesic){
  coords <- sp::coordinates(density_support)
  colnames(coords) <- c("x", "y")
  map_limits <- c(range(coords[,1]), range(coords[,2]))
  map_resolution <- raster::res(raster::raster(density_support))
  geomap <- t(1*raster::as.matrix(density_support))
  }
  
  # ............................................................
  # Calculate distances travelled and tally regions visited
  # ............................................................
  
  dist_euclid <- function(df){
    obj.ind <- df[, c("easting", "northing")]
    obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
    dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
    dd <- dd[!is.na(dd)]
    unname(dd)
  }
  
  dist_euclid <- function(df){
    obj.ind <- df[, c("easting", "northing")]
    obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
    dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
    dd <- dd[!is.na(dd)]
    unname(dd)
  }
  
  out <- purrr::set_names(cohort.ab) |> 
    purrr::map2(.y = cohort.names, .f = ~{
      
      # ............................................................
      # Calculate Euclidean distances between pairs of sequential points
      # ............................................................
      
      dist.euclid <- sapply(seq(n.ind), function(x) dist_euclid(locations[[.x]][,,x]), simplify = FALSE)
      
      dist.euclid.region <- lapply(seq(n.ind), function(x) cbind(locations[[.x]][,,x], c(dist.euclid[[x]], NA))) |> 
        do.call(what = rbind) |> 
        data.frame() |> 
        tibble::as_tibble() |> 
        dplyr::rename_at(vars(starts_with("V")), function(x) "step") |> 
        dplyr::mutate(region = sort(regions$region)[region]) |> 
        dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |> 
        dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region))
      # dist.euclid.region <- sort(regions$region)[dist.euclid.region$region])
      
      # dist.euclid.region <- dist.euclid.region |> dplyr::mutate(region = ifelse(sort(regions$region)[dist.euclid.region$region])

      # dist.euclid <- lapply(seq(n.ind), function(x){
      #   obj.ind <- locations[[.x]][,c("easting", "northing"),x]
      #   obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
      #   dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
      #   dd <- dd[!is.na(dd)]
      #   unname(dd)
      # })
      
      # ............................................................
      # Calculate geodesic distances between pairs of sequential points
      # ............................................................
      
      if(geodesic){
        dist.geo <- summary_geo(obj, map_limits, map_resolution, geomap)
      } else {
        dist.geo <- NA
      }
      
      # ............................................................
      # Total migration distances
      # ............................................................
      
      tot.dist.euclid <- purrr::map_dbl(.x = dist.euclid, .f = ~sum(.x))
      if(geodesic) tot.dist.geo <- purrr::map_dbl(.x = dist.geo, .f = ~sum(.x)) else tot.dist.geo <- NULL
      
      # ............................................................
      # Tally points by region
      # ............................................................
      
      obj.sp <- lapply(seq(dim(locations[[.x]])[3]), function(x){
        sp::SpatialPoints(coords = locations[[.x]][ , c("easting", "northing"), x], proj4string = narw_crs())
      })
      
      locations.per.region <- sp::over(regions, do.call(rbind, obj.sp), returnList = TRUE) |> 
        purrr::set_names(nm = regions$region)
      
      pr.ind <- purrr::map(.x = obj.sp, .f = ~sp::over(regions, .x, returnList = TRUE)) |> 
        purrr::map(.f = ~purrr::set_names(x = .x, nm = regions$region))
      
      regions.per.ind <- purrr::map_dbl(.x = pr.ind, .f = ~length(purrr::discard(.x = .x, .p = function(x) length(x)==0)))
      
      days.per.regions <- purrr::map(sort(regions$region), .f = ~{
        sapply(X = pr.ind, FUN = function(x) length(x[[.x]]))
      }) |> purrr::set_names(sort(regions$region))
      
      ind.per.regions <- purrr::map(.x = days.per.regions, .f = ~as.numeric(sum(.x>0)))
      
      locations.per.country <- 
        list("Canada" = length(unlist(locations.per.region[c("GSL", "SCOS", "CABOT", "BOF_lower", "BOF_upper")])),
             "U.S" = length(unlist(locations.per.region[c("CCB", "MIDA", "SNE", "SEUS", "GOM")])))
      
      locations.per.region <- purrr::map(.x = locations.per.region, .f = ~length(.x))
      
      # Return outputs
      
      list(dist = list(euclid = dist.euclid, euclid.region = dist.euclid.region, geo = dist.geo),
           totdist = list(euclid = tot.dist.euclid, geo = tot.dist.geo),
           n = list(locreg = locations.per.region,
                    locctry = locations.per.country,
                    dayreg = days.per.regions,
                    indreg = ind.per.regions,
                    regind = regions.per.ind))
    }) # End purrr loop
  
  # ............................................................
  # OUTPUTS
  # ............................................................
  
  par(mfrow = c(1,2))
  
  purrr::walk(.x = seq_along(out), .f = ~{
    
    deucl <- unlist(out[[cohort.ab[.x]]][["dist"]][["euclid"]])
    dreg <- out[[cohort.ab[.x]]][["dist"]][["euclid.region"]]
    
    # Step lengths -- overall
    hist(deucl, freq = FALSE, main = cohort.names[.x], xlab = "Step length (km)", breaks = 20)
    lines(density(deucl, adjust = 2), col = "orange", lwd = 1.5)
    
    # Step lengths -- by region
    boxplot(step~region, data = dreg, ylab = "Step lengths (km)", xlab = "")
    
  })
  
  par(mfrow = c(1,1))
  
  locs.by.region <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::enframe(.x$n$locreg) |> 
      dplyr::mutate(value = unlist(value)) |> 
      dplyr::rename(region = name, n = value) |>
      dplyr::mutate(cohort = .y) |> 
      dplyr::relocate(cohort, .after = region) |> 
      dplyr::mutate(percent = round(100 * n / sum(n), 1))
  }) |> do.call(what = rbind)
  
  if(percent){

    locs.by.region <- locs.by.region |> 
      dplyr::select(-n) |> 
      dplyr::arrange(region) |>
      tidyr::pivot_wider(names_from = cohort, values_from = percent) |> 
      janitor::adorn_totals() |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
    
    } else {
      
      locs.by.region <- locs.by.region |> 
        dplyr::select(-percent) |> 
        tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
        dplyr::arrange(region) |>
        janitor::adorn_totals()
      
    }
  
  locs.by.region <- format_table(locs.by.region)
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================")
  
  print(knitr::kable(locs.by.region, format = "simple"))
  
  locs.by.country <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::enframe(.x$n$locctry) |> 
      dplyr::mutate(value = unlist(value)) |> 
      dplyr::rename(country = name, n = value) |>
      dplyr::mutate(cohort = .y) |> 
      dplyr::mutate(percent = round(100 * n / sum(n), 1))
  }) |> do.call(what = rbind) 
  
  if(percent){
    locs.by.country <- locs.by.country |> 
      dplyr::select(-n) |> 
      dplyr::arrange(country) |>
      tidyr::pivot_wider(names_from = cohort, values_from = percent) |> 
      janitor::adorn_totals() |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  } else {
    locs.by.country <- locs.by.country |> 
      dplyr::select(-percent) |> 
      tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
      dplyr::arrange(country) |>
      janitor::adorn_totals()
  }
  
  locs.by.country <- format_table(locs.by.country)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(locs.by.country, format = "simple"))
  cat("\n")
  
  cat("=============================================================\n")
  cat("MOVEMENTS (km)\n")
  cat("=============================================================\n")
  
  # cat("\n+++++++++++ Daily steps +++++++++++", sep = "")
  
  dist.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    fivenum(unlist(.x$dist$euclid))[c(1,3,5)] |> 
      tibble::as_tibble(.name_repair = "minimal") |> 
      dplyr::mutate(param = c("min", "mean", "max")) |> 
      dplyr::bind_rows(tibble::tibble(value = sd(unlist(.x$dist$euclid)), param = "sd", cohort = .y)) |> 
      dplyr::rename(distance = value) |> 
      dplyr::mutate(distance = round(distance, 1)) |> 
      dplyr::mutate(cohort = .y)
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = param, values_from = distance)
  
  geodesic.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    fivenum(unlist(.x$dist$geo))[c(1,3,5)] |> 
      tibble::as_tibble(.name_repair = "minimal") |> 
      dplyr::mutate(param = c("min", "mean", "max")) |> 
      dplyr::bind_rows(tibble::tibble(value = sd(unlist(.x$dist$geo)), param = "sd", cohort = .y)) |> 
      dplyr::rename(distance = value) |> 
      dplyr::mutate(distance = round(distance, 1)) |> 
      dplyr::mutate(cohort = .y)
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = param, values_from = distance)
  
  dist.type <- "euclidean"
  
  if(geodesic){
    dist.type <-  c(dist.type, "geodesic")
    dist.df <- rbind(dist.df, geodesic.df)
  }
  
  dist.df <- dist.df |> 
    dplyr::mutate(type = dist.type) |> 
    dplyr::relocate(type, .before = cohort) |> 
    dplyr::mutate(step = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
    dplyr::select(-min, -max, -mean, -sd)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  # print(knitr::kable(dist.df, format = "simple"))
  # cat("\n")
  # 
  # cat("\n+++++++++++ Migration +++++++++++", sep = "")
  
  tot.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::tibble(d = c(.x$totdist$euclid, .x$totdist$geo), 
                   type = rep(dist.type, each = length(.x$totdist$euclid))) |> 
      dplyr::mutate(cohort = .y) |> 
      dplyr::group_by(type, cohort) |> 
      dplyr::summarise(min = format(round(min(d),0), big.mark = ","), 
                       max = format(round(max(d),0), big.mark = ","), 
                       mean = format(round(mean(d),0), big.mark = ","),
                       sd = format(round(sd(d),0), big.mark = ","), .groups = 'drop')
  }) |> do.call(what = rbind) |> 
    dplyr::mutate(migration = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
    dplyr::select(-min, -max, -mean, -sd)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  move.df <- dplyr::left_join(dist.df, tot.df, by = c("type", "cohort"))
  print(knitr::kable(move.df, format = "simple"))
  cat("\n")
  
  cat("=============================================================\n")
  cat("HABITAT USE\n")
  cat("=============================================================\n\n")
  
  # Number of individuals visiting each region 
  
  indreg.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::enframe(.x$n$indreg) |> 
      dplyr::mutate(value = unlist(value)) |> 
      dplyr::rename(n = value, region = name) |> 
      dplyr::mutate(cohort = .y)
  }) |> do.call(what = rbind)
  
  if(percent){
    indreg.df <- indreg.df |> 
    dplyr::mutate(perc = round(100 * n / n.ind)) |> 
    dplyr::select(-n) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = perc) |> 
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  } else {
    indreg.df <- indreg.df |> 
      tidyr::pivot_wider(names_from = cohort, values_from = n) 
  }
  
  # Days spent in each region
  
  days.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    
    purrr::map(.x = .x$n$dayreg, .f = ~{
      if(all(.x == 0)) 0 else .x[.x >0]
    }) |> tibble::enframe() |> 
      dplyr::rename(region = name, days = value) |> 
      dplyr::group_by(region) |> 
      dplyr::mutate(cohort = .y) |> 
      dplyr::group_by(region, cohort) |> 
      dplyr::summarise(min = min(unlist(days)),
                       max = max(unlist(days)),
                       mean = round(mean(unlist(days)),1),
                       sd = round(sd(unlist(days)), 1), .groups = 'drop')
  }) |> do.call(what = rbind) |> 
    dplyr::mutate(days = paste0(mean, " (± ", sd, ") [", min, " – ", max, "]")) |> 
    dplyr::select(-min, -max, -mean, -sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = days) 

  
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
  # 
  # Number of regions per individual
  
  regind.df <- suppressWarnings(purrr::map2(.x = out, .y = cohort.ab, .f = ~{
      tmp <- .x$n$regind |> 
      janitor::tabyl() |> dplyr::mutate(percent = 100 * percent)
      names(tmp)[1] <- "No.regions"
      tmp |>  dplyr::mutate(cohort = .y)
      })
      ) |> do.call(what = rbind)
  
  if(percent){
    regind.df <- regind.df |> 
      dplyr::select(-n) |> 
      tidyr::pivot_wider(names_from = cohort, values_from = percent) |> 
      janitor::adorn_totals() |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ paste(.x, "%")))
  } else {
    regind.df <- regind.df |> 
      dplyr::select(-percent) |> 
      tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
      janitor::adorn_totals()
  }
  
  regind.df <- format_table(regind.df)
  

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
  print(knitr::kable(indreg.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Days spent in each region +++++++++++")
  print(knitr::kable(days.df, format = "simple"))
  cat("\n")
  
  cat("+++++++++++ Total number of regions visited +++++++++++")
  print(knitr::kable(regind.df, format = "simple"))
  
 
  cat("\n\n")
  cat("=============================================================\n")
  cat("ACTIVITY BUDGETS\n")
  cat("=============================================================\n")
  
  activ <- purrr::set_names(cohort.ab) |> 
    purrr::map(.f = ~{
      lapply(X = seq_len(n.ind), FUN = function(i){
        if(cohort.id %in% c(4,5)){
          tmp <- cbind(obj[["locs"]][[.x]][, "region", i, drop = FALSE], 
                       obj[["sim"]][[.x]][[1]][, c("alive", "t_travel", "t_feed", "t_rest", "t_nurse", "feed"), i])
        } else {
          tmp <- cbind(obj[["locs"]][[.x]][, "region", i, drop = FALSE], 
                       obj[["sim"]][[.x]][, c("alive", "t_travel", "t_feed", "t_rest", "t_nurse", "feed"), i])
        }
        colnames(tmp)[1] <- "region"
        tmp
      }) |> do.call(what = rbind)
    })
  
  activ.df <- tibble::enframe(activ) |> 
    dplyr::rename(cohort = name, data = value) |> 
    dplyr::mutate(data = purrr::map(.x = data, .f = ~{
      tibble::as_tibble(.x)
    })) |> tidyr::unnest(cols = c(data)) |> 
    dplyr::mutate(region = sort(regions$region)[region]) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(total = sum(t_travel, t_feed, t_nurse, t_rest)) |> 
      dplyr::ungroup() |> 
    dplyr::filter(alive == 1)
  
  cat("Total:", mean(activ.df$total), "hrs\n") 
  
  # activ.df |> dplyr::mutate(ind = 1:nrow(activ.df)) |> dplyr::filter(total < 24)
  
  # which(activ.df$total < 24)
  
  if(mean(activ.df$total) < 24) warning("Inconsistent time allocation in activity budgets")

  
  travel.df <- activ.df |> dplyr::group_by(cohort, region) |> 
    dplyr::summarise(mean = mean(t_travel), sd = sd(t_travel), .groups = "drop") |> 
    dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |> 
    dplyr::select(-mean, -sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = value) |> 
    dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |> 
    dplyr::arrange(region)
  
  feed.df <- activ.df |> dplyr::filter(feed == 1) |> dplyr::select(-feed) |> dplyr::group_by(cohort, region) |> 
    dplyr::summarise(mean = mean(t_feed), sd = sd(t_feed), .groups = "drop") |> 
    dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |> 
    dplyr::select(-mean, -sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = value) |> 
    dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |> 
    dplyr::arrange(region)
  
  if(ncol(feed.df) == 1) feed.df[, cohort.ab] <- 0
  
  rest.df <- activ.df |> dplyr::group_by(cohort, region) |> 
    dplyr::summarise(mean = mean(t_rest), sd = sd(t_rest), .groups = "drop") |> 
    dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |> 
    dplyr::select(-mean, -sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = value) |> 
    dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |> 
    dplyr::arrange(region)
  
  nurse.df <- activ.df |> dplyr::group_by(cohort, region) |> 
    dplyr::summarise(mean = mean(t_nurse), sd = sd(t_nurse), .groups = "drop") |> 
    dplyr::mutate(value = paste0(round(mean,2), " (± ", round(sd,2), ")")) |> 
    dplyr::select(-mean, -sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = value) |> 
    dplyr::right_join(y = tibble::tibble(region = regions$region), by = "region") |> 
    dplyr::arrange(region)
  
  activ.list <- purrr::set_names(cohort.ab) |>  purrr::map(.f = ~{
    tmp <- dplyr::bind_cols(travel.df[, c("region", .x)], feed.df[, .x], rest.df[,.x], nurse.df[,.x], .name_repair = "minimal")
    names(tmp)[2:5] <- c("travel (hrs)", "feed (hrs)", "rest (hrs)", "nurse (hrs)")
    tmp
  })
  
  purrr::walk2(.x = cohort.ab, .y = cohort.names, .f = ~{
    cat("\n+++++++++++ ", .y, " +++++++++++")
    print(knitr::kable(activ.list[[.x]], format = "simple"))
  })
  
  
  activ.1 <- activ.df |>
    dplyr::filter(feed == 1) |> 
    dplyr::select(-feed) |> 
    tidyr::pivot_longer(!c("region", "cohort"), names_to = "behav", values_to = "hrs") |> 
    dplyr::mutate(behav = gsub("t_", "", behav)) |> 
    dplyr::right_join(y = tibble::as_tibble(expand.grid(cohort = cohort.ab, region = regions$region, behav = c("travel", "feed", "rest", "nurse"))), by = c("region","cohort", "behav")) |> 
    dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |> 
    dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region)) |> 
    dplyr::mutate(behav = stringr::str_to_upper(behav)) |> 
    dplyr::filter(behav == "FEED")
  
  activ.2 <- activ.df |> 
    dplyr::select(-feed) |> 
    tidyr::pivot_longer(!c("region", "cohort"), names_to = "behav", values_to = "hrs") |> 
    dplyr::mutate(behav = gsub("t_", "", behav)) |> 
    dplyr::right_join(y = tibble::as_tibble(expand.grid(cohort = cohort.ab, region = regions$region, behav = c("travel", "feed", "rest", "nurse"))), by = c("region","cohort", "behav")) |> 
    dplyr::mutate(region = ifelse(grepl("BOF", region), "BOF", region)) |> 
    dplyr::mutate(region = ifelse(grepl("CABOT", region), "SCOS", region)) |> 
    dplyr::mutate(behav = stringr::str_to_upper(behav)) |> 
    dplyr::filter(!behav == "FEED")
  
  activ.out <- rbind(activ.1, activ.2)
  
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
