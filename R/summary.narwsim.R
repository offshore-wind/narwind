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
                            do.plot = FALSE,
                            show.legend = FALSE,
                            geodesic = FALSE){
  
  # obj = m
  # do.plot = FALSE
  # show.legend = FALSE
  # geodesic = FALSE
  
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
  
  cohorts <- obj$param$cohorts
  cohort.names <- obj$param$cohort.names
  cohort.ab <- unname(obj$param$cohort.ab)
  
  if(show.legend){
    cat("BOF: Bay of Fundy\n")
    cat("CABOT: Cabot Strait\n")
    cat("CCB: Cape Cod Bay\n")
    cat("GOM: Gulf of Maine and Georges Bank\n")
    cat("GSL: Gulf of St Lawrence\n")
    cat("MIDA: Mid-Atlantic\n")
    cat("SCOS: Scotian Shelf\n")
    cat("SEUS: South-east United States\n")
    cat("SNE: Southern New England\n\n")
  }
  
  # Dynamic scoping to make sure geodesic calculations are fast
  environment(summary_geo) <- .GlobalEnv
  
  # ............................................................
  # Load spatial support
  # ............................................................
  
  coords <- sp::coordinates(density_support)
  colnames(coords) <- c("x", "y")
  map_limits <- c(range(coords[,1]), range(coords[,2]))
  map_resolution <- raster::res(raster::raster(density_support))
  geomap <- t(1*raster::as.matrix(density_support))
  
  # ............................................................
  # Calculate distances travelled and tally regions visited
  # ............................................................
  
  out <- purrr::set_names(cohort.ab) |> 
    purrr::map2(.y = cohort.names, .f = ~{
      
      # ............................................................
      # Calculate Euclidean distances between pairs of sequential points
      # ............................................................
      
      dist.euclid <- lapply(seq(n.ind), function(x){
        obj.ind <- locations[[.x]][,,x]
        obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
        dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
        dd <- dd[!is.na(dd)]
        unname(dd)
      }) 
      
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
      
      list(dist = list(euclid = dist.euclid, geo = dist.geo),
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
  
  if(do.plot){
    par(mfrow = c(1,2))
    boxplot(unlist(dist.euclid), main = "Euclidean distances (daily, km)")
    boxplot(unlist(dist.geo), main = "Approx. geodesic distances (daily, km)")
    par(mfrow = c(1,1))
  }
  
  locs.by.region <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::enframe(.x$n$locreg) |> 
      dplyr::mutate(value = unlist(value)) |> 
      dplyr::rename(region = name, n = value) |>
      dplyr::mutate(cohort = .y) |> 
      dplyr::relocate(cohort, .after = region)
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    dplyr::arrange(region) |>
    janitor::adorn_totals()
  
  locs.by.region <- format_table(locs.by.region)
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================\n")
  
  print(knitr::kable(locs.by.region, format = "simple"))
  
  locs.by.country <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::enframe(.x$n$locctry) |> 
      dplyr::mutate(value = unlist(value)) |> 
      dplyr::rename(country = name, n = value) |>
      dplyr::mutate(cohort = .y)
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) |> 
    dplyr::arrange(country) |>
    janitor::adorn_totals()
  
  locs.by.country <- format_table(locs.by.country)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(locs.by.country, format = "simple"))
  cat("\n")
  
  cat("=============================================================\n")
  cat("DAILY MOVEMENTS (km)\n")
  cat("=============================================================")
  
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
    dplyr::relocate(type, .before = cohort)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(dist.df, format = "simple"))
  cat("\n")
  
  cat("=============================================================\n")
  cat("MIGRATORY MOVEMENTS (km)\n")
  cat("=============================================================\n")
  
  tot.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    tibble::tibble(d = c(.x$totdist$euclid, .x$totdist$geo), 
                   type = rep(dist.type, each = length(.x$totdist$euclid))) |> 
      dplyr::mutate(cohort = .y) |> 
      dplyr::group_by(type, cohort) |> 
      dplyr::summarise(min = format(round(min(d),0), big.mark = ","), 
                       max = format(round(max(d),0), big.mark = ","), 
                       mean = format(round(mean(d),0), big.mark = ","),
                       sd = format(round(sd(d),0), big.mark = ","), .groups = 'drop')
  }) |> do.call(what = rbind)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(tot.df, format = "simple"))
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
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = n) 
  
  # Days spent in each region
  
  days.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    
    purrr::map(.x = .x$n$dayreg , .f = ~{
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
  }) |> do.call(what = rbind)
  
  days.min <- days.df |> dplyr::select(region, cohort, min) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = min)
  
  days.max <- days.df |> dplyr::select(region, cohort, max) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = max) 
  
  days.mean <- days.df |> dplyr::select(region, cohort, mean) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = mean) 
  
  days.sd <- days.df |> dplyr::select(region, cohort, sd) |> 
    tidyr::pivot_wider(names_from = cohort, values_from = sd) 
  
  # Number of regions per individual
  
  regind.df <- purrr::map2(.x = out, .y = cohort.ab, .f = ~{
    fivenum(.x$n$regind)[c(1,3,5)] |> 
      tibble::as_tibble(.name_repair = "minimal") |> 
      dplyr::mutate(param = c("min", "mean", "max")) |> 
      dplyr::rename(n = value) |> 
      dplyr::mutate(cohort = .y)
  }) |> do.call(what = rbind) |> 
    tidyr::pivot_wider(names_from = param, values_from = n) 
  
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  cat("** Number of animals visiting each region -------------------")
  print(knitr::kable(indreg.df, format = "simple"))
  cat("\n")
  
  cat("** Days spent in each region -------------------")
  cat("\n\n++ Minimum")
  print(knitr::kable(days.min, format = "simple"))
  cat("\n++ Maximum")
  print(knitr::kable(days.max, format = "simple"))
  cat("\n++ Mean")
  print(knitr::kable(days.mean, format = "simple"))
  cat("\n++ Standard deviation")
  print(knitr::kable(days.sd, format = "simple"))
  cat("\n")
  
  cat("** Number of regions visited by each animal -------------------")
  print(knitr::kable(regind.df, format = "simple"))
  
}
