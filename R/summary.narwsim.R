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

summary.narwsim <- function(obj, do.plot = FALSE, show.legend = FALSE, geodesic = FALSE){
  
  if(!"narwsim" %in% class(obj)) stop("Input must be of class <narwsim>")
  
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n")
  cat("\n")
  cat("     NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("                --- MODEL SUMMARY ---\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n\n")
  
  n.ind <- dim(obj)[3]
  
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
  
  geomap <- targets::tar_read(density_support)
  coords <- sp::coordinates(geomap)
  colnames(coords) <- c("x", "y")
  map_limits <- c(range(coords[,1]), range(coords[,2]))
  map_resolution <- raster::res(raster::raster(geomap))
  
  geomap <- geomap |> raster::as.matrix()
  geomap <- t(1*geomap)
  
  # ............................................................
  # Load regions
  # ............................................................
  
  regions <- targets::tar_read(regions)

  # ............................................................
  # Calculate Euclidean distances between pairs of sequential points
  # ............................................................
  
  dist.euclid <- lapply(seq(n.ind), function(x){
    obj.ind <- obj[,,x]
    obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
    dd <- sqrt((obj.ind[, 3]-obj.ind[, 1])^2+(obj.ind[, 4]-obj.ind[, 2])^2)
    dd <- dd[!is.na(dd)]
    unname(dd)
  }) 
  
  # D <- lapply(seq(dim(obj)[3]), function(x){
  #   obj.ind <- obj[,,x]
  #   cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))})
  
  # ............................................................
  # Calculate geodesic distances between pairs of sequential points
  # ............................................................

  if(geodesic){
    dist.geo <- summary_geo(obj, map_limits, map_resolution, geomap)
  } else {
    dist.geo <- list(NA)
  }

  # ............................................................
  # Total migration distances
  # ............................................................
  
  tot.dist.euclid <- purrr::map_dbl(.x = dist.euclid, .f = ~sum(.x))
  if(geodesic) tot.dist.geo <- purrr::map_dbl(.x = dist.geo, .f = ~sum(.x)) else tot.dist.geo <- NULL
  
  # ............................................................
  # Tally points by region
  # ............................................................
 
  #  obj.sf <- lapply(seq(dim(obj)[3]), function(x){
  #   tmp <- sf::st_as_sf(x = data.frame(obj[ , c("easting", "northing"), x]), coords = c("easting", "northing"), crs = narw_crs())
  #   sf::st_agr(tmp) <- "constant"
  #   tmp
  # })
   
   obj.sp <- lapply(seq(dim(obj)[3]), function(x){
     sp::SpatialPoints(coords = obj[ , c("easting", "northing"), x], proj4string = narw_crs())
   })
  
  locations.per.region <- sp::over(regions, do.call(rbind, obj.sp), returnList = TRUE) |> 
    purrr::set_names(nm = regions$region)
   
  # sf::st_agr(regions) <- "constant"
  
  pr.ind <- purrr::map(.x = obj.sp, .f = ~sp::over(regions, .x, returnList = TRUE)) |> 
    purrr::map(.f = ~purrr::set_names(x = .x, nm = regions$region))
  regions.per.ind <- purrr::map_dbl(.x = pr.ind, .f = ~length(purrr::discard(.x = .x, .p = function(x) length(x)==0)))
  # pr.ind <- purrr::map(.x = obj.sf, .f = ~sf::st_intersection(regions, .x)$region)
 
  # regions.per.ind <- purrr::map_dbl(.x = pr.ind, .f = ~length(unique(.x)))
  
  days.per.regions <- purrr::map(sort(regions$region), .f = ~{
    sapply(X = pr.ind, FUN = function(x) length(x[[.x]]))
  }) |> purrr::set_names(sort(regions$region))
  
  # days.per.regions <- purrr::set_names(regions$region) |> 
  #   purrr::map(.f = ~{sapply(X = pr.ind, function(i) sum(grepl(pattern = .x, x = i)))})
  
  ind.per.regions <- purrr::map(.x = days.per.regions, .f = ~as.numeric(sum(.x>0)))
    
  # locations.per.region <- purrr::set_names(sort(regions$region)) |> 
  #   purrr::map(.f = ~sf::st_intersection(regions[regions$region == .x,], do.call(what = rbind, obj.sf))$region)
  
  locations.per.country <- list("Canada" = length(unlist(locations.per.region[c("GSL", "SCOS", "CABOT", "BOF_lower", "BOF_upper")])),
                                "U.S" = length(unlist(locations.per.region[c("CCB", "MIDA", "SNE", "SEUS", "GOM")])))
  
  # locations.per.country <- 
  #   purrr::map(.x = locations.per.region, 
  # .f = ~unname(sapply(X = .x, FUN = function(x) ifelse(x %in% c("GSL", "SCOS", "CABOT", "BOF_lower", "BOF_upper"), "Canada", "U.S."))))
  
  # locations.per.region <- sp::over(regions, do.call(rbind, obj.sp), returnList = TRUE) |> purrr::map(.f = ~length(.x)) |> 
  #   purrr::set_names(nm = regions$region)
  
  locations.per.region <- purrr::map(.x = locations.per.region, .f = ~length(.x))
  # locations.per.country <- unname(unlist(locations.per.country))
  
  
  # ............................................................
  # OUTPUTS
  # ............................................................
  
  if(do.plot){
    par(mfrow = c(1,2))
    boxplot(unlist(dist.euclid), main = "Euclidean distances (daily, km)")
    boxplot(unlist(dist.geo), main = "Approx. geodesic distances (daily, km)")
    par(mfrow = c(1,1))
  }
  
  locs.by.region <- tibble::enframe(locations.per.region) |> 
    dplyr::mutate(value = unlist(value)) |> 
    dplyr::rename(region = name, n = value) |>
    dplyr::mutate(percent = 100*n/sum(n)) |> 
    janitor::adorn_totals() |>
    dplyr::mutate(n = formatC(n, big.mark = ","), percent = formatC(percent, digits = 1, format = "f")) 

  locs.by.region <- format_table(locs.by.region)
  
  cat("=============================================================\n")
  cat("LOCATIONS\n")
  cat("=============================================================\n")
  
  print(knitr::kable(locs.by.region, format = "simple"))

  locs.by.country <- tibble::enframe(locations.per.country) |> 
    dplyr::mutate(value = unlist(value)) |> 
    dplyr::rename(country = name, n = value) |>
    dplyr::mutate(percent = 100 * n/sum(n)) |> 
    janitor::adorn_totals() |>
    dplyr::mutate(n = formatC(n, big.mark = ","), percent = formatC(percent, digits = 1, format = "f")) 
  
  locs.by.country <- format_table(locs.by.country)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(locs.by.country, format = "simple"))
  # print.data.frame(locs.by.country, right = F, row.names = F)
  cat("\n")
  
  cat("=============================================================\n")
  cat("DAILY MOVEMENTS (km)\n")
  cat("=============================================================")
  
  dist.df <- fivenum(unlist(dist.euclid))[c(1,3,5)] |> 
    tibble::as_tibble(.name_repair = "minimal") |> 
    dplyr::mutate(param = c("min", "mean", "max")) |> 
    dplyr::rename(distance = value) |> 
    dplyr::mutate(distance = round(distance, 1)) |> 
    tidyr::pivot_wider(names_from = param, values_from = distance) |> 
    dplyr::relocate(mean, .after = max)
  
  geodesic.df <- fivenum(unlist(dist.geo))[c(1,3,5)] |> 
    tibble::as_tibble(.name_repair = "minimal") |> 
    dplyr::mutate(param = c("min", "mean", "max")) |> 
    dplyr::rename(distance = value) |> 
    dplyr::mutate(distance = round(distance, 1)) |> 
    tidyr::pivot_wider(names_from = param, values_from = distance) |> 
    dplyr::relocate(mean, .after = max)
  
  dist.type <- "euclidean"
  if(geodesic){
    dist.type <-  c(dist.type, "geodesic")
    dist.df <- rbind(dist.df, geodesic.df)
  }
    
  dist.df <- dist.df |> 
    dplyr::mutate(type = dist.type) |> 
    dplyr::relocate(type, .before = min)
  
  # dist.df <- format_table(dist.df, bottom = FALSE)

  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(dist.df, format = "simple"))
  # print.data.frame(dist.df, right = F, row.names = F)
  cat("\n")
  
  cat("=============================================================\n")
  cat("MIGRATORY MOVEMENTS (km)\n")
  cat("=============================================================\n")
  
  tot.df <- tibble::tibble(d = c(tot.dist.euclid, tot.dist.geo), 
                           type = rep(dist.type, each = length(tot.dist.euclid))) |> dplyr::group_by(type) |> 
    dplyr::summarise(min = format(round(min(d),0), big.mark = ","), 
                     max = format(round(max(d),0), big.mark = ","), 
                     mean = format(round(mean(d),0), big.mark = ","),
                     sd = format(round(sd(d),0), big.mark = ",")) |> 
    as.data.frame()
  
  # tot.df <- format_table(tot.df, bottom = FALSE)
  
  # ............................................................
  # Print outputs
  # ............................................................
  
  print(knitr::kable(tot.df, format = "simple"))
  # print.data.frame(tot.df, right = F, row.names = F)
  cat("\n")
  
  cat("=============================================================\n")
  cat("HABITAT USE\n")
  cat("=============================================================\n\n")
  
  # Number of individuals visiting each region 
  
  indreg.df <- tibble::enframe(ind.per.regions) |> 
    dplyr::mutate(value = unlist(value)) |> 
    dplyr::rename(n = value, region = name) |> 
    dplyr::mutate(percent = round(100 * n / n.ind), 1)
  
  # indreg.df <- format_table(indreg.df)
  
  # Days spent in each region
  
  days.df <- purrr::map(.x = days.per.regions, .f = ~{
    if(all(.x == 0)){
      0
    } else {.x[.x >0]}
  }) |> tibble::enframe() |> 
    dplyr::rename(region = name, days = value) |> 
    dplyr::group_by(region) |> 
    dplyr::summarise(min = min(unlist(days)),
                     max = max(unlist(days)),
                     mean = round(mean(unlist(days)),1),
                     sd = round(sd(unlist(days)), 1))
  
  days.df <- format_table(days.df)
  
  # Number of regions per individual
  
  regind.df <- fivenum(regions.per.ind)[c(1,3,5)] |> 
    tibble::as_tibble(.name_repair = "minimal") |> 
    dplyr::mutate(param = c("min", "mean", "max")) |> 
    dplyr::rename(n = value) |> 
    tidyr::pivot_wider(names_from = param, values_from = n) |> 
    dplyr::relocate(mean, .after = max)
  
  # regind.df <- format_table(regind.df, bottom = FALSE)

  # ............................................................
  # Print outputs
  # ............................................................
  
  cat("** Number of animals visiting each region -------------------")
  print(knitr::kable(indreg.df, format = "simple"))
  cat("\n")
  
  cat("** Days spent in each region -------------------")
  print(knitr::kable(days.df, format = "simple"))
  cat("\n")
  
  cat("** Number of regions visited by each animal -------------------")
  print(knitr::kable(regind.df, format = "simple"))
  
}
