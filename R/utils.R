# MODEL ------------------------------------------------------

#' Initialize model
#'
#' @export
initialize_model <- function(){
  Rcpp::sourceCpp("src/simtools.cpp")
  source("R/run_model.R")
  assign("init.model", value = TRUE, envir = .GlobalEnv)
}

# SPATIAL DATA ------------------------------------------------------
  
#' Projected coordinate system
#'
#' @return An object of class \code{CRS}.
#' 
narw_crs <- function(){
  sp::CRS("+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
}

predict_hgam <- function(object, newdata, type, n.chunks = 5){
  
  N <- nrow(newdata)
  n <- ceiling(N/n.chunks)
  new.data <- split(1:N, ceiling(seq_along(1:N)/n)) |> 
    purrr::map(.f = ~newdata[.x,])
  
  purrr::map(.x = new.data, .f = ~predict(object, newdata = .x, type = type)) |>
    do.call(what = c)
  
}

support_as_polygon <- function(){
  support_grd <- targets::tar_read(density_support) |> raster::raster()
  support_grd[support_grd == 0, ] <- NA
  support_grd <- raster::rasterToPolygons(support_grd)
  support_grd <- rgeos::gUnionCascaded(support_grd)
  support_grd
}

spatial_support <- function(min_density = 0.001){
  
  # ..........................................
  # Density surfaces
  # ..........................................
  
  map_rasters <- targets::tar_read(density_merge) |> 
    purrr::map(.f = ~as(.x, "SpatialGridDataFrame"))
  
  # ..........................................
  # Spatial support
  # ..........................................
  
  # Identify locations that exceed a minimum density value in some raster
  filtered_support <- map_rasters[[1]]
  filtered_support$Nhat <- FALSE
  for(m in map_rasters) {
    m$Nhat[is.na(m$Nhat)] <- 0
    filtered_support$Nhat <- filtered_support$Nhat | (m$Nhat > min_density)
  }
  
  # Identify connected components within the filtered support
  coords <- sp::coordinates(filtered_support)
  colnames(coords) <- c("x", "y")
  rowseq <- sort(unique(coords[,"x"]))
  colseq <- sort(unique(coords[,"y"]))
  fs <- spatstat.geom::im(mat = as.matrix(filtered_support), xcol = colseq, yrow = rowseq)
  components <- spatstat.geom::connected(fs)
  
  # Use the largest component to define the spatial mask
  filtered_support$Nhat <- as.numeric(components$v)
  largest_component <- which.max(table(filtered_support$Nhat))
  filtered_support$Nhat <- filtered_support$Nhat == largest_component
  filtered_support$Nhat[!is.finite(filtered_support$Nhat)] <- FALSE
  
  # ..........................................
  # Manually add support for Canada
  # ..........................................
  
  # Load region boundaries
  regions.support <- targets::tar_read(regions)
  
  # Filter out Canada
  regions.support <- regions.support[regions.support@data$region %in% c("GSL", "SCOS"),] 
  
  # Overlay regions onto current spatial support
  regions.raster <- raster::rasterize(
    x = regions.support,
    y = raster::raster(x = regions.support,
                       resolution = filtered_support@grid@cellsize,
                       origin = raster::origin(raster::raster(filtered_support))))
  
  regions.raster[!is.na(regions.raster)] <- 1
  
  combined_support <- raster::merge(regions.raster, raster::raster(filtered_support))
  
  # ..........................................
  # Mask out land
  # ..........................................
  
  wrld <- targets::tar_read(world)
  
  outlines.raster <- raster::rasterize(
    x = wrld, 
    y = raster::raster(x = wrld, 
                       resolution = filtered_support@grid@cellsize, 
                       origin = raster::origin(raster::raster(combined_support))))
  
  # Mark land regions with a mask value
  outlines.raster@data@values[outlines.raster@data@values > 0] <- -1
  outlines.raster <- raster::crop(outlines.raster, combined_support)
  
  # Combine data to form mask
  # Remove/mask land regions and non-labeled regions from spatial support
  cleaned_support <- merge(outlines.raster, combined_support)
  cleaned_support <- raster::crop(cleaned_support, raster::extent(filtered_support))
  cleaned_support[cleaned_support <= 0] <- 0
  
  # Convert to logical
  cleaned_support@data@values <- cleaned_support@data@values > 0
  names(cleaned_support) <- "Nhat"
  
  # Trim
  cleaned_support <- raster::trim(cleaned_support)
  
  # Convert to SGDF
  as(cleaned_support, "SpatialGridDataFrame")
  
  
  # sightings <- targets::tar_read(sightings)
  # 
  # # Project sighting coordinates
  # sightings.pp = sp::SpatialPoints(
  #   coords = sightings[,c("Longitude","Latitude")], 
  #   proj4string = sp::CRS("+proj=longlat"))
  # 
  # sightings.proj = sp::spTransform(x = sightings.pp, CRSobj = narw_crs())
  # 
  # # Calculate proportion of sightings that lie in support
  # perc.sightings = sp::over(sightings.proj, cleaned_support)
  # cat("Sightings in support:", round(100*mean(is.finite(unlist(perc.sightings))),1), "%")
  # 
  # # Make simulation support plottable
  # ddf = data.frame(cleaned_support)
  # 
  # # Build spatial support plot
  # pl <- ggplot2::ggplot() + 
  #   # Simulation support
  #   ggplot2::geom_raster(
  #     mapping = ggplot2::aes(x = s1, y = s2, fill = layer),
  #     data = ddf
  #   ) + 
  #   # Sightings
  #   ggplot2::geom_point(
  #     mapping = ggplot2::aes(x = Longitude, y = Latitude),
  #     data = as.data.frame(sightings.proj),
  #     size = .001
  #   ) + 
  #   # Land
  #   ggplot2::geom_sf(data = world_sf, fill = "lightgrey", color = "black", size = 0.25) +
  #   ggplot2::guides(fill = "none") +
  #   coord_sf(xlim = raster::extent(cleaned_support)[1:2], 
  #            ylim = raster::extent(cleaned_support)[3:4], expand = TRUE) +
  #   xlab("") + ylab("")
  # 
  # print(pl)
  
}

sum_Nhat <- function(df, reg = "GSL", Ntot = 336){
  
  # Ntot is the latest population abundance estimate
  # See NARWC report card 2021.
  # https://www.narwc.org/uploads/1/1/6/6/116623219/2021report_cardfinal.pdf
  
  if(!"month" %in% names(df)) stop("Cannot find <month> variable")
  
  # Assign grid cells to regions
  ppts <- sp::SpatialPoints(coords = df[, c("x", "y")], proj4string = narw_crs())
  myregions <- sp::spTransform(regions, CRSobj = narw_crs())
  df$region <- regions$region[sp::over(sp::geometry(ppts), sp::geometry(myregions), returnList = FALSE)]
  
  # Ntot <- df |> dplyr::group_by(month) |> dplyr::summarise(Ntot = round(sum(Nhat),0))
  
  df.out <- df |> 
    dplyr::group_by(region,month) |> 
    dplyr::summarise(Nhat = round(sum(Nhat),0)) |> 
    dplyr::arrange(region) |> 
    dplyr::filter(!is.na(region)) |> 
    dplyr::mutate(perc = round(100* Nhat / Ntot, 1)) |> 
    dplyr::filter(region == reg) |> 
    dplyr::ungroup()
  
  # df.out <- dplyr::left_join(df, Ntot, by = "month") |> 
  #   dplyr::mutate(perc = round(100* Nhat / Ntot, 1)) |> 
  #   dplyr::filter(region == reg) |> 
  #   dplyr::ungroup()
  
  print(df.out)
  cat(paste0("Average: ", round(mean(df.out[df.out$perc>0,]$perc),1), "%"))
  
}



# geodesic_dist <- function(grid.res = 85){
#   
#   # Density support
#   d.support <- targets::tar_read(density_support) |> raster::raster()
#   d.support[d.support < 1] <- NA
#   d.poly <- raster::rasterToPolygons(d.support) |> rgeos::gUnionCascaded()
#   d.support[is.na(d.support)] <- 3
#   
#   # Create grid
#   r.grid <- sp::makegrid(d.poly, cellsize = grid.res) |>
#     dplyr::mutate(cell = NA) |> 
#     dplyr::rename(x = x1, y = x2)
#   
#   r.grid <- sp::SpatialPointsDataFrame(coords = r.grid[, c("x", "y")], 
#                                        data = r.grid,
#                                        proj4string = narw_crs())
#   
#   r.grid <- r.grid[which(sp::over(r.grid, d.poly) == 1),]
#   
#   r.grid$cell <- purrr::map_dbl(.x = 1:length(r.grid), 
#                                 .f = ~raster::cellFromXY(d.support, r.grid[.x,]@coords))
#   
#   d <- targets::tar_read(density_narw)[[1]] |> 
#     raster::raster() |> raster::mask(mask = d.poly)
#   vor <- dismo::voronoi(r.grid)
#   vr <- raster::intersect(vor, d.poly)
#   vr <- raster::rasterize(vr, d, 'cell')
#   
#   # Calculate geodesic distances from each grid point -- return a raster
#   future::plan(future::multisession)
#   
#   suppressMessages(gd <- furrr::future_map(.x = r.grid$cell,
#                                            .f = ~{
#                                              library(raster)
#                                              r <- d.support
#                                              r[.x] <- 2
#                                              terra::gridDistance(r, origin = 2, omit = 3)}, 
#                                            .options = furrr::furrr_options(seed = TRUE)))
#   
#   # gd <- raster::stack(gd)
#   names(gd) <- r.grid$cell
#   
#   return(list(grid = r.grid, cellno = vr, raster = gd))
#   
# }

plot_preds <- function(dat, month = NULL, do.facet = TRUE, plot.sightings = TRUE, hide.lgd = FALSE, alpha = 1){
  
  if(!is.null(month)){
    do.facet <- FALSE
    dat <- dat[dat$month == month,]
  }
  
  colour.breaks <-  colour_breaks(dat)
  legend.title <- expression(atop(atop(textstyle("Individuals"),
                                       atop(textstyle("per 25 km"^2)))))
  
  dat <- dat |> 
    dplyr::mutate(Ncol = cut(Nhat, breaks = colour.breaks, include.lowest = TRUE)) |> 
    dplyr::mutate(Ncol = factor(Ncol))
  
  levels(dat$Ncol) <-  gsub(pattern = ",", replacement = " – ", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\[|\\]", replacement = "", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\(|\\)", replacement = "", levels(dat$Ncol))
  
  r.dat <- raster::rasterFromXYZ(dat[, c("x", "y", "Nhat")], res = d.res, crs = narw_crs())
  x.lim <- raster::extent(r.dat)[1:2]
  y.lim <- raster::extent(r.dat)[3:4]
  
  # world.df <- ggplot2::fortify(world) |> dplyr::rename(x = long, y = lat)
  
  p <- ggplot2::ggplot(data = dat) +  
    ggplot2::geom_raster(aes(x,y, fill = Ncol)) +
    # ggplot2::geom_polygon(data = world.df, aes(x,y,group = group),
    # fill = "lightgrey", color = "black", size = 0.25) +
    ggplot2::geom_sf(data = world_sf, fill = "lightgrey", color = "black", size = 0.25) +
    ylab("") + xlab("") +
    coord_sf(xlim = x.lim, ylim = y.lim, expand = FALSE) +
    scale_fill_manual(name = legend.title,
                      values = pals::viridis(length(levels(dat$Ncol))),
                      guide = guide_legend(reverse = TRUE))
  
  if(hide.lgd) p <- p + theme(legend.position = "none")
  
  if(plot.sightings){
    
    narwc_df <- narwc.sp@data
    
    if(!is.null(month)) narwc_df <- narwc_df[narwc_df$month == month,]
    
    p <- p + geom_point(data = narwc_df, aes(x, y), alpha = alpha)
  }
  
  if(do.facet) p <- p + facet_wrap(~month)
  return(p)
}

add_country <- function(df){
  ppts <- sp::SpatialPoints(coords = df[, c("x", "y")], proj4string = narw_crs())
  df$region <- "U.S."
  df$region[sp::over(sp::geometry(canada_ocean), sp::geometry(ppts), returnList = TRUE)[[1]]] <- "Canada"
  return(df)
}

add_xy <- function(dat, CRS.obj = narw_crs()){
  
  tmp <- sp::SpatialPoints(
    coords = dat[, c("longitude", "latitude")],
    proj4string = sp::CRS("+proj=longlat")) |>
    sp::spTransform(CRSobj = CRS.obj)
  
  xy <- sp::coordinates(tmp) |> as.data.frame()
  names(xy) <- c("x", "y")
  
  cbind(dat, xy)
  
}

add_latlon <- function(dat, CRS.obj = narw_crs()){
  
  tmp <- sp::SpatialPoints(
    coords = dat[, c("easting", "northing")],
    proj4string = narw_crs()) |>
    sp::spTransform(CRSobj = sp::CRS("+proj=longlat"))
  
  xy <- sp::coordinates(tmp) |> as.data.frame()
  names(xy) <- c("long", "lat")
  
  cbind(dat, xy)
  
}

pointcount = function(r, pts, mask = NULL){
  # Make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  
  # Get the cell index for each point and make a table:
  counts = table(raster::cellFromXY(r, pts))
  
  # Fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  
  if(!is.null(mask)){
    r2 <- raster::crop(r2, raster::extent(mask))
    r2 <- raster::mask(r2, mask)
  }
  
  return(r2)
}

toRaster_area <- function(poly, bufferwidth = NULL){
  
  if(!is.null(bufferwidth)) poly <- terra::buffer(poly, width = bufferwidth)
  
  r <- raster::raster(ext = raster::extent(poly), resolution = d.res, crs = narw_crs())
  r <- suppressWarnings(raster::projectRaster(from = r, to = d[[1]], alignOnly = TRUE))
  r <- raster::setValues(r, 0)
  
  # Crop the raster to Canadian waters
  r <- raster::crop(r, raster::extent(poly))
  r <- raster::mask(r, poly)
  
  # Convert raster cells to polygons to calculate areas
  r.poly <- raster::rasterToPolygons(r)
  r.overlap <- sp::over(r.poly, world, returnList = TRUE)
  r.overlap <- purrr::map_dbl(.x = r.overlap, .f = ~nrow(.x))
  r.overlap <- r.overlap == 1
  r.poly.overlap <- r.poly[r.overlap,]
  
  r.poly.overlap <- 
    rgeos::gIntersection(spgeom1 = poly, spgeom2 = r.poly.overlap, byid = TRUE, drop_lower_td = TRUE)
  
  r.area <- rgeos::gArea(r.poly, byid = TRUE)
  r.area[r.overlap] <- rgeos::gArea(r.poly.overlap, byid = TRUE)
  r.poly$area <- r.area
  r.poly$layer <- NULL
  
  r <- raster::rasterize(r.poly, r, field = "area")
  names(r) <- c("area")
  r$d <- 0
  r$d[is.na(r$area)] <- NA
  
  r
}

#' Clip and fill input raster
#'
#' Fill NA values within a mask region
#' 
#' @param x Source raster to fill/pad.
#' @param mask Extent of the spatial region to fill.
#' @param fill Value used for filling.
#' @param offset Additional value to add to all locations within mask.
#' @param clamp Logical. Set to TRUE to give all locations outside mask region an NA value.

fill_mask <- function(x, mask, fill, offset = 0, clamp = FALSE) {

  if(!is.logical(mask$Nhat)){
    stop('mask raster must have logical entries')
  }
  to_fill <- is.na(x$Nhat[mask$Nhat])
  x$Nhat[mask$Nhat][to_fill] <- fill
  x$Nhat[mask$Nhat] <- x$Nhat[mask$Nhat] + offset
  if(clamp) {
    x$Nhat[!mask$Nhat] = NA
  }
  return(x)
}

#' Title
#'
#' @export

clip_density <- function(){
  
  # Load density surfaces
  d <- targets::tar_read(density_merge)
  density.support <- targets::tar_read(density_support)
  
  # Apply spatial mask

  map_spgrd <- lapply(d, function(x) {
    dens <- as(x, "SpatialGridDataFrame")
    fill_mask(x = dens, mask = density.support, fill = 1e-6, offset = 1e-6, clamp = TRUE)
  })
  
  return(map_spgrd)
}

regions_matrix <- function(){
  regions <- targets::tar_read(regions)
  support <- targets::tar_read(density_support)
  regions.m <- raster::rasterize(regions, raster::raster(support)) |> as("SpatialGridDataFrame")
  regions.m$regionID <- as.numeric(factor(regions.m$region))
  regions.m$region <- NULL
  regions.m
}

#' Download density maps
#'
#' Download density maps from the Duke University website and unzip files
#'
#' @param version Version of the North Atlantic right whale model
#'
#' @export

get_density <- function(region = NULL, version = 12, start = 2010, end = 2019) {
  
  if(is.null(region)) stop("Missing <region> argument.")
  if(!region %in% c("US", "CA")) stop("Unrecognized <region> argument.")
  
  if(region == "US"){
    
    # Avoid time out issues
    options(timeout = 1000)
    version <- as.character(version)
    
    # Density maps published at:
    # https://seamap.env.duke.edu/models/Duke/EC/EC_North_Atlantic_right_whale_history.html
    
    # Direct link to dataset
    src <- paste0("https://seamap.env.duke.edu/seamap-models-files/Duke/EC/North_Atlantic_right_whale/v", version, "/EC_North_Atlantic_right_whale_v", version, ".zip")
    
    # Create data directory
    f <- file.path("data", "densitymaps")
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    # Download and unzip data
    dst <- file.path(f, basename(src))
    utils::download.file(url = src, destfile = dst)
    unzip(zipfile = file.path(f, basename(src)), exdir = f)
    
    # Remove zip file
    unlink(dst)
    
    # Delete unnecessary files and re-arrange folder
    unlink("data/densitymaps/Animations", recursive = TRUE)
    unlink("data/densitymaps/ArcGIS_Symbology", recursive = TRUE)
    unlink("data/densitymaps/Maps", recursive = TRUE)
    unlink(list.files(paste0("data/densitymaps/Rasters/", paste0(start, "-", end)), full.names = TRUE, pattern = "aux.xml"))
    invisible(file.copy(list.files(paste0("data/densitymaps/Rasters/", paste0(start, "-", end)), full.names = TRUE, pattern = "density_month"), "data/densitymaps/", recursive = TRUE))
    unlink(list.files("data/densitymaps", pattern = ".pdf", full.names = TRUE))
    unlink(list.files("data/densitymaps", pattern = ".txt", full.names = TRUE))
    unlink("data/densitymaps/Rasters/", recursive = TRUE)
  
  } else if(region == "CA"){
    
    # ........................................................
    # Spatial data ---
    # ........................................................
    
    # Import one of Jason's density surfaces
    d <- list.files("data/densitymaps", pattern = ".img")
    d <- purrr::map(.x = d, .f = ~{
      d.file <- raster::raster(file.path("data/densitymaps", .x))
      d.file <- raster::projectRaster(from = d.file, crs = narw_crs())
      25*d.file/100
    }) |> purrr::set_names(nm = month.abb)
    
    # Import and filter NARWC sightings
    # NARW only started using the GoSL as of 2015.
    # See Simard et al. (2019) and Meyer-Gutbrod et al. (2022)
    # There are only three sightings of NARW in 2021 in the current version of the database,
    # all of which are in U.S. waters.
    
    narwc <- targets::tar_read(sightings) |>
      janitor::clean_names() |>
      dplyr::filter((longitude > -85 & longitude < -44) & (latitude < 58 | latitude > 55)) |>
      dplyr::filter(sighting_year >= 2015 & sighting_year < 2021) |>
      dplyr::filter(!is.na(eg_no)) |> 
      dplyr::mutate(sighting_time = stringr::str_pad(sighting_time, 4, pad = "0")) |> 
      dplyr::mutate(sighting_time = paste0(stringr::str_sub(sighting_time, 1, 2), ":", stringr::str_sub(sighting_time, 3, 4))) |> 
      dplyr::mutate(date = lubridate::dmy_hm(
        paste0(sighting_day, "/", sighting_month,"/", sighting_year, " ", sighting_time))) |> 
      dplyr::mutate(month = factor(month.abb[sighting_month], levels = month.abb))
    
    # Add projected coordinates
    narwc <- add_xy(narwc)
    
    # We also remove potential duplicate sightings to avoid getting unrealistically high
    # densities in some grid cells.
    
    # Split the data by eg_no, year, and month
    narwc <- narwc |> dplyr::mutate(comb = glue::glue("{sighting_year}-{sighting_month}-{eg_no}"))
    
    narwc.thin <- split(x = narwc, f = factor(narwc$comb)) |> 
      purrr::map(.f = ~colMeans(.x[,c("x", "y")])) |> 
      do.call(what = rbind) |> 
      data.frame() |> 
      tibble::rownames_to_column(var = "comb")
    
    narwc <- narwc |> dplyr::select(-longitude, -latitude, -x, -y) |> 
      dplyr::left_join(y = narwc.thin, by = "comb") |> 
      dplyr::distinct(comb, .keep_all = TRUE) |> 
      tibble::as_tibble() |> 
      dplyr::select(-comb)
    
    narwc.sp <- sp::SpatialPointsDataFrame(
      coords = narwc[, c("x", "y")],
      data = narwc,
      proj4string = narw_crs())
    
    # Import regions and project to CRS
    regions <- targets::tar_read(regions) |> 
      sp::spTransform(CRSobj = narw_crs())
    
    # Using a zero-width buffer cleans up many topology problems
    regions <- rgeos::gBuffer(regions, byid = TRUE, width = 0)
    
    # Import shapefile for Canada and project to CRS
    canada <- targets::tar_read(canada) |> 
      sp::spTransform(CRSobj = narw_crs())
    
    # Import and crop coastline shapefile
    world <- targets::tar_read(world)
    
    # Extent of Canadian waters
    canada_ocean <- raster::erase(canada, world)
    
    # Identify sightings made in Canadian waters
    narwc <- add_country(narwc)
    narwc.sp@data <- narwc
    
    # Remove sightings on land and sightings outside the study area
    onland <- rgeos::gIntersects(world, narwc.sp, byid = T) |> rowSums()
    outarea <- rgeos::gIntersects(regions, narwc.sp, byid = T) |> rowSums()

    narwc.sp <- narwc.sp[-unique(c(which(outarea==0), which(onland==1))),]
    narwc <- narwc[-unique(c(which(outarea==0), which(onland==1))),]

    # Create raster with the same resolution and origin as the existing density surfaces
    # Fill the raster with zeroes
    d.res <- raster::res(d$Jan)
    r.canada <- toRaster_area(canada_ocean)
    
    # ........................................................
    # Point density ---
    # ........................................................
  
    min.n <- 10
    
    df.count <- month.abb |> 
      purrr::set_names() |> 
      purrr::map(.f = ~{narwc |> dplyr::filter(region == "Canada", month == .x) |> nrow()}) |> 
      tibble::enframe() |> dplyr::rename(month = name) |> 
      dplyr::mutate(count = unlist(value)) |> dplyr::select(-value) |> 
      dplyr::mutate(month = factor(month, levels = month.abb)) |> 
      dplyr::mutate(n = ifelse(count >= min.n, count, 0))
    
    narw.months <- df.count[df.count$n >0,]$month
    
    narwc.bymonth <- narwc |> 
      dplyr::filter(region == "Canada") |> 
      (function(x) split(x = x, f = x$month))()
    
    r.dens <- purrr::map(
      .x = month.abb,
      .f = ~ {
        if(df.count[df.count$month == .x,]$n == 0){
          r.canada$d
        } else {
          uyr <- unique(narwc.bymonth[[.x]]$sighting_year)
          r.out <- lapply(X = uyr,
                          FUN = function(yr) {
                            rpts <- narwc.bymonth[[.x]] |> dplyr::filter(sighting_year == yr)
                            rpts <- split(rpts, f = factor(lubridate::day(rpts$date)))
                            
                            rpts.count <- lapply(X = rpts, FUN = function(y){
                              rr <- sp::SpatialPointsDataFrame(
                                coords = y[, c("x", "y")],
                                data = y, proj4string = narw_crs())
                              pointcount(r.canada$d, rr, canada_ocean)
                            })
                            
                            raster::stack(rpts.count)
                            
                          }) |> purrr::set_names(nm = uyr)
          raster::stack(r.out)
        }}) |> purrr::set_names(nm = month.abb)
    
    r.narw <- purrr::map(.x = r.dens, .f = ~raster::calc(x = .x, fun = mean))
    
    df.narw <- purrr::map(.x = r.narw, .f = ~{raster::as.data.frame(.x, xy = TRUE) |>
        dplyr::filter(!is.na(layer))})
    
    df.narw <- purrr::map(.x = month.abb, .f = ~dplyr::mutate(.data = df.narw[[.x]], month = .x)) |> 
      do.call(what = rbind) |> 
      dplyr::rename(count = layer)
    
    # ........................................................
    # Cell areas ---
    # ........................................................
    
    # Convert raster cells to polygons to calculate areas
    r.poly <- raster::rasterToPolygons(r.canada)
    df.narw$area <- rep(r.poly$area, length(month.abb))
    df.narw <- add_country(df.narw)
    
    # ........................................................
    # Dataset ---
    # ........................................................
    
    # Set random seed
    set.seed(20221017)
    
    # Create regular grid across the study area to avoid over inflation of zeroes
    # Use the buffer function from <terra> to ensure no points fall on land
    # including islands
    r.grid <- sp::makegrid(canada_ocean, cellsize = 75) |>
      sp::SpatialPoints(proj4string = narw_crs())
    
    r.grid <- rgeos::gIntersection(spgeom1 = terra::buffer(canada_ocean, width = -10), spgeom2 = r.grid)
    row.names(r.grid) <- as.character(seq_len(length(r.grid)))
    
    df.month <- month.abb |> 
      purrr::set_names() |> 
      purrr::map(.f = ~{
        
        if(df.count[df.count$month == .x,]$n == 0){
          
          NULL
          
        } else {
          
          # Extract the data for the month of interest
          df.filter <- df.narw |> dplyr::filter(month == .x)
          
          # Extract cells that have sightings
          pos <- df.filter |> dplyr::filter(count > 0)
          
          # Convert these locations into a polygon and buffer the resulting area
          pos.poly <- raster::rasterFromXYZ(pos[, c("x", "y", "count")], res = d.res, crs = narw_crs()) |>
            raster::rasterToPolygons()
          
          r.dissolve <- rgeos::gUnionCascaded(pos.poly)
          r.buffer <- rgeos::gBuffer(r.dissolve, width = 25)
          
          # Remove grid points that fall within the buffer
          r.grid.zero <- rgeos::gDifference(r.grid, r.buffer)
          
          df.zero <- data.frame(sp::coordinates(r.grid.zero), 
                                count = 0, 
                                month = factor(unique(pos$month)),
                                area = raster::extract(r.canada$area, r.grid.zero))
          
          df.zero <- add_country(df.zero)
          df.out <- rbind(pos, df.zero)
          row.names(df.out) <- NULL
          
          df.out
        }
      }) 
    
    # ........................................................
    # Soap film ---
    # ........................................................
    
    ocean.xy <- broom::tidy(canada_ocean) |> dplyr:: rename(x = long, y = lat)
    N <- floor(abs(raster::extent(canada_ocean)[2]-raster::extent(canada_ocean)[1])/85)
    gx <- seq(raster::extent(canada_ocean)[1], raster::extent(canada_ocean)[2], length.out = N)
    gy <- seq(raster::extent(canada_ocean)[3], raster::extent(canada_ocean)[4], length.out = N)
    gp <- expand.grid(gx, gy)
    names(gp) <- c("x","y")
    
    # Returns an error if last element (25th) is included
    # Also returns an error Error in place.knots(x, nk) : more knots than unique data values is not allowed
    b.ind <- c(1:7)
    
    # The GAM needs the border coordinates as a list of lists,
    # where each list describes one border segment or island:
    oceancoords <- ocean.xy |> dplyr::select(x,y,piece)
    borderlist <- split(oceancoords, oceancoords$piece)
    border.narw <- lapply(borderlist, `[`, c(1,2))
  
    borderlist <- borderlist[b.ind]
    border.narw <- border.narw[b.ind]
    border.narw <- lapply(seq_along(b.ind), function(n) as.list.data.frame(border.narw[[n]]))
    
    # We can now use the inSide function from mgcv to select knots that are inside the border
    knots <- gp[with(gp, mgcv::inSide(bnd = border.narw, x, y)), ]
    
    # This passes the soap_check test but still returns Error in crunch.knots
    # Manually remove knots that are on or outside boundary
    
    text.counter <- NULL
    while (is.null(text.counter)){
      
      test.dsm <- evaluate::evaluate(
        "mgcv::gam(formula = count ~ offset(log(area)) + s(x, y, bs = \"so\", xt = list(bnd = border.narw)),
                         data = df.month[[\"Jul\"]][1:20,],
                         method = \"REML\",
                         family = tw(link = \"log\"),
                         knots = knots)")
      
      if(stringr::str_sub(as.character(test.dsm[[2]]), 1, 21) == "Error in crunch.knots"){
        
        error.msg <- as.character(test.dsm[[2]])
        error.msg <- stringr::str_sub(error.msg, gregexpr(":", error.msg)[[1]][1] + 2, nchar(error.msg))
        problem.knot <- readr::parse_number(error.msg)
        knots <- knots[-problem.knot, ]
        # cat("Removing knot:", problem.knot, "\n")
        text.counter <- NULL
        
      } else {
        
        text.counter <- "Stop"
      }
      
    } # End while loop
    
    # Just islands
    border.islands <- border.narw
    border.islands[2:length(border.islands)] <- 
      purrr::map(.x = border.islands[2:length(border.islands)], .f = ~{
        tmp <- .x
        tmp$f = rep(0, times = length(tmp$x))
        tmp
      })
    
    
    # Set constraint on boundary conditions around all islands but not coastline
    border.narw <- purrr::map(.x = border.narw, .f = ~{
      tmp <- .x
      tmp$f = rep(0, times = length(tmp$x))
      tmp
    })
    
    # ........................................................
    # hGAM ---
    # ........................................................
    
    model.hgam.isl <- "count ~ offset(log(area)) + s(x, y, m = 2, bs = \"so\", xt = list(bnd = border.islands)) + s(x, y, by = month, bs = \"so\", m = 1, xt = list(bnd = border.islands)) + s(month, bs = \"re\", k = length(narw.months))"
    
    model.hgam.main <- "count ~ offset(log(area)) + s(x, y, m = 2, bs = \"so\", xt = list(bnd = border.narw)) + s(x, y, by = month, bs = \"so\", m = 1, xt = list(bnd = border.narw)) + s(month, bs = \"re\", k = length(narw.months))"
    
    df.hgam <- do.call(rbind, df.month) |> dplyr::mutate(month = factor(month))
    
    hgam.isl <- mgcv::bam(
      formula = as.formula(model.hgam.isl),
      data = df.hgam,
      method = "fREML",
      discrete = TRUE,
      family = tw(),
      knots = knots)
    
    hgam.main <- mgcv::bam(
      formula = as.formula(model.hgam.main),
      data = df.hgam,
      method = "fREML",
      discrete = TRUE,
      family = tw(),
      knots = knots)
    
    # Prediction df
    preds.df.hgam <- month.abb |> 
      purrr::set_names() |> 
      purrr::map(.f = ~{
        as.data.frame(r.canada, xy = TRUE, na.rm = TRUE) |> 
          dplyr::select(-d) |> 
          dplyr::mutate(month = factor(.x))
      }) |> do.call(what = rbind) |> 
      dplyr::mutate(Nhat_isl = 0, Nhat_main = 0, Nhat = 0)
    
    future::plan(multisession)
    
    # Model predictions
    preds.hgam.isl <-
      furrr::future_map(.x = narw.months, 
                        .f = ~predict_hgam(object = hgam.isl, 
                                           newdata = preds.df.hgam[preds.df.hgam$month %in% .x,], 
                                           type = "response")) |> do.call(what = c)
    
    preds.hgam.main <-
      furrr::future_map(.x = narw.months, 
                        .f = ~predict_hgam(object = hgam.main, 
                                           newdata = preds.df.hgam[preds.df.hgam$month %in% .x,], 
                                           type = "response")) |> do.call(what = c)
    
    preds.df.hgam$Nhat_isl[preds.df.hgam$month %in% narw.months] <- preds.hgam.isl
    preds.df.hgam$Nhat_main[preds.df.hgam$month %in% narw.months] <- preds.hgam.main
    
    preds.df.hgam <- preds.df.hgam |> 
      dplyr::rowwise() |> 
      dplyr::mutate(Nhat = min(Nhat_isl, Nhat_main)) |> 
      dplyr::ungroup()
    
    # Rasters - Canada only
    d.narwc.canada <- month.abb |>
      purrr::set_names() |> 
      purrr::map(.f = ~{
        preds.list.hgam <- split(preds.df.hgam, f = factor(preds.df.hgam$month))
        raster::rasterFromXYZ(preds.list.hgam[[.x]][, c("x", "y", "Nhat")], res = d.res, crs = narw_crs())})
  
    d.out <- gsub(pattern = "density", replacement = "density_canada", x = list.files("data/densitymaps"))
    
    # Save outputs
    purrr::walk(.x = 1:12, .f = ~{
      raster::writeRaster(d.narwc.canada[[.x]], 
                          filename = file.path("data/densitymaps/", d.out[.x]), 
                          overwrite = TRUE, 
                          format = "HFA")
    })
    
    unlink(list.files("data/densitymaps", pattern = ".aux", full.names = TRUE))
    
  }
  
}

merge_density <- function(){

  d_usa <- list.files("data/densitymaps", pattern = "density_month", full.names = TRUE)
  d_canada <- list.files("data/densitymaps", pattern = "canada", full.names = TRUE)
  
  # Whole range
  d.narwc <- purrr::map2(.x = d_usa,
                         .y = d_canada,
                         .f = ~{

      r_usa <- raster::raster(.x) |> raster::projectRaster(crs = narw_crs())
      r_usa <- 25*r_usa/100
      r_canada <- raster::raster(.y)
      r_narw <- raster::merge(r_usa, r_canada)
      names(r_narw) <- "Nhat"
      r_narw}) |> purrr::set_names(nm = month.abb)

  d.narwc  
  
}


#' Download density maps
#'
#' Download density maps from the Duke University website and unzip files
#'
#' @param version Version of the North Atlantic right whale model
#'
#' @export

get_ais <- function(year = 2021) {
  
  # Avoid time out issues
  options(timeout = 4000)
  
  # AIS data published at:
  # https://coast.noaa.gov/digitalcoast/data/
  
  # Direct link to dataset
  src <- paste0("https://marinecadastre.gov/downloads/data/ais/ais", year, "/AISVesselTracks", year, ".zip")
  
  # Create data directory
  f <- file.path("data", "ais")
  dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
  
  # Download and unzip data
  dst <- file.path(f, basename(src))
  utils::download.file(url = src, destfile = dst)
  unzip(zipfile = file.path(f, basename(src)), exdir = f)
  
  # Remove zip file
  unlink(dst)
  
  # Delete unnecessary files and re-arrange folder
  unlink("data/densitymaps/Animations", recursive = TRUE)
  unlink("data/densitymaps/ArcGIS_Symbology", recursive = TRUE)
  unlink("data/densitymaps/Maps", recursive = TRUE)
  unlink(list.files(paste0("data/densitymaps/Rasters/", paste0(start, "-", end)), full.names = TRUE, pattern = "aux.xml"))
  invisible(file.copy(list.files(paste0("data/densitymaps/Rasters/", paste0(start, "-", end)), full.names = TRUE, pattern = "density_month"), "data/densitymaps/", recursive = TRUE))
  unlink(list.files("data/densitymaps", pattern = ".pdf", full.names = TRUE))
  unlink(list.files("data/densitymaps", pattern = ".txt", full.names = TRUE))
  unlink("data/densitymaps/Rasters/", recursive = TRUE)
  
}


#' Load wind farm shapefiles
#'
#' @export

get_farms <- function(){

  windfarms <- lapply(X = list.files("data/windfarms", pattern = "windfarm_.+shp", full.names = TRUE), FUN = function(x) raster::shapefile(x))
  windfarms <- do.call(raster::bind, windfarms)
  windfarms <- sp::spTransform(windfarms, narw_crs())
  return(windfarms)
}

TL <- function(r, a = 0.185){
  loss <- 20*log10(r*1000)
  loss <- 15*log10(r*1000)
  loss[loss < 0] <- 0
  loss <- loss + a * r
  return(loss)
}

summary_geo <- function(obj, lims, res, geomap) {
  lapply(seq(dim(obj)[3]), function(x){
    obj.ind <- obj[,,x]
    obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
    dd <- sapply(1:nrow(obj.ind), FUN = function(r) geoD(geomap, obj.ind[r,1], obj.ind[r,2], obj.ind[r,3], obj.ind[r,4], lims, res))
    unname(dd)
  })
}

piling_noise <- function(source.lvl = 180, x, y){
  
  # ambient <- targets::tar_read(density_narw)[[1]] |> raster::raster()
  # ambient[!is.na(ambient)] <- 60
  
  d.support <- targets::tar_read(density_support) |> raster::raster()
  d.support[d.support < 1] <- NA
  d.support[is.na(d.support)] <- 3
  
  d.support[raster::cellFromXY(d.support, cbind(x,y))] <- 2

  rdist <- terra::gridDistance(d.support, origin = 2, omit = 3)
  
  dB <- source.lvl - TL(rdist)
  dB[dB < 60] <- 60
  plot(terra::rast(dB), col = pals::viridis(100), main = "Noise levels (dB)", xlab = "Easting (km)")
}

geodesic_dist <- function(){

  
  r.grid <- r.grid[which(sp::over(r.grid, d.poly) == 1),]
  
  r.grid$cell <- purrr::map_dbl(.x = 1:length(r.grid), 
                                .f = ~raster::cellFromXY(d.support, r.grid[.x,]@coords))
  
  d <- targets::tar_read(density_narw)[[1]] |> 
    raster::raster() |> raster::mask(mask = d.poly)
  vor <- dismo::voronoi(r.grid)
  vr <- raster::intersect(vor, d.poly)
  vr <- raster::rasterize(vr, d, 'cell')
  
  # Calculate geodesic distances from each grid point -- return a raster
  future::plan(future::multisession)
  
  suppressMessages(gd <- furrr::future_map(.x = r.grid$cell,
                                           .f = ~{
                                             library(raster)
                                             r <- d.support
                                             r[.x] <- 2
                                             terra::gridDistance(r, origin = 2, omit = 3)}, 
                                           .options = furrr::furrr_options(seed = TRUE)))
  
  # gd <- raster::stack(gd)
  names(gd) <- r.grid$cell
  
  return(list(grid = r.grid, cellno = vr, raster = gd))
  
}








get_turbines <- function(){
  
  turbines <- readr::read_csv("data/windfarms/turbine_locs.csv", col_types = "fdd") |> 
    janitor::clean_names()
  
  if(!all(names(turbines) %in% c("farm", "longitude", "latitude"))) 
    stop("Cannot find all required fields")
  
  turbines <- add_xy(dat = turbines)
  
  split(x = turbines, f = turbines$farm)
    
}

get_regions <- function(eez = FALSE) {
  
  if(eez){
    
    regions <- raster::shapefile(dir(path = "data/regions",
                                     pattern = "regions_eez.shp$", 
                                     full.names = TRUE)) |> 
      sp::spTransform(CRSobj = narw_crs())
    
  } else {
  
  # Load shapefile of spatial regions (CCB, MIDA, SEUS, etc.)
  regions <- raster::shapefile(dir(path = "data/regions", 
                                   pattern = "regions.shp$", 
                                   full.names = TRUE))
  regions$Id <- NULL

  sp::spTransform(x = regions, CRSobj = narw_crs())
  }
}

get_world <- function(sc = "medium"){
  
  wrld <- rnaturalearth::ne_countries(scale = sc, returnclass = "sp")
  wrld <- sp::spTransform(wrld, CRSobj = narw_crs()) |> 
    raster::crop(y = raster::extent(c(-574.1478, 1934.829, -1534.189, 2309.078)))
  wrld <- wrld[wrld$admin %in% c("Canada", "United States of America"),]
  return(wrld)
}

# BIOENERGETICS ------------------------------------------------------

#' Decrease in milk feeding efficiency as a function of calf age
#'
#' @param a Calf age in days.
#' @param Tstop Time at which milk stops nursing and becomes nutritionally independent (in days).
#' @param Tdecrease Time at which milk comsuption starts to decrease (in days).
#' @param E Steepness of the decline

milk_assimilation <- function(a = seq(0, 500), Tstop = 365, Tdecrease = 100, E = 0.9){
  p <- (a-Tdecrease)/(Tstop-Tdecrease)
  res <- (1 - p)/ (1 - (E*p))
  res[res<0] <- 0
  res[res>1 & a>Tdecrease] <- 0
  res[res>1 & a<Tdecrease] <- 1
  return(res)
}

#' Increase in milk provisioning by female as a function of her body condition.
#'
#' @param E Steepness of the non-linear relationship.
#' @param blubber_mass Mass of blubber (kg).
#' @param maintenance_mass Mass of structural tissues (kg).
#' @param target_condition Target body condition, expressed as the ratio of reserve to maintenance mass.
#' @param starvation Starvation threshold

milk_provisioning <- function(E = -2, blubber_mass, maintenance_mass, target_condition = 0.3, starvation = 0.15){
  p1 <- 1-E
  p2 <- blubber_mass - starvation * maintenance_mass
  p3 <- maintenance_mass * (target_condition - starvation)
  p4 <- E * (blubber_mass - (starvation*maintenance_mass))
  res <- p1*p2/(p3-p4)
  res[blubber_mass/maintenance_mass < starvation] <- 0
  res[blubber_mass/maintenance_mass >= target_condition] <- 1
  return(res)
}


# Feeding effort as a function of body condition
#' Title
#'
#' @param eta Steepness of the non-linear relationship 
#' @param target_condition 
#' @param maintenance_mass 
#' @param blubber_mass 

feeding_effort <- function(eta, target_condition = 0.3, maintenance_mass, blubber_mass){
  1/(1+exp(-eta * ((target_condition*maintenance_mass/blubber_mass)-1)))
}

#' Feeding effort as function of copepod density
#'
#' @param gamma Feeding threshold
#' @param D Coepepod density

cop_threshold <- function(gamma, D){
  1/(1+exp(gamma-D))
}

# PLOTTING ------------------------------------------------------

# Jason's colour scale
colour_breaks <- function(dat){
  colour.breaks <- 25*c(0,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3,10)/100
  colour.breaks <- c(colour.breaks, seq(max(colour.breaks), ceiling(max(dat$Nhat, na.rm = TRUE)), 
                                        length.out = 5))
  colour.breaks <- round(colour.breaks[!duplicated(colour.breaks)],3)
  return(colour.breaks)
}

# UTILITY ------------------------------------------------------

format_table <- function(df, top = TRUE, bottom = TRUE, sign = "-"){
  # dashes <- purrr::map_dbl(.x = names(df), .f = ~nchar(.x))
  # dashes <- purrr::map(.x = dashes, .f = ~paste0(rep(sign, .x), collapse = ""))
  rbind(df[1:nrow(df) - 1,], c("---", rep("",ncol(df)-1)), df[nrow(df),])
}
