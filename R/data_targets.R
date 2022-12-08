data_targets = list(
  
  # ........................................................................
  # Import sightings from NARWC
  # ........................................................................
  
  targets::tar_target(name = sightings, command = read.csv("data/sightings/narwc_sightings.csv")),
  
  # ........................................................................
  # Import shapefiles and other spatial data
  # ........................................................................
  
  targets::tar_target(name = windfarms, command = get_farms()),  
  targets::tar_target(name = regions, command = get_regions()),
  targets::tar_target(name = regions_m, command = regions_matrix()),
  targets::tar_target(name = regions_eez, command = get_regions(eez = TRUE)),
  targets::tar_target(name = canada, command = get_regions(canada = TRUE)),
  targets::tar_target(name = world, command = get_world()),
  targets::tar_target(name = turbines, command = get_turbines()),
  
  # ........................................................................
  # Density maps 
  # ........................................................................
   
  targets::tar_target(name = density_usa, command = get_density(version = 12, region = "US")),
  targets::tar_target(name = density_canada, command = get_density(region = "CA")),
  targets::tar_target(name = density_merge, command = merge_density()),
  
  # ........................................................................
  # Compute spatial support and clip density rasters
  # ........................................................................
  
  targets::tar_target(name = density_support, command = spatial_support()),
  targets::tar_target(name = density_narw, command = spatial_support()),
  targets::tar_target(name = support_poly, command = support_as_polygon())
  
  # ........................................................................
  # Calculate geodesic distances
  # ........................................................................
  
  # targets::tar_target(name = geodesic, command = geodesic_dist(grid.res = 100))
  
  # ........................................................................
  # Download AIS data
  # ........................................................................
  # targets::tar_target(AIS_download, get_ais()),
  


)