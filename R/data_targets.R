data_targets = list(
  
  # ........................................................................
  # Import sightings from NARWC
  # ........................................................................
  
  targets::tar_target(name = sightings, command = read.csv("data/sightings/narwc_sightings.csv")),
  
  # ........................................................................
  # Import and process model parameter spreadsheet
  # ........................................................................
  
  targets::tar_target(name = params, command = {
    readr::read_csv("data/parameters/BOEM_140M0121C0008_ModelParameters.csv", na = c("-", "NA"), skip = 1, 
                    col_types = readr::cols()) |> 
      janitor::clean_names() |> 
      dplyr::select(-row) |> 
      dplyr::mutate(min = as.numeric(min), 
                    max = as.numeric(max),
                    mean_median = as.numeric(mean_median),
                    sd_se = as.numeric(sd_se),
                    sample_size = as.numeric(sample_size))}),
  
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
  targets::tar_target(name = support_poly, command = support_as_polygon()),
  targets::tar_target(name = density_narw, command = clip_density()),
  targets::tar_target(name = density_weighted_seus, command = w_density(target = "SEUS", option = 2)),
  targets::tar_target(name = density_weighted_gsl, command = w_density(target = "GSL", option = 2)),
  
  # ........................................................................
  # Surfaces
  # ........................................................................
  
  targets::tar_target(name = dummy_prey, command = proxy_prey()),
  # targets::tar_target(name = dummy_noise, command = proxy_noise()),
  targets::tar_target(name = dummy_vessels, command = proxy_vessels()),
  targets::tar_target(name = fishing_layer, command = entgl_surface()),
  targets::tar_target(name = noise_layer, command = noise_surface()),
  targets::tar_target(name = daylight, command = get_daylight()),
  targets::tar_target(name = entgl_d, command = get_entglD()),
  
  # ........................................................................
  # Dose-response
  # ........................................................................
  
  targets::tar_target(name = doseresponse, command = get_doseresponse()),
  
  # ........................................................................
  # Minimum body condition required to bring fetus to term 
  # ........................................................................
  
  targets::tar_target(name = gam_gest, command = gestation_threshold()),
  
  # ........................................................................
  # Save data to package
  # ........................................................................
  
  targets::tar_target(name = save_data, command = save_objects()),

  # ........................................................................
  # Offshore wind scenarios
  # ........................................................................
  
  targets::tar_target(name = scenarios, command = get_scenarios())
  
)