data_targets = list(
  
  # Download density maps and unzip files
  tar_target(
    name = density_map_download, 
    command = {
      
      # Density maps (v11.1) published at:
      # https://seamap.env.duke.edu/models/Duke/EC/EC_North_Atlantic_right_whale_history.html
      
      # direct link to dataset
      src = 'https://seamap.env.duke.edu/seamap-models-files/Duke/EC/North_Atlantic_right_whale/v11.1/EC_North_Atlantic_right_whale_v11.1.zip'
      
      # create data directory
      f = file.path('workflows', 'simulation', 'data', 'densitymaps')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      
      # download and unzip data
      dst = file.path(f, basename(src))
      download.file(url = src, destfile = dst)
      unzip(zipfile = file.path(f, basename(src)), exdir = f)
      
      # remove zip file
      unlink(dst)
      
    }
  ),
  
  # List paths to density map files
  tar_target(
    name = density_maps,
    command = {
      f = file.path('workflows', 'simulation', 'data', 'densitymaps', 'Rasters',
                    '2010-2018')
      maps = dir(path = f, pattern = 'density\\.img$', full.names = TRUE)
      maps
    }
  ),
  
  tar_target(
    name = windfarms,
    command = {
      windfarm1 = readOGR(
        dsn = file.path('data', 'WindFarms_USA', 'WindFarm1.shp'), 
        p4s = '+proj=lonlat'
      )
      
      windfarm2 = readOGR(
        dsn = file.path('data', 'WindFarms_USA', 'WindFarm2.shp'), 
        p4s = '+proj=lonlat'
      )
      
      windfarm3 = readOGR(
        dsn = file.path('data', 'WindFarms_USA', 'WindFarm3.shp'), 
        p4s = '+proj=lonlat'
      )
      windfarms <- union(windfarm1, windfarm2)
      windfarms <- union(windfarms, windfarm3)
      windfarms
    }
  ),
  
  tar_target(
    name = maximum_density_support,
    command = {
      
      # load and label density maps
      f = targets::tar_read(density_maps)
      map_rasters = lapply(f, function(f) readGDAL(fname = f, silent = TRUE))
      map_month = as.numeric(
        str_match(string = f, pattern = 'month_([0-9]{2})')[,2]
      )
      
      map_rasters = lapply(f, function(f) readGDAL(fname = f, silent = TRUE))
      
      # # spatial support for all maps
      # map_supports = lapply(map_rasters, function(m) {
      #   is.finite(m$band1)
      # })
      # 
      # # build spatial mask for raster with largest spatial support
      # max_ind = which.max(sapply(map_supports, sum))
      # max_support = map_rasters[[max_ind]]
      # max_support$band1 = is.finite(max_support$band1)
      
      # identify locations that exceed a minimum density value in some raster
      filtered_support = map_rasters[[1]]
      filtered_support$band1 = FALSE
      for(m in map_rasters) {
        m$band1[is.na(m$band1)] = 0
        filtered_support$band1 = filtered_support$band1 | (m$band1 > 0.001)
      }
      
      # identify connected components within the filtered support
      coords = coordinates(filtered_support)
      rowseq = sort(unique(coords[,'x']))
      colseq = sort(unique(coords[,'y']))
      fs = im(mat = as.matrix(filtered_support), xcol = colseq, yrow = rowseq)
      components = connected(fs)
      
      # use the largest component to define the spatial mask
      filtered_support$band1 = as.numeric(components$v)
      largest_component = which.max(table(filtered_support$band1))
      filtered_support$band1 = filtered_support$band1 == largest_component
      filtered_support$band1[!is.finite(filtered_support$band1)] = FALSE
      
      
      #
      # manually add spatial support for CA
      #
      
      # load region boundaries
      regions = readOGR(
        dsn = file.path('workflows', 'simulation', 'data', 'regions', 'NARW_regions.shp'), 
        p4s = '+proj=lonlat'
      )
      
      # project boundaries
      regions.proj = spTransform(
        x = regions, 
        CRSobj = filtered_support@proj4string
      )
      
      # save a copy for plotting
      regions.proj.full = regions.proj
      
      # filter regions
      regions.proj = regions.proj[
        regions.proj@data$region %in% c('GSL', 'ELSW'),
      ]
      
      # overlay regions onto current spatial support
      regions.raster = rasterize(
        x = regions.proj, 
        y = raster(x = regions.proj, 
                   resolution = filtered_support@grid@cellsize, 
                   origin = origin(raster(filtered_support)))
      )
      combined_support = merge(regions.raster, raster(filtered_support))
      
      # get raw country outlines
      outlines = map_data(map = 'world', region = c('usa', 'canada'))
      
      # filter country outlines
      outlines = outlines[
        !(outlines$subregion %in% c('Alaska', 'Hawaii')) & 
          !(outlines$lat > 61),
      ]
      
      # project country outlines
      outlines.sp = SpatialPoints(
        coords = outlines[,c('long','lat')], 
        proj4string = CRS('+proj=longlat')
      )
      outlines.proj = spTransform(
        x = outlines.sp, 
        CRSobj = filtered_support@proj4string
      )
      outlines[,c('long','lat')] = outlines.proj@coords
      
      # convert to polygons
      outlines.poly = df_to_SpatialPolygons(
        df = outlines, 
        keys = 'group', 
        coords = c('long', 'lat'),
        proj = outlines.proj@proj4string
      )
      
      # overlay regions onto spatial support
      outlines.raster = rasterize(
        x = outlines.poly, 
        y = raster(x = outlines.poly, 
                   resolution = filtered_support@grid@cellsize, 
                   origin = origin(raster(combined_support)))
      )
      
      # mark land regions with a mask value
      outlines.raster@data@values[outlines.raster@data@values > 0] = -1
      
      # combine data to form mask
      cleaned_support = merge(outlines.raster, combined_support)
       
      # remove/mask land regions and non-labeled regions from spatial support
      cleaned_support@data@values[cleaned_support@data@values <= 0] = NA
      
      # convert to logical
      cleaned_support@data@values = cleaned_support@data@values > 0
      
      # trim spatial support
      cleaned_support = trim(cleaned_support)
      
      # convert to SGDF
      cleaned_support <- as(cleaned_support, 'SpatialGridDataFrame')
      
      
      #
      # check fit!
      #
      
      # project sighting coordinates
      sightings.pp = SpatialPoints(
        coords = sightings[,c('Longitude','Latitude')], 
        proj4string = CRS('+proj=longlat')
      )
      sightings.proj = spTransform(
        x = sightings.pp, 
        CRSobj = cleaned_support@proj4string
      )
      
      # calculate proportion of sightings that lie in support
      d = over(sightings.proj, cleaned_support)
      mean(is.finite(unlist(d)))
      
      
      #
      # plot simulation region alongside coast line
      #
      
      # make simulation support plottable
      df = data.frame(cleaned_support)

      # build spatial support plot
      pl = ggplot() + 
        # country outlines
        geom_polygon(
          mapping = aes(x = long, y = lat, group = group), 
          data = outlines
        ) + 
        # simulation support
        geom_raster(
          mapping = aes(x = s1, y = s2, fill = layer),
          data = df
        ) + 
        # overlay sightings
        geom_point(
          mapping = aes(x = Longitude, y = Latitude),
          data = data.frame(sightings.proj),
          size = .001
        ) + 
        # formatting
        scale_x_continuous(
          name = 'Easting (m)',
          limits = c(-.9e6, 1.75e6), 
          oob = scales::oob_keep
        ) + 
        scale_y_continuous(
          name = 'Northing (m)',
          limits = c(-1e6, 2.3e6), 
          oob = scales::oob_keep
        ) + 
        scale_fill_manual(values = c('TRUE' = 'cornflowerblue')) + 
        guides(fill = 'none') +
        theme_few() +
        coord_equal()
      
      # save plot
      f = file.path('workflows', 'simulation', 'output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, 'simulation_support.png'), 
             dpi = 'print')
      
      #
      # plot regions alongside coastline
      #
      
      df = fortify(regions.proj.full)
      df$id = factor(df$id, levels = 0:5, labels = regions.proj.full$region)
      
      # build spatial support plot
      pl = ggplot() + 
        # country outlines
        geom_polygon(
          mapping = aes(x = long, y = lat, group = group), 
          data = outlines
        ) + 
        # region outlines
        geom_polygon(
          mapping = aes(x = long, y = lat, group = group, col = id, fill = id), 
          data = df,
          alpha = .2
        ) + 
        # formatting
        scale_x_continuous(
          name = 'Easting (m)',
          limits = c(-.9e6, 1.75e6), 
          oob = scales::oob_keep
        ) + 
        scale_y_continuous(
          name = 'Northing (m)',
          limits = c(-1e6, 2.3e6), 
          oob = scales::oob_keep
        ) + 
        scale_fill_brewer('Region', type = 'qual', palette = 'Dark2') + 
        scale_color_brewer('Region', type = 'qual', palette = 'Dark2') + 
        theme_few() +
        coord_equal()
      
      # save plot
      f = file.path('workflows', 'simulation', 'output', 'figures')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, 'simulation_regions.png'), 
             dpi = 'print')
      
      # export the mask!
      # max_support
      filtered_support
    }
  ),
  
  tar_target(
    name = spatial_regions,
    command = {
      # load shapefile of spatial regions CCB, MIDA, SEUS, etc.
      d = file.path('workflows', 'simulation', 'data', 'regions')
      f = dir(path = d, pattern = '.shp$', full.names = TRUE)
      regions = readOGR(dsn = f, p4s = '+proj=longlat +datum=WGS84')
      # remove null field
      regions$Id = NULL
      # project coordinates to raster grid
      m = readGDAL(density_maps[1])
      regions_proj = spTransform(x = regions, CRSobj = m@proj4string)
      # export projected regions
      regions_proj
    }
  ),
  
  tar_target(
    name = sightings,
    command = {
      read.csv(file.path('workflows', 'simulation', 'data', 'sightings',
                         'Sightings_IDs.csv'))
    }
  )
  
)