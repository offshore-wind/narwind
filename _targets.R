# library(targets)
# library(tarchetypes) #  As an extension to 'targets', the 'tarchetypes' package provides convenient user-side functions to make 'targets' easier to use. By establishing reusable archetypes for common kinds of targets and pipelines, these functions help express complicated reproducible workflows concisely and compactly. 
# library(future)
# library(future.batchtools) # A Future API for Parallel and Distributed Processing using ‘batchtools’

future::plan(future::multisession)

# Set packages to load
targets::tar_option_set(packages = c('rgdal', 'stringr', 'spatstat.geom', 'ggplot2', 'raster',), deployment = 'main')

# Load R files and workflows
lapply(list.files(file.path("R"), full.names = TRUE, recursive = TRUE, pattern = '\\.R$'), source)

# Assemble workflow
list(data_targets)
