library(targets)

# make workflow accessible from interactive R sessions whose working directory 
# is the overall project/repository's root directory
workflow = file.path('workflows', 'simulation')
writeLines(
  text = c(
    paste('script:', file.path(workflow, '_targets.R')),
    paste('store:', file.path(workflow, '_targets'))
  ), 
  con = '_targets.yaml'
)

tar_make(density_map_download)