# MODEL ------------------------------------------------------

write <- function(x, ...) UseMethod("write")

find_lean <- function(target.max, maxBC = 0.6, age = 69){
  opt.fun <- function(x, age, target.max, maxBC){
    return((length2mass(age2length(age, agL(age)), mL(), x)/(1-maxBC)-target.max)^2)}
  out <- optimize(f = opt.fun, interval = c(0,1), target.max = target.max, maxBC = maxBC, age = age)
  return(out$minimum)
}

meta <- function(seed = 215513, 
                 fill.NA = FALSE,
                 min.yr = NULL,
                 max.yr = NULL,
                 comb = NULL,
                 input = list(),
                 theoretic.bounds = NULL,
                 remove.study = NULL,
                 verbose = FALSE){
  
  set.seed(seed)
  options(pillar.sigfig = 7)
  
  if(length(input) > 0){
    param.input <- input$param.input
    bycols <- input$bycols
    fill.NA <- input$fill.NA
    fillers <- input$fillers
    min.yr <- input$min.yr
    max.yr <- input$max.yr
    comb <- input$comb
    remove.study <- input$remove.study
    if(is.null(bycols)) bycols <- ""
  }

  # See https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
  # for a detailed explanation
  
  # dat <- readr::read_csv("data/parameters/BOEM_140M0121C0008_ModelParameters.csv", na = c("-", "NA"), col_types = readr::cols()) |> 
  #   janitor::clean_names()
  dat <- params
  dat <- dplyr::select(dat, 1:which(names(dat) == "first_author")) |> 
    dplyr::select(-note_comments, -description_definition, -symbol)
  
  allp <- sort(unique(dat$parameter))
  allp <- allp[!allp == ""]
  
  if(length(input) == 0){
    print(data.frame(Parameter = allp), right = FALSE)
    cat("\n")
    param.input <- readline(prompt = "Select parameter: ")
  }
  
  isnumeric <- suppressWarnings(!is.na(as.numeric(param.input)))
  if(isnumeric){
      param <- allp[as.numeric(param.input)]
    } else {
      param <- param.input
    }

  dat.f <- dat |> dplyr::filter(parameter == param)
  if(!is.null(min.yr)) dat.f <- dat.f |> dplyr::filter(year > min.yr)
  if(!is.null(max.yr)) dat.f <- dat.f |> dplyr::filter(year < max.yr)
  if(!is.null(remove.study)) dat.f <- dat.f |> dplyr::filter(!first_author %in% remove.study)
  
  # Separate rows where needed
  dat.f <- tidyr::separate_rows(dat.f, "age_class", sep = ",") |> 
    dplyr::mutate(age_class = trimws(age_class))
  
  dat.f <- tidyr::separate_rows(dat.f, "sex", sep = ",") |> 
    dplyr::mutate(sex = trimws(sex))
  
  dat.f <- tidyr::separate_rows(dat.f, "reproductive_state", sep = ",") |> 
    dplyr::mutate(reproductive_state = trimws(reproductive_state))
  
  cat("\n")
  if(verbose) print(dat.f, n = 100)
  cat("\n")
  
  if(length(input) == 0){
    print(data.frame(Column = names(dat.f)), right = FALSE)
    bycols <- readline(prompt = "Select columns: ")
  } else {
    isnumeric <- suppressWarnings(!all(is.na(as.numeric(bycols))))
    if(isnumeric){
      bycols <- names(dat.f)[as.numeric(bycols)]
    } else {
      bycols <- bycols
    }
  }
  
  if(all(bycols == "")){
    dat.f$dummy <- "All"
    bycols <- "dummy"
  } else {
    if(length(input) == 0) bycols <- names(dat.f)[as.numeric(stringr::str_trim(unlist(strsplit(bycols, ","))))]
      # bycols <- names(dat.f)[bycols] else 
  }
  
  # Fill NA values if necessary
  if(fill.NA){
    if(length(input) > 0){
      for(k in 1:length(bycols)){
        dat.f[[bycols[k]]] <- ifelse(is.na(dat.f[[bycols[k]]]), fillers[k], dat.f[[bycols[k]]])
      }
    } else {
      fillers <- vector(length = length(bycols))
      for(k in bycols){
        word <- readline(prompt = paste0("Fill NA for <", k, "> with: "))
        dat.f[[k]] <- ifelse(is.na(dat.f[[k]]), word, dat.f[[k]])
      }
    }
  } else {
      for(k in bycols){
        dat.f[[k]] <- ifelse(is.na(dat.f[[k]]), "<NA>", dat.f[[k]])
      }
  }
  
  # If sample size is not given, assume 1
  dat.f <- dat.f |> 
    dplyr::rowwise() |> 
    dplyr::mutate(sample_size = ifelse(is.na(sample_size), 1, sample_size)) |> 
    dplyr::ungroup()
  
  # Combine categories if needed
  if(!is.null(comb)){
    for(k in 1:length(bycols)){
      dat.f[[bycols[k]]] <- ifelse(dat.f[[bycols[k]]] %in% comb[[1]], paste0(comb[[1]], collapse = ", "), dat.f[[bycols[k]]])
    }
  }
  
  dat.list <- split(dat.f, purrr::map(.x = bycols, .f = ~dplyr::pull(dat.f, .x)), drop = FALSE)
  
  dat.out <- purrr::map(.x = dat.list, .f = ~{
    
    if(nrow(.x) > 0){
    
    for(j in 1:nrow(.x)){

      myvars <- .x[j, c("min", "max", "mean_median", "sd_se")]
      myvars <- which(!is.na(myvars))
      
      varmin <- ifelse(all(is.na(.x$min)), NA, min(.x$min, na.rm = TRUE))
      varmax <- ifelse(all(is.na(.x$max)), NA, max(.x$max, na.rm = TRUE))
      
      if(is.na(varmin) & !is.null(theoretic.bounds)) varmin <- theoretic.bounds[1] else varmin <- min(dat.f$min, na.rm = TRUE)
      if(is.na(varmax) & !is.null(theoretic.bounds)) varmax <- theoretic.bounds[2] else varmax <- max(dat.f$max, na.rm = TRUE)
      
      if(length(myvars) == 1){
        
        if(myvars == 3){ # Only mean
          
          # If sample size > 1, estimate SD based on methods of Wan et al. (2014) - https://doi.org/10.1186/1471-2288-14-135
          # else: take weighted mean of other SD, weighted by sample size
          
          if(.x[j,]$sample_size > 1){

          .x[j,]$sd_se <- (varmax-varmin)/(2*(qnorm((max(.x[j,]$sample_size, 1, na.rm = TRUE)-0.375)/(max(.x[j,]$sample_size, 1, na.rm = TRUE)+0.25),0,1)))
            
          } else {

          w <- .x[, c("sd_se", "sample_size")]
          wm <- weighted.mean(w$sd_se, w = w$sample_size, na.rm = TRUE)
          if(nrow(.x) > 1) .x[j,]$sd_se <- wm
          }

        } else if(myvars == 1){ # Only min
          
          # Only min -- take mean/SD of random draws from Uniform bounded by min value and max across records
          
          n <- runif(10000, .x[j,]$min, varmax)
          .x[j,]$mean_median <- mean(n)
          .x[j,]$sd_se <- sd(n)
          
        } else if(myvars == 2){ # Only max
          
          # Only max -- take mean/SD of random draws from Uniform bounded by min across records and max value
          
          n <- runif(10000, varmin, .x[j,]$max)
          .x[j,]$mean_median <- mean(n)
          .x[j,]$sd_se <- sd(n)
          
        }
        
      } else {
        
        if(sum(myvars) == 1){
          
          # Only min
          
          n <- runif(10000, .x[j,]$min, varmax)
          .x[j,]$mean_median <- mean(n)
          .x[j,]$sd_se <- sd(n)
          
        } else if(sum(myvars) == 2){
          
          # Only max
          
          n <- runif(10000, varmin, .x[j,]$max)
          .x[j,]$mean_median <- mean(n)
          .x[j,]$sd_se <- sd(n)
          
        } else if(sum(myvars) == 3){
          
          # Only min and max
          
          n <- runif(10000, varmin, varmax)
          .x[j,]$mean_median <- mean(n)
          .x[j,]$sd_se <- sd(n)
          
        } else if(sum(myvars) == 4 | sum(myvars) == 6){ 
          
          # Only min and mean
          
          if(.x[j,]$sample_size > 1){
          
          .x[j,]$sd_se <- (varmax-.x[j,]$min)/(2*(qnorm((max(.x[j,]$sample_size, 1, na.rm = TRUE)-0.375)/(max(.x[j,]$sample_size, 1, na.rm = TRUE)+0.25),0,1)))
          
          } else {
            
            w <- .x[, c("sd_se", "sample_size")]
            w <- w[complete.cases(w),]
            wm <- weighted.mean(w$sd_se, w = w$sample_size, na.rm = TRUE)
            if(nrow(.x) > 1) .x[j,]$sd_se <- wm
            
          }
        }
      }
    }

    if(nrow(.x) > 1) {
      
      # Apply the inverse-variance method to calculate a weighted mean, using adjusted random-effects weights
# 
#       if(all(is.na(.x$sd_se)) | sum(!is.na(.x$sd_se)) == 1){
#         
#         combined.mean <- mean(.x$mean_median)
#         combined.se <- NA
#         
#       } else {

      if(all(is.na(.x$sd_se))){
        
        combined.mean <- mean(.x$mean_median)
        combined.se <- NA
        
      } else {
      
      weights <- 1 / ((.x$sd_se^2) + sd(.x$mean_median)^2)
      combined.mean <- sum(weights * .x$mean_median) / sum(weights)
      combined.se <- sqrt (sum(weights ^ 2 * .x$sd_se ^ 2) / sum(weights ^ 2))
      # }

      }
      tibble::tibble(mean = combined.mean, se = combined.se)

    } else if (nrow(.x) <= 1){
      
      tibble::tibble(mean = .x$mean_median, se = .x$sd_se)
     
    } 
    }
  }) |> tibble::enframe() |> 
    tidyr::unnest(cols = c(value)) |> 
    dplyr::arrange(name)

  # return(list(dat.f, dat.out))
  return(dat.out)
}

#' #' Initialize model
#' #'
#' #' @export
#' load <- function(){
#'   Rcpp::sourceCpp("src/simtools.cpp")
#'   # Rcpp::sourceCpp("src/bioenergetics_functions.cpp")
#'   source("R/run_model.R")
#'   assign("init.model", value = TRUE, envir = .GlobalEnv)
#' }

w_density <- function(target = "SEUS", option = 1){
  
  # Create raster from regions polygon
  reg <- targets::tar_read(regions)
  supportp <- targets::tar_read(support_poly)
  dens <- targets::tar_read(density_narw)
  densr <- purrr::map(.x = dens, .f = ~raster::raster(.x) |> raster::crop(dens[[1]]))
  reg$rank <- rank(reg$region)
  r <- raster::raster()
  raster::extent(r) <- raster::extent(reg)
  rp <- raster::rasterize(reg, r, 'rank')
  rp <- raster::resample(rp, densr[[1]], method = "ngb")
  rp <- raster::mask(rp, supportp)
  
  rp.SEUS <- rp.GSL <- rp
  
  # cabot.xy <- sp::SpatialPoints(coords = cbind(1400, 1600), proj4string = narw_crs())
  # 
  # circular.d <- raster::distanceFromPoints(object = densr[[1]], xy = cabot.xy) |> 
  #   raster::mask(mask = densr[[1]])
  
  if(target == "SEUS"){
    
    # Set value to 1 in SEUS
    rp.SEUS[!rp.SEUS$layer == which(sort(reg$region) == "SEUS")] <- 0
    rp.SEUS[rp.SEUS>0] <- 1
    
  } else if(target == "GSL"){
    
    # Set value to 1 in GSL
    rp.GSL[!rp.GSL$layer == which(sort(reg$region) == "GSL")] <- 0
    rp.GSL[rp.GSL>0] <- 1 
    
  }
  
  # # Set value to 1 in CABOT
  # rp.CA[rp.CA$layer == which(sort(reg$region) == "GSL")] <- 99
  # rp.CA[rp.CA$layer == which(sort(reg$region) == "CABOT")] <- 99
  # rp.CA[rp.CA$layer == which(sort(reg$region) == "SCOS")] <- 99
  # rp.CA[rp.CA<99] <- 0
  # rp.CA[rp.CA==99] <- 1
  
  # Create latitudinal gradient
  r_y <- sp::coordinates(rp.GSL, na.rm = TRUE)
  r_y <- raster::rasterFromXYZ(cbind(r_y, r_y[, 2, drop = FALSE]))
  r_y <- raster::mask(r_y, densr[[1]])
  
  r_north <- 1/(1+exp(-0.004*(r_y-1380)))
  r_south <- 1/(1+exp(0.004*(r_y-12)))
  
  r_southbound <- 1/(1+exp(0.004*(r_y-500)))
  r_northbound <- 1/(1+exp(-0.004*(r_y-800)))
  
  # Breeding season in the SEUS (Nov through to February - Krystan et al., 2018)
  # Foraging season in the GSL (June though to Oct – Crowe et al., 2021)
  
  wdens <- purrr::map(.x = seq_along(densr), .f = ~{
    
    if(.x == 11){ # Migrating south in November
      
      if(target == "SEUS"){
        
        # OPTION 1
        
        if(option == 1){
          
          # # Determine mass outside of SEUS
          mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.SEUS)), na.rm = TRUE)
          # Redistribute mass across range according to north-south gradient
          d1 <- densr[[.x]] * rp.SEUS
          d2 <- mout * r_south/sum(raster::getValues(r_south), na.rm = TRUE) # Sums to mout
          d.out <- d1 + d2
          
        } else if (option == 2) {
          
          # OPTION 2
          
          # Clip density to SEUS
          d1 <- densr[[.x]] * rp.SEUS
          # Re-scale gradient outside of SEUS
          d1[d1==0] <- NA
          dv <- na.omit(apply(X = raster::as.matrix(d1), MARGIN = 1, FUN = mean, na.rm = TRUE))
          d2 <- rescale_raster(x = r_south * (1-rp.SEUS), new.max = dv[1])
          d2[d2==0] <- NA
          d.out <- raster::merge(d1, d2)
          
        } else if (option == 3){
          
          # OPTION 3
          pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs())
          pts <- sp::spsample(pathGSL, n = 1000, type = "regular")
          distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |>
            raster::mask(mask = densr[[1]])
          dist.r <- (1/distances)^(1/4) # Fourth root
          d.out <- rescale_raster(dist.r, new.min = raster::cellStats(dist.r, "min"),
                                  new.max = raster::cellStats(dist.r*(1-rp.CA), "max")) + r_southbound
          
        }  
        
        as(d.out, "SpatialGridDataFrame")
        
      } else {
        
        dens[[.x]]  
        
      }
      
    } else if(.x == 6) { # Migrating north in June
      
      if(target == "GSL"){
        
        if(option == 1){
          
          # OPTION 1
          
          # Determine mass outside of GSL
          mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.GSL)), na.rm = TRUE)
          # Redistribute mass across range according to north-south gradient
          d1 <- densr[[.x]] * rp.GSL
          d2 <- mout * r_north/sum(raster::getValues(r_north), na.rm = TRUE) # sums to mout
          d.out <- d1 + d2
          
        } else if (option == 2){
          
          # OPTION 2
          
          # Clip density to GSL
          d1 <- densr[[.x]] * rp.GSL
          # Re-scale gradient outside of GSL
          d1[d1==0] <- NA
          dv <- na.omit(apply(X = raster::as.matrix(d1), MARGIN = 1, FUN = mean, na.rm = TRUE))
          d2 <- rescale_raster(x = r_north * (1-rp.GSL), new.max = dv[1])
          d2[d2==0] <- NA
          d.out <- raster::merge(d1, d2)
          
        } else if (option == 3){
          
          # OPTION 3
          pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs()) |> rgeos::gBuffer(width = 50)
          pts <- sp::spsample(pathGSL, n = 10000, type = "regular")
          distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |> raster::mask(mask = densr[[1]])
          dist.r <- (1/distances)^(1/4) # Fourth root
          d.out <- rescale_raster(dist.r, new.min = raster::cellStats(dist.r, "min"), 
                                  new.max = raster::cellStats(dist.r*rp.CA, "max")) + r_northbound
          
        }
        
        as(d.out, "SpatialGridDataFrame")
        
      } else {
        
        dens[[.x]]  
        
      }
      
    } else if(.x %in% c(1:2, 12)) { # Breeding season in the SEUS
      
      if(target == "SEUS"){
        
        # Mass outside of SEUS
        mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.SEUS)), na.rm = TRUE)
        
        # Redistribute mass
        d <- densr[[.x]] * rp.SEUS
        d <- d + d*mout/sum(raster::getValues(d), na.rm = TRUE) + 1e-06
        as(d, "SpatialGridDataFrame")
        
      } else {
        
        dens[[.x]]  
        
      }
      
    } else if(.x %in% c(7:10)) { # Foraging season in the GSL
      
      if(target == "GSL"){
        
        # Mass outside of GSL
        mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.GSL)), na.rm = TRUE)
        
        # Redistribute mass
        d <- densr[[.x]] * rp.GSL
        d <- d + d*mout/sum(raster::getValues(d), na.rm = TRUE) + 1e-06
        as(d, "SpatialGridDataFrame")
        
      } else {
        
        dens[[.x]]  
        
      }
      
    } else { # March through to May
      
      dens[[.x]]
      
    }
    
  }) |> purrr::set_names(nm = month.abb)
  
  return(wdens)
  
}

warp <- function(p, coords, lower = 700, upper = 850, left = 500, right = 750){
  
  inside <- which(coords[,2] > lower & coords[,2] < upper & coords[,1] > left & coords[,1] < right)
  outside <- which(!1:nrow(coords) %in% inside)
  N.outside <- sum(p[outside])
  p[outside] <- 0
  p[inside][p[inside] > 0] <- p[inside][p[inside] > 0] + (N.outside / length(p[inside][p[inside] > 0]))
  return(p)
}

initiate_xy <- function(maplist,
                         coords,
                         migrate,
                         init.month,
                         nsim){
  
  out <- do.call(rbind, lapply(1:12, function(mo) {
    
    sampled.xy <- purrr::map(.x = maplist, .f = ~{
      sample(x = 1:283077, size = nsim, replace = TRUE, prob = .x[[mo]])
    })
    
    purrr::map_int(.x = 1:nsim, .f = ~ {
      if (mo %in% c(1:2, 11:12) & migrate[.x, 1] == 1) {
        sampled.xy[[1]][.x]
      } else if (mo %in% c(6:9) & migrate[.x, 2] == 1) {
        sampled.xy[[2]][.x]
      } else if(mo == init.month){
        sampled.xy[[3]][.x]
      } else {
        sampled.xy[[4]][.x]
      }
    })
    
  }))
  row.names(out) <- month.abb
  return(out)
}
# 
# 
# init_xy <- function(maps, 
#                     coords, 
#                     cohort.id, 
#                     nsim){
#   
#   out <- do.call(rbind, lapply(seq_along(maps), function(mo) {
#     
#     p <- as.numeric(maps[[mo]])
#     p[!is.finite(p)] <- 0
#     
#     # Lactating females
#     if(cohort.id == 5){
#       
#       # Breeding season (Nov through to February - Krystan et al. 2018)
#       if(mo %in% c(1:2, 11:12)){
#         p <- warp(p, coords, -5000, -300, -1000, 3000)
#       } else {
#         p <- warp(p, coords)
#       }
#       
#     } else {
#       
#       if(mo %in% c(6:9)){
#         # Initiate on feeding grounds in GSL
#         p <- warp(p, coords, 1400, 5000, -1000, 3000)
#       } else {
#         p <- warp(p, coords)
#       }
#       
#     }
# 
#     sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
#   }))
#   row.names(out) <- month.abb
#   return(out)
# }

# init_xy <- function(maps, 
#                     # maps.weighted, 
#                     coords, 
#                     cohort.id, 
#                     nsim, 
#                     northSEUS = -12, 
#                     southGSL = 1400){
#   
#   m <- maps
#   # if(cohort.id == 5) m <- maps.weighted else m <- maps
#   # if(cohort.id > 0) m <- maps.weighted else m <- maps
#   
#   out <- do.call(rbind, lapply(seq_along(m), function(mo) {
#     
#     p <- as.numeric(m[[mo]])
#     p[!is.finite(p)] <- 0
#     
#     # Foraging season in the GSL (June though to Oct – Crowe et al., 2021)
#     if(mo %in% c(7:10)){
#       notGSL <- which(coords[,2] < southGSL)
#       GSL <- which(coords[,2] >= southGSL)
#       N.notGSL <- sum(p[notGSL])
#       p[notGSL] <- 0
#       p[GSL][p[GSL] > 0] <- p[GSL][p[GSL] > 0] + (N.notGSL / length(p[GSL][p[GSL] > 0]))
#       
#       # Breeding season (Nov through to February - Krystan et al. 2018)
#     } else if (mo %in% c(1:2, 11:12)){
#       notSEUS <- which(coords[,2] > northSEUS)
#       SEUS <- which(coords[,2] <= northSEUS)
#       N.notSEUS <- sum(p[notSEUS])
#       p[notSEUS] <- 0
#       p[SEUS][p[SEUS] > 0] <- p[SEUS][p[SEUS] > 0] + (N.notSEUS / length(p[SEUS][p[SEUS] > 0]))
#       
#     } else if (mo == 5){
#       
#       area <- which(coords[,2] <= 1200)
#       notarea <- setdiff(seq_along(p), area)
#       N.notarea <- sum(p[notarea])
#       p[notarea] <- 0
#       p[area][p[area] > 0] <- p[area][p[area] > 0] + (N.notarea / length(p[area][p[area] > 0]))
#       
#     } else if (mo == 6){
#       
#       SCOS <- which(coords[,1] >= 1400 & coords[,2] <=1500)
#       notSCOS <- setdiff(seq_along(p), SCOS)
#       N.notSCOS <- sum(p[notSCOS])
#       p[notSCOS] <- 0
#       p[SCOS][p[SCOS] > 0] <- p[SCOS][p[SCOS] > 0] + (N.notSCOS / length(p[SCOS][p[SCOS] > 0]))
#       
#     } else {
#       p[p <= 0.01] <- 0
#     }
#     sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
#   }))
#   
#   return(out)
# }

# 
# init_inds <- function(x, nsim, northSEUS = -12, southGSL = 1500, cohort.id){
#   # -12 northing marks the northermost boundary of the calving grounds
#   coords <- sp::coordinates(density_narw[[1]])
#   out <- do.call(rbind, lapply(seq_along(x), function(m) {
#     
#     p <- as.numeric(x[[m]])
#     p[!is.finite(p)] <- 0
#     
#     # Pregnant females
#     if(cohort.id == 5){
#       
#       # During the breeding season (Nov through to February - Krystan et al. 2018)
#       if(m %in% c(1:2, 11:12)){
#         notSEUS <- which(coords[,2] > northSEUS)
#         SEUS <- which(coords[,2] <= northSEUS)
#         N.notSEUS <- sum(p[notSEUS])
#         p[notSEUS] <- 0
#         p[SEUS][p[SEUS] > 0] <- p[SEUS][p[SEUS] > 0] + (N.notSEUS / length(p[SEUS][p[SEUS] > 0]))
#       }
#       
#       # Foraging season in the GSL (June though to Oct – Crowe et al., 2021)
#       if(m %in% c(6:10)){
#         notGSL <- which(coords[,2] < southGSL)
#         GSL <- which(coords[,2] >= southGSL)
#         N.notGSL <- sum(p[notGSL])
#         p[notGSL] <- 0
#         p[GSL][p[GSL] > 0] <- p[GSL][p[GSL] > 0] + (N.notGSL / length(p[GSL][p[GSL] > 0]))
#       }
#       
#     }
#     sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
#   }))
#   row.names(out) <- month.abb
#   # plot(world, axes = TRUE, col = "grey")
#   # points(raster::xyFromCell(raster::raster(density_narw$Jan), out[,1]), pch = 16)
#   return(out)
# 
# }

transpose_array <- function(input, cohortID, dates) {
  
  out.array <- purrr::map2(.x = input, .y = cohortID, .f = ~ {
    
    tp <- lapply(X = .x, FUN = function(l) {
      tmp <- aperm(l, c(3, 1, 2))
      if(dim(tmp)[3] > 1){
      dimnames(tmp)[[1]] <- format(dates)
      dimnames(tmp)[[3]] <- paste0("whale.", seq_len(dim(tmp)[3]))
      }
      tmp
    })
    
    outl <- list()
    outl[["locs"]] <- tp[["locs"]]
    outl[["stress"]] <- tp[["stressors"]]
    outl[["activ"]] <- tp[["activity"]]
    
    if(.y == 5){ # Lactating females / calves
      
      outl[["E"]] <- list(adjuv = tp[["E"]], calves = tp[["E_calves"]])
      
      tp[grepl("fetus", names(tp))] <- NULL  
      outl[["inits"]] <- list(adjuv = tp[["inits"]], calves = tp[["inits_calves"]])
      outl[["attrib"]] <- list(adjuv = tp[["attrib"]], calves = tp[["attrib_calves"]])
      outl[["kj"]] <- list(adjuv = do.call(abind::abind, list(tp[["in.kj"]], tp[["out.kj"]], along = 2)),
                           calves = do.call(abind::abind, list(tp[["in.kj_calves"]], tp[["out.kj_calves"]], along = 2)))
      
    } else if(.y == 4){ # Pregnant females
      
      outl[["E"]] <- tp[["E"]]
      outl[["inits"]] <- tp[["inits"]]
      tp[grepl("calves", names(tp))] <- NULL  
      outl[["attrib"]] <- list(adjuv = tp[["attrib"]], fetus = tp[["attrib_fetus"]])
      outl[["kj"]] <- do.call(abind::abind, list(tp[grepl("kj", names(tp))], along = 2))
      

    } else {
      
      outl[["E"]] <- tp[["E"]]
      outl[["inits"]] <- tp[["inits"]]
      tp[grepl("fetus", names(tp))] <- NULL  
      tp[grepl("calves", names(tp))] <- NULL  
      outl[["attrib"]] <- tp[["attrib"]]
      outl[["kj"]] <- do.call(abind::abind, list(tp[grepl("kj", names(tp))], along = 2))
    }
    
    outl
    
  })
  
  return(out.array)
}

consolidate <- function(dtl, nsim, cnames, dates, months){
  purrr::map2(.x = dtl,
              .y = cnames, 
              .f = ~{
                dt <- data.table::data.table(day = rep(0:365, times = nsim), 
                                             date = rep(dates, times = nsim),
                                             month = rep(months, times = nsim),
                                             whale = rep(1:nsim, each = 366))
                a <- cbind(array2dt(.x), dt)
                a$region <- sort(regions$region)[a$region]
                a$cohort_name <- .y
                data.table::setcolorder(a, c((ncol(a)-4):(ncol(a)-1), 1:(ncol(a)-5))) 
              })
}


# transpose_array2 <- function(input, array.name, dates){
#   purrr::map(.x = input, .f = ~{ # different cohorts
# 
#     dat <- .x[[array.name]]
# 
#     if(array.name == "attrib" & is.list(dat)){
#       tmp <- purrr::map(.x = dat, .f = ~{
# 
#         tmp.nested <- aperm(.x, c(3,1,2))
#         dimnames(tmp.nested)[[1]] <- format(dates)
#         dimnames(tmp.nested)[[3]] <- paste0("whale.", seq_len(dim(tmp.nested)[3]))
#         tmp.nested
# 
#       })
# 
#     } else if(array.name == "kj" & is.list(dat)) {
# 
#       if(class(dat[["in"]]) == "array"){
# 
#         tmp.in <- aperm(dat[["in"]], c(3,1,2))
#         dimnames(tmp.in)[[1]] <- format(dates)
#         dimnames(tmp.in)[[3]] <- paste0("whale.", seq_len(dim(tmp.in)[3]))
# 
#         tmp.out <- aperm(dat[["out"]], c(3,1,2))
#         dimnames(tmp.out)[[1]] <- format(dates)
#         dimnames(tmp.out)[[3]] <- paste0("whale.", seq_len(dim(tmp.out)[3]))
# 
#       } else {
# 
#         tmp.in <- purrr::map(.x = dat[["in"]], .f = ~{
#           tmp <- aperm(.x, c(3,1,2))
#           dimnames(tmp)[[1]] <- format(dates)
#           dimnames(tmp)[[3]] <- paste0("whale.", seq_len(dim(tmp)[3]))
#           tmp
#         })
# 
#         tmp.out <- purrr::map(.x = dat[["out"]], .f = ~{
#           tmp <- aperm(.x, c(3,1,2))
#           dimnames(tmp)[[1]] <- format(dates)
#           dimnames(tmp)[[3]] <- paste0("whale.", seq_len(dim(tmp)[3]))
#           tmp
#         })
#       }
# 
#       tmp <- list(`in` = tmp.in, out = tmp.out)
# 
#     } else {
#       tmp <- aperm(dat, c(3,1,2))
#       dimnames(tmp)[[1]] <- format(dates)
#       dimnames(tmp)[[3]] <- paste0("whale.", seq_len(dim(tmp)[3]))
#     }
# 
#     .x[[array.name]] <- tmp
#     .x
#   })
# }

optim_feeding <- function(bounds = list(c(0,0.2)), nm = NULL, linecol = "black", verbose = TRUE){
  
  # Assuming B > 0 And S > 0:
  # - A is the value of the horizontal asymptote when x tends to -Inf
  # - D is the value of the horizontal asymptote when x tends to Inf
  # - B describes how rapidly the curve makes its transition between the two asymptotes;
  # - C is a location parameter, which does not have a nice interpretation (except if S = 1)
  # - S describes the asymmetry of the curve (the curve is symmetric when S = 1)
  
  # Define a 5-parameter logistic curve
  f <- function(x, A = 1, D = 0, B, C, S) A + (D-A) / (1 + exp(B*(C-x)))^S
  
  # Define function to optimize
  opt.f <- function(param, bounds, err = 0.001){
    first <- (1 - f(bounds[1], B = param[1], C = param[2], S = param[3]) - err)^2
    second <- (0 - f(bounds[2], B = param[1], C = param[2], S = param[3]) +  err)^2
    sum(first, second)
  }
  
  # Run optimization
  res <- purrr::map(.x = seq_along(bounds), 
             .f = ~ {
               out <- optim(par = c(15, 1, 0.1, 0.01), fn = opt.f, bounds = bounds[[.x]])
               
               # Print results
               if(verbose){
               cat("\nEstimated parameter values:\n")
               cat("B",out$par[1],"\n")
               cat("C",out$par[2],"\n")
               cat("s",out$par[3],"\n\n")
               print(data.frame(BC = bounds[[.x]], 
                     effort = c(f(bounds[[.x]][1], 
                                  B = out$par[1], 
                                  C = out$par[2], 
                                  S = out$par[3]), 
                                f(bounds[[.x]][2], 
                                  B = out$par[1],
                                  C = out$par[2],
                                  S = out$par[3]))))
               }
               return(out)
             })

  # Plot
  x <- seq(0, 1, length.out = 100)
  plot(x, f(x, B = res[[1]]$par[1], C = res[[1]]$par[2], S = res[[1]]$par[3]), 
       type = "n", ylab = "Feeding effort", xlab = "Body condition (relative blubber mass, %)")
  purrr::walk(.x = seq_along(bounds), .f = ~{
    lines(x, f(x, B = res[[.x]]$par[1], C = res[[.x]]$par[2], S = res[[.x]]$par[3]), lty = .x, col = linecol)
  })
  if(!is.null(nm)) legend("topright", legend = nm, lty = seq_along(bounds), col = linecol)
  
}

# Function to add new individual
new_calf <- function(arr, year, attrib){
  
  if(!ncol(arr) == length(mat.attribs)) stop("Cannot add new individual")
  
  newind <- array(data = NA, dim = c(dim(arr)[1:2], 1), dimnames = list(
    paste0("yr ", 0:yrs), mat.attribs, paste0("whale ", dim(arr)[3] + 1)))
  
  newind[year,"alive", 1] <- 1
  newind[year,"cohort", 1] <- 0
  newind[year,"female", 1] <- rbinom(n = 1, size = 1, prob = 0.5) # 1:1 sex ratio at birth
  newind[year, "age", 1] <- 0
  
  l.params <- agL(0)
  newind[year,"length", 1] <- age2length(0, l.params)
  newind[year:dim(arr)[1],"length_a", 1] <- l.params[,1]
  newind[year:dim(arr)[1],"length_b", 1] <- l.params[,2]
  newind[year:dim(arr)[1],"length_c", 1] <- l.params[,3]
  
  m.params <- mL(1)
  mass <- length2mass(newind[year,"length", 1], m.params, FALSE)
  newind[year:dim(arr)[1],"mass_a", 1] <- m.params[, 1]
  newind[year:dim(arr)[1],"mass_b", 1] <- m.params[, 2]
  newind[year,"tot_mass", 1] <- mass
  
  newind[year,"bc",] <- start_bcondition(0)
  
  newind[year,"lean_mass", 1] <- mass - (newind[year,"bc",] * mass)

  newind[year,"trest",] <- 0
  newind[year,"t2calf",] <- 0
  newind[year,"min_bc",] <- 0
  newind[year,"birth",] <- 0
  newind[year,"p_surv",] <- 1
  
  #' ---------------------------
  
  return(abind::abind(arr, newind, along = 3))
  
}

predict_backup <- function(obj,
         yrs = 35,
         n = 100,
         popr = 1,
         do.plot = FALSE,
         seed = 125897,
         ...) {
  
  set.seed(seed)
  
  if(sum(suppressWarnings(purrr::map_lgl(.x = obj$gam$pred, .f = ~any(is.na(.x)))) > 0)) 
    stop("Insufficient sample size. Cannot make predictions.") 
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values
  if(length(args) == 0){
    spline <- TRUE
    progress <- TRUE
  } else {
    if("spline" %in% names(args)) spline <- args[["spline"]] else spline <- TRUE
    if("progress" %in% names(args)) progress <- args[["progress"]] else progress <- TRUE
  }
  
  cat("Initializing ...\n")
  
  # if(is.null(obj$gam)) stop("Insufficient data available. Cannot proceed with population projections.")
  # if(!identical(cohortID, 1:6)) stop("Missing cohorts in input <obj>. Cannot proceed with population projections.")
  # if(length(obj$gam$fit$surv$xlevels[[1]]) < 6 | length(obj$gam$fit$bc$xlevels[[1]]) < 6) stop("Missing factor levels in input <obj>. Cannot proceed with population projections.")
  
  # plogis("link" predictions + error)
  
  # Adapted from original code by Scott Creel
  # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
  # 
  # Prediction intervals
  # https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html
  
  # Population estimate as per 2022 NARW report card is 340 (+/- 7).
  
  # test <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = seq(0,0.6, 0.01)), type = "link", se.fit = TRUE)
  # plot(seq(0,0.6, 0.01), plogis(test$fit), type = "l")
  # lines(seq(0,0.6, 0.01), plogis(Reduce("+", test)), lty = 2)
  # lines(seq(0,0.6, 0.01), plogis(Reduce("-", test)), lty = 2)
  # test2 <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = 0.2), type = "response")
  # abline(v = 0.2)
  # abline(h = tesct2)
  
  cohortID <- obj$param$cohortID
  cohorts <- obj$param$cohorts |> dplyr::slice(1) |> 
    dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |> 
    dplyr::bind_rows(obj$param$cohorts) |> 
    dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
  cohorts <- cohorts[c(1,3,5,2,4,6,7,8)]
  
  # Attributes to monitor during projection
  mat.attribs <- c("alive", "cohort", "female", "age", "length", "length_a", "length_b", "length_c",
                   "tot_mass", "lean_mass", "bc", "mass_a", "mass_b", "p_surv", "min_bc", "trest", "t2calf", "birth")
  
  # Current year
  current.yr <- lubridate::year(lubridate::now())
  
  # Extract terminal functions
  mod <- obj$gam$fit
  mod[["gest"]] <- gam_gest
  
  #'------------------------------------------------------
  # GAM PARAMETERS
  #'------------------------------------------------------
  
  mbc_preds <- obj$gam$pred$bc_gest
  bc_preds <- obj$gam$pred$bc
  surv_preds <- obj$gam$pred$surv
  
  #'------------------------------------------------------
  # INITIALIZATION
  #'------------------------------------------------------
  
  # Define initial population vector
  # c(F1 = 0, F2 = 0, F3 = 2, F4 = 6, F5 = 4, F6 = 0, F7 = 6, F8 = 0, 
  # F9 = 5, "F10+" = 48, FC = 7, FR = 0, FB = 61, 
  # M1 = 0, M2 = 0, 
  # M3 = 2, M4 = 5, 
  # "M5+" = 212)
  N0 <- c(0, 7, 212, 0, 71, 1, 7, 61)
  names(N0) <- cohorts[, name]
  totn <- sum(N0)* (1 + popr)
  
  cat("Setting up ...\n")
  
  # Create matrices and initialize them
  # rows: years <yrs>
  # columns: <attributes>
  # layers: individuals <n>
  # 4th dimension: replicate projection -> later converted to list
  narw.indiv <- array(data = NA, c(yrs + 1, length(mat.attribs), totn, n), 
                      dimnames = list(paste0("yr ", 0:yrs), 
                                      mat.attribs,
                                      paste0("whale ", 1:totn),
                                      paste0("prj ", 1:n)))
  
  cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts), FUN = function(x) rep(cohorts$id[x], each = N0[x])))
  
  animals <- 1:sum(N0)
  
  # Alive and population cohort
  narw.indiv[1, "alive", animals, ] <- rep(1, )
  narw.indiv[1, "cohort", animals, ] <- rep(cohort.vec, n)
  
  # Sex
  #  -- Calves (male)
  narw.indiv[, "female", 1:N0[1], ] <- 0
  #  -- Calves (female)
  fem <- which(cohort.vec == 0)
  narw.indiv[, "female", fem[fem > N0[1]], ] <- 1
  #  -- Juveniles and adults
  narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6)), ] <- 1
  narw.indiv[, "female", which(cohort.vec %in% c(1, 3)), ] <- 0
  
  cat("Age ...\n")
  
  # Age
  ages <- start_age_vec(rep(cohort.vec, n))
  narw.indiv[1, "age", animals, ] <- ages
  
  cat("Length ...\n")
  
  # Total body length
  l.params <- agL_vec(ages)
  lengths <- age2length_vec(ages, l.params)
  narw.indiv[1, "length", animals, ] <- lengths
  
  narw.indiv[, "length_a", animals, ] <- rep(l.params[, 1], each = yrs + 1)
  narw.indiv[, "length_b", animals, ] <- rep(l.params[, 2], each = yrs + 1)
  narw.indiv[, "length_c", animals, ] <- rep(l.params[, 3], each = yrs + 1)
  
  cat("Mass ...\n")
  
  # Lean mass
  m.params <- mL(n * sum(N0))
  narw.indiv[, "mass_a", animals, ] <- rep(m.params[, 1], each = yrs + 1)
  narw.indiv[, "mass_b", animals, ] <- rep(m.params[, 2], each = yrs + 1)
  mass <- length2mass_vec(lengths, m.params)
  narw.indiv[1, "lean_mass", animals, ] <- mass
  
  # Body conditon
  bc <- start_bcondition_vec(rep(cohort.vec, n))
  narw.indiv[1, "bc", animals, ] <- bc
  
  # Total mass
  narw.indiv[1, "tot_mass", animals, ] <- mass / (1-bc)
  
  # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
  # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
  # This corresponds to a Normal (7.7, 1.45)
  # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
  # 
  # Stewart et al. 2022 -- 
  # The degree to which the energetic reserves of females are depleted during lactation may govern
  # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).
  
  cat("Reproduction ...\n")
  
  t2calf <- (rep(cohort.vec, n) == 6) * random_int(sum(N0) * n)
  narw.indiv[1, "t2calf", animals, ] <- t2calf
  narw.indiv[1, "trest", animals, ] <- as.numeric(ifelse(t2calf == 0, 13, 1)) * (as.numeric(narw.indiv[1, "cohort", animals, ]) == 6)
  
  if (!spline) {
    narw.indiv[1, "min_bc", animals, ] <- predict_m(
      model = mod,
      values = as.vector(narw.indiv[1, "tot_mass", animals, ]),
      prediction = "gest"
    ) * as.vector(narw.indiv[1, "cohort", animals, ] == 6)
  } else {
    narw.indiv[1, "min_bc", animals, ] <- mbc_preds(narw.indiv[1, "tot_mass", animals, ]) * (narw.indiv[1, "cohort", animals, ] == 6)
  }
  
  narw.indiv[1, "birth", animals, ] <- ifelse(narw.indiv[1, "trest", animals, ] == 13 & narw.indiv[1, "t2calf", animals, ] == 0, 1, 0)
  narw.indiv[1, "p_surv", animals, ] <- 1
  
  cat("List ...\n")
  
  #' ---------------------------
  # IMPORTANT 
  #' ---------------------------
  # Turn array into a list
  narw.indiv <- purrr::array_branch(narw.indiv, 4) # bottleneck
  
  cat("totpop ...\n")
  
  # Number of individuals in each cohort
  narw.pop <- array(
    data = NA, c(n, yrs + 1, nrow(cohorts)),
    dimnames = list(
      paste0("prj ", 1:n),
      paste0("yr ", 0:yrs),
      cohorts$name
    )
  )
  narw.pop[, 1, ] <- rep(N0, each = n)
  
  cat("totpopinit ...\n")
  
  # Total population size
  tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj", 1:n), paste0("yr ", 0:yrs)))
  tot.pop[, 1] <- sum(N0)
  
  #'------------------------------------------------------
  # RUN PROJECTIONS
  #'------------------------------------------------------
  
  # This uses nested loops. 
  # The prj loop (outermost loop) replicates the projection <n> times.
  # The i loop is next, and steps across all years of projection from an initial population vector.
  
  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n, clear = FALSE, width = 80
  )
  
  cat("Running projections ...\n")
  start.time <- Sys.time()
  
  for(prj in 1:n){
    
    if(progress) pb$tick() # Update progress bar
    
    animals <- 1:sum(N0)
    
    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      current.dat <- as.matrix(narw.indiv[[prj]][i-1, , animals])
      alive <- current.dat["alive", animals] * (current.dat["age", animals] <=69)
      
      #' ----------------------------
      # SURVIVAL
      #' ----------------------------
      # Predict survival probability based on body condition
      
      if(!spline){
        
        ps <- alive * predict_m(model = mod, cohort = current.dat["cohort",animals], 
                                values = current.dat["bc",animals], prediction = "surv")
        
      } else {
        
        ps <- alive * (surv_preds[["0"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 0) +
                         surv_preds[["1"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 1) +
                         surv_preds[["2"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 2) +
                         surv_preds[["3"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 3) +
                         surv_preds[["4"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 4) +
                         surv_preds[["5"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 5) +
                         surv_preds[["6"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 6))
      }
      
      # ps <- 1
      narw.indiv[[prj]][i, "p_surv", animals] <- ps
      
      # Determine whether the animal survived
      alive <- rbinom(n = animals, size = 1, prob = ps) * (current.dat["age", animals] <=69)
      narw.indiv[[prj]][i, "alive", animals] <- alive
      
      # Sex remains the same
      narw.indiv[[prj]][i, "female", animals] <- current.dat["female", animals]
      
      #' ----------------------------
      # GROWTH
      #' ----------------------------
      
      # Increment age
      narw.indiv[[prj]][i, "age", animals] <- alive * (current.dat["age", animals] + 1)
      
      # Increment length
      newLp <- agL_vec(animals)
      
      narw.indiv[[prj]][i,"length_a", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,1], current.dat["length_a", animals])
      
      narw.indiv[[prj]][i,"length_b", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,2], current.dat["length_b", animals])
      
      narw.indiv[[prj]][i,"length_c", animals] <- 
        ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,3], current.dat["length_c", animals])
      
      narw.indiv[[prj]][i, "length", animals] <-
        alive * age2length_vec(
          narw.indiv[[prj]][i, "age", animals],
          t(narw.indiv[[prj]][i, c("length_a", "length_b", "length_c"), animals])
        )
      
      # Increment lean mass
      narw.indiv[[prj]][i, "lean_mass", animals] <- alive * length2mass_vec(narw.indiv[[prj]][i, "length", animals],
                                                                            t(narw.indiv[[prj]][i, c("mass_a", "mass_b"), animals]), lean = TRUE)
      
      # Predict new body condition from current body condition
      if (!spline) {
        narw.indiv[[prj]][i, "bc", animals] <- alive * predict_m(
          model = mod, cohort = current.dat["cohort", animals],
          values = current.dat["bc", animals], prediction = "bc")
      } else {
        narw.indiv[[prj]][i, "bc", animals] <-
          alive * (bc_preds[["0"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 0) +
                     bc_preds[["1"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 1) +
                     bc_preds[["2"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 2) +
                     bc_preds[["3"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 3) +
                     bc_preds[["4"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 4) +
                     bc_preds[["5"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 5) +
                     bc_preds[["6"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 6))
      }
      
      # Increment total mass
      narw.indiv[[prj]][i, "tot_mass", animals] <-
        alive * narw.indiv[[prj]][i, "lean_mass", animals] / (1 - narw.indiv[[prj]][i, "bc", animals])
      
      #' ----------------------------
      # REPRODUCTION
      #' ----------------------------
      
      # Which animals are resting females?
      rest.females <- (current.dat["cohort", animals] == 6)
      
      # Which animals are juvenile females that are ready to start reproducing
      juvenile.females.ofage <- (current.dat["cohort", animals] == 2) * (current.dat["age", animals] >= 9)
      
      # Which animals calved in previous step?
      prev.births <- current.dat["birth", animals]
      
      newt2calf <- ifelse(prev.births == 1, random_int(sum(prev.births), lwr = 1), 
                          ifelse(current.dat["t2calf", animals] == 0, 0, 
                                 current.dat["t2calf", animals] - 1))
      
      # Time spent in resting state - only incremented if calving event hasn't occurred, otherwise reset
      narw.indiv[[prj]][i, "t2calf", animals] <- alive * rest.females * newt2calf
      
      # Years until next calving event
      narw.indiv[[prj]][i, "trest", animals] <- 
        alive * (rest.females | juvenile.females.ofage) * 
        ifelse(narw.indiv[[prj]][i - 1, "trest", animals] == 13, 1, current.dat["trest", animals] + 1)
      
      # Minimum body condition needed to successfully bring fetus to term without starving
      # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
      
      if (!spline) {
        narw.indiv[[prj]][i, "min_bc", animals] <-
          alive * predict_m(model = mod, values = narw.indiv[[prj]][i, "tot_mass", animals], prediction = "gest") * rest.females
      } else {
        narw.indiv[[prj]][i, "min_bc", animals] <-
          alive * mbc_preds(narw.indiv[[prj]][i, "tot_mass", animals]) * rest.females
      }
      
      # Birth of new calf, conditional on the mother being alive, in pregnant state
      narw.indiv[[prj]][i, "birth", animals] <- alive * (current.dat["cohort", animals] == 4)
      
      # Maturity - transitions between cohorts
      narw.indiv[[prj]][i, "cohort", animals] <-
        alive * increment_cohort(
          cohort = narw.indiv[[prj]][i - 1, "cohort", animals],
          age = narw.indiv[[prj]][i, "age", animals],
          female = narw.indiv[[prj]][i, "female", animals],
          bc = narw.indiv[[prj]][i, "bc", animals],
          min_bc = narw.indiv[[prj]][i, "min_bc", animals],
          trest = narw.indiv[[prj]][i, "trest", animals],
          t2calf = narw.indiv[[prj]][i, "t2calf", animals])
      
      new.births <- sum(narw.indiv[[prj]][i, "birth", animals])
      
      if(new.births > 0){
        narw.indiv[[prj]][i, , (max(animals)+1):(max(animals)+new.births)] <- add_calf(n = new.births, attr = mat.attribs)
        animals <- 1:(length(animals) + new.births)
      }
      
      # Number of animals in each cohort
      # Calves (male)
      narw.pop[prj, i, 1] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
                                   (narw.indiv[[prj]][i, "female", animals] == 0) *
                                   (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Juveniles and adults (male)
      narw.pop[prj, i, 2] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 1) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 3] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 3) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Calves (female)
      narw.pop[prj, i, 4] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
                                   (narw.indiv[[prj]][i, "female", animals] == 1) *
                                   (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Juvenile and reproductive adults (female)
      narw.pop[prj, i, 5] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 2) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 6] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 4) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 7] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 5) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      narw.pop[prj, i, 8] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 6) * (narw.indiv[[prj]][i, "alive", animals] == 1))
      
      # Total population size
      tot.pop[prj, i] <- sum(narw.indiv[[prj]][i, "alive", animals], na.rm = TRUE)
      
    } # End years
  } # End projections
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  cat("Processing outputs ...\n")
  
  # narw.out <- purrr::map(.x = 1:n, .f = ~{
  #   reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
  #     dplyr::mutate(attr = mat.attribs[attr]) |> 
  #     tidyr::pivot_wider(names_from = attr, values_from = value) |> 
  #     dplyr::mutate(prj = .x) |> 
  #     dplyr::relocate(prj, .before = yr)
  # }) |> do.call(what = rbind) |> 
  #   data.table::data.table()
  
  # narw.out <- narw.out[is.finite(rowSums(narw.out)),]
  
  narw.df <- purrr::map(.x = cohorts$name, .f = ~{
    narw.pop[,,.x] |> 
      tibble::as_tibble() |> 
      tibble::rownames_to_column(var = "prj") |> 
      tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
      dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
      dplyr::mutate(cohort = stringr::str_to_sentence(.x))
  }) |> do.call(what = rbind) |> data.table::data.table()
  
  tot.df <- tibble::as_tibble(tot.pop) |> 
    tibble::rownames_to_column(var = "prj") |> 
    tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
    dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
    dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
  
  # births.df <- purrr::map(.x = 1:n, .f = ~{
  #   m <- matrix(rowSums(narw.indiv[[.x]][2:(yrs+1),"birth",], na.rm = TRUE), ncol = 1)
  #   colnames(m) <- .x
  #   m
  # }) |> do.call(what = cbind) |> 
  #   tibble::as_tibble() |> 
  #   tibble::rownames_to_column(var = "year") |> 
  #   dplyr::mutate(year = as.numeric(year)) |> 
  #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |> 
  #   dplyr::select(prj, year, birth) |> 
  #   dplyr::arrange(prj, year) |> 
  #   data.table::data.table()
  # 
  # deaths.df <- purrr::map(.x = 1:n, .f = ~{
  #   m <- matrix(apply(X = narw.indiv[[.x]][2:(yrs+1),"alive",],
  #                     MARGIN = 1,
  #                     FUN = function(x) {
  #     r <- x[!is.na(x)]
  #     r <- sum(r == 0)
  #     r
  #     }), ncol = 1)
  #   colnames(m) <- .x
  #   m
  # }) |> do.call(what = cbind) |>
  #   tibble::as_tibble() |>
  #   tibble::rownames_to_column(var = "year") |>
  #   dplyr::mutate(year = as.numeric(year)) |>
  #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
  #   dplyr::select(prj, year, death) |>
  #   dplyr::arrange(prj, year) |>
  #   data.table::data.table()
  
  narw.conf <- narw.df[
    , list(
      mean = mean(N),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>
    dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)"
    )))
  
  tot.conf <- tot.df[
    , list(
      mean = mean(N),
      lwr = quantile(N, 0.025),
      uppr = quantile(N, 0.975)
    ),
    list(year, cohort)
  ] |>
    dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)",
      "North Atlantic right whales"
    )))
  
  p1 <- plot_projection(narw.df, narw.conf)
  p2 <- plot_projection(tot.df, tot.conf)
  
  if(do.plot){
    print(p1)
    print(p2)
  }
  
  # Find 95% confidence intervals on final population size
  cat("Final population size:\n")
  final.pop <- unname(tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.5, 0.025, 0.975))])
  cat("N = ", round(final.pop[1],0), " (95% CI: ", round(final.pop[2],0), "–", round(final.pop[3],0), ")\n", sep = "")
  
  cat(paste0("Time elapsed: ", run_time))
  cat("\n")
  
  return(list(
    # dat = narw.out,
    out = list(df = rbind(narw.df, tot.df), 
               ci = rbind(narw.conf, tot.conf),
               plot = list(p1, p2))))
  
}

# Function to predict survival or body condition from body condition
predict_m <- function(model, cohort = NULL, values, prediction = "surv"){
  if(prediction %in% c("surv", "bc")) newdat <- data.frame(start_bc = values, cohort = cohort)
  if(prediction == "gest") newdat <- data.frame(mass = values)
  mgcv::predict.gam(object = model[[prediction]], newdata = newdat, newdata.guaranteed = TRUE, type = "response")
}

#' predictnarwsim <- function(obj,
#'                             yrs = 35,
#'                             n = 100) {
#'   
#'   # if(is.null(obj$gam)) stop("Insufficient data available. Cannot proceed with population projections.")
#'   # if(!identical(cohortID, 1:6)) stop("Missing cohorts in input <obj>. Cannot proceed with population projections.")
#'   # if(length(obj$gam$fit$surv$xlevels[[1]]) < 6 | length(obj$gam$fit$bc$xlevels[[1]]) < 6) stop("Missing factor levels in input <obj>. Cannot proceed with population projections.")
#'   
#'   # plogis("link" predictions + error)
#'   
#'   # Adapted from original code by Scott Creel
#'   # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
#'   # 
#'   # Prediction intervals
#'   # https://stat.ethz.ch/pipermail/r-help/2011-April/275632.html
#'   
#'   # Population estimate as per 2022 NARW report card is 340 (+/- 7).
#'   
#'   # test <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = seq(0,0.6, 0.01)), type = "link", se.fit = TRUE)
#'   # plot(seq(0,0.6, 0.01), plogis(test$fit), type = "l")
#'   # lines(seq(0,0.6, 0.01), plogis(Reduce("+", test)), lty = 2)
#'   # lines(seq(0,0.6, 0.01), plogis(Reduce("-", test)), lty = 2)
#'   # test2 <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = 0.2), type = "response")
#'   # abline(v = 0.2)
#'   # abline(h = tesct2)
#'   
#'   cohortID <- obj$param$cohortID
#'   cohorts <- obj$param$cohorts |> dplyr::slice(1) |> 
#'     dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |> 
#'     dplyr::bind_rows(obj$param$cohorts) |> 
#'     dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
#'   cohorts <- cohorts[c(1,3,5,2,4,6,7,8)]
#'   
#'   # Attributes to monitor during projection
#'   mat.attribs <- c("alive", "cohort", "female", "age", "length", "length_a", "length_b", "length_c",
#'                    "tot_mass", "lean_mass", "bc", "mass_a", "mass_b", "p_surv", "min_bc", "calving_int", "t2calf", "birth")
#'   
#'   # Current year
#'   current.yr <- lubridate::year(lubridate::now())
#'   
#'   # Extract terminal functions
#'   mod <- obj$gam$fit
#'   mod[["gest"]] <- gam_gest
#'   
#'   #'------------------------------------------------------
#'   # INITIALIZATION
#'   #'------------------------------------------------------
#'   
#'   # Define initial population vector
#'   N0 <- c(2, 5, 212, 2, 69, 1, 7, 61)
#'   names(N0) <- cohorts[, name]
#'   
#'   # Create matrices and initialize them
#'   # rows: years <yrs>
#'   # columns: <attributes>
#'   # layers: individuals <n>
#'   # 4th dimension: replicate projection -> later converted to list
#'   narw.indiv <- array(data = NA, c(yrs + 1, length(mat.attribs), sum(N0), n), 
#'                       dimnames = list(paste0("yr ", 0:yrs), 
#'                                       mat.attribs,
#'                                       paste0("whale ", 1:sum(N0)),
#'                                       paste0("prj ", 1:n)))
#'   
#'   cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts), FUN = function(x) rep(cohorts$id[x], each = N0[x])))
#'   
#'   # Alive and population cohort
#'   narw.indiv[1,"alive",,] <- 1
#'   narw.indiv[1,"cohort",,] <- rep(cohort.vec, n)
#'   
#'   # Sex
#'   #  -- Calves (male)
#'   narw.indiv[,"female", 1:N0[1], ] <- 0  
#'   #  -- Calves (female)
#'   fem <- which(cohort.vec == 0)
#'   narw.indiv[,"female", fem[fem > N0[1]], ] <- 1   
#'   #  -- Juveniles and adults
#'   narw.indiv[,"female", which(cohort.vec %in% c(2,4,5,6)), ] <- 1
#'   narw.indiv[,"female", which(cohort.vec %in% c(1,3)), ] <- 0
#'   
#'   # Age
#'   ages <- start_age_vec(rep(cohort.vec, n))
#'   narw.indiv[1,"age", , ] <- ages
#'   
#'   # Total body length
#'   l.params <- agL_vec(ages)
#'   lengths <- age2length_vec(ages, l.params)
#'   narw.indiv[1,"length", , ] <- lengths
#'   narw.indiv[,"length_a",,] <- rep(l.params[,1], each = yrs + 1)
#'   narw.indiv[,"length_b",,] <- rep(l.params[,2], each = yrs + 1)
#'   narw.indiv[,"length_c",,] <- rep(l.params[,3], each = yrs + 1)
#'   
#'   # Total mass
#'   m.params <- mL(n*sum(N0))
#'   mass <- length2mass_vec(lengths, m.params, FALSE)
#'   narw.indiv[,"mass_a",,] <- rep(m.params[,1], each = yrs + 1)
#'   narw.indiv[,"mass_b",,] <- rep(m.params[,2], each = yrs + 1)
#'   narw.indiv[1,"tot_mass", , ] <- mass
#'   
#'   # Body conditon
#'   narw.indiv[1,"bc",,] <- start_bcondition_vec(rep(cohort.vec, n))
#'   
#'   # lean mass
#'   narw.indiv[1,"lean_mass", , ] <- mass - (narw.indiv[1,"bc",,] * mass)
#'   
#'   # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
#'   # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
#'   # This corresponds to a Normal (7.7, 1.45)
#'   # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
#'   # 
#'   # Stewart et al. 2022 -- 
#'   # The degree to which the energetic reserves of females are depleted during lactation may govern
#'   # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).
#'   narw.indiv[1,"calving_int",,] <- (rep(cohort.vec, n) == 6) * rnorm(sum(N0)*n, mean = 7.7, sd = 1.45)
#'   narw.indiv[1,"t2calf",,] <- round(narw.indiv[1,"calving_int",,],0)
#'   
#'   narw.indiv[1,"min_bc",,] <- predict_m(model = mod,
#'                                         values = as.vector(narw.indiv[1,"tot_mass",,]), 
#'                                         prediction = "gest") * as.vector(narw.indiv[1, "cohort",, ] == 6)
#'   narw.indiv[1,"birth",,] <- 0
#'   narw.indiv[1,"p_surv",,] <- 1
#'   
#'   #' ---------------------------
#'   # IMPORTANT 
#'   #' ---------------------------
#'   # Turn array into a list - otherwise cannot add new calves within projections
#'   narw.indiv <- purrr::array_branch(narw.indiv, 4) 
#'   
#'   # Number of individuals in each cohort
#'   narw.pop <- array(data = NA, c(n, yrs + 1, length(N0)), 
#'                     dimnames = list(paste0("prj ", 1:n),
#'                                     paste0("yr ", 0:yrs), 
#'                                     cohorts$name))
#'   narw.pop[,1,] <- rep(N0, each = n)
#'   
#'   # Total population size
#'   tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj",1:n), paste0("yr ", 0:yrs)))
#'   tot.pop[,1] <- sum(N0)
#'   
#'   #'------------------------------------------------------
#'   # RUN PROJECTIONS
#'   #'------------------------------------------------------
#'   
#'   # This uses nested loops. 
#'   # The prj loop (outermost loop) replicates the projection <n> times.
#'   # The i loop is next, and steps across all years of projection from an initial population vector.
#'   
#'   # Set up progress bar
#'   pb <- progress::progress_bar$new(
#'     format = "[:bar] :percent eta: :eta",
#'     total = n, clear = FALSE, width = 80
#'   )
#'   
#'   for(prj in 1:n){
#'     
#'     pb$tick() # Update progress bar
#'     
#'     #'------------------------------------------------------
#'     # Loop over years
#'     #'------------------------------------------------------
#'     
#'     for(i in 2:(yrs+1)){
#'       
#'       current.dat <- as.matrix(narw.indiv[[prj]][i-1, , ])
#'       alive <- narw.indiv[[prj]][i-1, "alive", ]
#'       
#'       #' ----------------------------
#'       # SURVIVAL
#'       #' ----------------------------
#'       # Predict survival probability based on body condition
#'       ps <- alive * predict_m(model = mod, cohort = current.dat["cohort",], values = current.dat["bc",], prediction = "surv")
#'       narw.indiv[[prj]][i, "p_surv", ] <- ps
#'       
#'       # Determine whether the animal survived
#'       # Maximum longevity of 69 years
#'       alive <- rbinom(n = dim(narw.indiv[[prj]])[3], size = 1, prob = ps) * (narw.indiv[[prj]][i-1, "age", ] <=69)
#'       narw.indiv[[prj]][i, "alive", ] <- alive
#'       
#'       # Sex remains the same
#'       narw.indiv[[prj]][i, "female", ] <- narw.indiv[[prj]][i-1, "female", ]
#'       
#'       #' ----------------------------
#'       # GROWTH
#'       #' ----------------------------
#'       
#'       # Increment age
#'       narw.indiv[[prj]][i, "age", ] <- alive * narw.indiv[[prj]][i-1, "age", ] + 1
#'       
#'       # Increment length
#'       narw.indiv[[prj]][i,"length", ] <- alive * age2length_vec(narw.indiv[[prj]][i, "age", ],
#'                                                                 t(narw.indiv[[prj]][i, c("length_a", "length_b", "length_c"),]))
#'       
#'       # Increment lean mass
#'       narw.indiv[[prj]][i,"lean_mass", ] <- alive * length2mass_vec(narw.indiv[[prj]][i, "length", ],
#'                                                                     t(narw.indiv[[prj]][i, c("mass_a", "mass_b"),]), lean = TRUE)
#'       
#'       # Predict new body condition from current body condition
#'       narw.indiv[[prj]][i, "bc", ] <- alive * predict_m(model = mod, cohort = current.dat["cohort",],
#'                                                         values = current.dat["bc",], prediction = "bc")
#'       
#'       # Increment total mass
#'       narw.indiv[[prj]][i,"tot_mass", ] <- alive * narw.indiv[[prj]][i,"lean_mass", ] / (1 - narw.indiv[[prj]][i,"bc", ])
#'       
#'       #' ----------------------------
#'       # REPRODUCTION
#'       #' ----------------------------
#'       
#'       narw.indiv[[prj]][i, "calving_int", ] <- alive * narw.indiv[[prj]][i-1, "calving_int", ]
#'       narw.indiv[[prj]][i, "t2calf", ] <- ifelse(narw.indiv[[prj]][i-1, "t2calf", ] - 1 < 0, 0, narw.indiv[[prj]][i-1, "t2calf", ] - 1)
#'       
#'       # Minimum body condition needed to successfully bring fetus to term without starving
#'       # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
#'       narw.indiv[[prj]][i, "min_bc", ] <- alive * predict_m(model = mod, values = narw.indiv[[prj]][i, "tot_mass", ], 
#'                                                             prediction = "gest") * (narw.indiv[[prj]][i-1, "cohort", ] == 6)
#'       
#'       # Birth of new calf
#'       narw.indiv[[prj]][i, "birth", ] <- alive * (narw.indiv[[prj]][i, "bc", ] >= narw.indiv[[prj]][i, "min_bc", ]) * 
#'         (narw.indiv[[prj]][i-1, "cohort", ] == 6) * (narw.indiv[[prj]][i-1, "birth", ] == 0)
#'       
#'       # Maturity - transitions between cohorts
#'       narw.indiv[[prj]][i, "cohort", ] <- alive * increment_cohort(cohort = narw.indiv[[prj]][i-1, "cohort", ],
#'                                                                    age = narw.indiv[[prj]][i, "age", ],
#'                                                                    female = narw.indiv[[prj]][i, "female", ],
#'                                                                    t2calf = narw.indiv[[prj]][i, "t2calf", ])
#'       
#'       # New births
#'       # narw.indiv[[prj]] <- new_calf(narw.indiv[[prj]][,,], year = i, attrib = mat.attribs)
#'       
#'       # Number of animals in each cohort
#'       # Calves (male)
#'       narw.pop[prj, i, 1] <- sum((narw.indiv[[prj]][i, "cohort", ] == 0) *
#'                                    (narw.indiv[[prj]][i, "female", ] == 0) *
#'                                    (narw.indiv[[prj]][i, "alive", ] == 1))
#'       
#'       # Juveniles and adults (male)
#'       narw.pop[prj, i, 2] <- sum((narw.indiv[[prj]][i, "cohort", ] == 1) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       narw.pop[prj, i, 3] <- sum((narw.indiv[[prj]][i, "cohort", ] == 3) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       
#'       # Calves (female)
#'       narw.pop[prj, i, 4] <- sum((narw.indiv[[prj]][i, "cohort", ] == 0) *
#'                                    (narw.indiv[[prj]][i, "female", ] == 1) *
#'                                    (narw.indiv[[prj]][i, "alive", ] == 1))
#'       
#'       # Juvenile and reproductive adults (female)
#'       narw.pop[prj, i, 5] <- sum((narw.indiv[[prj]][i, "cohort", ] == 2) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       narw.pop[prj, i, 6] <- sum((narw.indiv[[prj]][i, "cohort", ] == 4) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       narw.pop[prj, i, 7] <- sum((narw.indiv[[prj]][i, "cohort", ] == 5) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       narw.pop[prj, i, 8] <- sum((narw.indiv[[prj]][i, "cohort", ] == 6) * (narw.indiv[[prj]][i, "alive", ] == 1))
#'       
#'       # Total population size
#'       tot.pop[prj, i] <- sum(narw.indiv[[prj]][i, "alive", ], na.rm = TRUE)
#'       
#'     } # End years
#'   } # End projections
#'   
#'   # Compile outputs
#'   narw.df <- purrr::map(.x = cohorts$name, .f = ~{
#'     narw.pop[,,.x] |> 
#'       tibble::as_tibble() |> 
#'       tibble::rownames_to_column(var = "prj") |> 
#'       tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
#'       dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
#'       dplyr::mutate(cohort = stringr::str_to_sentence(.x))
#'   }) |> do.call(what = rbind) |> data.table::data.table()
#'   
#'   tot.df <- tibble::as_tibble(tot.pop) |> 
#'     tibble::rownames_to_column(var = "prj") |> 
#'     tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
#'     dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
#'     dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
#'   
#'   narw.conf <- narw.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)] |> 
#'     dplyr::mutate(cohort = factor(cohort, levels = c("Calves (male)",
#'                                                      "Calves (female)",
#'                                                      "Juveniles (male)", 
#'                                                      "Juveniles (female)", 
#'                                                      "Adults (male)",
#'                                                      "Adults (female, pregnant)",
#'                                                      "Adults (female, lactating)",
#'                                                      "Adults (female, resting)")))
#'   
#'   tot.conf <- tot.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)]|> 
#'     dplyr::mutate(cohort = factor(cohort, levels = c("Calves (male)",
#'                                                      "Calves (female)",
#'                                                      "Juveniles (male)", 
#'                                                      "Juveniles (female)", 
#'                                                      "Adults (male)",
#'                                                      "Adults (female, pregnant)",
#'                                                      "Adults (female, lactating)",
#'                                                      "Adults (female, resting)",
#'                                                      "North Atlantic right whales")))
#'   
#'   p1 <- plot_projection(narw.df, narw.conf)
#'   p2 <- plot_projection(tot.df, tot.conf)
#'   
#'   print(p1)
#'   print(p2)
#'   
#'   # Find 95% confidence intervals on final population size
#'   cat("Final population size:\n")
#'   final.pop <- unname(tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.5, 0.025, 0.975))])
#'   cat("N = ", final.pop[1], " (95% CI: ", final.pop[2], "–", final.pop[3], ")\n", sep = "")
#'   
#'   return(list(ind = narw.indiv, df = rbind(narw.df, tot.df), ci = rbind(narw.conf, tot.conf)))
#'   
#' }

# Stripped down version of predict function for TP splines (from <mgcv>)
# predict_thinplate <- function(object, dat) {
#   x <- array(0, 0)
#   for (i in 1:object$dim)
#   {
#     xx <- dat[[object$term[i]]]
#     xx <- xx - object$shift[i]
#     if (i == 1) {
#       n <- length(xx)
#     } else
#     if (length(xx) != n) stop("arguments of smooth not same dimension")
#     if (length(xx) < 1) stop("no data to predict at")
#     x <- c(x, xx)
#   }
# 
#   by <- 0
#   by.exists <- FALSE
#   ## following used to be object$null.space.dim, but this is now *post constraint*
#   M <- mgcv:::null.space.dimension(object$dim, object$p.order[1])
# 
#   ind <- 1:object$bs.dim
# 
#   if (is.null(object$drop.null)) object$drop.null <- 0 ## pre 1.7_19 compatibility
# 
#   if (object$drop.null > 0) object$bs.dim <- object$bs.dim + M
# 
#   X <- matrix(0, n, object$bs.dim)
#   oo <- .C(mgcv:::C_predict_tprs, as.double(x), as.integer(object$dim), as.integer(n), as.integer(object$p.order[1]),
#     as.integer(object$bs.dim), as.integer(M), as.double(object$Xu),
#     as.integer(nrow(object$Xu)), as.double(object$UZ), as.double(by), as.integer(by.exists),
#     X = as.double(X)
#   )
#   X <- matrix(oo$X, n, object$bs.dim)
#   if (object$drop.null > 0) {
#     if (FALSE) { ## not param
#       X <- (X %*% object$P)[, ind] ## drop null space
#     } else { ## original
#       X <- X[, ind]
#       X <- sweep(X, 2, object$cmX)
#     }
#   }
#   X
# } ## Predict.matrix.tprs.smooth

# dh <- list(mass = 30000)
# j <- attr(gam_gest$smooth[[1]], "nCons")
# coefs <- coef(gam_gest)

# Function to make predictions from GAM model linking body mass to the minimum body condition
# required to successfully complete a pregnancy - stripped down version of mgcv::predict.gam 
# predict_gest <- function(object, data.list, qrc) 
# {
#   # object = gam_gest$smooth[[1]]
#   # pm <- mgcv:::Predict.matrix3(gam_gest$smooth[[1]], data)[["X"]]
#   # dk <- mgcv:::ExtractData(object, data, NULL)
#   
#   pm <- predict_thinplate(object, data.list)
#   # k <- ncol(pm)
#   # t(qr.qty(qrc, t(pm))[(j + 1):k, , drop = FALSE])
#   out <- t(qr.qty(qrc, t(pm)))
#   out[,1] <- 1
#   out
# }

predict_leslie <- function(obj,
                            n = 100,
                            yrs = 35) {
  
  gg.opts <- ggplot2::theme(axis.text = element_text(size = 10, color = "black"),
                            axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
                            axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                            axis.title = element_text(size = 12),
                            axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
                            axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
                            strip.background = element_rect(fill = "grey20"),
                            strip.text = element_text(colour = 'white', size = 12))
  
  
  # Adapted from original code by Scott Creel
  # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
  
  # Population estimate as per 2022 NARW report card is 340 (+/- 7).
  # The percentage of the population estimated 
  
  # cohort.id <- obj$param$cohort.id
  # cohort.ab <- obj$param$cohort.ab
  current.yr <- lubridate::year(lubridate::now())
  
  cohorts <- c("Calves (male)",
               "Juveniles (male)", 
               "Adults (male)",
               "Calves (female)",
               "Juveniles (female)", 
               "Adults (female, pregnant)",
               "Adults (female, lactating)",
               "Adults (female, resting)") |> 
    tolower() |> abbreviate(minlength = 6) |> tibble::enframe() |> 
    dplyr::rename(cohort = name, abbr = value)
  
  # Vector to store total population size for N replicate projections
  narwpop <- purrr::set_names(x = cohorts$cohort) |> 
    purrr::map(.f = ~matrix(0, n, yrs, dimnames = list(paste0("proj",1:n), paste0("yr ", 1:yrs))))
  
  totpop <- matrix(0, n, yrs, dimnames = list(paste0("proj",1:n), paste0("yr ", 1:yrs)))
  
  # 1: Calves (male)
  # 2: Juveniles (male)
  # 3: Adults (male)
  # 4: Calves (female)
  # 5: Juveniles (female)
  # 6: Adults (female, pregnant)
  # 7: Adults (female, lactating)
  # 8: Adults (female, resting)
  
  p.fem <- 0.5 # Fraction of females at birth
  p.birth <- 0.5 # Probability that a pregnant female will give birth to a calf (fecundity)
  p.survival <- rep(0.8, nrow(cohorts)) # Probability of survival
  names(p.survival) <- cohorts
  p.recruit <- c(0.7, 0.7) # Male, female
  tau_rp <- 0.7 # Transition prob from resting to pregnant
  tau_pl <- 0.7 # Transition prob from pregnant to lactating
  tau_lr <- 0.8 # Transition prob from lactating to resting
  
  popmat <- matrix(0, nrow(cohorts), nrow(cohorts))
  popmat[1, 6] <- p.survival[6] * p.birth * (1 - p.fem)
  
  popmat[2, 1] <- p.survival[1]
  popmat[2, 2] <- p.survival[2]
  
  popmat[3, 2] <- p.survival[2] * p.recruit[1]
  popmat[3, 3] <- p.survival[3]
  
  popmat[4, 4] <- p.survival[6] * p.birth * p.fem
  
  popmat[5, 5] <- p.survival[4]
  popmat[5, 6] <- p.survival[5] * (1 - p.recruit[2])
  
  popmat[6, 5] <- p.survival[5] * p.recruit[2]
  popmat[6, 6] <- p.survival[6] * (1 - tau_pl)
  popmat[6, 8] <- p.survival[8] * tau_rp
  
  popmat[7, 6] <- p.survival[6] * tau_pl
  popmat[7, 7] <- p.survival[7] * (1 - tau_lr)
  
  popmat[8, 7] <- p.survival[7] * tau_lr
  popmat[8, 8] <- p.survival[8] * (1 - tau_rp)
  
  
  popmat <- matrix(popmat, nrow = nrow(cohorts), byrow = FALSE, dimnames = list(cohorts$abbr, cohorts$abbr))
  popmat <- round(popmat, 5)
  
  # STOCHASTIC LESLIE PROJECTION
  #' ----------------------------
  # This uses 4 nested loops. 
  # The p loop (outermost loop) replicates the projection <n> times.
  # The y loop is next, and steps across all years of projection from an initial population vector.
  # The i and j loops are innermost and draw stochastic parameters (body condition and survival)
  
  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n, clear = FALSE, width = 80
  )
  
  for(p in 1:n){
    
    pb$tick() # Update progress bar
    
    # Matrix of age-class values, with one row per time step. 
    # The columns represent the numbers of individuals in each age class.
    N.mat <- matrix(NA, nrow = yrs, ncol = nrow(cohorts), dimnames = list(paste("yr", 1:yrs), cohorts$abbr))
    
    # Initial population vector
    # Sums to 340 (abundance estimate as of 2022).
    N.mat[1,] <- c(7, 10, 163, 8, 10, 50, 15, 77)
    
    # Begin loop for projections
    for(i in 2:yrs){ 		
      
      # Set up a temporary Leslie matrix for each iteration 
      #  with the correct mean fecundities and survival rates
      leslie.mat <- popmat 				
      
      # Randomly draw fecundities for each class
      # Mean of Poisson is fecundity value
      # for(j in 1:n.stages){
      #   leslie.mat[1,j] <- ifelse(popmat[1,j] == 0, 0, popmat[1,j] + rnorm(1, 0, 0.01))
      # }
      
      for(j in 1:nrow(cohorts)){
        for(k in 1:nrow(cohorts)){
          leslie.mat[k,j] <- ifelse(popmat[k,j] == 0, 0, popmat[k,j] + rnorm(1, 0, 0.01))
        }
      }
      
      # Randomly draw survival probabilities for each class.
      # n.ind is number of individuals to calculate live/die for
      # for(k in 1:(n.stages-1)){
      #   n.ind <- N.mat[(i-1),k]
      #   # Need ifelse statement to deal with the possibility that there are no individuals in that age class
      #   leslie.mat[(k+1),k] <- ifelse(n.ind > 1, rbinom(1,size = round(n.ind), p = les.mat[(k+1),k])/n.ind,0)
      # }
      
      # Matrix multiplication for next time-step.
      N.mat[i,] <- leslie.mat %*% N.mat[(i-1),]          
      
    } # End i loop over time
    
    for(k in seq_along(narwpop)){
      narwpop[[k]][p,] <- N.mat[,k]
    }
    totpop[p,] <- rowSums(N.mat)
    
  } # Ends p loop over replicate projections
  
  # Compile outputs
  narw.df <- purrr::imap(.x = narwpop, .f = ~{
    tibble::as_tibble(.x) |> 
      tibble::rownames_to_column(var = "proj") |> 
      tidyr::pivot_longer(!proj, names_to = "year", values_to = "N") |> 
      dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
      dplyr::mutate(cohort = stringr::str_to_sentence(.y))
  }) |> do.call(what = rbind) |> data.table::data.table()
  
  tot.df <- tibble::as_tibble(totpop) |> 
    tibble::rownames_to_column(var = "proj") |> 
    tidyr::pivot_longer(!proj, names_to = "year", values_to = "N") |> 
    dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
    dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
  
  narw.conf <- narw.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)]
  tot.conf <- tot.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)]
  
  # Generate plots
  make_plot <- function(df, conf){
    ggplot2::ggplot() +
      # ggplot2::geom_path(data = df, aes(x = year, y = N, group = factor(proj)), colour = "lightgrey", linewidth = 0.2) +
      ggplot2::geom_path(data = conf, aes(x = year, y = mean), colour = "#1565C0") +
      ggplot2::geom_ribbon(data = conf, aes(x = year, y = mean, ymin = lwr, ymax = uppr), alpha = 0.25, fill = "#1565C0") +
      ggplot2::facet_wrap(~cohort) + 
      ggplot2::scale_x_continuous(breaks = pretty(seq(current.yr, current.yr + yrs))) +
      ggplot2::scale_y_continuous(breaks = pretty(df$N), labels = scales::comma) +
      xlab("") + ylab("Abundance") +
      gg.opts
  }
  
  p1 <- make_plot(tot.df, tot.conf)
  p2 <- make_plot(narw.df, narw.conf)
  
  print(p1)
  print(p2)
  
  # Find 95% confidence intervals on final population size
  ci <- tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.025, 0.975))]
  
  return(ci)
  
}


reshape_array <- function(a, value, ..., .id = NULL) 
{
  qs <- rlang::quos(...)
  if (missing(value)) {
    evalue <- rlang::sym("var")
  }
  else {
    evalue <- rlang::enquo(value)
  }
  len <- length(qs)
  d <- dim(a)
  if (len > 0) {
    dimnames <- purrr::map(qs, rlang::quo_name) |> purrr::as_vector()
  }
  else {
    dimnames <- paste0("dim_", 1:length(d))
  }
  l <- list()
  for (i in 1:length(d)) {
    l[[i]] <- 1:d[i]
  }
  names(l) <- dimnames
  tidy <- expand.grid(l)
  tidy[[rlang::quo_name(evalue)]] <- a[as.matrix(tidy)]
  if (!is.null(.id)) 
    tidy[[.id]] <- rlang::expr_name(a)
  return(tidy)
}

scale_growth <- function(target.max = 45000, maxBC = 0.6, age = 69){
  # When lean = 1, length2mass returns the total mass of an animal given its length(age), as per Fortune et al. (2021)
  # By setting lean < 1, we can scale down this relationship to obtain a growth curve for lean mass.
  # We use an optimisation to find esstimate a scalar that will prevent unrealistic values of total mass (>45t) 
  # given fluctuating values of reserve mass and a maximum body condition of 60%.
  opt.fun <- function(x, age, target.max, maxBC){
    return((length2mass(age2length(age, agL(age)), mL(), x)/(1-maxBC)-target.max)^2)}
  out <- optimize(f = opt.fun, interval = c(0,1), target.max = target.max, maxBC = maxBC, age = age)
  return(out$minimum)
}

gestation_threshold <- function(){
  
  set.seed(20230412)
  
  # Simulate pregnant females - set stressors and growth to 0 to ensure that animals do not die
  # Note: Fetus growth is not tied to growth argument
  mg <- narw(nsim = 1000, cohortID = 4, progress = 1, stressors = 0, growth = 0)
  mdat <- mg$sim$`ad(f,p)`
  
  # Energy density of lipids -- (kJ / kg - converted to GJ/kg)
  ED_lipids = 39300/1000000
  
  # Percent lipid breakdown during catabolism
  lipid_cat = runif(1000, 0.53, 0.6)
  
  # Compile data
  m_kj <- mdat[, list(mass = unique(mass), E_gest = sum(E_gest)/1000), whale]
  m_kj[, kg:= E_gest/(ED_lipids*lipid_cat) + 0.05*mass]
  m_kj[, min_bc:=kg/mass]
  
  gam_gest <- mgcv::gam(min_bc ~ s(mass), data = m_kj, method = "REML", family = Gamma(link = "log"))
  
  # Residual checks
  gratia::appraise(gam_gest)
  
  # Extract fitted relationship
  b0 <- coef(gam_gest)[1] # Intercept
  gest.df <- gratia::smooth_estimates(gam_gest) |>
    gratia::add_confint() |>
    dplyr::mutate(across(tidyselect::all_of(c("est", "lower_ci", "upper_ci")), .fns = \(x) x + b0)) |>
    gratia::transform_fun(fun = \(x) exp(x))
  
  # Create plot
  p <- ggplot2::ggplot(data = gest.df) + 
    geom_ribbon(aes(x = mass, ymin = lower_ci, ymax = upper_ci), alpha = 0.25) +
    geom_line(aes(mass, est)) + 
    xlab("Total mass (kg)") + ylab("Minimum body condition (%)") +
    theme_narw() +
    scale_x_continuous(breaks = pretty(range(m_kj$mass))) +
    scale_y_continuous(breaks = pretty(range(m_kj$min_bc), n = 10)) +
    geom_rug(data = m_kj, aes(x = mass), sides = "b")

  print(p)
  
  return(gam_gest)
  
}

# NOISE ------------------------------------------------------

response_km <- function(n = 100, max.d = 200, SL = 200, logfac = 15, a = 1.175, main = ""){
  
  # max.d is the maximum distance from the source to be considered
  # n gives the number of curves to extract and plot
  
  set.seed(23579)
  
  x <- seq(0, max.d, length.out = 1000) # Distances from source
  db <- km2dB(x, SL = SL, logfac = logfac, a = a) # Received levels
  db[db<85] <- 85 # So that splinefun returns 0 below lower bound
  
  # Dose-response function
  dose.range <- seq(85, 200, length.out = 1000)
  f <- splinefun(dose.range, doseresponse[,1])
  
  # Apply dose-response function to received levels
  y <- f(db)
  y[y > 1] <- 1 # In case SL is > 200
  
  plot(x, y, type = "l", col = "grey", ylab = "p(response)", xlab = "Distance from source (km)", main = main)
  for(k in sample(x = 2:ncol(doseresponse), size = n, replace = FALSE)){
    f <- splinefun(dose.range, doseresponse[,k])
    y <- f(db)
    y[y>1] <- 1
    lines(x, y, col = "grey")
  }
}


TL <- function(r, logfac = 15, a = 1.175){
  # r in km
  # alpha is in dB/km
  loss <- logfac * log10(r*1000)
  loss[loss < 0] <- 0
  loss <- loss + a * r
  return(loss)
}

dB2km <- function(target.dB, SL = 200, logfac = 15, a = 1.175){
  opt.fun <- function(r, SL, target.L, logfac){
    return((SL-TL(r, logfac = logfac, a = a)-target.L)^2)}
  out <- optimize(f = opt.fun, interval = c(0,30000), SL = SL, target.L = target.dB, logfac = logfac)
  return(out$minimum)
}

km2dB <- function(r, SL = 200, logfac = 15, a = 1.175){
  return(SL-TL(r, logfac = logfac, a = a))
}

get_doseresponse <- function(input.ee){
  data.table::fread("data/elicitation/EE_results.csv", select = 2:3) |> as.matrix()
  # data.table::fread("data/elicitation/EE_results.csv", select = 4:5003) |> as.matrix()
}

proxy_noise <- function(ambient = 60, source.lvl = 220, mitigation = 10, x, y){
  
  # Spatial support
  support.poly <- targets::tar_read(support_poly)
  
  # Load turbine locations for each wind farm
  turbine.locs <- targets::tar_read(turbines)
  
  # Calculate central locations
  farm.locs <- purrr::map(.x = turbine.locs, .f = ~{
    x.avg <- mean(.x$x); y.avg <- mean(.x$y)
    .x[, c("longitude", "latitude", "x", "y")] <- NULL
    .x$easting <- x.avg; .x$northing <- y.avg
    .x <- add_latlon(.x)
    dplyr::distinct(.x)
  }) |> do.call(what = rbind)
  
  # Calculate distances to source
  r <- ambient.r <- targets::tar_read(density_support) |> raster::raster()
  raster::values(r) <- NA
  ambient.r <- ambient.r - 1
  ambient.r[ambient.r < 0] <- NA
  ambient.r[ambient.r == 0] <- ambient
  
  rdist <- purrr::map(.x = seq_along(turbine.locs), .f = ~{
    tmp <- r
    tmp[raster::cellFromXY(tmp, farm.locs[.x, c("easting", "northing")])] <- 1
    terra::distance(tmp) |> terra::mask(mask = support.poly)}) |> 
    purrr::set_names(nm = names(turbine.locs))

  # Compute sound levels
  db <- purrr::map(.x = rdist, .f = ~{
    db.out <- source.lvl - mitigation - TL(.x)
    db.out[db.out < ambient] <- NA
    db.out})
  
  db.all <- raster::stack(db) |> raster::calc(fun = max)
  db.all <- raster::stack(ambient.r, db.all) |> raster::calc(fun = sum, na.rm = TRUE)
  db.all[db.all == 0] <- NA
  
  # plot(terra::rast(db.all), col = pals::viridis(100), main = "Noise levels (dB)", xlab = "Easting (km)")
  
  db.all <- rescale_raster(db.all, new.min = ambient, new.max = source.lvl)
  db.list <- purrr::map(.x = month.abb, .f = ~db.all)
  
  out <- purrr::map(.x = db.list, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}

# FISHING GEAR ------------------------------------------------------

entgl_surface <- function(scalar = 1){
  
  d <- targets::tar_read(density_narw)
  regions <- targets::tar_read(regions)
  rd <- raster::raster(d[[1]])
  GSL_scalar <- 1.114327 # Gulf of St Lawrence
  ELSWC_scalar <- 0.03025689 # Elsewhere Canada
  
  purrr::map(.x = 1:12, .f = ~{
    
    print(month.abb[.x])
    
    # United States
    lines.r <- raster::raster(x = paste0("data/gear/LineNum_PostAction/DST_LineNum_PostAction_m", 
                                         stringr::str_pad(.x, 2, pad = "0"), ".tif")) |> 
      raster::projectRaster(to = rd)
    
    
    regions.raster <- raster::rasterize(regions[, "region"], rd)
    US.waters <- sum(lines.r[], na.rm = TRUE)
    
    # US.waters <- 
    # sum(lines.r[which(regions.raster[] %in% which(regions$region %in% c("CCB", "SNE", "MIDA", "SEUS", "GOM")))], na.rm = TRUE)
    
    # Apply scalars
    GSL.waters <- GSL_scalar * US.waters
    ELSWC.waters <- ELSWC_scalar * US.waters
    
    GSL.cells <- which(regions.raster[] == which(regions$region == "GSL"))
    ELSWC.cells <- which(regions.raster[] == which(regions$region %in% c("CABOT, SCOS, BOF_lower, BOF_upper")))
    
    lines.all <- regions.raster
    lines.all[!is.na(lines.all)] <- 0
    lines.all <- raster::merge(lines.r, lines.all)
    lines.all[GSL.cells] <- scalar * (GSL.waters / length(GSL.cells))
    lines.all[ELSWC.cells] <- scalar * (ELSWC.waters / length(ELSWC.cells))
    
    as(lines.all / raster::maxValue(lines.all), "SpatialGridDataFrame")
    
  }) |> purrr::set_names(nm = month.abb)
}

get_entglD <- function(){
  readr::read_csv("data/gear/duration_entanglements_BOEM.csv") |> 
    janitor::clean_names()
}

entgl_durations <- function(){
  
  # Duration of entanglement events
  entgl <- targets::tar_read(entgl_d)
  severity <- unique(entgl$severity)
  
  # Negative binomial fit
  fnb <- purrr::set_names(severity) |> 
    purrr::map(.f = ~fitdistrplus::fitdist(entgl[entgl$severity == .x, ]$duration_days, "nbinom"))
  
  par(mfrow = c(3, 2))
  purrr::walk(
    .x = severity,
    .f = ~ {
      hist(entgl[entgl$severity == .x, ]$duration_days,
           freq = FALSE,
           main = .x,
           xlab = "Entanglement duration (days)"
      )
      x <- seq(0, max(entgl[entgl$severity == .x, ]$duration_days), by = 1)
      lines(x, dnbinom(x, size = fnb[[.x]]$estimate["size"], mu = fnb[[.x]]$estimate["mu"]), lwd = 1.5)
      fitdistrplus::cdfcomp(fnb[[.x]])
    }
  )
  print(fnb)
  
  cat("Truncation at 365 days:\n\n")
  purrr::map(.x = fnb, .f = ~1-pnbinom(365, size = .x$estimate["size"], mu = .x$estimate["mu"]))
}

proxy_fishing <- function(){
  
  # Same as dummy prey layer, but rescale to 0–1 range
  
  support_poly <- targets::tar_read(support_poly)
  density_narw <- targets::tar_read(density_narw)
  
  # Create a correlated random field as a dummy prey surface
  grid <- sp::makegrid(support_poly, cellsize = 15) |>
    dplyr::rename(x = x1, y = x2)
  area.ex <- raster::extent(support_poly)
  grid <- expand.grid(seq(area.ex[1], area.ex[2], length.out = 200),
                      seq(area.ex[3]-500, area.ex[4], length.out = 200))
  names(grid) <- c("x", "y")
  
  target.r <- raster::raster(density_narw[[1]])
  
  out <- purrr::map(.x = month.abb, .f = ~{
    g.dummy <- gstat::gstat(formula = z ~ 1, 
                            locations = ~ x + y, 
                            dummy = TRUE, 
                            beta = 1, 
                            model = gstat::vgm(psill = 1000, range = 150, model = 'Exp'), 
                            nmax = 20)
    
    yy <- predict(g.dummy, newdata = grid, nsim = 4)
    sp::gridded(yy) = ~x+y
    yy <- raster::raster(yy)
    
    # Clip to area of interest
    raster::resample(yy, target.r) |> 
      raster::crop(target.r) |> 
      raster::mask(target.r) |> 
      rescale_raster(new.max = 0.4)
  })
  
  out <- purrr::map(.x = out, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}


# VESSEL TRAFFIC ------------------------------------------------------

proxy_vessels <- function(pmax = 0.05){
  
  # Data from https://globalmaritimetraffic.org/ - monthly rasters between Jan and Dec 2022
  
  poly <- targets::tar_read(support_poly)
  poly.v <- as(poly, "SpatVector")
  
  A.files <- list.files(path = "data/vessels/", pattern = "A.tif")
  B.files <- list.files(path = "data/vessels/", pattern = "B.tif")
  
  # Prepare spatial rasters
  A.maps <- purrr::map(.x = A.files, .f = ~{
    raster::raster(file.path("data/vessels", .x)) |> 
    as("SpatRaster") |> 
    terra::project(y = as.character(narw_crs())) |>  
    terra::mask(mask = poly.v)
    }) |> purrr::set_names(nm = month.abb)
    
  B.maps <- purrr::map(.x = B.files, .f = ~{
    raster::raster(file.path("data/vessels", .x)) |> 
      as("SpatRaster") |> 
      terra::project(y = as.character(narw_crs())) |>  
      terra::mask(mask = poly.v)
  }) |> purrr::set_names(nm = month.abb)  
  
  # Resample to desired resolution
  target.r <- targets::tar_read(density_support) |> raster::raster() |> as("SpatRaster")
  A.maps <- purrr::map(.x = A.maps, .f = ~terra::resample(.x, target.r))
  B.maps <- purrr::map(.x = B.maps, .f = ~terra::resample(.x, target.r))
  
  # Merge rasters
  vessels.r <- purrr::map2(.x = A.maps, .y = B.maps, 
                           .f = ~{
                             out.r <- raster::merge(raster::raster(.x), raster::raster(.y))
                             out.r <- sqrt(out.r)
                             out.r <- rescale_raster(out.r, new.max = pmax)
                             out.r[is.na(out.r)] <- 0
                             raster::mask(out.r, mask = poly)})
  
  out <- purrr::map(.x = vessels.r, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}

# PREY ------------------------------------------------------
  
proxy_prey <- function(){

  
  # Sorochan et al. (2022) --
  # In Grand Manan Basin, where Calanus spp. is almost exclusively composed of C. finmarchicus 
  # (e.g. Michaud and Taggart, 2007), Baumgartner and Mate (2003) reported NARW foraging dive
  #  depths ranging from ∼80 to 175 m and maximum concentrations of stage CV of C. finmarchicus 
  #  near foraging NARW ranging from 2000 to 20000 ind m−3.

  #  Gavrilchuk et al. (2021) --
  #  Estimated an abundance threshold of 4476–6700 ind m−3 stage
  #  CV C. finmarchicus required to meet energetic demands of a resting female NARW and higher
  #  values for pregnant and lactating females. 
  
  # Grieve et al. (2017)
  # Figure 5 shows predicted/sampled values on scale of up to 12e5 ind / 100 m3 => 12,000 ind/m3
  
  support_poly <- targets::tar_read(support_poly)
  density_narw <- targets::tar_read(density_narw)
  regions.support <- targets::tar_read(regions)
  
  target.r <- raster::raster(density_narw[[1]])

  # regions.raster <- raster::rasterize(
  #   x = regions.support,
  #   y = raster::raster(x = regions.support,
  #                      resolution = raster::res(target.r),
  #                      origin = raster::origin(target.r),
  #                      ext = raster::extent(target.r),
  #                      crs = narw_crs())) |> 
  #   raster::resample(target.r, method = "ngb") |> 
  #   raster::mask(target.r)
  # 
  # # To zero out the prey layer in the SEUS and MIDA
  # cells.OUT <- which(regions$region %in% c("SEUS", "MIDA"))
  # for(i in seq_along(cells.OUT)) regions.raster[regions.raster$layer == cells.OUT[i]] <- -1
  # regions.raster[regions.raster >= 0] <- 1
  # regions.raster[regions.raster <0] <- 0
  
  gradient.r <- raster::as.data.frame(target.r, xy = TRUE)
  gradient.r$val <- 1/(1+exp(-0.02*(gradient.r$x-500)))
  gradient.r$val[which(is.na(gradient.r$Nhat))] <- NA
  gradient.r <- raster::rasterFromXYZ(xyz = gradient.r[,c("x", "y", "val")], res = raster::res(target.r), crs = narw_crs())
  # plot(world, col = "grey")
  # plot(gradient.r, add = T)
  
  # Create a correlated random field as a dummy prey surface
  grid <- sp::makegrid(support_poly, cellsize = 15) |>
    dplyr::rename(x = x1, y = x2)
  area.ex <- raster::extent(support_poly)
  grid <- expand.grid(seq(area.ex[1], area.ex[2], length.out = 200),
                      seq(area.ex[3]-500, area.ex[4], length.out = 200))
  names(grid) <- c("x", "y")
  
  out <- purrr::map(.x = 1:12, .f = ~{

    seed <- as.numeric(paste0("250", .x))
    set.seed(seed)
    
    g.dummy <- gstat::gstat(formula = 
                            z ~ 1 + y, 
                            locations = ~ x + y, 
                            dummy = TRUE, 
                            beta = 12000, 
                            model = gstat::vgm(psill = 2000, range = 500, model = 'Exp'), 
                            nmax = 20)
    
    yy <- predict(g.dummy, newdata = grid, nsim = 4)
    sp::gridded(yy) = ~x+y
    yy <- raster::raster(yy)

    # Clip to area of interest
    rout <- raster::resample(yy, target.r) |> 
      raster::crop(target.r) |> 
      raster::mask(target.r)
    
    rout * gradient.r
    # r <- rout * gradient.r
    # plot_raster(r, zero = TRUE)
  })
  
  out <- purrr::map(.x = out, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}

# SPATIAL DATA ------------------------------------------------------


#' Projected coordinate system
#'
#' @return An object of class \code{CRS}.
#' 
narw_crs <- function(){
  sp::CRS("+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
}

# rescale_raster <- function(x, new.min = 0, new.max = 1, x.min = NULL, x.max = NULL) {
#   if(is.null(x.min)) x.min = raster::cellStats(x, "min")
#   if(is.null(x.max)) x.max = raster::cellStats(x, "max")
#   new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
# }

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

  rpts <- raster::shapefile(file.path(getwd(), "data/regions/r_pts.shp")) |> sp::spTransform(CRSobj = narw_crs())
  
  # ..........................................
  # Manually add support for Canada
  # ..........................................
  
  # Load region boundaries
  regions.support <- targets::tar_read(regions)
  
  # Filter out Canada
  regions.support <- regions.support[regions.support@data$region %in% c("GSL", "SCOS", "CABOT"),] 
  
  # Overlay regions onto current spatial support
  regions.raster <- raster::rasterize(
    x = regions.support,
    y = raster::raster(x = regions.support,
                       resolution = filtered_support@grid@cellsize,
                       origin = raster::origin(raster::raster(filtered_support))))
  
  regions.raster[!is.na(regions.raster)] <- 1
  # regions.raster <- raster::mask(regions.raster, world, inverse = TRUE)
  
  # world.raster <- raster::rasterize(
  #   x = world,
  #   y = raster::raster(filtered_support))
  # world.raster[!is.na(world.raster)] <- 1
  
  fr <- raster::raster(filtered_support)
  fcn <- raster::cellFromXY(fr, coordinates(fr))
  fr <- raster::setValues(fr, fcn)
  which.cells <- raster::extract(fr, rpts)
  
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
  
  cleaned_support[which.cells] <- 0
  
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
    dplyr::summarise(Nhat = round(sum(Nhat),0), .groups = "keep") |> 
    dplyr::arrange(region) |> 
    dplyr::filter(!is.na(region)) |> 
    dplyr::mutate(perc = round(100* Nhat / Ntot, 1)) |> 
    dplyr::filter(region == reg) |> 
    dplyr::ungroup()
  
  # df.out <- dplyr::left_join(df, Ntot, by = "month") |> 
  #   dplyr::mutate(perc = round(100* Nhat / Ntot, 1)) |> 
  #   dplyr::filter(region == reg) |> 
  #   dplyr::ungroup()
  

  cat(paste0("Average: ", round(mean(df.out[df.out$perc>0,]$perc),1), "%"))
  return(df.out)
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
  regions.m <- raster::rasterize(regions[, "region"], raster::raster(support)) |> as("SpatialGridDataFrame")
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

summary_geo <- function(obj, lims, res, geomap) {
  lapply(seq(dim(obj)[3]), function(x){
    obj.ind <- obj[,,x]
    obj.ind <- cbind(obj.ind, c(obj.ind[-1, 1], NA), c(obj.ind[-1, 2], NA))
    dd <- sapply(1:nrow(obj.ind), FUN = function(r) geoD(geomap, obj.ind[r,1], obj.ind[r,2], obj.ind[r,3], obj.ind[r,4], lims, res))
    unname(dd)
  })
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
    
    regions <- raster::shapefile(dir(path = "data/regions", pattern = "regions_eez.shp$", full.names = TRUE)) |> 
      sp::spTransform(CRSobj = narw_crs())
    
  } else {
  
  # Load shapefile of spatial regions (CCB, MIDA, SEUS, etc.)
  regions <- raster::shapefile(dir(path = "data/regions", pattern = "regions.shp$", full.names = TRUE))
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
# 
# body_condition <- function(n, cohorts, alpha = 6, beta = 20){
#   
#   sapply(X = cohorts, FUN = function(x) rbeta(n = n, shape1 = alpha, shape2 = beta), simplify = FALSE, USE.NAMES = TRUE)
#   
#   # curve(expr = dbeta(x = x, shape1 = 6, shape2 = 20),
#   #       xlab = "", ylab = "", main = "Beta PDF with alpha = 6 and beta = 20",
#   #       lwd = 2, col = 1, xlim = c(0, 1))
#   
# }


#' Decrease in milk feeding efficiency as a function of calf age
#'
#' @param a Calf age in days.
#' @param Tstop Time at which milk stops nursing and becomes nutritionally independent (in days).
#' @param Tdecrease Time at which milk comsuption starts to decrease (in days).
#' @param E Steepness of the decline

# milk_assimilation_R <- function(a = seq(0, 500), Tstop = 365, Tdecrease = 100, E = 0.9){
#   p <- (a-Tdecrease)/(Tstop-Tdecrease)
#   res <- (1 - p)/ (1 - (E*p))
#   res[res<0] <- 0
#   res[res>1 & a>Tdecrease] <- 0
#   res[res>1 & a<Tdecrease] <- 1
#   return(res)
# }

#' Increase in milk provisioning by female as a function of her body condition.
#'
#' @param E Steepness of the non-linear relationship.
#' @param blubber_mass Mass of blubber (kg).
#' @param maintenance_mass Mass of structural tissues (kg).
#' @param target_condition Target body condition, expressed as the ratio of reserve to maintenance mass.
#' @param starvation Starvation threshold

# milk_provisioning_R <- function(E = -2, blubber_mass, maintenance_mass, target_condition = 0.3, starvation = 0.15){
#   p1 <- 1-E
#   p2 <- blubber_mass - starvation * maintenance_mass
#   p3 <- maintenance_mass * (target_condition - starvation)
#   p4 <- E * (blubber_mass - (starvation*maintenance_mass))
#   res <- p1*p2/(p3-p4)
#   res[blubber_mass/maintenance_mass < starvation] <- 0
#   res[blubber_mass/maintenance_mass >= target_condition] <- 1
#   return(res)
# }


# Feeding effort as a function of body condition
#' Title
#'
#' @param eta Steepness of the non-linear relationship 
#' @param target_condition 
#' @param maintenance_mass 
#' @param blubber_mass 

# feeding_effort_R <- function(eta, target_condition = 0.3, maintenance_mass, blubber_mass){
#   1/(1+exp(-eta * ((target_condition*maintenance_mass/blubber_mass)-1)))
# }

#' Feeding effort as function of copepod density
#' 
#' @param gamma Feeding threshold
#' @param D Coepepod density

# feeding_threshold_R <- function(gamma, D){
#   1/(1+exp(gamma-D))
# }

# mammary_mass_R <- function(M){
#   10^(0.902*log(M, 10)-1.965)
# }

# milk_production_R <- function(mu){
#   0.0835*mu^1.965
# }

# LC <- function(M, sr, sm, phi, psi, delta){
#   
#   (sum(phi*delta) * sr * (1.46+0.0005*M) + sum(psi*delta) * sm * (5.17+0.0002 *M))*M
#   
# }

# PLOTTING ------------------------------------------------------

plot_gam <- function(obj){
  
  x <- seq(0, 1, 0.01)
  
  mbc_preds <- obj$gam$pred$bc_gest
  bc_preds <- obj$gam$pred$bc
  surv_preds <- obj$gam$pred$surv
  
  par(mfrow = c(3,3))
  for(i in obj$param$cohorts$id){
    plot(x, surv_preds[[as.character(i)]](x), type = 'l', xlab = "", ylab = "p(surv)", main = obj$param$cohorts[i, name])
  }
  
  par(mfrow = c(3,3))
  for(i in obj$param$cohorts$id){
    plot(x, bc_preds[[as.character(i)]](x), type = 'l', xlab = "", ylab = "BC", main = obj$param$cohorts[i, name])
  }
}

growth_curve <- function(param, obj, cohortID, ylabel){

  cohorts <- obj$param$cohorts
  dat <- data.table(obj$sim[[cohorts[id==cohortID, abb]]])
  ind.survived <- dat[day == 365 & alive == 1, whale,]
  dat <- dat[whale %in% ind.survived]
  n.ind <- length(ind.survived)
  # n.ind <- obj$param$nsim
  if(nrow(dat) > 0){

  param_dat <- matrix(unlist(dat[day > 0 & cohort == cohortID, param, with = FALSE]),
                      ncol = 365, 
                      nrow = n.ind, byrow = TRUE)
  param_dat[param_dat == 0] <- NA

  param_median <- apply(param_dat, MARGIN = 2, FUN = function(x) median(x, na.rm = TRUE))
  param_lower <- apply(param_dat, MARGIN = 2, FUN = function(x) quantile(x, 0.025, na.rm = TRUE))
  param_upper <- apply(param_dat, MARGIN = 2, FUN = function(x) quantile(x, 0.975, na.rm = TRUE))
  
  param_df <- data.table::data.table(day = 1:365, median = param_median, lwr = param_lower, uppr = param_upper)
  
  y.range <- c(min(param_df$lwr, na.rm = TRUE), max(param_df$uppr, na.rm = TRUE))
  
  p <- ggplot2::ggplot(data = param_df) +
    ggplot2::geom_ribbon(aes(x = day, ymin = lwr, ymax = uppr), alpha = 0.25, fill = "deepskyblue4") +
    ggplot2::geom_path(aes(x = day, y = median), col = "deepskyblue4") + 
    theme_narw() +
    ylab(ylabel) + xlab("Day of the year") +
    ggtitle(label = ifelse(grepl(pattern = "calf", param), cohorts[1,name], cohorts[id==cohortID,name])) +
    ggplot2::scale_x_continuous(breaks = seq(5, 365, by = 40)) +
    ggplot2::scale_y_continuous(
      limits = ~ y.range,
      breaks = ~ pretty(y.range, 10))
  
  } else {
    p <- list(NULL)
  }
  return(p)
}


theme_narw <- function(vertical = FALSE){
  
  # font <- "Georgia"   # Assign font family up front
  
  theme_grey() %+replace%    # Replace elements we want to change
    
    ggplot2::theme(
      
      # Axis elements
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.text.y = ggplot2::element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), 
                                          angle = ifelse(vertical, 90, 0),
                                          vjust = ifelse(vertical, 0.5, 0),
                                          hjust = ifelse(vertical, 0.5, 0)),
      axis.text.x = ggplot2::element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
      axis.title = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), angle = 90),
      axis.title.x = ggplot2::element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),

      # Panel elements
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = ifelse(vertical, 0.5, 0)),

      # Facet elements
      strip.background = ggplot2::element_rect(fill = "grey20"),
      strip.text = ggplot2::element_text(colour = 'white', size = 12),
      strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "cm")),
      strip.text.y = element_text(margin = margin(0,0.2,0,0.2, "cm"), angle = -90)
      
    )
}

plot_raster <- function(r, duke = FALSE, zero = FALSE, breaks = NULL){
  
  if(class(r) == "SpatialGridDataFrame") r <- raster::raster(r)
  
  dat <- raster::as.data.frame(r, xy = TRUE)
  names(dat)[3] <- "Nhat"
  dat <- dat[complete.cases(dat),]
  
  if(is.null(breaks)){
  if(zero){
    dat.c <- dat |> dplyr::filter(Nhat > 0)
    colour.breaks <- unique(c(0, quantile(dat.c$Nhat, seq(0,1, 0.1))))
  }else {
    if(duke){
      colour.breaks <- colour_breaks(dat)
    } else {
      colour.breaks <- c(0, rgeoda::natural_breaks(10, dat[, "Nhat", drop = FALSE]), max(dat$Nhat))
      # colour.breaks <- c(seq(0, 0.01, length.out = 5), seq(0.02, 0.1, length.out = 5), 0.2, 0.5, 0.8, 1)
    }
  }
  } else {
    colour.breaks <- breaks
  }
  
  dat <- dat |> 
    dplyr::mutate(Ncol = cut(Nhat, breaks = colour.breaks, include.lowest = TRUE)) |> 
    dplyr::mutate(Ncol = factor(Ncol))
  
  levels(dat$Ncol) <-  gsub(pattern = ",", replacement = " – ", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\[|\\]", replacement = "", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\(|\\)", replacement = "", levels(dat$Ncol))
  
  r.dat <- raster::rasterFromXYZ(dat[, c("x", "y", "Nhat")], res = raster::res(r), crs = narw_crs())
  x.lim <- raster::extent(r.dat)[1:2]
  y.lim <- raster::extent(r.dat)[3:4]
  
  world.sf <- sf::st_as_sf(world)
  
  ggplot2::ggplot(data = dat) +  
    ggplot2::geom_raster(aes(x,y, fill = Ncol)) +
    ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    ylab("") + xlab("") +
    theme_narw() + 
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::coord_sf(xlim = x.lim, ylim = y.lim, expand = FALSE) +
    ggplot2::scale_fill_manual(values = pals::viridis(length(levels(dat$Ncol))),
                      guide = guide_legend(reverse = TRUE))
  
}

# Jason's colour scale
colour_breaks <- function(dat){
  colour.breaks <- 25*c(0,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3,10)/100
  colour.breaks <- c(colour.breaks, seq(max(colour.breaks), ceiling(max(dat$Nhat, na.rm = TRUE)), length.out = 5))
  colour.breaks <- round(colour.breaks[!duplicated(colour.breaks)],3)
  return(colour.breaks)
}

deg2radians_R <- function(angle){
  angle*pi/180;
}

gape_size_R <- function(L, omega, alpha){
  (deg2radians_R(alpha) * ((omega^2)/4 + (0.2077*L - 1.095)^2))/2
}


# UTILITIES ------------------------------------------------------

inspect <- function(obj, cohort = "ad(f,l)", whaleID = 1, calf = FALSE){
  pdf(file = "allplots.pdf")
  par(mfrow = c(3,3))
  if(calf){
    nn <- names(m$sim[[1]])[7:112][grepl(pattern = "calf", names(m$sim[[1]])[7:112])]
  } else {
    nn <- names(m$sim[[1]])[7:112]
  }
  for(i in nn){
    dat <- obj[["sim"]][[cohort]][whale == whaleID, c("day", i), with = FALSE]
    plot(dat, xlab = "Day", ylab = "", main = i, type = "l")
  }
  par(mfrow = c(1,1))
  dev.off()
}

format_dt <- function(dt, direction = "col"){
  dt |>  janitor::adorn_percentages(denominator = direction) |>
    janitor::adorn_pct_formatting() |> 
    janitor::adorn_ns()
}

array2dt <- function(a){
  y <- aperm(a, c(1, 3, 2))
  dim(y) <- c(prod(dim(a)[-2]), dim(a)[2])
  y <- data.table::data.table(y)
  names(y) <- colnames(a)
  return(y)
}

add_whale <- function(y, n.ind){
  y$row_id <- 1:nrow(y)
  y$day <- rep(1:365, times = n.ind)
  y$whale <- rep(1:n.ind, each = 365)
  return(y)
}

extract <- function(obj){
  if(!"narwsim" %in% class(obj)) stop("Input object must be of class <narwsim>")
  coh <- obj$param$cohort.ab
  out <- lapply(X = coh, FUN = function(x){
      abind::abind(obj[["locs"]][[x]], obj[["sim"]][[x]], along = 2)
  }) |> purrr::set_names(nm = coh)
  return(out)
}

get_daylight <- function(){
  
  out <- list()
  
  # For each month
  for(m in 1:12){
    
    # Retrieve density surface
    d <- targets::tar_read(density_narw)[[m]]
    
    # Convert to data.frame and add lat/lon coordinates
    r <- d |> 
      data.frame() |> 
      dplyr::select(-Nhat) |> 
      dplyr::rename(easting = s1, northing = s2) |> 
      add_latlon()
    
    target.date <- lubridate::as_date(paste0("2023-", m, "-15"))
    cat("\n")
    print(target.date)
    cat("\n")
    
    # Calculate duration of daylight hours on the 15th day of each month
    pb <- progress::progress_bar$new(total = nrow(r))
    for(u in 1:nrow(r)){
      pb$tick()
      daylight <- suncalc::getSunlightTimes(date = target.date, lat = r[u,"lat"], lon = r[u,"long"], tz = "EST")
      r$daylight[u] <- daylight$sunset - daylight$sunrise
    }
    
    # Return output raster
    rdaylight <- raster::rasterFromXYZ(xyz = r[, c("easting", "northing", "daylight")],
                                       res = raster::res(raster::raster(d)), crs = narw_crs())
    rdaylight <- raster::resample(x = rdaylight, y = raster::raster(d))
    rdaylight <- as(rdaylight, "SpatialGridDataFrame")
    
    out[[month.abb[m]]] <- rdaylight
  }
  
  return(out)
}

estBetaParams <- function(mu, std) {
  var <- std * std
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estGammaParams <- function(mu, std) {
  v <- std * std
  shape <- (mu^2)/v
  scale <- v/mu
  return(params = list(shape = shape, scale = scale))
}

estRayleighParams <- function(target.mean, upper, e = 0.001){
  
  opt.f <- function(param, p, e){
    set.seed(20230126)
    n <- extraDistr::rrayleigh(n = 100000, sigma = param[1])
    m <- (mean(n)-p[1])^2+ (max(n)-p[2])^2
    return(m)
  }
  
  out <- suppressWarnings(optim(par = 1, fn = opt.f, p = c(target.mean, upper), e = e))
  cat("Lambda:", out$par, "\n")
  cat("Sigma (movement):", out$par/sqrt(2), "\n")
}

estTruncNorm <- function(target.mean, L, U){
  
  opt.f <- function(param, bounds, target.mean){
    set.seed(20230126)
    n <- sapply(1:10000, FUN = function(x) rtnorm(bounds[1], param[1], bounds[1], bounds[2]))
    return((mean(n)-target.mean)^2)
  }
  
  out <- optim(par = 2, fn = opt.f, bounds = c(L,U), target.mean = target.mean)
  return(out$par)
}

draw <- function(distrib = "norm", mu, std = NULL, L = -Inf, U = Inf, seed = 215513, x.lab = "",
                 y.lab = "", title = NULL, verbose = TRUE, linecol = "black"){
  
  set.seed(seed)
  
  if(!distrib %in% c("tnorm", "norm", "gamma", "beta", "halfnorm")) stop("Unrecognized distribution")
  
  if(distrib == "halfnorm"){
    
    truncN.param <- estTruncNorm(target.mean = mu, L = L, U = U)
    if(verbose) print(truncN.param)
    n <- truncnorm::rtruncnorm(10000, L, U, L, truncN.param)
    
    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = truncnorm::dtruncnorm(x = x, L, U, L, truncN.param), from = min(n), to = max(n), 
          main = ifelse(is.null(title), "Half Normal distribution", title), ylab = y.lab, xlab = x.lab, lwd = 1.5, col = linecol)
    
  } else if(distrib == "tnorm"){
    
    n <- truncnorm::rtruncnorm(10000, L, U, mu, std)

    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = truncnorm::dtruncnorm(x = x, a = L, b = U, mean = mu, sd = std), from = min(n), to = max(n),
          main = ifelse(is.null(title), "Truncated Normal distribution", title), ylab = y.lab, xlab = x.lab, lwd = 1.5, col = linecol)
    
  } else if(distrib == "norm"){
    
    n <- rnorm(10000, mu, std)
    
    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = dnorm(x = x, mu, std), from = min(n), to = max(n),
          main = ifelse(is.null(title), " Normal distribution", title), ylab = "", xlab = "", lwd = 1.5, col = linecol)
    
  } else if(distrib == "gamma") {
    
    gamma.params <- estGammaParams(mu, std)
    if(verbose) print(gamma.params)
    n <- rgamma(10000, shape = gamma.params$shape, scale = gamma.params$scale)
    
    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = dgamma(x = x, shape = gamma.params$shape, scale = gamma.params$scale), from = min(n), to = max(n), 
          main = ifelse(is.null(title), "Gamma distribution", title), ylab = y.lab, xlab = x.lab, lwd = 1.5, col = linecol)
    
  } else if(distrib == "beta") {
    
    beta.params <- estBetaParams(mu, std)
    if(verbose) print(beta.params)
    n <- rbeta(10000, beta.params$alpha, beta.params$beta)

    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = dbeta(x = x, beta.params$alpha, beta.params$beta), from = 0, to = 1, 
          main = ifelse(is.null(title), "Beta distribution", title), ylab = y.lab, xlab = x.lab, lwd = 1.5, col = linecol)
  }
}

rescale_raster <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min <- raster::cellStats(x, "min")
  if(is.null(x.max)) x.max <- raster::cellStats(x, "max")
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

save_object <- function(obj, tg = TRUE, internal = FALSE) {
  if (internal) {
    tmp_env <- new.env(hash = FALSE)
    load("R/sysdata.rda", envir = tmp_env)
    if(tg) dat <- suppressWarnings(targets::tar_read_raw(obj)) else dat <- get(obj, envir = .GlobalEnv)
    tmp_env[[obj]] <- dat
    save(list = names(tmp_env), file = "R/sysdata.rda", envir = tmp_env)
    # usethis::use_data(daylight, density_support, doseresponse, entgl_d, params, regions, regions_m, support_poly, wL, world, internal = TRUE, overwrite = TRUE)
  } else {
    suppressWarnings(targets::tar_load(obj))
    for (i in obj) {
      file.name <- paste0("data/", i, ".rda")
      if (file.exists(file.name)) file.remove(file.name)
      do.call(save, c(lapply(i, as.name), file = paste0("data/", i, ".rda")))
    }
  }
}

format_table <- function(df, top = TRUE, bottom = TRUE, sign = "-"){
  # dashes <- purrr::map_dbl(.x = names(df), .f = ~nchar(.x))
  # dashes <- purrr::map(.x = dashes, .f = ~paste0(rep(sign, .x), collapse = ""))
  rbind(df[1:nrow(df) - 1,], c("---", rep("",ncol(df)-1)), df[nrow(df),])
}
