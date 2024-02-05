augment <- function (x, ...) UseMethod("augment", x)
write <- function(x, ...) UseMethod("write")

# MODEL ------------------------------------------------------

get_scenarios <- function(scenario = 1){
  
  if(scenario > 0){
  
  turbines <- targets::tar_read(turbines)
  vessel_routes <- targets::tar_read(vessel_routes)
  vessel_transits <- targets::tar_read(vessel_transits)
  
  # vparams <- vessel_transits[, list(speed_knt = unique(speed_knt), Nvessels = sum(Nvessels)), list(category, windfarm, routeID, phase)]
  
  vparams <- vessel_transits[phase == ifelse(scenario < 3, "Construction", "Operation and maintenance")]
  
  if(scenario == 1){
    start.month <- c(1,2,1)
    start.day = c(15,1,1)
  } else if(scenario == 2){
    start.month <- c(7,9,5)
    start.day = c(15,15,1)
  } else {
    start.month <- c(1,1,1)
    start.day = c(1,1,1)
  }
  
  phase <- ifelse(scenario <=2, 1, 2)
  
  out <-
    list(phase = phase,
      locs = data.table::as.data.table(turbines[[paste0("scenario_0", scenario)]]),
      routes = vessel_routes,
      vessels = vparams,
      start.month = start.month,
      start.day = start.day,
      piles.per.day = 1,
      ambient = 80,
      sourceLvL = 220,
      lowerdB = c(0,10,0)[scenario],
      logrange = 15,
      absorb = 1.175
    )
  
  } else {
    
    out <-
      list(phase = 0,
           locs = NULL,
           routes = NULL,
           vessels = NULL,
           start.month = NULL,
           start.day = NULL,
           piles.per.day = NULL,
           ambient = NULL,
           sourceLvL = NULL,
           lowerdB = NULL,
           logrange = NULL,
           absorb = NULL
      )
    
  }
  
  class(out) <- c("narwscenario", class(out))
  return(out)
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
  
  # if(verbose){
  #   cat("\n*** SINGLE summary stat ***\n")
  #   cat("mean –– N > 1 –– SD estimated based on methods of Wan et al. (2014)\n")
  #   cat("mean –– N == 1 –– SD equals weighted mean of other SD values, weighted by sample size\n")
  #   cat("min –– mean and SD derived from random draws of a Uniform bounded by min and max across records\n")
  #   cat("max –– mean and SD derived from random draws of a Uniform bounded by min and max across records\n")
  #   cat("\n")
  #   cat("*** MULTIPLE summary stats ***\n")
  #   cat("min & max –– mean and SD derived from random draws of a Uniform bounded by min and max across records\n")
  #   cat("min & mean –– N > 1 –– SD estimated based on methods of Wan et al. (2014)\n")
  #   cat("min & mean –– N == 1 –– SD equals weighted mean of other SD values, weighted by sample size\n")
  #   cat("max & mean –– N > 1 –– SD estimated based on methods of Wan et al. (2014)\n")
  #   cat("max & mean –– N == 1 –– SD equals weighted mean of other SD values, weighted by sample size\n")
  #   cat("\n")
  # }
  
  if(verbose){
    cat("Original data --------\n\n")
    print(dat.f, n = 100)
    cat("\n")
  } 
  
  impute_meta <- function(N = NULL,
                          minval = NULL,
                          maxval = NULL,
                          q1 = NULL,
                          q3 = NULL){

      if(all(is.null(c(q1,q3)))){

        sd.out <- (maxval-minval)/(2*(qnorm((N-0.375)/(N+0.25),0,1)))

        } else{

        sd.out <- ((maxval-minval)/(4*(qnorm((N-0.375)/(N+0.25),0,1)))) +
     ((q3-q1)/(4*(qnorm((0.75*N-0.125)/(N+0.25),0,1))))

      }

    sd.out[is.infinite(sd.out)] <- NA
    return(sd.out)
  }
  
  dat.out <- purrr::map(.x = dat.list, .f = ~{
    
    if(nrow(.x) > 0){
      
      varmin <- ifelse(all(is.na(.x$min)), NA, min(.x$min, na.rm = TRUE))
      varmax <- ifelse(all(is.na(.x$max)), NA, max(.x$max, na.rm = TRUE))
      
      if(is.na(varmin) & !is.null(theoretic.bounds))
        varmin <- theoretic.bounds[1] else varmin <- min(dat.f$min, na.rm = TRUE)
      
      if(is.na(varmax) & !is.null(theoretic.bounds)) 
        varmax <- theoretic.bounds[2] else varmax <- max(dat.f$max, na.rm = TRUE)
      
      # Impute min and max
      
      .x$min <- tidyr::replace_na(.x$min, replace = varmin)
      .x$max <- tidyr::replace_na(.x$max, replace = varmax)
      
      # Impute means
      
      .x$mean_median <- sapply(1:nrow(.x), FUN = function(x){
        if(is.na(.x$mean_median[x])) mean(c(.x$min[x], .x$max[x])) else .x$mean_median[x]
      })
      
      # Impute sample sizes
      
      .x$sample_size <- tidyr::replace_na(.x$sample_size, replace = 1)
      
      # Impute SD (if N>1, based on methods of Wan et al. 2014)
      
      .x$sd_se[is.na(.x$sd_se)] <- impute_meta(N = .x$sample_size, minval = .x$min, maxval = .x$max)[is.na(.x$sd_se)]
      .x$sd_se[is.na(.x$sd_se)] <- weighted.mean(x = .x$sd_se[!is.na(.x$sd_se)], w = .x$sample_size[!is.na(.x$sd_se)]) 
      
    
    # for(j in 1:nrow(.x)){
    # 
    #   myvars <- .x[j, c("min", "max", "mean_median", "sd_se")]
    #   myvars <- which(!is.na(myvars))
    # 
    #   varmin <- ifelse(all(is.na(.x$min)), NA, min(.x$min, na.rm = TRUE))
    #   varmax <- ifelse(all(is.na(.x$max)), NA, max(.x$max, na.rm = TRUE))
    # 
    #   if(is.na(varmin) & !is.null(theoretic.bounds))
    #     varmin <- theoretic.bounds[1] else varmin <- min(dat.f$min, na.rm = TRUE)
    # 
    #   if(is.na(varmax) & !is.null(theoretic.bounds))
    #     varmax <- theoretic.bounds[2] else varmax <- max(dat.f$max, na.rm = TRUE)
    # 
    #   # If sample size > 1, estimate SD based on methods of Wan et al. (2014)
    #   # https://doi.org/10.1186/1471-2288-14-135
    #   # else: take weighted mean of other SD, weighted by sample size
    # 
    #   if(length(myvars) == 1){
    # 
    #     if(myvars == 3){ # Only mean
    # 
    #       if(.x[j,]$sample_size > 1){
    # 
    #       .x[j,]$sd_se <- (varmax-varmin)/(2*(qnorm((max(.x[j,]$sample_size, 1, na.rm = TRUE)-0.375)/(max(.x[j,]$sample_size, 1, na.rm = TRUE)+0.25),0,1)))
    # 
    #       } else {
    # 
    #       w <- .x[, c("sd_se", "sample_size")]
    #       wm <- weighted.mean(w$sd_se, w = w$sample_size, na.rm = TRUE)
    #       if(nrow(.x) > 1) .x[j,]$sd_se <- wm
    #       }
    # 
    #     } else if(myvars == 1){ # Only min
    # 
    #       # Only min -- take mean/SD of random draws from Uniform bounded by min value and max across records
    # 
    #       n <- runif(10000, .x[j,]$min, varmax)
    #       .x[j,]$mean_median <- mean(n)
    #       .x[j,]$sd_se <- sd(n)
    # 
    #     } else if(myvars == 2){ # Only max
    # 
    #       # Only max -- take mean/SD of random draws from Uniform bounded by min across records and max value
    # 
    #       n <- runif(10000, varmin, .x[j,]$max)
    #       .x[j,]$mean_median <- mean(n)
    #       .x[j,]$sd_se <- sd(n)
    # 
    #     }
    # 
    #   } else {
    # 
    #     # if(sum(myvars) == 1){
    #     #
    #     #   # Only min
    #     #
    #     #   n <- runif(10000, .x[j,]$min, varmax)
    #     #   .x[j,]$mean_median <- mean(n)
    #     #   .x[j,]$sd_se <- sd(n)
    #     #
    #     # } else if(sum(myvars) == 2){
    #     #
    #     #   # Only max
    #     #
    #     #   n <- runif(10000, varmin, .x[j,]$max)
    #     #   .x[j,]$mean_median <- mean(n)
    #     #   .x[j,]$sd_se <- sd(n)
    # 
    #     if(sum(myvars) == 3){
    # 
    #       # Only min and max
    # 
    #       n <- runif(10000, varmin, varmax)
    #       .x[j,]$mean_median <- mean(n)
    #       .x[j,]$sd_se <- sd(n)
    # 
    #     } else if(sum(myvars) == 4 | sum(myvars) == 6){
    # 
    #       # Only min and mean
    # 
    #       if(.x[j,]$sample_size > 1){
    # 
    #       .x[j,]$sd_se <- (varmax-.x[j,]$min)/(2*(qnorm((max(.x[j,]$sample_size, 1, na.rm = TRUE)-0.375)/(max(.x[j,]$sample_size, 1, na.rm = TRUE)+0.25),0,1)))
    # 
    #       } else {
    # 
    #         w <- .x[, c("sd_se", "sample_size")]
    #         w <- w[complete.cases(w),]
    #         wm <- weighted.mean(w$sd_se, w = w$sample_size, na.rm = TRUE)
    #         if(nrow(.x) > 1) .x[j,]$sd_se <- wm
    # 
    #       }
    #     }
    #   }
    # }

      if(verbose){
      cat("Imputed data --------\n\n")
      print(.x)
      cat("\n---------------------\n\n")
      }
      
    if(nrow(.x) > 1) {
      
      # Apply the inverse-variance method to calculate a weighted mean, using adjusted random-effects weights

      if(any(is.na(.x$sd_se))){
        
        combined.mean <- mean(.x$mean_median)
        combined.se <- NA
        
      } else {
      
      weights <- 1 / ((.x$sd_se^2) + sd(.x$mean_median)^2)
      combined.mean <- sum(weights * .x$mean_median) / sum(weights)
      combined.se <- sqrt(sum(weights ^ 2 * .x$sd_se ^ 2) / sum(weights ^ 2))

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
          d2 <- rescale(x = r_south * (1-rp.SEUS), new.max = dv[1])
          d2[d2==0] <- NA
          d.out <- raster::merge(d1, d2)
          
        } else if (option == 3){
          
          # OPTION 3
          pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs())
          pts <- sp::spsample(pathGSL, n = 1000, type = "regular")
          distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |>
            raster::mask(mask = densr[[1]])
          dist.r <- (1/distances)^(1/4) # Fourth root
          d.out <- rescale(dist.r, new.min = raster::cellStats(dist.r, "min"),
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
          d2 <- rescale(x = r_north * (1-rp.GSL), new.max = dv[1])
          d2[d2==0] <- NA
          d.out <- raster::merge(d1, d2)
          
        } else if (option == 3){
          
          # OPTION 3
          pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs()) |> rgeos::gBuffer(width = 50)
          pts <- sp::spsample(pathGSL, n = 10000, type = "regular")
          distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |> raster::mask(mask = densr[[1]])
          dist.r <- (1/distances)^(1/4) # Fourth root
          d.out <- rescale(dist.r, new.min = raster::cellStats(dist.r, "min"), 
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
  
  out <- do.call(rbind, lapply(seq_along(maplist[[1]]), function(mo) {
    
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
  row.names(out) <- names(maplist[[1]])
  return(out)
}

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
                dt <- data.table::data.table(day = rep(0:(length(dates)-1), times = nsim), 
                                             date = rep(dates, times = nsim),
                                             month = rep(months, times = nsim),
                                             whale = rep(1:nsim, each = length(dates)))
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

optim_feeding <- function(bounds = list(), 
                           err = 0.05,
                           nm = NULL, 
                           cex = 1.15,
                           linecol = "black",
                           verbose = TRUE){
  
  # Assuming B > 0 And S > 0:
  # - A is the value of the horizontal asymptote when x tends to -Inf
  # - D is the value of the horizontal asymptote when x tends to Inf
  # - B describes how rapidly the curve makes its transition between the two asymptotes;
  # - C is a location parameter, which does not have a nice interpretation (except if S = 1)
  # - S describes the asymmetry of the curve (the curve is symmetric when S = 1)
  
  # Define a 5-parameter logistic curve
  f <- function(x, A = 1, D = 0, B, C, S) A + (D-A) / (1 + exp(B*(C-x)))^S
  
  # Define function to optimize
  opt.f <- function(param, bounds, e = err){
    first <- (1 - f(bounds[1], B = param[1], C = param[2], S = param[3]))^2
    second <- (0 - f(bounds[2], B = param[1], C = param[2], S = param[3]) + e)^2
    third <- (1-param[3])^2
    sum(first, second, third)
  }
  
  # Run optimization
  res <- purrr::map(.x = seq_along(bounds), 
                    .f = ~ {
                      out <- optim(par = c(1, 0.1, 0.01), fn = opt.f, bounds = bounds[[.x]])
                      
                      # Print results
                      if(verbose){
                        cat("\n", nm[.x], "\n")
                        cat("-----------------\n")
                        cat("\nEstimated parameter values:\n")
                        cat("B",out$par[1],"\n")
                        cat("C",out$par[2],"\n")
                        cat("S",out$par[3],"\n\n")
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
  tg <- paste("Target BC: \n", nm)
  ypos <- c(0.7, 0.75)
  plot(x, f(x, B = res[[1]]$par[1], C = res[[1]]$par[2], S = res[[1]]$par[3]), 
       type = "n", ylab = "Feeding effort", xlab = "Body condition (relative blubber mass, %)",
       cex.axis = cex, cex.main = cex, cex.lab = cex)
  purrr::walk(.x = seq_along(bounds), .f = ~{
    abline(v = bounds[[.x]][2], lty = 3)
    text(bounds[[.x]][2]+ 0.01, ypos[.x], tg[[.x]], cex = 0.9, adj = 0) 
    lines(x, f(x, B = res[[.x]]$par[1], C = res[[.x]]$par[2], S = res[[.x]]$par[3]), lty = .x, col = linecol)
  })
  if(!is.null(nm)) legend("topright", legend = nm, lty = seq_along(bounds), col = linecol, cex = cex)
  
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

#' predict_backup <- function(obj,
#'          yrs = 35,
#'          n = 100,
#'          popr = 1,
#'          do.plot = FALSE,
#'          seed = 125897,
#'          ...) {
#'   
#'   set.seed(seed)
#'   
#'   if(sum(suppressWarnings(purrr::map_lgl(.x = obj$gam$pred, .f = ~any(is.na(.x)))) > 0)) 
#'     stop("Insufficient sample size. Cannot make predictions.") 
#'   
#'   # Function ellipsis –– optional arguments
#'   args <- list(...)
#'   
#'   # Default values
#'   if(length(args) == 0){
#'     spline <- TRUE
#'     progress <- TRUE
#'   } else {
#'     if("spline" %in% names(args)) spline <- args[["spline"]] else spline <- TRUE
#'     if("progress" %in% names(args)) progress <- args[["progress"]] else progress <- TRUE
#'   }
#'   
#'   cat("Initializing ...\n")
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
#'   # test <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = seq(0, find_maxBC(), 0.01)), type = "link", se.fit = TRUE)
#'   # plot(seq(0,find_maxBC(), 0.01), plogis(test$fit), type = "l")
#'   # lines(seq(0, find_maxBC(), 0.01), plogis(Reduce("+", test)), lty = 2)
#'   # lines(seq(0, find_maxBC(), 0.01), plogis(Reduce("-", test)), lty = 2)
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
#'                    "tot_mass", "lean_mass", "bc", "mass_a", "mass_b", "p_surv", "min_bc", "trest", "t2calf", "birth")
#'   
#'   # Current year
#'   current.yr <- lubridate::year(lubridate::now())
#'   
#'   # Extract terminal functions
#'   mod <- obj$gam$fit
#'   mod[["gest"]] <- gam_gest
#'   
#'   #'------------------------------------------------------
#'   # GAM PARAMETERS
#'   #'------------------------------------------------------
#'   
#'   mbc_preds <- obj$gam$pred$bc_gest
#'   bc_preds <- obj$gam$pred$bc
#'   surv_preds <- obj$gam$pred$surv
#'   
#'   #'------------------------------------------------------
#'   # INITIALIZATION
#'   #'------------------------------------------------------
#'   
#'   # Define initial population vector
#'   # c(F1 = 0, F2 = 0, F3 = 2, F4 = 6, F5 = 4, F6 = 0, F7 = 6, F8 = 0, 
#'   # F9 = 5, "F10+" = 48, FC = 7, FR = 0, FB = 61, 
#'   # M1 = 0, M2 = 0, 
#'   # M3 = 2, M4 = 5, 
#'   # "M5+" = 212)
#'   N0 <- c(0, 7, 212, 0, 71, 1, 7, 61)
#'   names(N0) <- cohorts[, name]
#'   totn <- sum(N0)* (1 + popr)
#'   
#'   cat("Setting up ...\n")
#'   
#'   # Create matrices and initialize them
#'   # rows: years <yrs>
#'   # columns: <attributes>
#'   # layers: individuals <n>
#'   # 4th dimension: replicate projection -> later converted to list
#'   narw.indiv <- array(data = NA, c(yrs + 1, length(mat.attribs), totn, n), 
#'                       dimnames = list(paste0("yr ", 0:yrs), 
#'                                       mat.attribs,
#'                                       paste0("whale ", 1:totn),
#'                                       paste0("prj ", 1:n)))
#'   
#'   cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts), FUN = function(x) rep(cohorts$id[x], each = N0[x])))
#'   
#'   animals <- 1:sum(N0)
#'   
#'   # Alive and population cohort
#'   narw.indiv[1, "alive", animals, ] <- rep(1, )
#'   narw.indiv[1, "cohort", animals, ] <- rep(cohort.vec, n)
#'   
#'   # Sex
#'   #  -- Calves (male)
#'   narw.indiv[, "female", 1:N0[1], ] <- 0
#'   #  -- Calves (female)
#'   fem <- which(cohort.vec == 0)
#'   narw.indiv[, "female", fem[fem > N0[1]], ] <- 1
#'   #  -- Juveniles and adults
#'   narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6)), ] <- 1
#'   narw.indiv[, "female", which(cohort.vec %in% c(1, 3)), ] <- 0
#'   
#'   cat("Age ...\n")
#'   
#'   # Age
#'   ages <- start_age_vec(rep(cohort.vec, n))
#'   narw.indiv[1, "age", animals, ] <- ages
#'   
#'   cat("Length ...\n")
#'   
#'   # Total body length
#'   l.params <- agL_vec(ages)
#'   lengths <- age2length_vec(ages, l.params)
#'   narw.indiv[1, "length", animals, ] <- lengths
#'   
#'   narw.indiv[, "length_a", animals, ] <- rep(l.params[, 1], each = yrs + 1)
#'   narw.indiv[, "length_b", animals, ] <- rep(l.params[, 2], each = yrs + 1)
#'   narw.indiv[, "length_c", animals, ] <- rep(l.params[, 3], each = yrs + 1)
#'   
#'   cat("Mass ...\n")
#'   
#'   # Lean mass
#'   m.params <- mL(n * sum(N0))
#'   narw.indiv[, "mass_a", animals, ] <- rep(m.params[, 1], each = yrs + 1)
#'   narw.indiv[, "mass_b", animals, ] <- rep(m.params[, 2], each = yrs + 1)
#'   mass <- length2mass_vec(lengths, m.params)
#'   narw.indiv[1, "lean_mass", animals, ] <- mass
#'   
#'   # Body conditon
#'   bc <- start_bcondition_vec(rep(cohort.vec, n))
#'   narw.indiv[1, "bc", animals, ] <- bc
#'   
#'   # Total mass
#'   narw.indiv[1, "tot_mass", animals, ] <- mass / (1-bc)
#'   
#'   # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
#'   # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
#'   # This corresponds to a Normal (7.7, 1.45)
#'   # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
#'   # 
#'   # Stewart et al. 2022 -- 
#'   # The degree to which the energetic reserves of females are depleted during lactation may govern
#'   # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).
#'   
#'   cat("Reproduction ...\n")
#'   
#'   t2calf <- (rep(cohort.vec, n) == 6) * random_int(sum(N0) * n)
#'   narw.indiv[1, "t2calf", animals, ] <- t2calf
#'   narw.indiv[1, "trest", animals, ] <- as.numeric(ifelse(t2calf == 0, 13, 1)) * (as.numeric(narw.indiv[1, "cohort", animals, ]) == 6)
#'   
#'   if (!spline) {
#'     narw.indiv[1, "min_bc", animals, ] <- predict_m(
#'       model = mod,
#'       values = as.vector(narw.indiv[1, "tot_mass", animals, ]),
#'       prediction = "gest"
#'     ) * as.vector(narw.indiv[1, "cohort", animals, ] == 6)
#'   } else {
#'     narw.indiv[1, "min_bc", animals, ] <- mbc_preds(narw.indiv[1, "tot_mass", animals, ]) * (narw.indiv[1, "cohort", animals, ] == 6)
#'   }
#'   
#'   narw.indiv[1, "birth", animals, ] <- ifelse(narw.indiv[1, "trest", animals, ] == 13 & narw.indiv[1, "t2calf", animals, ] == 0, 1, 0)
#'   narw.indiv[1, "p_surv", animals, ] <- 1
#'   
#'   cat("List ...\n")
#'   
#'   #' ---------------------------
#'   # IMPORTANT 
#'   #' ---------------------------
#'   # Turn array into a list
#'   narw.indiv <- purrr::array_branch(narw.indiv, 4) # bottleneck
#'   
#'   cat("totpop ...\n")
#'   
#'   # Number of individuals in each cohort
#'   narw.pop <- array(
#'     data = NA, c(n, yrs + 1, nrow(cohorts)),
#'     dimnames = list(
#'       paste0("prj ", 1:n),
#'       paste0("yr ", 0:yrs),
#'       cohorts$name
#'     )
#'   )
#'   narw.pop[, 1, ] <- rep(N0, each = n)
#'   
#'   cat("totpopinit ...\n")
#'   
#'   # Total population size
#'   tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj", 1:n), paste0("yr ", 0:yrs)))
#'   tot.pop[, 1] <- sum(N0)
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
#'   cat("Running projections ...\n")
#'   start.time <- Sys.time()
#'   
#'   for(prj in 1:n){
#'     
#'     if(progress) pb$tick() # Update progress bar
#'     
#'     animals <- 1:sum(N0)
#'     
#'     #'------------------------------------------------------
#'     # Loop over years
#'     #'------------------------------------------------------
#'     
#'     for(i in 2:(yrs+1)){
#'       
#'       current.dat <- as.matrix(narw.indiv[[prj]][i-1, , animals])
#'       alive <- current.dat["alive", animals] * (current.dat["age", animals] <=69)
#'       
#'       #' ----------------------------
#'       # SURVIVAL
#'       #' ----------------------------
#'       # Predict survival probability based on body condition
#'       
#'       if(!spline){
#'         
#'         ps <- alive * predict_m(model = mod, cohort = current.dat["cohort",animals], 
#'                                 values = current.dat["bc",animals], prediction = "surv")
#'         
#'       } else {
#'         
#'         ps <- alive * (surv_preds[["0"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 0) +
#'                          surv_preds[["1"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 1) +
#'                          surv_preds[["2"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 2) +
#'                          surv_preds[["3"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 3) +
#'                          surv_preds[["4"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 4) +
#'                          surv_preds[["5"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 5) +
#'                          surv_preds[["6"]](current.dat["bc",animals]) * (current.dat["cohort",animals] == 6))
#'       }
#'       
#'       # ps <- 1
#'       narw.indiv[[prj]][i, "p_surv", animals] <- ps
#'       
#'       # Determine whether the animal survived
#'       alive <- rbinom(n = animals, size = 1, prob = ps) * (current.dat["age", animals] <=69)
#'       narw.indiv[[prj]][i, "alive", animals] <- alive
#'       
#'       # Sex remains the same
#'       narw.indiv[[prj]][i, "female", animals] <- current.dat["female", animals]
#'       
#'       #' ----------------------------
#'       # GROWTH
#'       #' ----------------------------
#'       
#'       # Increment age
#'       narw.indiv[[prj]][i, "age", animals] <- alive * (current.dat["age", animals] + 1)
#'       
#'       # Increment length
#'       newLp <- agL_vec(animals)
#'       
#'       narw.indiv[[prj]][i,"length_a", animals] <- 
#'         ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,1], current.dat["length_a", animals])
#'       
#'       narw.indiv[[prj]][i,"length_b", animals] <- 
#'         ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,2], current.dat["length_b", animals])
#'       
#'       narw.indiv[[prj]][i,"length_c", animals] <- 
#'         ifelse(narw.indiv[[prj]][i, "age", animals] == 1, newLp[,3], current.dat["length_c", animals])
#'       
#'       narw.indiv[[prj]][i, "length", animals] <-
#'         alive * age2length_vec(
#'           narw.indiv[[prj]][i, "age", animals],
#'           t(narw.indiv[[prj]][i, c("length_a", "length_b", "length_c"), animals])
#'         )
#'       
#'       # Increment lean mass
#'       narw.indiv[[prj]][i, "lean_mass", animals] <- alive * length2mass_vec(narw.indiv[[prj]][i, "length", animals],
#'                                                                             t(narw.indiv[[prj]][i, c("mass_a", "mass_b"), animals]), lean = TRUE)
#'       
#'       # Predict new body condition from current body condition
#'       if (!spline) {
#'         narw.indiv[[prj]][i, "bc", animals] <- alive * predict_m(
#'           model = mod, cohort = current.dat["cohort", animals],
#'           values = current.dat["bc", animals], prediction = "bc")
#'       } else {
#'         narw.indiv[[prj]][i, "bc", animals] <-
#'           alive * (bc_preds[["0"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 0) +
#'                      bc_preds[["1"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 1) +
#'                      bc_preds[["2"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 2) +
#'                      bc_preds[["3"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 3) +
#'                      bc_preds[["4"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 4) +
#'                      bc_preds[["5"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 5) +
#'                      bc_preds[["6"]](current.dat["bc", animals]) * (current.dat["cohort", animals] == 6))
#'       }
#'       
#'       # Increment total mass
#'       narw.indiv[[prj]][i, "tot_mass", animals] <-
#'         alive * narw.indiv[[prj]][i, "lean_mass", animals] / (1 - narw.indiv[[prj]][i, "bc", animals])
#'       
#'       #' ----------------------------
#'       # REPRODUCTION
#'       #' ----------------------------
#'       
#'       # Which animals are resting females?
#'       rest.females <- (current.dat["cohort", animals] == 6)
#'       
#'       # Which animals are juvenile females that are ready to start reproducing
#'       juvenile.females.ofage <- (current.dat["cohort", animals] == 2) * (current.dat["age", animals] >= 9)
#'       
#'       # Which animals calved in previous step?
#'       prev.births <- current.dat["birth", animals]
#'       
#'       newt2calf <- ifelse(prev.births == 1, random_int(sum(prev.births), lwr = 1), 
#'                           ifelse(current.dat["t2calf", animals] == 0, 0, 
#'                                  current.dat["t2calf", animals] - 1))
#'       
#'       # Time spent in resting state - only incremented if calving event hasn't occurred, otherwise reset
#'       narw.indiv[[prj]][i, "t2calf", animals] <- alive * rest.females * newt2calf
#'       
#'       # Years until next calving event
#'       narw.indiv[[prj]][i, "trest", animals] <- 
#'         alive * (rest.females | juvenile.females.ofage) * 
#'         ifelse(narw.indiv[[prj]][i - 1, "trest", animals] == 13, 1, current.dat["trest", animals] + 1)
#'       
#'       # Minimum body condition needed to successfully bring fetus to term without starving
#'       # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
#'       
#'       if (!spline) {
#'         narw.indiv[[prj]][i, "min_bc", animals] <-
#'           alive * predict_m(model = mod, values = narw.indiv[[prj]][i, "tot_mass", animals], prediction = "gest") * rest.females
#'       } else {
#'         narw.indiv[[prj]][i, "min_bc", animals] <-
#'           alive * mbc_preds(narw.indiv[[prj]][i, "tot_mass", animals]) * rest.females
#'       }
#'       
#'       # Birth of new calf, conditional on the mother being alive, in pregnant state
#'       narw.indiv[[prj]][i, "birth", animals] <- alive * (current.dat["cohort", animals] == 4)
#'       
#'       # Maturity - transitions between cohorts
#'       narw.indiv[[prj]][i, "cohort", animals] <-
#'         alive * increment_cohort(
#'           cohort = narw.indiv[[prj]][i - 1, "cohort", animals],
#'           age = narw.indiv[[prj]][i, "age", animals],
#'           female = narw.indiv[[prj]][i, "female", animals],
#'           bc = narw.indiv[[prj]][i, "bc", animals],
#'           min_bc = narw.indiv[[prj]][i, "min_bc", animals],
#'           trest = narw.indiv[[prj]][i, "trest", animals],
#'           t2calf = narw.indiv[[prj]][i, "t2calf", animals])
#'       
#'       new.births <- sum(narw.indiv[[prj]][i, "birth", animals])
#'       
#'       if(new.births > 0){
#'         narw.indiv[[prj]][i, , (max(animals)+1):(max(animals)+new.births)] <- add_calf(n = new.births, attr = mat.attribs)
#'         animals <- 1:(length(animals) + new.births)
#'       }
#'       
#'       # Number of animals in each cohort
#'       # Calves (male)
#'       narw.pop[prj, i, 1] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
#'                                    (narw.indiv[[prj]][i, "female", animals] == 0) *
#'                                    (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       
#'       # Juveniles and adults (male)
#'       narw.pop[prj, i, 2] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 1) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       narw.pop[prj, i, 3] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 3) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       
#'       # Calves (female)
#'       narw.pop[prj, i, 4] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 0) *
#'                                    (narw.indiv[[prj]][i, "female", animals] == 1) *
#'                                    (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       
#'       # Juvenile and reproductive adults (female)
#'       narw.pop[prj, i, 5] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 2) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       narw.pop[prj, i, 6] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 4) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       narw.pop[prj, i, 7] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 5) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       narw.pop[prj, i, 8] <- sum((narw.indiv[[prj]][i, "cohort", animals] == 6) * (narw.indiv[[prj]][i, "alive", animals] == 1))
#'       
#'       # Total population size
#'       tot.pop[prj, i] <- sum(narw.indiv[[prj]][i, "alive", animals], na.rm = TRUE)
#'       
#'     } # End years
#'   } # End projections
#'   
#'   end.time <- Sys.time()
#'   run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
#'   
#'   cat("Processing outputs ...\n")
#'   
#'   # narw.out <- purrr::map(.x = 1:n, .f = ~{
#'   #   reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
#'   #     dplyr::mutate(attr = mat.attribs[attr]) |> 
#'   #     tidyr::pivot_wider(names_from = attr, values_from = value) |> 
#'   #     dplyr::mutate(prj = .x) |> 
#'   #     dplyr::relocate(prj, .before = yr)
#'   # }) |> do.call(what = rbind) |> 
#'   #   data.table::data.table()
#'   
#'   # narw.out <- narw.out[is.finite(rowSums(narw.out)),]
#'   
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
#'   # births.df <- purrr::map(.x = 1:n, .f = ~{
#'   #   m <- matrix(rowSums(narw.indiv[[.x]][2:(yrs+1),"birth",], na.rm = TRUE), ncol = 1)
#'   #   colnames(m) <- .x
#'   #   m
#'   # }) |> do.call(what = cbind) |> 
#'   #   tibble::as_tibble() |> 
#'   #   tibble::rownames_to_column(var = "year") |> 
#'   #   dplyr::mutate(year = as.numeric(year)) |> 
#'   #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |> 
#'   #   dplyr::select(prj, year, birth) |> 
#'   #   dplyr::arrange(prj, year) |> 
#'   #   data.table::data.table()
#'   # 
#'   # deaths.df <- purrr::map(.x = 1:n, .f = ~{
#'   #   m <- matrix(apply(X = narw.indiv[[.x]][2:(yrs+1),"alive",],
#'   #                     MARGIN = 1,
#'   #                     FUN = function(x) {
#'   #     r <- x[!is.na(x)]
#'   #     r <- sum(r == 0)
#'   #     r
#'   #     }), ncol = 1)
#'   #   colnames(m) <- .x
#'   #   m
#'   # }) |> do.call(what = cbind) |>
#'   #   tibble::as_tibble() |>
#'   #   tibble::rownames_to_column(var = "year") |>
#'   #   dplyr::mutate(year = as.numeric(year)) |>
#'   #   tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
#'   #   dplyr::select(prj, year, death) |>
#'   #   dplyr::arrange(prj, year) |>
#'   #   data.table::data.table()
#'   
#'   narw.conf <- narw.df[
#'     , list(
#'       mean = mean(N),
#'       lwr = quantile(N, 0.025, na.rm = TRUE),
#'       uppr = quantile(N, 0.975, na.rm = TRUE)
#'     ),
#'     list(year, cohort)
#'   ] |>
#'     dplyr::mutate(cohort = factor(cohort, levels = c(
#'       "Calves (male)",
#'       "Calves (female)",
#'       "Juveniles (male)",
#'       "Juveniles (female)",
#'       "Adults (male)",
#'       "Adults (female, pregnant)",
#'       "Adults (female, lactating)",
#'       "Adults (female, resting)"
#'     )))
#'   
#'   tot.conf <- tot.df[
#'     , list(
#'       mean = mean(N),
#'       lwr = quantile(N, 0.025),
#'       uppr = quantile(N, 0.975)
#'     ),
#'     list(year, cohort)
#'   ] |>
#'     dplyr::mutate(cohort = factor(cohort, levels = c(
#'       "Calves (male)",
#'       "Calves (female)",
#'       "Juveniles (male)",
#'       "Juveniles (female)",
#'       "Adults (male)",
#'       "Adults (female, pregnant)",
#'       "Adults (female, lactating)",
#'       "Adults (female, resting)",
#'       "North Atlantic right whales"
#'     )))
#'   
#'   p1 <- plot_projection(narw.df, narw.conf)
#'   p2 <- plot_projection(tot.df, tot.conf)
#'   
#'   if(do.plot){
#'     print(p1)
#'     print(p2)
#'   }
#'   
#'   # Find 95% confidence intervals on final population size
#'   cat("Final population size:\n")
#'   final.pop <- unname(tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.5, 0.025, 0.975))])
#'   cat("N = ", round(final.pop[1],0), " (95% CI: ", round(final.pop[2],0), "–", round(final.pop[3],0), ")\n", sep = "")
#'   
#'   cat(paste0("Time elapsed: ", run_time))
#'   cat("\n")
#'   
#'   return(list(
#'     # dat = narw.out,
#'     out = list(df = rbind(narw.df, tot.df), 
#'                ci = rbind(narw.conf, tot.conf),
#'                plot = list(p1, p2))))
#'   
#' }

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
#'   # test <- mgcv::predict.gam(mod$gam[[1]]$surv, newdata = data.frame(start_bc = seq(0,find_maxBC(), 0.01)), type = "link", se.fit = TRUE)
#'   # plot(seq(0,find_maxBC(), 0.01), plogis(test$fit), type = "l")
#'   # lines(seq(0,find_maxBC(), 0.01), plogis(Reduce("+", test)), lty = 2)
#'   # lines(seq(0,find_maxBC(), 0.01), plogis(Reduce("-", test)), lty = 2)
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

#' predict_leslie <- function(obj,
#'                             n = 100,
#'                             yrs = 35) {
#'   
#'   gg.opts <- ggplot2::theme(axis.text = element_text(size = 10, color = "black"),
#'                             axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
#'                             axis.text.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
#'                             axis.title = element_text(size = 12),
#'                             axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
#'                             axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
#'                             strip.background = element_rect(fill = "grey20"),
#'                             strip.text = element_text(colour = 'white', size = 12))
#'   
#'   
#'   # Adapted from original code by Scott Creel
#'   # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
#'   
#'   # Population estimate as per 2022 NARW report card is 340 (+/- 7).
#'   # The percentage of the population estimated 
#'   
#'   # cohort.id <- obj$param$cohort.id
#'   # cohort.ab <- obj$param$cohort.ab
#'   current.yr <- lubridate::year(lubridate::now())
#'   
#'   cohorts <- c("Calves (male)",
#'                "Juveniles (male)", 
#'                "Adults (male)",
#'                "Calves (female)",
#'                "Juveniles (female)", 
#'                "Adults (female, pregnant)",
#'                "Adults (female, lactating)",
#'                "Adults (female, resting)") |> 
#'     tolower() |> abbreviate(minlength = 6) |> tibble::enframe() |> 
#'     dplyr::rename(cohort = name, abbr = value)
#'   
#'   # Vector to store total population size for N replicate projections
#'   narwpop <- purrr::set_names(x = cohorts$cohort) |> 
#'     purrr::map(.f = ~matrix(0, n, yrs, dimnames = list(paste0("proj",1:n), paste0("yr ", 1:yrs))))
#'   
#'   totpop <- matrix(0, n, yrs, dimnames = list(paste0("proj",1:n), paste0("yr ", 1:yrs)))
#'   
#'   # 1: Calves (male)
#'   # 2: Juveniles (male)
#'   # 3: Adults (male)
#'   # 4: Calves (female)
#'   # 5: Juveniles (female)
#'   # 6: Adults (female, pregnant)
#'   # 7: Adults (female, lactating)
#'   # 8: Adults (female, resting)
#'   
#'   p.fem <- 0.5 # Fraction of females at birth
#'   p.birth <- 0.5 # Probability that a pregnant female will give birth to a calf (fecundity)
#'   p.survival <- rep(0.8, nrow(cohorts)) # Probability of survival
#'   names(p.survival) <- cohorts
#'   p.recruit <- c(0.7, 0.7) # Male, female
#'   tau_rp <- 0.7 # Transition prob from resting to pregnant
#'   tau_pl <- 0.7 # Transition prob from pregnant to lactating
#'   tau_lr <- 0.8 # Transition prob from lactating to resting
#'   
#'   popmat <- matrix(0, nrow(cohorts), nrow(cohorts))
#'   popmat[1, 6] <- p.survival[6] * p.birth * (1 - p.fem)
#'   
#'   popmat[2, 1] <- p.survival[1]
#'   popmat[2, 2] <- p.survival[2]
#'   
#'   popmat[3, 2] <- p.survival[2] * p.recruit[1]
#'   popmat[3, 3] <- p.survival[3]
#'   
#'   popmat[4, 4] <- p.survival[6] * p.birth * p.fem
#'   
#'   popmat[5, 5] <- p.survival[4]
#'   popmat[5, 6] <- p.survival[5] * (1 - p.recruit[2])
#'   
#'   popmat[6, 5] <- p.survival[5] * p.recruit[2]
#'   popmat[6, 6] <- p.survival[6] * (1 - tau_pl)
#'   popmat[6, 8] <- p.survival[8] * tau_rp
#'   
#'   popmat[7, 6] <- p.survival[6] * tau_pl
#'   popmat[7, 7] <- p.survival[7] * (1 - tau_lr)
#'   
#'   popmat[8, 7] <- p.survival[7] * tau_lr
#'   popmat[8, 8] <- p.survival[8] * (1 - tau_rp)
#'   
#'   
#'   popmat <- matrix(popmat, nrow = nrow(cohorts), byrow = FALSE, dimnames = list(cohorts$abbr, cohorts$abbr))
#'   popmat <- round(popmat, 5)
#'   
#'   # STOCHASTIC LESLIE PROJECTION
#'   #' ----------------------------
#'   # This uses 4 nested loops. 
#'   # The p loop (outermost loop) replicates the projection <n> times.
#'   # The y loop is next, and steps across all years of projection from an initial population vector.
#'   # The i and j loops are innermost and draw stochastic parameters (body condition and survival)
#'   
#'   # Set up progress bar
#'   pb <- progress::progress_bar$new(
#'     format = "[:bar] :percent eta: :eta",
#'     total = n, clear = FALSE, width = 80
#'   )
#'   
#'   for(p in 1:n){
#'     
#'     pb$tick() # Update progress bar
#'     
#'     # Matrix of age-class values, with one row per time step. 
#'     # The columns represent the numbers of individuals in each age class.
#'     N.mat <- matrix(NA, nrow = yrs, ncol = nrow(cohorts), dimnames = list(paste("yr", 1:yrs), cohorts$abbr))
#'     
#'     # Initial population vector
#'     # Sums to 340 (abundance estimate as of 2022).
#'     N.mat[1,] <- c(7, 10, 163, 8, 10, 50, 15, 77)
#'     
#'     # Begin loop for projections
#'     for(i in 2:yrs){ 		
#'       
#'       # Set up a temporary Leslie matrix for each iteration 
#'       #  with the correct mean fecundities and survival rates
#'       leslie.mat <- popmat 				
#'       
#'       # Randomly draw fecundities for each class
#'       # Mean of Poisson is fecundity value
#'       # for(j in 1:n.stages){
#'       #   leslie.mat[1,j] <- ifelse(popmat[1,j] == 0, 0, popmat[1,j] + rnorm(1, 0, 0.01))
#'       # }
#'       
#'       for(j in 1:nrow(cohorts)){
#'         for(k in 1:nrow(cohorts)){
#'           leslie.mat[k,j] <- ifelse(popmat[k,j] == 0, 0, popmat[k,j] + rnorm(1, 0, 0.01))
#'         }
#'       }
#'       
#'       # Randomly draw survival probabilities for each class.
#'       # n.ind is number of individuals to calculate live/die for
#'       # for(k in 1:(n.stages-1)){
#'       #   n.ind <- N.mat[(i-1),k]
#'       #   # Need ifelse statement to deal with the possibility that there are no individuals in that age class
#'       #   leslie.mat[(k+1),k] <- ifelse(n.ind > 1, rbinom(1,size = round(n.ind), p = les.mat[(k+1),k])/n.ind,0)
#'       # }
#'       
#'       # Matrix multiplication for next time-step.
#'       N.mat[i,] <- leslie.mat %*% N.mat[(i-1),]          
#'       
#'     } # End i loop over time
#'     
#'     for(k in seq_along(narwpop)){
#'       narwpop[[k]][p,] <- N.mat[,k]
#'     }
#'     totpop[p,] <- rowSums(N.mat)
#'     
#'   } # Ends p loop over replicate projections
#'   
#'   # Compile outputs
#'   narw.df <- purrr::imap(.x = narwpop, .f = ~{
#'     tibble::as_tibble(.x) |> 
#'       tibble::rownames_to_column(var = "proj") |> 
#'       tidyr::pivot_longer(!proj, names_to = "year", values_to = "N") |> 
#'       dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
#'       dplyr::mutate(cohort = stringr::str_to_sentence(.y))
#'   }) |> do.call(what = rbind) |> data.table::data.table()
#'   
#'   tot.df <- tibble::as_tibble(totpop) |> 
#'     tibble::rownames_to_column(var = "proj") |> 
#'     tidyr::pivot_longer(!proj, names_to = "year", values_to = "N") |> 
#'     dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
#'     dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
#'   
#'   narw.conf <- narw.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)]
#'   tot.conf <- tot.df[, list(mean = mean(N), lwr = quantile(N, 0.025), uppr = quantile(N, 0.975)), list(year, cohort)]
#'   
#'   # Generate plots
#'   make_plot <- function(df, conf){
#'     ggplot2::ggplot() +
#'       # ggplot2::geom_path(data = df, aes(x = year, y = N, group = factor(proj)), colour = "lightgrey", linewidth = 0.2) +
#'       ggplot2::geom_path(data = conf, aes(x = year, y = mean), colour = "#1565C0") +
#'       ggplot2::geom_ribbon(data = conf, aes(x = year, y = mean, ymin = lwr, ymax = uppr), alpha = 0.25, fill = "#1565C0") +
#'       ggplot2::facet_wrap(~cohort) + 
#'       ggplot2::scale_x_continuous(breaks = pretty(seq(current.yr, current.yr + yrs))) +
#'       ggplot2::scale_y_continuous(breaks = pretty(df$N), labels = scales::comma) +
#'       xlab("") + ylab("Abundance") +
#'       gg.opts
#'   }
#'   
#'   p1 <- make_plot(tot.df, tot.conf)
#'   p2 <- make_plot(narw.df, narw.conf)
#'   
#'   print(p1)
#'   print(p2)
#'   
#'   # Find 95% confidence intervals on final population size
#'   ci <- tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.025, 0.975))]
#'   
#'   return(ci)
#'   
#' }


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



#' Estimate daily p(dead) to match a target annual mortality rate
#' 
#' @param annual.rate Double. Target annual mortality rate (in %)
#' @note If the daily probability of a mortality event occurring is p, then the daily probability of
#' survival is (1-p). The probability of survival over a year is then given by (1-p)^365.
#' Therefore the probability of a mortality event occurring at least once during the year
#' (annual mortality rate) is p_year = 1-((1-p)^365). It follows that p_day = 1-((1-p_year)^(1/365))
#' @return Daily mortality probability.
find_mortality <- function(annual.rate = 0.02){
  p_daily = 1-((1-annual.rate)^(1/365))
  return(p_daily)
}

# find_mortality <- function(annual.rate = 0.02){
#   
#   opt.fun <- function(x, r, verbose = FALSE){
#     
#     r0 <- rbinom(n = 100, size = 365, prob = x) 
#     r0 <- length(r0[r0>0])/100
#     
#     r1 <- rbinom(n = 1000, size = 365, prob = x)
#     r1 <- length(r1[r1>0])/1000
#     
#     r2 <- rbinom(n = 10000, size = 365, prob = x)
#     r2 <- length(r2[r2>0])/10000
#     
#     r3 <- rbinom(n = 100000, size = 365, prob = x)
#     r3 <- length(r3[r3>0])/100000
#     
#     out <- (r0-r)^2 + (r1-r)^2 + (r2-r)^2 + (r3-r)^2
#     
#     if(verbose){
#       cat("Target p:", rate, "\n")
#       cat("-------------------\n")
#       cat("N = 100:", r0, "\n")
#       cat("N = 1,000:", r1, "\n")
#       cat("N = 10,000:", r2, "\n")
#       cat("N = 100,000:", r3, "\n\n")
#     }
#     
#     return(out)
#   }
#   
#   est <- optim(par = 0.001, fn = opt.fun, r = rate)
#   opt.fun(x = est$par, r = rate, verbose = TRUE)
#   return(est$par)
# }

scale_growth <- function(target.max = 45000, maxBC = find_maxBC(cohort = 2), age = 69){
  # maxBC (previous) = 0.6586636
  # When lean = 1, length2mass returns the total mass of an animal given its length(age), as per Fortune et al. (2021)
  # By setting lean < 1, we can scale down this relationship to obtain a growth curve for lean mass.
  # We use an optimisation to estimate a scalar that will prevent unrealistic values of total mass (>45t) 
  # given fluctuating values of reserve mass and a maximum body condition of ca. 60%.
  opt.fun <- function(x, age, target.max, maxBC){
    return((length2mass(age2length(age, agL(age)), mL(), x)/(1-maxBC)-target.max)^2)}
  out <- optimize(f = opt.fun, interval = c(0,1), target.max = target.max, maxBC = maxBC, age = age)
  return(out$minimum)
}

find_lean <- function(max.mass = 45000, maxBC = find_maxBC(cohort = 6), age = 69){
  
  # Mean total mass from published growth curve 
  mt <- length2mass(age2length(69, agL(69)), mL(), 1)
  
  # Maximum lean mass given max total mass and max BC
  mlm <- (1-maxBC)*max.mass
  
  # Resulting scalar - we basically
  return(mlm/mt)
  
  # opt.fun <- function(x, age, target.max, maxBC){
  #   return(((length2mass(age2length(age, agL(age)), mL(), x)/(1-maxBC))-target.max)^2)}
  # out <- optimize(f = opt.fun, interval = c(0,1), target.max = max.mass, maxBC = maxBC, age = age)
  # return(out$minimum)
}

find_maxBC <- function(cohort = 4, threshold = 0.05){
  
  if(cohort == 4){
  myfunc <- function(par){
   ( threshold - feeding_effort(30, 0.6, par))^2
  }
  } else {
    myfunc <- function(par){
      ( threshold - feeding_effort_vec(10, 0.41718, par))^2
    }
  }
  optimize(myfunc, c(0,1), tol = 0.0001)$minimum
}


gestation_threshold <- function(){
  
  set.seed(20230412)
  
  # Simulate pregnant females - set stressors and growth to 0 to ensure that animals do not die
  # Note: Fetus growth is not tied to growth argument
  mg <- narw(nsim = 1000, cohortID = 4, progress = 1, stressors = 0, growth = 0)
  mdat <- mg$sim$`ad(f,p)`
  
  # Energy density of lipids -- (kJ / kg - converted to GJ/kg)
  ED_lipids = 39539/1000000
  
  # Percent lipid breakdown during catabolism
  lipid_cat = 0.80
  
  # Compile data
  m_kj <- mdat[, list(mass = unique(mass), E_gest = sum(E_gest)/1000), whale]
  m_kj[, kg:= E_gest/(ED_lipids*lipid_cat) + 0.05*mass]
  m_kj[, min_bc:=kg/mass]
  
  gam_gest <- mgcv::gam(min_bc ~ s(mass), data = m_kj, method = "REML", family = betar(link = "logit"))
  
  # Residual checks
  gratia::appraise(gam_gest)
  
  # Extract fitted relationship
  b0 <- coef(gam_gest)[1] # Intercept
  invlink <- gam_gest$family$linkinv
  gest.df <- gratia::smooth_estimates(gam_gest) |>
    gratia::add_confint() |>
    dplyr::mutate(across(tidyselect::all_of(c("est", "lower_ci", "upper_ci")), .fns = \(x) x + b0)) |>
    gratia::transform_fun(fun = \(x) invlink(x))
  
  # Create plot
  p <- ggplot2::ggplot(data = gest.df) + 
    geom_point(data = m_kj, aes(x = mass, y = min_bc), alpha = 0.2) +
    geom_ribbon(aes(x = mass, ymin = lower_ci, ymax = upper_ci), alpha = 0.25, fill = "#104E8B") +
    geom_line(aes(mass, est), colour = "#104E8B") + 
    xlab("Total mass (kg)") + ylab("Minimum body condition (%)") +
    theme_narw() +
    scale_x_continuous(breaks = pretty(range(gam_gest_model$model$mass))) +
    scale_y_continuous(breaks = pretty(c(min(gest.df$lower_ci), max(gest.df$upper_ci)), n = 10)) +
    geom_rug(data = gam_gest_model$model, aes(x = mass), sides = "b")
  
  print(p)
  
  return(gam_gest)
  
}

get_pva <- function(){
  
  median_minpop <- readr::read_csv(file = "data/popviability/Linden2023_figure14.csv", 
                                col_types = "dd") |>  
    dplyr::mutate(year = 2023 + year) |> 
    dplyr::mutate(label = "NOAA")
    
  median_pop <- readr::read_csv(file = "data/popviability/Linden2023_figure11A.csv", 
                                 col_types = "dd") |> 
    dplyr::mutate(label = "NOAA")

  return(list(pop = median_pop, minpop = median_minpop))
  
}

get_bc <- function(){
  fred <- readr::read_csv("data/condition/christiansen2022.csv") |> 
    dplyr::filter(tissue == "blubber")
  out.spline <- splinefun(x = fred$index, y = fred$bc)
  return(out.spline)
}

init_bc <- function(calves = c(-0.61, 14.64272),
                    juveniles = c(-13.1, 2.9),
                    adults = c(-16.7, 2),
                    lactating = c(9.4, 4.8)){
  
  # From Christiansen et al. (2022) J Physiol
  # The birth condition value is -0.0061 (or -0.61%) with a SE of 0.0049 (or -0.49%)
  # The sample size is: n = 992
  # This gives an SD = 0.1543306
  # 
  # From Christiansen et al. (2020) MEPS
  # lactating NARW females (mean = −9.4%, SE = 4.8)
  # adults (mean = −16.7%, SE = 2.0)
  # immature NARWs (mean = −13.1%, SE = 2.9)
  
  out <- list()
  
  # Calculate relative blubber mass from body condition indices (mean +/- 2SE)
  # using a spline function <bc_index> fitted to data from Christiansen et al. (2022) -- Figure 7D
  tmp <- c(bc_index(calves[1]/100),
    bc_index(calves[1]/100-2*(calves[2]/100)),
    bc_index(calves[1]/100+2*(calves[2]/100)))
  
  # Estimate associated SD in relative blubber mass
  out[["calves"]] <- c(tmp[1], mean(c((tmp[3] - tmp[1])/2, (tmp[1] - tmp[2])/2)))
  
  tmp <- c(bc_index(juveniles[1]/100),
                bc_index(juveniles[1]/100-2*(juveniles[2]/100)),
                bc_index(juveniles[1]/100+2*(juveniles[2]/100)))
  
  out[["juveniles"]] <- c(tmp[1], mean(c((tmp[3] - tmp[1])/2, (tmp[1] - tmp[2])/2)))
  
  tmp <- c(bc_index(adults[1]/100),
                bc_index(adults[1]/100-2*(adults[2]/100)),
                bc_index(adults[1]/100+2*(adults[2]/100)))
  
  out[["adults"]] <- c(tmp[1], mean(c((tmp[3] - tmp[1])/2, (tmp[1] - tmp[2])/2)))
  
  tmp <- c(bc_index(lactating[1]/100),
                bc_index(lactating[1]/100-2*(lactating[2]/100)),
                bc_index(lactating[1]/100+2*(lactating[2]/100)))
  
  out[["lactating"]] <- c(tmp[1], mean(c((tmp[3] - tmp[1])/2, (tmp[1] - tmp[2])/2)))
  
  return(out)
  
}

# NOISE ------------------------------------------------------

add_SPL <- function(spl, raster = FALSE){
  if(raster){
    10*log(Reduce("+", lapply(X = spl, FUN = function(x) 10^(x/10))), base = 10)
  } else {
    10*log(sum(10^(spl/10)), base = 10)
  }
}

get_turbines <- function(){
  
  # Import CSV of turbine locations
  turbines <- readr::read_csv("data/windfarms/turbine_locs.csv", col_types = "fdd") |> 
    janitor::clean_names() |> 
    dplyr::rename(windfarm = farm)
  
  if(!all(names(turbines) %in% c("windfarm", "longitude", "latitude"))) 
    stop("Cannot find all required fields")
  
  # Add projected coordinates and split into a list
  turbines <- add_xy(dat = turbines)
  turbines <- split(x = turbines, f = turbines$windfarm)
  
  # Month of simulation start
  init.month <- 10
  
  # Load turbine locations for each wind farm
  no.piles <- purrr::map_dbl(.x = turbines, .f = ~nrow(.x))
  
  # Assume construction occurs N-S and W-E
  turbines <- purrr::map(.x = turbines, .f = ~dplyr::arrange(.x, -latitude, longitude))
  turbines <- purrr::map(.x = turbines, .f = ~dplyr::mutate(.x, sortorder = 1:nrow(.x)))
  
  sortorder <- purrr::map(.x = no.piles, .f = ~split(seq_len(.x), ceiling(seq_len(.x)/10))) |> 
    purrr::map_depth(.depth = 2, .f = ~rev(.x)) |> 
    purrr::map(.f = ~unlist(.x)) 
  
  turbines <- purrr::map2(.x = turbines, .y = sortorder, .f = ~dplyr::mutate(.x, sortorder = .y))
  turbines <- purrr::map(.x = turbines, .f = ~dplyr::arrange(.x, sortorder))
  
  # Define start dates for piling activities in each scenario
  current.year <- lubridate::year(lubridate::now())
  
  start.dates <- list(
    
    # Scenario 1 – Jan 15th start at Farm 1, Feb 1st start at Farm 2, and Jan 1st start at Farm 3
    c(paste0(current.year+1, stringr::str_pad(1, width = 2, pad = "0"), "15"),
      paste0(current.year+1, stringr::str_pad(2, width = 2, pad = "0"), "01"),
      paste0(current.year+1, stringr::str_pad(1, width = 2, pad = "0"), "01")),
    
    # Scenario 2 – Jul 15th start at Farm 1, Sep 15th start at Farm 2, and May 1st start at Farm 3
    c(paste0(current.year+1, stringr::str_pad(7, width = 2, pad = "0"), "15"),
      paste0(current.year+1, stringr::str_pad(9, width = 2, pad = "0"), "15"),
      paste0(current.year+1, stringr::str_pad(5, width = 2, pad = "0"), "01"))
  )
  
  piling.dates <- lapply(X = 1:2, FUN = function(sc) {
    purrr::map2(
      .x = start.dates[[sc]],
      .y = no.piles,
      .f = ~ {
        out <- seq(from = lubridate::ymd(.x), by = "day", length.out = .y)
        if(any(grepl(pattern = "02-29", x = out))){
          out <- out[which(!grepl(pattern = "02-29", x = out))]
          out <- c(out, out[length(out)] + 1)
        } 
        out
      }
    ) |> purrr::set_names(nm = paste0("windfarm_0", names(no.piles)))
  }) |> purrr::set_names(nm = c(1,2)) |> 
    tibble::enframe() |> 
    tidyr::unnest(cols = c("value")) |> 
    dplyr::mutate(windfarm = rep(1:3, 2)) |> 
    tidyr::unnest(cols = c("value")) |>
    dplyr::rename(scenario = name, date = value) |> 
    dplyr::mutate(scenario = as.numeric(scenario)) |> 
    data.table::as.data.table()
  
  # Compile data
  scenarios <- lapply(X = 1:2, FUN = function(sc) {
    purrr::map(
      .x = 1:3,
      .f = ~ cbind(turbines[[.x]], piling.dates[scenario == sc & windfarm == .x, "date"])
    ) |> do.call(what = "rbind") |> tibble::as_tibble() |> dplyr::select(-sortorder)
  }) |> purrr::set_names(nm = paste0("scenario_0", 1:2))
  
  # Note: We only allow N=60 turbines to be installed at site 1 under scenario 2 to
  # align with Southall et al. (2021)
  scenarios[[2]] <- split(scenarios[[2]], f = factor(scenarios[[2]]$windfarm))
  scenarios[[2]][[1]] <- dplyr::slice(scenarios[[2]][[1]], 1:60)
  scenarios[[2]] <- do.call(rbind, scenarios[[2]])
  
  # Add third scenario - only maintenance operations
  scenarios["scenario_03"] <- list(NULL)
  
  return(scenarios)
  
}

dummy_raster <- function(value = NULL){
  dummy.r <- targets::tar_read(density_support) |> raster::raster()
  dummy.r <- dummy.r - 1
  dummy.r[dummy.r < 0] <- NA
  if(!is.null(value)) dummy.r[dummy.r == 0] <- value
  return(dummy.r)
}

# noise_surface <- function(ambient.db = 80, 
#                           source.lvl = 220,
#                           mitigation = c(0, 10)){
# 
#   # Spatial support
#   support.poly <- targets::tar_read(support_poly)
#   
#   # Turbines
#   turbines <- targets::tar_read(turbines)
#   # data("turbines", envir = environment())
#   
#   # Strip years from dates
#   turbines <- purrr::map(.x = turbines[1:2], .f = ~dplyr::mutate(.data = .x, date = format(as.Date(date), "%d-%m")))
#   
#   # Calculate distances to each source (turbine location)
#   r <- ambient.r <- targets::tar_read(density_support) |> raster::raster()
#   raster::values(r) <- NA
#   ambient.r <- ambient.r - 1
#   ambient.r[ambient.r < 0] <- NA
#   ambient.r[ambient.r == 0] <- ambient.db
#   
#   # Compute distances to turbine locations
#   dist2turbines <- purrr::map(
#     .x = turbines,
#     .f = ~ {
#       .x$dsource <- lapply(X = 1:nrow(.x), FUN = function(turb) {
#         tmp <- r
#         tmp[raster::cellFromXY(tmp, data.frame(.x[turb, c("x", "y")]))] <- 1
#         # If y is missing from the distance function, computes the distance for all cells that are NA to the nearest
#         # cell that is not NA
#         terra::distance(tmp) |> 
#           terra::mask(mask = support.poly)
#       })
#       .x
#     }
#   )
#   
#   # Piling dates & days of the year
#   piling_dates <- purrr::map(.x = turbines, .f = ~unique(.x$date))
#   date_seq <- get_dates()
#   
#   # Compute sound levels
#   soundLvls <- purrr::map2(.x = dist2turbines,
#                       .y = mitigation,
#                       .f = ~ {
#      .x$db <- lapply(X = 1:nrow(.x), FUN = function(d) {
#       db.out <- (source.lvl - .y) - TL(.x$dsource[[d]])
#       db.out[db.out < ambient.db] <- NA
#       db.out
#     })
#     .x
#   })
#   
#   # Split by date
#   soundLvls_bydate <- purrr::map(.x = soundLvls, .f = ~{
#     split(.x, f = factor(.x$date, levels = unique(.x$date)))
#   })
#   
#   # Sum sound fields
#   sound.layer <- purrr::map_depth(.x = soundLvls_bydate, .depth = 2, .f = ~ {
#     db.lyr <- raster::stack(.x$db) |>
#       raster::calc(fun = max) |>
#       raster::stack(y = ambient.r) |>
#       raster::calc(fun = sum, na.rm = TRUE)
#     db.lyr[db.lyr == 0] <- NA
#     rescale(db.lyr, new.min = ambient.db, new.max = source.lvl)
#   })
#   
#   sound.out <- lapply(X = names(sound.layer), FUN = function(sc) {
#     purrr::map(.x = date_seq[-1], .f = ~ {
#       if (.x %in% names(sound.layer[[sc]])){
#         sound.layer[[sc]][[.x]]
#       } else {
#         ambient.r
#       }
#     })
#   })
#   
#   names(sound.out) <- names(turbines)
#   names(sound.out[[1]]) <- date_seq[-1]
#   names(sound.out[[2]]) <- date_seq[-1]
#   
#   out <- purrr::map_depth(.x = sound.out, .depth = 2, .f = ~as(.x, "SpatialGridDataFrame"))
#   
#   out["scenario_03"] <- list(NULL)
#   return(out)
#   
# }

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

dB2km <- function(target.dB, SL = 220, logfac = 15, a = 1.175, mitigation = 0){
  opt.fun <- function(r, SL, target.L, logfac){
    return(((SL-mitigation)-TL(r, logfac = logfac, a = a)-target.L)^2)}
  out <- optimize(f = opt.fun, interval = c(0,30000), SL = SL, target.L = target.dB, logfac = logfac)
  return(out$minimum)
}

km2dB <- function(r, SL = 220, logfac = 15, a = 1.175, mitigation = 0){
  return((SL-mitigation)-TL(r, logfac = logfac, a = a))
}

get_doseresponse <- function(input.ee){
  data.table::fread("data/elicitation/EE_results.csv", select = 2:3) |> as.matrix()
  # data.table::fread("data/elicitation/EE_results.csv", select = 4:5003) |> as.matrix()
}

proxy_noise <- function(ambient = 80, source.lvl = 220, mitigation = 10, x, y){
  
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
  
  db.all <- rescale(db.all, new.min = ambient, new.max = source.lvl)
  db.list <- purrr::map(.x = month.abb, .f = ~db.all)
  
  out <- purrr::map(.x = db.list, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}

# FISHING GEAR ------------------------------------------------------

entgl_surface_risk <- function(scalar = 1){
  
  # Retrieve density layers & regions
  d <- targets::tar_read(density_narw)
  regions <- targets::tar_read(regions)
  support.poly <- targets::tar_read(support_poly)
  
  GSL_scalar <- 1.114327 # Gulf of St Lawrence
  ELSWC_scalar <- 0.03025689 # Elsewhere Canada
  
  # Rasterize regions polygons
  regions.raster <- raster::rasterize(regions[, "region"], raster::raster(d[[1]]))
  
  # Identify cells within regions of interest
  GSL.cells <- which(regions.raster[] == which(regions$region == "GSL"))
  ELSWC.cells <- which(regions.raster[] %in% which(regions$region %in% c("CABOT", "SCOS", "BOF_lower", "BOF_upper")))
  
  purrr::map(.x = seq_along(month.abb), .f = ~{
    
    cat("Calculating risk for ...", month.abb[.x], "\n")
    
    rd <- raster::raster(d[[.x]])
    
    # Import risk raster for U.S. waters
    risk.r <- raster::raster(x = paste0(
      "data/gear/Risk_PostAction/DST_Risk_PostAction_m",
      stringr::str_pad(.x, 2, pad = "0"), ".tif"
    )) |> raster::projectRaster(to = rd)
    
    # Divide the DST risk values by whale densities in each cell
    risk.scaled <- risk.r / rd
    
    # Scale all values by the maximum total risk across regions
    risk.max <- risk.scaled / 6435.2916
    
    # Multiply by PCoMS scalar to estimate entanglement probability
    # Posterior median and 95% credible interval: 0.45457 [0.433815 - 0.47574]
    risk.prob <- risk.max * 0.45457
    
    # Sum of line numbers in U.S. waters
    US.waters <- sum(risk.prob[], na.rm = TRUE)
    
    # Apply scalars
    GSL.waters <- GSL_scalar * US.waters
    ELSWC.waters <- ELSWC_scalar * US.waters
    
    risk.all <- regions.raster
    risk.all[!is.na(risk.all)] <- 0
    risk.all <- raster::merge(risk.prob, risk.all)
    
    risk.all[GSL.cells] <- scalar * (GSL.waters / length(GSL.cells))
    risk.all[ELSWC.cells] <- scalar * (ELSWC.waters / length(ELSWC.cells))
    
    risk.all <- raster::mask(risk.all, support.poly)
    
    as(risk.all, "SpatialGridDataFrame")
    
  }) |> purrr::set_names(nm = month.abb)
}

# entgl_surface_lines <- function(scalar = 1){
#   
#   # Retrieve density layers
#   d <- targets::tar_read(density_narw)
#   regions <- targets::tar_read(regions)
#   rd <- raster::raster(d[[1]])
#   
#   GSL_scalar <- 1.114327 # Gulf of St Lawrence
#   ELSWC_scalar <- 0.03025689 # Elsewhere Canada
#   
#   purrr::map(.x = 1:12, .f = ~{
#     
#     print(month.abb[.x])
# 
#     # Import raster of line numbers for U.S. waters
#     lines.r <- raster::raster(x = paste0("data/gear/LineNum_PostAction/DST_LineNum_PostAction_m", 
#                                          stringr::str_pad(.x, 2, pad = "0"), ".tif")) |> 
#       raster::projectRaster(to = rd)
#     
#     # lines.r <- raster::raster(x = paste0("data/gear/LineNum_PostAction/DST_LineNum_PostAction_m", 
#     #                                      stringr::str_pad(.x, 2, pad = "0"), ".tif")) 
#     
#     # Rasterize regions polygons
#     regions.raster <- raster::rasterize(regions[, "region"], rd)
#     # regions.raster <- raster::rasterize(sp::spTransform(regions[, "region"], sp::proj4string(lines.r)), lines.r)
#     
#     # Sum of line numbers in U.S. waters
#     US.waters <- sum(lines.r[], na.rm = TRUE)
#     
#     # US.waters <- 
#     # sum(lines.r[which(regions.raster[] %in% which(regions$region %in% c("CCB", "SNE", "MIDA", "SEUS", "GOM")))], na.rm = TRUE)
#     
#     # Apply scalars
#     GSL.waters <- GSL_scalar * US.waters
#     ELSWC.waters <- ELSWC_scalar * US.waters
#     
#     # Identify cells within regions of interest
#     GSL.cells <- which(regions.raster[] == which(regions$region == "GSL"))
#     ELSWC.cells <- which(regions.raster[] %in% which(regions$region %in% c("CABOT", "SCOS", "BOF_lower", "BOF_upper")))
#     
#     lines.all <- regions.raster
#     lines.all[!is.na(lines.all)] <- 0
#     lines.all <- raster::merge(lines.r, lines.all)
#     lines.all[GSL.cells] <- scalar * (GSL.waters / length(GSL.cells))
#     lines.all[ELSWC.cells] <- scalar * (ELSWC.waters / length(ELSWC.cells))
#     
#     as(lines.all / raster::maxValue(lines.all), "SpatialGridDataFrame")
#     
#   }) |> purrr::set_names(nm = month.abb)
# }

get_entglD <- function(){
  readr::read_csv("data/gear/duration_entanglements_BOEM.csv") |> 
    janitor::clean_names()
}

rzinbinom <- function(n, mu, theta, size, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta <- size
  rval <- rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi] <- 0
  rval
}

entgl_durations <- function(){
  
  # Duration of entanglement events
  entgl <- targets::tar_read(entgl_d)
  severity <- unique(entgl$severity)
  
  # Zero-truncated negative binomial fit - as entanglement duration must be a minimum of 1 day
  fnb <- purrr::set_names(severity) |> 
    purrr::map(.f = ~{
      tab <- table(entgl[entgl$severity == .x, "duration_days"])
      mat <- data.frame(j = as.numeric(names(tab)), n_j = as.numeric(tab))
      unlist(suppressWarnings(preseqR::preseqR.ztnb.em(mat))[1:2])
    })
  
  par(mfrow = c(3, 2))
  purrr::walk(
    .x = severity,
    .f = ~ {
      hist(entgl[entgl$severity == .x, ]$duration_days,
           freq = FALSE,
           main = .x,
           breaks = 30,
           cex.main = 1.2,
           cex.axis = 1.1,
           cex.lab = 1.1,
           xlab = "Entanglement duration (days)")
      x <- seq(0, max(entgl[entgl$severity == .x, ]$duration_days), by = 1)
      lines(x, countreg::dztnbinom(x, size = fnb[[.x]]["size"], mu = fnb[[.x]]["mu"]), lwd = 1.5)
      plot(ecdf(entgl[entgl$severity == .x, ]$duration_days), ylab = "CDF", main = "Empirical and theoretical CDFs")
      lines(countreg::pztnbinom(x, size = fnb[[.x]]["size"], mu = fnb[[.x]]["mu"]), col = "firebrick", lwd = 1.5)
    }
  )
  print(fnb)
  
  cat("% Truncation at 365 days:\n\n")
  purrr::map(.x = fnb, .f = ~round(100*(1-countreg::pztnbinom(365, size = .x["size"], mu = .x["mu"])), 3))
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
      rescale(new.max = 0.4)
  })
  
  out <- purrr::map(.x = out, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}


# VESSEL TRAFFIC ------------------------------------------------------

# csv2sql <- function(csv.path, dbName, tblName, delim = ","){
# 
#   current_wd <- getwd()
#   target_dir <- dirname(csv.path)
#   setwd(target_dir)
#   csv <- basename(csv.path)
#   
#   if(!grepl(pattern = ".sqlite", x = dbName)) dbName <- paste0(dbName, ".sqlite")
#   
#   inborutils::csv_to_sqlite(csv_file = csv.path, 
#                             sqlite_file = dbName, 
#                             table_name = tblName, 
#                             pre_process_size = 1000, 
#                             chunk_size = 50000, 
#                             delim = delim,
#                             show_progress_bar = TRUE)
#   
#   setwd(current_wd)
#   
# }

get_transits <- function(){
  vt <- readr::read_csv("data/vessels/vessel_transits.csv", col_types = c("fffiiiid")) |>
    data.table::as.data.table()
  vt[, routeID:=paste0(windfarm, "-", route)]
  return(vt)
}

get_routes <- function(){
  ro <- purrr::map(.x = 1:3, .f = ~{
    out <- raster::shapefile(paste0("/Users/philbouchet/OneDrive - University of St Andrews/BOEM/project/gis/BOEM_Renewable_Energy_Areas/WindFarm", .x, "_Route.shp")) |> 
      sp::spTransform(CRSobj = narw_crs())
    out$windfarm <- .x
    out$site <- ifelse(out$windfarm < 3, "Southern New England", "Virginia")
    out
  }) |> do.call(what = rbind)
  return(ro)
}

shp2ais <- function(shp){
  
  # Input shapefile
  shp <- raster::shapefile("/Users/philbouchet/OneDrive - University of St Andrews/BOEM/project/gis/BOEM_Renewable_Energy_Areas/WindFarm1_Route.shp")
  
  # Checks for spatial overlap
  shp_sf <- sf::st_as_sf(shp) |> sf::st_transform(crs = narw_crs())
  inside <- suppressWarnings(sf::st_intersection(vessel_grid$sf, shp_sf))
  
  speed <- 10
  b <- 11
  lg <- 40
  
  # If all vessel positions are outside of the vessel_grid, then can safely discard
  if(nrow(inside) > 0){
    
    if(nrow(shp_sf) > 1){
      
      # Intersect with grid
      ais_out <- suppressWarnings(sf::st_intersection(shp_sf, vessel_grid$sf))
      ais_out$cellID <- ais_out$area_km2 <- NULL
      
      # Join with underlying grid
      # Use the terra package as it makes things much quicker!
      ais_vect <- terra::vect(ais_out)
      ais_grid <- rbind(
        terra::disagg(terra::intersect(ais_vect, vessel_grid_vect)),
        terra::erase(ais_vect, vessel_grid_vect)
      ) |> sf::st_as_sf()
      
      # Update distances after grid join
      ais_grid$dist_km <- as.numeric(sf::st_length(ais_grid))
      ais_grid <- ais_grid[!is.na(ais_grid$cellID),]
      
      # Add speed, beam and length
      ais_grid$speed <- speed
      ais_grid$beam_m <- b
      ais_grid$length_m <- lg
      
      }
    }
  }

create_vesselgrid <- function(){
  
  # Load a density surface and rescale to 10 km resolution
  # Then convert to polygons
  d <- targets::tar_read(density_narw)
  world <- targets::tar_read(world)
  regions <- targets::tar_read(regions)

  vessel_grid_r <- raster::raster(d$Jan) |> 
    raster::aggregate(fact = 2)
  names(vessel_grid_r) <- "risk"
  vessel_grid_r[!is.na(vessel_grid_r)] <- 0
  
  # Convert to spatial polygons
  vessel_grid_sf <- raster::rasterToPolygons(vessel_grid_r)
  vessel_grid_sf$cellID <- seq_len(nrow(vessel_grid_sf))
  vessel_grid_sf$risk <- NULL
  
  # Erase land and calculate cell areas
  vessel_grid_sf_noland <- vessel_grid_sf |> 
    raster::erase(y = world) |> 
    sf::st_as_sf()

  # Merge back
  vessel_grid_sf$area_km2 <- as.numeric(sf::st_area(vessel_grid_sf_noland)) 
  vessel_grid_sf <- sf::st_as_sf(vessel_grid_sf)
  
  # Add region and country by performing a spatial join on cell centroids
  centroids <- sf::st_join(sf::st_centroid(vessel_grid_sf), sf::st_as_sf(regions), join = sf::st_nearest_feature)
  vessel_grid_sf <- dplyr::left_join(vessel_grid_sf, sf::st_drop_geometry(centroids)[, c(1,3,4)], by = "cellID")
  
  # Re-create the raster
  r.out <- raster::stack(
    fasterize::fasterize(
      sf = vessel_grid_sf[, "cellID"],
      raster = vessel_grid_r,
      field = "cellID"),
    fasterize::fasterize(
    sf = vessel_grid_sf[, "area_km2"],
    raster = vessel_grid_r,
    field = "area_km2"))
  names(r.out) <- c("cellID", "area_km2")
  
  return(list(raster = r.out, sf = vessel_grid_sf))

}

proxy_vessels <- function(pmax = 0.005){
  
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
                             out.r <- rescale(out.r, new.max = pmax)
                             out.r[is.na(out.r)] <- 0
                             raster::mask(out.r, mask = poly)})
  
  out <- purrr::map(.x = vessels.r, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}

# PREY ------------------------------------------------------
  
get_glorys <- function(){
  
  # Project regions shapefile to lat/lon coordinate system
  support <- targets::tar_read(support_poly)
  support.ll <- sp::spTransform(support, CRSobj = narw_crs(TRUE))
  
  # Load GLORYS data, then crop to spatial support
  # Source: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/files?subdataset=cmems_mod_glo_phy_my_0.083deg_static_202311--ext--coords
  glorys <- suppressWarnings(raster::stack("data/calanus/GLO-MFC_001_030_coordinates.nc")) |> 
    raster::crop(y = support.ll) |>
    raster::mask(mask = support.ll)
  
  # Extract grid coordinates
  glorys <- raster::as.data.frame(glorys, xy = TRUE, na.rm = TRUE) |> 
    dplyr::rename(longitude = x, latitude = y)
  glorys <- glorys[, c("longitude", "latitude")]
  
  # Convert to <data.table>
  glorys <- data.table::as.data.table(glorys)
  
  return(glorys)
  
}

redraw_prey <- function(){
  
  prey_layer <- targets::tar_read(prey_layer)
  
  # Maximum values of each raster
  rval <- purrr::map(.x = prey_layer, .f = ~raster::raster(.x) |> raster::getValues()) |> do.call(what = c)
  rval <- rval[!is.na(rval)]
  
  purrr::walk(.x = names(prey_layer), .f = ~{

    calmap <- plot_raster(prey_layer[[.x]],
    breaks = c(0, 7.5e-06, 5e-03, 0.015, 0.03, 0.05, 0.065, 0.0875, 0.125, 0.25, 0.5, 0.75, 1, 2.16),
    title = month.name[which(month.abb == .x)])
    
    ggsave(plot = calmap, 
           filename = paste0("out/calanus/calanus_month", 
           stringr::str_pad(which(month.abb == .x), 2, pad = "0"), ".jpeg"),
           dpi = 300, width = 1800, height = 2250, units = "px")
  })
  
}


ipdw_raster <- function (sf_ob, 
                         costras,
                         range, 
                         paramlist, 
                         overlapped = FALSE, 
                         yearmon = "default", removefile = TRUE, step = 16, dist_power = 1, 
                         trim_rstack = FALSE) 
{
  
  pathdists <- ipdw::pathdistGen(sf_ob, costras, range, step)
  pathdists@crs <- raster::crs(costras)
  
  final.ipdw <- ipdwInterp_raster(sf_ob, 
                                  pathdists, 
                                  paramlist, 
                                  yearmon, 
                                  removefile = TRUE, 
                                  overlapped = overlapped, 
                                  dist_power = dist_power, 
                                  trim_rstack = trim_rstack)
  
  return(final.ipdw)
}

ipdwInterp_raster <- function (sf_ob, rstack, paramlist, overlapped = FALSE, yearmon = "default", 
                               removefile = TRUE, dist_power = 1, trim_rstack = FALSE) 
{
  if (missing(paramlist)) {
    stop("Must pass a specific column name to the paramlist argument.")
  }
  if (any(!(paramlist %in% names(sf_ob)))) {
    stop(paste0("Variable(s) '", paste0(paramlist[!(paramlist %in% 
                                                      names(sf_ob))], collapse = "', '"), "' does not exist in sf_ob object."))
  }
  range <- slot(rstack, "range")
  
  if (trim_rstack) {
    rstack <- raster::mask(rstack, sf::st_convex_hull(sf_ob), 
                           inverse = FALSE)
  }
  
  
  for (k in seq_len(length(paramlist))) {
    
    points_layers <- ipdw::rm_na_pointslayers(param_name = paramlist[k], sf_ob = sf_ob, rstack = rstack)
    
    sf_ob <- points_layers$sf_ob
    rstack <- points_layers$rstack
    
    rstack.sum <- raster::calc(rstack, fun = function(x) {
      sum(x^dist_power, na.rm = TRUE)
    })
    
    rstack.sum <- raster::reclassify(rstack.sum, cbind(0, NA))
    
    pb <- progress::progress_bar$new(total = dim(rstack)[3])
    for (i in 1:dim(rstack)[3]) {
      pb$tick()
      ras.weight <- rstack[[i]]^dist_power/rstack.sum
      param.value <- data.frame(sf_ob[i, paramlist[k]])
      param.value2 <- as.vector(unlist(param.value[1]))
      ras.mult <- ras.weight * param.value2
      rf <- raster::writeRaster(ras.mult, 
                                filename = file.path(tempdir(),paste(paramlist[k], "A5raster", i, ".grd", sep = "")), 
                                overwrite = TRUE)
    }
    
    
    raster_data_full <- list.files(path = file.path(tempdir()), 
                                   pattern = paste(paramlist[k], "A5raster*", sep = ""), 
                                   full.names = TRUE)
    
    raster_data <- raster_data_full[grep(".grd", raster_data_full, fixed = TRUE)]
    
    fileNum <- as.numeric(gsub(".*A5raster([0123456789]*)\\.grd$", "\\1", raster_data))
    
    raster_data <- raster_data[order(fileNum)]
    
    rstack.mult <- raster::stack(raster_data)
    
    finalraster <- raster::calc(rstack.mult, fun = function(x) {
      sum(x, na.rm = TRUE)
    })
    
    if (overlapped == TRUE) {
      finalraster <- raster::reclassify(finalraster, cbind(0, 
                                                           NA))
    }
    
    r <- raster::rasterize(sf_ob, rstack[[1]], paramlist[k])
    
    finalraster <- raster::cover(r, finalraster)
    # finalraster <- new("ipdwResult", finalraster, range = range, dist_power = dist_power)
    
    
    file.remove(raster_data_full)
    return(finalraster)
  }
  
  
}


#' Prey surfaces
#'
#' Generates monthly prey abundance rasters from the DFO data
#'

get_prey <- function(months = 1:12, do.plot = TRUE){
  
  # load("/Volumes/GoogleDrive/My Drive/Documents/git/narwinddev/R/sysdata.rda")
  
  support_poly <- targets::tar_read(support_poly)
  prey_grid <- targets::tar_read(prey_grid)
  world <- targets::tar_read(world)
  
  # List input files
  prey.files <- list.files(file.path("data/calanus"), pattern = ".rds", full.names = TRUE)
  
  # Import prey data
  prey_layer <- purrr::map(.x = months, .f = ~{
    
    # File to import
    prey.file <- prey.files[grepl(pattern = paste0("month", stringr::str_pad(.x, 2, pad = "0")), x = prey.files)]
    
    # Read in the data, clean column names and convert to data.table
    calanus <- readRDS(prey.file) |> 
      janitor::clean_names() |> 
      data.table::as.data.table() |> 
      dplyr::rename(longitude = x, latitude = y)
    
    # Only retain cells within extent of spatial support for simulation
    # This prevents erroneous NA values when joining with Calanus data
    calanus <- add_xy(calanus)
    calanus$label <- as.numeric(calanus$label)
    
    calanus.sf <- sf::st_as_sf(unique(calanus[, c("label", "x", "y")]), coords = c("x", "y"), crs = narw_crs())
    calanus.sf <- sf::st_intersection(x = calanus.sf, y = sf::st_as_sf(support_poly))
    calanus <- calanus[label %in% calanus.sf$label]
    
    # Remove unnecessary columns
    calanus <- calanus[, .SD, .SDcols = names(calanus)[!grepl(pattern = "m2", x = names(calanus))]]
    calanus <- calanus[, .SD, .SDcols = names(calanus)[!grepl(pattern = "mg", x = names(calanus))]]
    calanus <- calanus[, .SD, .SDcols = names(calanus)[!grepl(pattern = "cfin_glac", x = names(calanus))]]
    
    # Need to perform a spatial join with the GLORYS grid as there are 
    # inconsistencies with floating point precision in the coordinates reported 
    # in the Calanus data, which prevent a traditional <left_join> in the first instance
    # However, the use of <st_join> results in data loss, so a <left_join> is needed to
    # retrieve the missing information
    calanusxy <- unique(calanus[, list(label, longitude, latitude)])
    latlon <- sf::st_join(x = sf::st_as_sf(x = prey_grid, 
                                           coords = c("longitude", "latitude"), 
                                           crs = narw_crs(latlon = TRUE)),
                          y = sf::st_as_sf(x = calanusxy,
                                           coords = c("longitude", "latitude"),
                                           crs = narw_crs(latlon = TRUE)),
                          join = nngeo::st_nn, 
                          maxdist = 5700,
                          k = 1, 
                          progress = TRUE)

    # Drop geometry and add coordinates to data
    latlon.df <- cbind(sf::st_drop_geometry(latlon), sf::st_coordinates(latlon)) |> 
      dplyr::rename(longitude = X, latitude = Y) |> 
      data.table::as.data.table()
    
    # Add new labels to grid cells outside original SDM prediction space
    latlon.df <- latlon.df[order(label)]
    latlon.df$label <- as.numeric(latlon.df$label)
    latlon.df[is.na(label)]$label <- seq(max(latlon.df$label, na.rm = TRUE) + 1, 
                                      by = 1, length.out = nrow(latlon.df[is.na(label)]))
    
    # Join data and add projected coordinates
    calanus <- calanus[, .SD, .SDcols = !names(calanus) %in% c("longitude", "latitude", "x", "y")]
    calanus <- calanus[latlon.df, on = "label"]
    calanus <- add_xy(calanus)

    # Make sure that data are in correct order
    calanus <- calanus[order(label, zlayer)]
    
    # Round up zlayer values to match concentration calculations from spreadsheet
    calanus[, zlayer:=round(zlayer,0)]
    
    # Determine shallow vs deep
    # calanus[, deep:= ifelse(zlayer <= 100, 0, 1)]
    # calanus$deep <- factor(calanus$deep, labels = c("0–100 m", ">100 m"))
    
    # Sum across species
    calanus[, gm3_tot := rowSums(.SD), .SDcols = names(calanus)[grepl(pattern = "gm3", x = names(calanus))]]
    
    # Calculate thickness of each depth layer
    # calanus[, zdiff:=dplyr::lag(zlayer), label]
    # calanus[, zdiff:=ifelse(is.na(zdiff), zlayer, zlayer - zdiff), label]
    
    # Depth-integrated prey concentration
    # calanus[, gm3_avg:= weighted.mean(gm3_tot, zdiff), list(label, deep)]
    
    # Maximum concentration
    # calanus[, gm3_max:= max(gm3_tot), list(label, deep)]
    
    # Maximum depth that right whales can dive to
    # Baumgartner et al. (2005)
    max.depth <- 306
    
    # Maximum concentration within diving range
    calanus[zlayer <= max.depth, gm3_max_dive:= max(gm3_tot), label]
    
    # Aggregate down to one value per grid cell
    calanus.out <- calanus[, list(label = unique(label), 
                   longitude = unique(longitude),
                   latitude = unique(latitude),
                   gm3_max_dive = unique(gm3_max_dive[!is.na(gm3_max_dive)])), 
                   label]
    
    # Gap-fill g/m3 values
    calanus.out[, gm3_max_dive:= ifelse(is.na(gm3_max_dive), 0, gm3_max_dive)]
    
    # There should only be one set of x,y coordinates per label
    check.labels <- calanus.out[, list(n = length(unique(paste0(longitude, "-", latitude)))), label]
    if(any(check.labels[!is.na(label),]$n > 1)) stop("More than 1 data point per grid cell")
   
    ### // ==== RASTERS ====
    
    calanusmax <-
      raster::rasterFromXYZ(
        xyz = unique(calanus.out[, list(longitude, latitude, gm3_max_dive)]),
        res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
        crs = narw_crs(TRUE)
      )
    
    #' ------------------------------------------
    # Remove grid cells where extrapolation occurred and do some gap-filling/interpolation
    
    northumberland <- raster::shapefile("data/calanus/Northumberland_BayFundy_remove.shp")
    
    # The R package <ipdw> provides functions for interpolation of geo-referenced point data
    # via Inverse Path Distance Weighting. It is useful for coastal marine applications
    # where barriers in the landscape preclude interpolation with Euclidean distances.
    # See https://cran.r-project.org/web/packages/ipdw/vignettes/ipdw2.html
    
    world_sf <- sf::st_as_sf(world) |> sf::st_transform(crs = narw_crs(TRUE))
    
    # Begin by converting the calanus raster to data.frame
    calanusmax_pts <- raster::as.data.frame(calanusmax, xy = TRUE) |> 
      sf::st_as_sf(coords = c("x", "y"), crs = narw_crs(TRUE)) |> 
      dplyr::mutate(cellID = dplyr::row_number()) |> 
      dplyr::filter(!is.na(gm3_max_dive))
    
    # **************************************
    # (1) Northumberland Strait 
    # **************************************
    
    ex_north <- raster::extent(c(-65, -61, 45.35, 47.5))
    
    calanusmax_gsl <- calanusmax_pts |> sf::st_crop(ex_north)

    # Next, use the world_sf object to create a cost raster defining travel
    # through land areas with a very high cost. As a result, interpolation
    # neighborhoods will be defined based on in-water rather than Euclidean distances.
    # Cost raster creation is accomplished with the ipdw function costrasterGen
    costras <- ipdw::costrasterGen(xymat = calanusmax_gsl,
                                   pols = world_sf,
                                   extent = "points", 
                                   projstr = narw_crs(TRUE),
                                   resolution = raster::res(calanusmax)[1]) |> 
      raster::projectRaster(to = calanusmax, method = "ngb")

    
    # Identify grid cells that fall within the Northumberland polygon
    calanus_north <- sp::over(sf::as_Spatial(calanusmax_gsl), northumberland, returnList = TRUE)
    calanus_north <- purrr::map_dbl(.x = calanus_north, .f = ~nrow(.x)) |> unname()
    
    # Cells to use for IPDW
    calanus_ipdw <- calanusmax_gsl[which(calanus_north == 0),]
    
    # Values to be replaced
    calanus_north <- calanusmax_gsl[which(calanus_north == 1),]

    # Find average nearest neighbour
    W <- spatstat.geom::owin(range(c(sf::st_bbox(calanus_ipdw)["xmin"],
                                     sf::st_bbox(calanus_ipdw)["xmax"])),
                           range(c(sf::st_bbox(calanus_ipdw)["ymin"],
                                   sf::st_bbox(calanus_ipdw)["ymax"])))

    kat.pp <- spatstat.geom::ppp(sf::st_coordinates(calanus_ipdw)[,1],
                                 sf::st_coordinates(calanus_ipdw)[,2], window = W)

    mean.neighdist <- mean(spatstat.geom::nndist(kat.pp))

    # ipdw functions rely on rasters stored in the temp directory
    # If these are not cleared, issues may arise
    tmpfiles.to.remove <- list.files(tempdir(), full.names = T, pattern = "gm3_max_dive")
    for(ff in tmpfiles.to.remove) file.remove(ff)
    
    # Perform inverse-distance weighted interpolation
    final.ipdw <- ipdw_raster(sf_ob = calanus_ipdw, 
                             costras = costras,
                             range = mean.neighdist * 10,
                             paramlist = "gm3_max_dive",
                             overlapped = TRUE)
    
    # Retrieve interpolated values
    calanus_north$gm3_max_dive <- NULL
    calanus_north$gm3_max_dive <- raster::extract(final.ipdw, calanus_north)
    
    # Replace values
    # plot(calanusmax, xlim = ex_north[1:2], ylim = ex_north[3:4], zlim = c(0,0.8))
    calanusmax[calanus_north$cellID] <- calanus_north$gm3_max_dive
    # plot(calanusmax, xlim = ex_north[1:2], ylim = ex_north[3:4], zlim = c(0,0.8))
    
    # **************************************
    # (2) Bay of Fundy
    # **************************************
    
    ex_fundy <- raster::extent(c(-68, -63, 44, 46))
    
    # Begin by converting the calanus raster to data.frame
    calanusmax_bof <- calanusmax_pts |> sf::st_crop(ex_fundy)
    
    # Next, use the world_sf object to create a cost raster defining travel
    # through land areas with a very high cost. As a result, interpolation
    # neighborhoods will be defined based on in-water rather than Euclidean distances.
    # Cost raster creation is accomplished with the ipdw function costrasterGen
    costras <- ipdw::costrasterGen(xymat = calanusmax_bof,
                                   pols = world_sf,
                                   extent = "points", 
                                   projstr = narw_crs(TRUE),
                                   resolution = raster::res(calanusmax)[1]) |> 
      raster::projectRaster(to = calanusmax, method = "ngb")
    
    
    # Identify grid cells that fall within the Northumberland polygon
    calanus_fundy <- sp::over(sf::as_Spatial(calanusmax_bof), northumberland, returnList = TRUE)
    calanus_fundy <- purrr::map_dbl(.x = calanus_fundy, .f = ~nrow(.x)) |> unname()
    
    # Cells to use for IPDW
    calanus_ipdw <- calanusmax_bof[which(calanus_fundy == 0),]
    
    # Values to be replaced
    calanus_fundy <- calanusmax_bof[which(calanus_fundy == 1),]
    
    # Find average nearest neighbour
    W <- spatstat.geom::owin(range(c(sf::st_bbox(calanus_ipdw)["xmin"],
                                     sf::st_bbox(calanus_ipdw)["xmax"])),
                             range(c(sf::st_bbox(calanus_ipdw)["ymin"],
                                     sf::st_bbox(calanus_ipdw)["ymax"])))
    
    kat.pp <- spatstat.geom::ppp(sf::st_coordinates(calanus_ipdw)[,1],
                                 sf::st_coordinates(calanus_ipdw)[,2], window = W)
    
    mean.neighdist <- mean(spatstat.geom::nndist(kat.pp))
    
    # ipdw functions rely on rasters stored in the temp directory
    # If these are not cleared, issues may arise
    tmpfiles.to.remove <- list.files(tempdir(), full.names = T, pattern = "gm3_max_dive")
    for(ff in tmpfiles.to.remove) file.remove(ff)
    
    # Perform inverse-distance weighted interpolation
    final.ipdw <- ipdw_raster(sf_ob = calanus_ipdw, 
                              costras = costras,
                              range = mean.neighdist * 10,
                              paramlist = "gm3_max_dive",
                              overlapped = TRUE)
    
    # Retrieve interpolated values
    calanus_fundy$gm3_max_dive <- NULL
    calanus_fundy$gm3_max_dive <- raster::extract(final.ipdw, calanus_fundy)
    
    # Replace values
    # plot(calanusmax, xlim = ex_fundy[1:2], ylim = ex_fundy[3:4], zlim = c(0,2.5))
    calanusmax[calanus_fundy$cellID] <- calanus_fundy$gm3_max_dive
    # plot(calanusmax, xlim = ex_fundy[1:2], ylim = ex_fundy[3:4], zlim = c(0,2.5))
    
    raster::projectRaster(from = calanusmax, crs = narw_crs())
    
    }) |> purrr::set_names(nm = month.abb[months])
  
  if (do.plot) {
    
    # Maximum values of each raster
    rval <- purrr::map(.x = prey_layer, .f = ~raster::getValues(.x)) |> do.call(what = c)
    rval <- rval[!is.na(rval)]
    
    purrr::walk(.x = names(prey_layer), .f = ~{
      
      calmap <- plot_raster(prey_layer[[.x]], breaks = c(0, 7.5e-06, 5e-03, 0.015, 0.03, 0.05, 0.065, 0.0875, 0.125, 0.25, 0.5, 0.75, 1, 2.16))
      # calmap <- plot_raster(prey_layer[[.x]], breaks = quantile(rval, seq(0,1,by = 0.1)))
      # breaks = pretty(c(0, max(rmax)), n = 10))
      
      ggsave(plot = calmap, filename = paste0("out/calanus_month", 
             stringr::str_pad(which(month.abb == .x), 2, pad = "0"), "_maxdepth.jpeg"),
             dpi = 300, width = 1800, height = 2250, units = "px")
    })
    
    # # Map depth of each grid cell
    # calanus.depth <- calanus[, list(zmax = max(bathymetry), 
    #                                 layer = as.numeric(sum(deep) > 0)), 
    #                          list(label, longitude, latitude)]
    # 
    # rdepth <- raster::rasterFromXYZ(xyz = calanus.depth[, list(longitude, latitude, zmax)],
    #                                 res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #                                 crs = narw_crs(TRUE))
    # 
    # # Map layers 0–100 m vs. 100+
    # rlayer <- raster::rasterFromXYZ(xyz = calanus.depth[, list(longitude, latitude, layer)],
    #                                 res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #                                 crs = narw_crs(TRUE))
    # 
    # # Maximum/average prey concentration
    # calanusr <- 
    #   list(max = raster::rasterFromXYZ(xyz = unique(calanus[, list(longitude, latitude, gm3_max)]), 
    #                                    res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #                                    crs = narw_crs(TRUE)),
    #        avg = raster::rasterFromXYZ(xyz = unique(calanus[, list(longitude, latitude, gm3_avg)]), 
    #                                    res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #                                    crs = narw_crs(TRUE)))
    # 
    # # Projects rasters
    # rdepth <- raster::projectRaster(from = rdepth, crs = narw_crs(), method = "bilinear")
    # rlayer <- raster::projectRaster(from = rlayer, crs = narw_crs(), method = "ngb")
    # calanusr <- purrr::map(.x = calanusr, .f = ~raster::projectRaster(from = .x, crs = narw_crs()))
    # 
    # 
    # ### // ==== PLOTS ====
    # 
    # # Difference between 0–100 m and 100+ m
    # calanusr.diff <- split(calanus, f = calanus$deep) |>
    #   purrr::map(.f = ~ {
    #     list(
    #       max = raster::rasterFromXYZ(
    #         xyz = .x[, list(longitude, latitude, gm3_max)],
    #         res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #         crs = narw_crs(TRUE)
    #       ),
    #       avg = raster::rasterFromXYZ(
    #         xyz = .x[, list(longitude, latitude, gm3_avg)],
    #         res = c(0.0833333, 0.0833333), # Resolution of GLORYS model
    #         crs = narw_crs(TRUE)
    #       )
    #     )
    #   })
    # 
    # calanusr.diff <-
    #   purrr::set_names(c("max", "avg")) |>
    #   purrr::map(.f = ~ calanusr.diff[[2]][[.x]] - raster::crop(calanusr.diff[[1]][[.x]], calanusr.diff[[2]][[.x]]))
    # 
    # calanusr.diff <- purrr::map(.x = calanusr.diff, .f = ~ raster::projectRaster(from = .x, crs = narw_crs()))
    # 
    # world.sf <- sf::st_as_sf(world) |>
    #   sf::st_transform(crs = narw_crs(TRUE))
    # 
    # ybrks <- c(seq(0, 0.9, 0.1), 0.95, 0.99, 0.999, 0.9999, 1)
    # 
    # calanus.df <- calanus[!is.na(gm3_max)] |>
    #   dplyr::mutate(
    #     Ncol.max = cut(gm3_max,
    #                    breaks = quantile(gm3_max, ybrks),
    #                    include.lowest = TRUE
    #     ),
    #     Ncol.avg = cut(gm3_avg,
    #                    breaks = quantile(gm3_avg, ybrks),
    #                    include.lowest = TRUE
    #     )
    #   ) |>
    #   dplyr::mutate(Ncol.max = factor(Ncol.max), Ncol.avg = factor(Ncol.avg))
    # 
    # levels(calanus.df$Ncol.max) <-
    #   gsub(pattern = ",", replacement = " – ", levels(calanus.df$Ncol.max))
    # levels(calanus.df$Ncol.max) <-
    #   gsub(pattern = "\\[|\\]", replacement = "", levels(calanus.df$Ncol.max))
    # levels(calanus.df$Ncol.max) <-
    #   gsub(pattern = "\\(|\\)", replacement = "", levels(calanus.df$Ncol.max))
    # 
    # levels(calanus.df$Ncol.avg) <-
    #   gsub(pattern = ",", replacement = " – ", levels(calanus.df$Ncol.avg))
    # levels(calanus.df$Ncol.avg) <-
    #   gsub(pattern = "\\[|\\]", replacement = "", levels(calanus.df$Ncol.avg))
    # levels(calanus.df$Ncol.avg) <-
    #   gsub(pattern = "\\(|\\)", replacement = "", levels(calanus.df$Ncol.avg))
    # 
    # # MAXIMUM CONCENTRATION
    # 
    # p.max <- ggplot2::ggplot(data = calanus.df) +
    #   ggplot2::geom_raster(aes(x = longitude, y = latitude, fill = Ncol.max)) +
    #   ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    #   ylab("") +
    #   xlab("") +
    #   theme_narw() +
    #   ggplot2::facet_wrap(~deep) +
    #   ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    #   ggplot2::coord_sf(
    #     xlim = range(calanus[!is.na(gm3_max), ]$longitude),
    #     ylim = range(calanus[!is.na(gm3_max), ]$latitude), expand = TRUE
    #   ) +
    #   ggplot2::scale_fill_manual(
    #     values = pals::viridis(length(levels(calanus.df$Ncol.max))),
    #     guide = guide_legend(reverse = TRUE),
    #     na.value = "transparent"
    #   )
    # p.max$labels$fill <- "g/m3"
    # 
    # # DEPTH-INTEGRATED AVERAGE CONCENTRATION
    # 
    # p.avg <- ggplot2::ggplot(data = calanus.df) +
    #   ggplot2::geom_raster(aes(x = longitude, y = latitude, fill = Ncol.avg)) +
    #   ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    #   ylab("") +
    #   xlab("") +
    #   theme_narw() +
    #   ggplot2::facet_wrap(~deep) +
    #   ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    #   ggplot2::coord_sf(
    #     xlim = range(calanus[!is.na(gm3_avg), ]$longitude),
    #     ylim = range(calanus[!is.na(gm3_avg), ]$latitude), expand = TRUE
    #   ) +
    #   ggplot2::scale_fill_manual(
    #     values = pals::viridis(length(levels(calanus.df$Ncol.avg))),
    #     guide = guide_legend(reverse = TRUE),
    #     na.value = "transparent"
    #   )
    # p.avg$labels$fill <- "g/m3"
    # 
    # # DIFFERENCE IN CONCENTRATIONS
    # 
    # # Move max by layer to columns
    # diff.max <- calanus |>
    #   dplyr::filter(!is.na(deep)) |>
    #   tidyr::pivot_wider(names_from = deep, values_from = gm3_max) |>
    #   data.table::as.data.table()
    # 
    # # Retrieve unique value per label
    # diff.max <- diff.max[, list(
    #   `0–100 m` = mean(`0–100 m`, na.rm = TRUE),
    #   `>100 m` = mean(`>100 m`, na.rm = TRUE)
    # ), label]
    # 
    # # Move avg by layer to columns
    # diff.avg <- calanus |>
    #   dplyr::filter(!is.na(deep)) |>
    #   tidyr::pivot_wider(names_from = deep, values_from = gm3_avg) |>
    #   data.table::as.data.table()
    # 
    # # Retrieve unique value per label
    # diff.avg <- diff.avg[, list(
    #   `0–100 m` = mean(`0–100 m`, na.rm = TRUE),
    #   `>100 m` = mean(`>100 m`, na.rm = TRUE)
    # ), label]
    # 
    # # Join back to data.table and calculate difference
    # diff.df <- calanus[!is.na(label)]
    # 
    # diff.df <- diff.df[diff.max, on = "label"]
    # diff.df[, diffmax := `>100 m` - `0–100 m`]
    # diff.df <- diff.df[, .SD, .SDcols = !names(diff.df) %in% c(">100 m", "0–100 m")]
    # 
    # diff.df <- diff.df[diff.avg, on = "label"]
    # diff.df[, diffavg := `>100 m` - `0–100 m`]
    # diff.df <- diff.df[, .SD, .SDcols = !names(diff.df) %in% c(">100 m", "0–100 m")]
    # 
    # # Determine which layer has higher concentration
    # diff.df[, diffmax.f := ifelse(diffmax > 0, ">100 m", "0–100 m")]
    # diff.df[, diffavg.f := ifelse(diffavg > 0, ">100 m", "0–100 m")]
    # 
    # pdiff.max <- ggplot2::ggplot(data = diff.df[!is.na(diffmax)]) +
    #   ggplot2::geom_tile(aes(x = longitude, y = latitude, fill = diffmax.f)) +
    #   ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    #   ylab("") +
    #   xlab("") +
    #   theme_narw() +
    #   ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    #   ggplot2::coord_sf(
    #     xlim = range(calanus[!is.na(gm3_max), ]$longitude),
    #     ylim = range(calanus[!is.na(gm3_max), ]$latitude), expand = TRUE
    #   ) +
    #   ggplot2::scale_fill_manual(
    #     name = "Higher max g/m3",
    #     values = c("royalblue4", "lightblue"),
    #     guide = guide_legend(reverse = TRUE),
    #     na.value = "transparent"
    #   )
    # 
    # pdiff.avg <- ggplot2::ggplot(data = diff.df[!is.na(diffavg)]) +
    #   ggplot2::geom_tile(aes(x = longitude, y = latitude, fill = diffavg.f)) +
    #   ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    #   ylab("") +
    #   xlab("") +
    #   theme_narw() +
    #   ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    #   ggplot2::coord_sf(
    #     xlim = range(calanus[!is.na(gm3_avg), ]$longitude),
    #     ylim = range(calanus[!is.na(gm3_avg), ]$latitude), expand = TRUE
    #   ) +
    #   ggplot2::scale_fill_manual(
    #     name = "Higher avg g/m3",
    #     values = c("royalblue4", "lightblue"),
    #     guide = guide_legend(reverse = TRUE),
    #     na.value = "transparent"
    #   )
    # 
    # # Assemble plots
    # p1 <- p.max + pdiff.max + patchwork::plot_layout(guides = "collect")
    # ggplot2::ggsave(p1, filename = "out/Calanus_August_MAXgm3.jpeg", dpi = 300)
    # p2 <- p.avg + pdiff.avg + patchwork::plot_layout(guides = "collect")
    # ggplot2::ggsave(p2, filename = "out/Calanus_August_AVGgm3.jpeg", dpi = 300)
  }
  
  # Convert to spgdf
  prey_layer <- purrr::map(.x = prey_layer, .f = ~as(.x, "SpatialGridDataFrame"))
  
  return(prey_layer)
  
}

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
  gradient.r$val <- 1/(1+exp(-0.025*(gradient.r$x-150)))
  gradient.r$val[which(is.na(gradient.r$Nhat))] <- NA
  gradient.r <- raster::rasterFromXYZ(xyz = gradient.r[,c("x", "y", "val")], res = raster::res(target.r), crs = narw_crs())
  # plot(world, col = "grey")
  # plot(gradient.r, add = T)
  
  # Create a correlated random field as a dummy prey surface
  grid <- sp::makegrid(support_poly, cellsize = 5) |>
    dplyr::rename(x = x1, y = x2)
  area.ex <- raster::extent(support_poly)
  grid <- expand.grid(seq(area.ex[1], area.ex[2], length.out = 200),
                      seq(area.ex[3]-500, area.ex[4], length.out = 200))
  names(grid) <- c("x", "y")
  
  out <- purrr::map(.x = 1:12, .f = ~{

    seed <- as.numeric(paste0("250", .x))
    set.seed(seed)
    
    xyz <- grid
    e <- rnorm(nrow(xyz), mean = 0, sd=100)
    xyz$z <- rlnorm(n = nrow(xyz), meanlog = log(-1*xyz$x + 1.2*xyz$y), sdlog = log(10)) + e
    xyz$z <- ifelse(xyz$z > 20000, 20000, xyz$z)
    xyz$z <- ifelse(is.na(xyz$z), 0, xyz$z)
    yy <- xyz[, c("x", "y", "z")]
    sp::gridded(yy) = ~x+y
    yy <- raster::raster(yy)
    
    # Clip to area of interest
    rout <- raster::resample(yy, target.r) |> 
      raster::crop(target.r) |> 
      raster::mask(target.r)
    
    # Add prey hotspots in foraging grounds of GOM and GSL
    
    N <- 5e5 # Number of random samples

    # Target parameters for univariate normal distributions
    rho <- 0
    mu1.gom <- 800; s1.gom <- 250
    mu2.gom <- 1000; s2.gom <- 250
    
    mu1.gsl <- 1100; s1.gsl <- 250
    mu2.gsl <- 1650; s2.gsl <- 250
    
    max.N <- c(1000, 1000, 1000, 5000, 5000, 15000, 15000, 15000, 15000, 10000, 5000, 1000)
    
    # Parameters for bivariate normal distribution
    mu.gom <- c(mu1.gom, mu2.gom) # Mean
    sigma.gom <- matrix(c(s1.gom^2, s1.gom*s2.gom*rho, s1.gom*s2.gom*rho, s2.gom^2), 2) # Covariance matrix
    
    mu.gsl <- c(mu1.gsl, mu2.gsl) # Mean
    sigma.gsl <- matrix(c(s1.gsl^2, s1.gsl*s2.gsl*rho, s1.gsl*s2.gsl*rho, s2.gsl^2), 2) # Covariance matrix
    
    bvn.gom <- mvtnorm::rmvnorm(N, mu.gom, sigma.gom, method = "svd")
    colnames(bvn.gom) <- c("x","y")
    
    bvn.gsl <- mvtnorm::rmvnorm(N, mu.gsl, sigma.gsl, method = "svd")
    colnames(bvn.gsl) <- c("x","y")
    
    nc.gom <- raster::rasterize(bvn.gom, raster::raster(density_narw$Jan), fun = 'count', background = NA)
    nc.gom <- rescale(nc.gom, new.min = 200, new.max = max.N[.x])
    nc.gom <- raster::mask(nc.gom, support_poly)
    
    nc.gsl <- raster::rasterize(bvn.gsl, raster::raster(density_narw$Jan), fun = 'count', background = NA)
    nc.gsl <- rescale(nc.gsl, new.min = 200, new.max = max.N[.x])
    nc.gsl <- raster::mask(nc.gsl, support_poly)

    nc.gom <- raster::merge(nc.gom, rout) |> 
      raster::mask(mask = regions.support[!regions.support$region == "GSL",]) 
    
    nc.gsl <- raster::merge(nc.gsl, rout) |> 
      raster::mask(mask = regions.support[regions.support$region == "GSL",]) 
    
    nc.out <- raster::merge(nc.gom, nc.gsl)
    # plot(nc.out)
    # plot(world, add = T)
    
    nc.out[!is.na(nc.out) & nc.out<0] <- 0
    
    nc.out
    
    # g.dummy <- gstat::gstat(formula = 
    #                         z ~ 1+y, 
    #                         locations = ~ x + y, 
    #                         dummy = TRUE, 
    #                         beta = 12, 
    #                         model = gstat::vgm(psill = 10, nugget = 2, range = 50, model = 'Sph'), 
    #                         nmax = 10)
    # 
    # yy <- predict(g.dummy, newdata = grid, nsim = 10)
    # yy$z <- apply(X = yy[, 3:ncol(yy)], MARGIN = 1, FUN = function(x) min(x, na.rm = TRUE))
    # yy <- yy[, c("x", "y", "z")]
    # sp::gridded(yy) = ~x+y
    # yy <- raster::raster(yy)
    # 
    # # Clip to area of interest
    # rout <- raster::resample(yy, target.r) |> 
    #   raster::crop(target.r) |> 
    #   raster::mask(target.r)
    # 
    # # plot(rout)
    # 
    # g.dummy01 <- gstat::gstat(formula = z ~ x + y, 
    #                         locations = ~ x + y, 
    #                         dummy = TRUE, 
    #                         beta = 12000, 
    #                         model = gstat::vgm(psill = 10, nugget = 2, range = 100, model = 'Sph'), 
    #                         nmax = 20)
    # 
    # yy01 <- predict(g.dummy01, newdata = grid, nsim = 4)
    # sp::gridded(yy01) = ~x+y
    # yy01 <- raster::raster(yy01) |> rescale(new.min = 0.5)
    # 
    # # Clip to area of interest
    # rout01 <- raster::resample(yy01, target.r) |> 
    #   raster::crop(target.r) |> 
    #   raster::mask(target.r)
    # 
    # # plot(rout01)
    # 
    # rout * gradient.r * rout01
    # plot_raster(rout * gradient.r * rout01, zero = TRUE)
  })
  
  out <- purrr::map(.x = out, .f = ~as(.x, "SpatialGridDataFrame"))
  names(out) <- month.abb
  
  return(out)
}


# IMPORTANCE SAMPLING -----------------------------------------------------

# included in this repo
# importance sampling based on code from mgcv
gam.imp <- function(b, ns = 10000){
  
  beta <- coef(b)
  Vb <- vcov(b)
  X <- model.matrix(b)
  
  # beta proposals
  bs <- rmvn(ns, beta, Vb)
  
  # log proposal density
  lfp <- dmvn(t(bs), beta, Vb)
  # log density for penalized MLE
  lfp1 <- dmvn(matrix(beta, ncol=1), beta, Vb)
  
  # storage
  wts <- rep(0, ns)
  
  # loglik for penalized MLE
  lpl0 <- lpl(b, beta, X)
  
  for (i in 1:ns) {
    # loglik of proposal...
    lpl1 <- lpl(b, bs[i,], X)
    # store weight
    wts[i] <- exp(lfp1-lfp[i]+lpl1-lpl0)
  }
  list(bs=bs, wts=wts)
}
# fixed gam.mh
## Simple post fit mcmc for mgcv.
## (c) Simon N. Wood (2020)

# mgcv::gam.mh with fix via Len Thomas

gam.mh <- function(b,ns=10000,burn=1000,t.df=40,rw.scale=.25,thin=1) {
  ## generate posterior samples for fitted gam using Metroplois Hastings sampler
  ## alternating fixed proposal and random walk proposal, both based on Gaussian
  ## approximation to posterior...
  if (inherits(b,"bam")) stop("not usable with bam fits")
  beta <- coef(b);Vb <- vcov(b)
  X <- model.matrix(b); burn <- max(0,burn)
  prog <- interactive();iprog <- 0
  di <- floor((ns+burn)/100)
  if (prog) prg <- txtProgressBar(min = 0, max = ns+burn, initial = 0,
                                  char = "=",width = NA, title="Progress", style = 3)
  bp <- rmvt(ns+burn,beta,Vb,df=t.df) ## beta proposals
  bp[1,] <- beta ## Don't change this after density step!!
  lfp <- dmvt(t(bp),beta,Vb,df=t.df) ## log proposal density
  
  rw <- is.finite(rw.scale)&&rw.scale>0
  if (rw) {
    R <- chol(Vb) 
    step <- rmvn(ns+burn,beta*0,Vb*rw.scale) ## random walk steps (mgcv::rmvn)
  }
  u <- runif(ns+burn);us <- runif(ns+burn) ## for acceptance check
  bs <- bp;j <- 1;accept <- rw.accept <- 0
  lpl0 <- lpl(b,bs[1,],X)
  for (i in 2:(ns+burn)) { ## MH loop
    ## first a static proposal...
    lpl1 <- lpl(b,bs[i,],X)
    if (u[i] < exp(lfp[j]-lfp[i]+lpl1-lpl0)) {
      lpl0 <- lpl1;accept <- accept + 1
      j <- i ## row of bs containing last accepted beta
    } else bs[i,] <- bs[i-1,]
    ## now a random walk proposal...
    if (rw) {
      lpl1 <- lpl(b,bs[i,]+step[i,],X)
      if (us[i] < exp(lpl1-lpl0)) { ## accept random walk step
        lpl0 <- lpl1;j <- i
        bs[i,] <- bs[i,] + step[i,]
        rw.accept <- rw.accept+1 
        #lfp[i] <- dmvt(bs[i,],beta,Vb,df=4,R=R) ## have to update static proposal density
        # FIX via LJT 5/10/20
        lfp[i] <- dmvt(bs[i,],beta,Vb,df=t.df,R=R) 
      }
    }  
    if (i==burn) accept <- rw.accept <- 0
    if (prog&&i%%di==0) setTxtProgressBar(prg, i)
  } ## MH loop
  if (burn>0) bs <- bs[-(1:burn),]
  if (thin>1) bs <- bs[seq(1,ns,by=thin),]
  if (prog) {
    setTxtProgressBar(prg, i);cat("\n")
    cat("fixed acceptance = ",accept/ns,"  RW acceptance = ",rw.accept/ns,"\n")
  }  
  list(bs=bs,rw.accept = rw.accept/ns,accept=accept/ns)
} ## gam.mh

# other parts of mcmc.r in mgcv that are not exported

## t-distribution stuff for mgcv.
## (c) Simon N. Wood (2020)

## some useful densities (require mgcv::rmvn)...

rmvt <- function(n,mu,V,df) {
  ## simulate multivariate t variates  
  y <- rmvn(n,mu*0,V)
  v <- rchisq(n,df=df)
  t(mu + t(sqrt(df/v)*y))
}

r.mvt <- function(n,mu,V,df) rmvt(n,mu,V,df)

dmvt <- function(x,mu,V,df,R=NULL) {
  ## multivariate t log density...
  p <- length(mu);
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  k <- - sum(log(diag(R))) - p*log(df*pi)/2 + lgamma((df+p)/2) - lgamma(df/2)
  k - if (is.matrix(z)) (df+p)*log1p(colSums(z^2)/df)/2 else (df+p)*log1p(sum(z^2)/df)/2
}

d.mvt <- function(x,mu,V,df,R=NULL) dmvt(x,mu,V,df,R)

# taken from mgcv source code

## Simple post fit mcmc for mgcv.
## (c) Simon N. Wood (2020)


## some functions to extract important components of joint density from
## fitted gam...

## evaluate penalty for fitted gam, possibly with new beta
# patched to include parapen
bSb <- function(b,beta=coef(b)) {
  bSb <- k <-  0
  sp <- if (is.null(b$full.sp)) b$sp else b$full.sp ## handling linked sp's
  
  
  # the parapen bits
  # need to do something clever with L at some point
  if(!is.null(b$paraPen)){
    for (i in 1:length(b$paraPen$S)) {
      k <- k + 1
      # get indices
      ii <- b$paraPen$off[i]
      ii <- ii:(ii+ncol(b$paraPen$S[[i]])-1)
      # add to penalty
      bSb <- bSb + sp[b$paraPen$full.sp.names[i]]*
        (t(beta[ii])%*%b$paraPen$S[[i]]%*%beta[ii])
    }
  }
  
  
  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + sp[k]*(t(beta[ii])%*%b$smooth[[i]]$S[[j]]%*%beta[ii])
      }
    }
  }
  
  
  bSb
} ## bSb

devg <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## evaluate the deviance of a fitted gam given possibly new coefs, beta
  ## for general families this is simply -2*log.lik
  if (inherits(b$family,"general.family")) {
    -2*b$family$ll(b$y,X,beta,b$prior.weights,b$family,offset=b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y,b$family$linkinv(X%*%beta+b$offset),b$prior.weights))
  }
} ## devg

lpl <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## log joint density for beta, to within uninteresting constants
  -(devg(b,beta,X)/b$sig2+bSb(b,beta)/b$sig2)/2
}


# SPATIAL DATA ------------------------------------------------------

#' Projected coordinate system
#'
#' @return An object of class \code{CRS}.
#' 
narw_crs <- function(latlon = FALSE){
  if(latlon) sp::CRS("+init=epsg:4326") else sp::CRS("+proj=aea +lat_0=34 +lon_0=-78 +lat_1=27.3333333333333 +lat_2=40.6666666666667 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs") 
}

get_limits <- function(spgdf){
  out <- sp::coordinates(spgdf)
  return(c(range(out[,1]), range(out[,2])))
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
  
  world <- targets::tar_read(world)
  
  
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
  regions <- targets::tar_read(regions)
  
  # Filter out Canada
  regions.support <- regions[regions@data$region %in% c("GSL", "SCOS", "CABOT"),] 
  
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

  # Need to add small area of GOM back in
  rastermask <- combined_support
  rastermask[rastermask == 0] <- NA
  rastermask <- raster::rasterToPolygons(x = rastermask, dissolve = TRUE)
  w1 <- raster::intersect(world, regions[regions$region == "GOM",])
  combined_area <- rgeos::gUnion(w1, rastermask)
  combined_area <- smoothr::fill_holes(combined_area, 100000)
  
  combined_area <- raster::rasterize(x = combined_area, y = combined_support)
  combined_support <- raster::merge(combined_area, combined_support)
  
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
    # The latest version (v12) of the Duke NARW density surface model works at 
    # a 5 km resolution; however, for consistency with previous versions, 
    # model predictions are expressed as numbers of individuals per 100 km2 (i.e., 10 x 10 km resolution). 
    # We therefore re-scale the monthly rasters to calculate density per 25 km2 (i.e., 5 x 5 km resolution).
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

# plot(world, xlim = c(500,650), ylim = c(720,850), axes = T)
# for(i in 1:nrow(turbines$scenario_01[turbines$scenario_01$farm == 1,])){
#   points(turbines$scenario_01[turbines$scenario_01$farm == 1,]$x[i],
#          turbines$scenario_01[turbines$scenario_01$farm == 1,]$y[i], pch = 16)
#   Sys.sleep(0.5)
# }


# for(i in 1:30){
#   points(turbines[[1]]$x[i], turbines[[1]]$y[i], pch = 16)
#   Sys.sleep(0.5)
# }
# 
# plot(turbines[[1]]$longitude, turbines[[1]]$latitude, pch = 16)
# for(i in 1:30){
#   points(turbines[[1]]$longitude[i], turbines[[1]]$latitude[i], pch = 16, col = "red")
#   Sys.sleep(0.5)
# }
# 
# plot(turbines[[2]]$longitude, turbines[[2]]$latitude, pch = 16)
# for(i in 1:30){
#   points(turbines[[2]]$longitude[i], turbines[[2]]$latitude[i], pch = 16, col = "red")
#   Sys.sleep(0.5)
# }
# 
# plot(turbines[[3]]$longitude, turbines[[3]]$latitude, pch = 16)
# for(i in 1:30){
#   points(turbines[[3]]$longitude[i], turbines[[3]]$latitude[i], pch = 16, col = "red")
#   Sys.sleep(0.5)
# }


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


# POSTERIOR SAMPLING ------------------------------------------------------

# https://github.com/dill/GAMsampling

# taken from mgcv source code

## Simple post fit mcmc for mgcv.
## (c) Simon N. Wood (2020)


## some functions to extract important components of joint density from
## fitted gam...

## evaluate penalty for fitted gam, possibly with new beta
# patched to include parapen
bSb <- function(b,beta=coef(b)) {
  bSb <- k <-  0
  sp <- if (is.null(b$full.sp)) b$sp else b$full.sp ## handling linked sp's
  
  
  # the parapen bits
  # need to do something clever with L at some point
  if(!is.null(b$paraPen)){
    for (i in 1:length(b$paraPen$S)) {
      k <- k + 1
      # get indices
      ii <- b$paraPen$off[i]
      ii <- ii:(ii+ncol(b$paraPen$S[[i]])-1)
      # add to penalty
      bSb <- bSb + sp[b$paraPen$full.sp.names[i]]*
        (t(beta[ii])%*%b$paraPen$S[[i]]%*%beta[ii])
    }
  }
  
  
  for (i in 1:length(b$smooth)) {
    m <- length(b$smooth[[i]]$S)
    if (m) {
      ii <- b$smooth[[i]]$first.para:b$smooth[[i]]$last.para
      for (j in 1:m) {
        k <- k + 1
        bSb <- bSb + sp[k]*(t(beta[ii])%*%b$smooth[[i]]$S[[j]]%*%beta[ii])
      }
    }
  }
  
  
  bSb
} ## bSb

devg <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## evaluate the deviance of a fitted gam given possibly new coefs, beta
  ## for general families this is simply -2*log.lik
  if (inherits(b$family,"general.family")) {
    -2*b$family$ll(b$y,X,beta,b$prior.weights,b$family,offset=b$offset)$l
  } else { ## exp or extended family
    sum(b$family$dev.resids(b$y,b$family$linkinv(X%*%beta+b$offset),b$prior.weights))
  }
} ## devg

lpl <- function(b,beta=coef(b),X=model.matrix(b)) {
  ## log joint density for beta, to within uninteresting constants
  -(devg(b,beta,X)/b$sig2+bSb(b,beta)/b$sig2)/2
}

## t-distribution stuff for mgcv.
## (c) Simon N. Wood (2020)

## some useful densities (require mgcv::rmvn)...

rmvt <- function(n,mu,V,df) {
  ## simulate multivariate t variates  
  y <- rmvn(n,mu*0,V)
  v <- rchisq(n,df=df)
  t(mu + t(sqrt(df/v)*y))
}

r.mvt <- function(n,mu,V,df) rmvt(n,mu,V,df)

dmvt <- function(x,mu,V,df,R=NULL) {
  ## multivariate t log density...
  p <- length(mu);
  if (is.null(R)) R <- chol(V)
  z <- forwardsolve(t(R),x-mu)
  k <- - sum(log(diag(R))) - p*log(df*pi)/2 + lgamma((df+p)/2) - lgamma(df/2)
  k - if (is.matrix(z)) (df+p)*log1p(colSums(z^2)/df)/2 else (df+p)*log1p(sum(z^2)/df)/2
}

d.mvt <- function(x,mu,V,df,R=NULL) dmvt(x,mu,V,df,R)


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

# ANIMATION ------------------------------------------------------

alignmove <- function (mobj, res = "min", digit = "min", unit = "days", spaceMethod = "greatcircle") 
{
  if (!inherits(mobj, c("Move", "MoveStack"))) out("Argument 'm' must be of class 'Move' or 'MoveStack'.", type = 3)
  m.length <- if (inherits(mobj, "MoveStack")) sapply(move::split(mobj), length) else length(mobj)
  if (any(m.length < 2))  out(paste0("Individual track(s) ", paste0(which(m.length < 2), collapse = ", "), " of 'm' consist(s) of less than 2 locations only. A minimum of 2 locations per indvidual track is required for alignment."), type = 3)
  
  ts <- move::timestamps(mobj)
  time.digits <- unique(as.numeric(format(ts, moveVis:::.convert_units(unit))))
  
  if (all(!c(inherits(res, "numeric"), inherits(res, "character"))))  out("Argument 'res' must either be numeric or one of c('min', 'max', 'mean').", type = 3)
  
  if (res == "min") res <- min(unique(unlist(move::timeLag(mobj, unit))))
  if (res == "max") res <- max(unique(unlist(timeLagtimeLag(mobj, unit))))
  if (res == "mean") res <- round(mean(unique(unlist(timeLagtimeLag(mobj, unit)))))
  
  res <- as.difftime(res, units = unit)
  
  if (all(!c(inherits(digit, "numeric"), inherits(digit, "character")))) out("Argument 'digit' must either be numeric or one of c('min', 'max', 'mean').", type = 3)
  if (digit == "min") digit <- min(time.digits)
  if (digit == "max") digit <- max(time.digits)
  if (digit == "mean") digit <- round(mean(time.digits))
  
  ts.shoulder <- list(min(ts), max(ts))
  
  # set.fun <- list(secs = `second<-`, mins = `minute<-`, hours = `hour<-`, days = `day<-`)
  # set.fun <- set.fun[[match(unit, names(set.fun))]]
  # ts.shoulder <- lapply(ts.shoulder, set.fun, value = digit)
  ts.target <- seq.POSIXt(ts.shoulder[[1]], ts.shoulder[[2]], by = res)
  m.indi <- if (inherits(mobj, "MoveStack")) move::split(mobj) else list(mobj)
  mobj <- move::moveStack(lapply(m.indi, function(x) {
    ts.m <- move::timestamps(x)
    ts.t <- ts.target[ts.target >= min(ts.m) & ts.target <= max(ts.m)]
    move::interpolateTime(x, ts.t, spaceMethod)
  }))
  mobj[, c("x", "y")] <- mobj@coords
  mobj[, "time"] <- move::timestamps(mobj)
  return(mobj)
}

getMap <- function (gg.ext, map_service, map_type, map_res, m.crs) 
{
  if (inherits(gg.ext, "bbox")) gg.ext <- list(gg.ext)
  
  r <- lapply(gg.ext, function(y) {
    
    gg.ext.ll <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(y), crs = sf::st_crs("+init=epsg:4326")))
    
    tg <- slippymath::bbox_to_tile_grid(gg.ext.ll, max_tiles = ceiling(map_res * 20))
    
    images <- apply(tg$tiles, MARGIN = 1, function(x) {
      file <- file.path("data/basemap", map_service, paste0(map_service, "_", map_type, "_", x[1], "_", x[2], ".png"))
      magick::image_write(magick::image_convert(magick::image_read(file), format = "PNG24"), file)
      return(file)
    })
    
    r <- quiet(slippymath::compose_tile_grid(tg, images))
    raster::crop(raster::projectRaster(r, crs = m.crs), raster::extent(y[1], y[3], y[2], y[4]))
    # raster::crop(raster::projectRaster(r, crs = m.crs), raster::extent(y[1], y[3], y[2], y[4]), snap = "moveVis:::out")
  })
  
  if (length(r) > 1) {
    ext.both <- list(east = extent(r$east), west = extent(r$west))
    rg <- c(east = diff(c(ext.both$east@xmin, ext.both$east@xmax)), 
            west = diff(c(ext.both$west@xmin, ext.both$west@xmax)))
    ext.both <- moveVis:::.expand_ext(ext.both, rg)
    extent(r$east) <- ext.both$east
    extent(r$west) <- ext.both$west
    ext.combi <- moveVis:::.combine_ext(ext.both)
    r[[which.min(rg)]] <- extend(r[[which.min(rg)]], ext.combi)
    r[[which.max(rg)]] <- resample(r[[which.max(rg)]], r[[which.min(rg)]])
    r <- list(merge(r[[1]], r[[2]]))
  }
  return(r)
}

addtext <- function (frames, labels, x, y, colour = "white", label.size = NA, size = 3, fill = "black", type = "label", verbose = TRUE, alpha = 1) 
{
  if (inherits(verbose, "logical")) 
    options(moveVis.verbose = verbose)
  if (!inherits(frames, "list")) 
    moveVis:::out("Argument 'frames' needs to be a list of ggplot objects. See frames_spatial()).", 
                  type = 3)
  if (!all(sapply(frames, function(x) inherits(x, "ggplot")))) 
    moveVis:::out("At least one element of argument 'frames' is not a ggplot object.", 
                  type = 3)
  if (!is.character(labels)) 
    moveVis:::out("Argument 'labels' must be of type 'character'.", 
                  type = 3)
  if (!is.character(colour)) 
    moveVis:::out("Argument 'colour' must be of type 'character'.", 
                  type = 3)
  if (!is.numeric(x)) 
    moveVis:::out("Argument 'x' must be of type 'numeric'.", type = 3)
  if (!is.numeric(y)) 
    moveVis:::out("Argument 'y' must be of type 'numeric'.", type = 3)
  if (!is.numeric(size)) 
    moveVis:::out("Argument 'size' must be of type 'numeric'.", type = 3)
  check <- list(labels = labels, x = x, y = y, colour = colour, size = size, fill = fill, alpha = alpha, label.size = label.size)
  data <- sapply(1:length(check), function(i) {
    if (length(check[[i]]) == 1) 
      v <- rep(check[[i]], length(frames))
    else v <- check[[i]]
    if (length(v) != length(frames)) 
      moveVis:::out(paste0("Length of argument ", names(check)[[i]], 
                           " must either be 1 or equal to the length of agrument 'frames'."), type = 3)
    return(v)
  }, simplify = F)
  data.classes <- sapply(data, class)
  data <- as.data.frame(do.call(cbind, data), stringsAsFactors = F)
  for (i in 1:ncol(data)) class(data[, i]) <- data.classes[i]
  data <- split(data, seq(nrow(data)))
  moveVis::add_gg(frames, gg = expr(ggplot2::annotate(type, x = data[[2]], y = data[[3]], label = data[[1]], colour = data[[4]], size = data[[5]], fill = data[[6]], alpha = data[[7]], label.size = data[[8]])), data = data, type = type)
}

addtimestamps <- function (frames, m = NULL, x = NULL, y = NULL, fill = "black", label.size = NA, alpha = 1, size = 3, type = "label", verbose = TRUE) 
{
  if (inherits(verbose, "logical")) 
    options(moveVis.verbose = verbose)
  if (!inherits(frames, "list")) 
    moveVis:::out("Argument 'frames' needs to be a list of ggplot objects. See frames_spatial()).", 
                  type = 3)
  if (!all(sapply(frames, function(x) inherits(x, "ggplot")))) 
    moveVis:::out("At least one element of argument 'frames' is not a ggplot object.", 
                  type = 3)
  if (!is.null(m)) {
    if (all(!c(inherits(m, "MoveStack"), inherits(m, "Move")))) 
      moveVis:::out("Argument 'm' must be of class 'Move' or 'MoveStack', if defined. Do not define 'm', if timestamps should be extracted from the attributes of 'frames' directly.", 
                    type = 3)
    ts <- as.character(sort(unique(move::timestamps(m))))
    ts <- lubridate::as_date(ts)
    ts <- paste0(month.abb[lubridate::month(ts)], "-", stringr::str_pad(lubridate::day(ts), 2, pad = "0"))
    
    if (length(ts) != length(frames)) 
      moveVis:::out("Lengths of unique timestamps of 'm' and 'frames' do not match. if 'm' is defined. Do only use the same move or moveStack object that you have used to create 'frames'. Otherwise, you can leave 'm' undefined so that frame times are extracted from the attributes of 'frames' directly.", 
                    type = 3)
  }
  else {
    ts <- as.character(moveVis::get_frametimes(frames))
  }
  if (is.null(x)) {
    gg.xy <- lapply(ggplot2::ggplot_build(frames[[1]])$data, function(x) cbind.data.frame(x = x$x, 
                                                                                          y = x$y))
    gg.xy <- do.call(rbind, gg.xy[!sapply(gg.xy, is.null)])
    x <- min(gg.xy$x) + ((max(gg.xy$x) - min(gg.xy$x))/2)
    y <- max(gg.xy$y) - ((max(gg.xy$y) - min(gg.xy$y)) * 
                           0.05)
  }
  addtext(frames, ts, x, y, size = size, alpha = alpha, fill = fill, type = type, label.size = label.size)
}


frames_sp <- function (mobj, 
                       r_list = NULL, 
                       r_times = NULL,
                       r_type = "gradient", 
                       fade_raster = FALSE, 
                       crop_raster = TRUE,
                       map_service = "osm", 
                       map_type = "streets", 
                       map_res = 1,
                       map_token = NULL,
                       margin_factor = 1.1,
                       equidistant = NULL, 
                       ext = NULL, 
                       path_size = 3, 
                       path_end = "round", 
                       path_join = "round",
                       path_mitre = 10, 
                       path_arrow = NULL, 
                       path_colours = NA, 
                       path_alpha = 1, path_fade = FALSE, 
                       path_legend = TRUE, 
                       path_legend_title = "Names", 
                       tail_length = 19, 
                       tail_size = 1, 
                       tail_colour = "white", 
                       trace_show = FALSE, 
                       trace_colour = "white", 
                       cross_dateline = FALSE,
                       ..., 
                       verbose = TRUE) 

# r_list = NULL 
# r_times = NULL
# r_type = "gradient" 
# fade_raster = FALSE 
# crop_raster = TRUE
# map_res = 1
# map_token = NULL
# margin_factor = 1.1
# equidistant = NULL 
# ext = NULL 
# path_size = 3 
# path_end = "round" 
# path_join = "round"
# path_mitre = 10 
# path_arrow = NULL 
# path_colours = NA 
# path_alpha = 1 
# path_fade = FALSE 
# path_legend = TRUE 
# path_legend_title = "Names" 
# tail_length = 19 
# tail_size = 1 
# tail_colour = "white" 
# trace_show = FALSE 
# trace_colour = "white" 
# cross_dateline = FALSE
# verbose = TRUE


{
  if (inherits(verbose, "logical")) options(moveVis.verbose = verbose)
  if (all(!c(inherits(mobj, "MoveStack"), inherits(mobj, "Move")))) moveVis:::out("Argument 'm' must be of class 'Move' or 'MoveStack'.", type = 3)
  if (inherits(mobj, "Move")) mobj <- moveStack(mobj)
  if (!is.null(r_list)) {
    if (all(!is.list(r_list), inherits(r_list, "Raster"))) r_list <- list(r_list)
    if (any(!sapply(r_list, compareCRS, y = mobj))) moveVis:::out("Projections of 'm' and 'r_list' differ.", type = 3)
    if (length(unique(sapply(r_list, nlayers))) > 1) moveVis:::out("Number of layers per raster object in list 'r' differ.", type = 3)
    if (!inherits(r_times, "POSIXct")) moveVis:::out("Argument 'r_times' must be of type 'POSIXct' if 'r_list' is defined.", type = 3)
    if (!isTRUE(r_type %in% c("gradient", "discrete", "RGB"))) moveVis:::out("Argument 'r_type' must eihter be 'gradient', 'discrete' or 'RGB'.", type = 3)
    if (!is.logical(fade_raster)) moveVis:::out("Argument 'fade_raster' has to be either TRUE or FALSE.", type = 3)
    if (!is.logical(crop_raster)) moveVis:::out("Argument 'crop_raster' has to be either TRUE or FALSE.", type = 3)
  } else {
    # if (isFALSE(tolower(map_service) %in% names(moveVis::get_maptypes()))) moveVis:::out(paste0("Argument 'map_service' must be ", paste0(names(moveVis::get_maptypes()), collapse = ", ")))
    # if (isFALSE(tolower(map_type) %in% moveVis::get_maptypes(map_service))) moveVis:::out("The defined map type is not supported for the selected service. Use get_maptypes() to get all available map types.", type = 3)
    if (!is.numeric(map_res)) moveVis:::out("Argument 'map_res' must be 'numeric'.", type = 3)
    if (any(map_res < 0, map_res > 1)) moveVis:::out("Argument 'map_res' must be a value between 0 and 1.", type = 3)
    if (all(!inherits(map_token, "character"), map_service == "mapbox")) moveVis:::out("Argument 'map_token' must be defined to access a basemap, if 'r_list' is not defined and 'map_service' is 'mapbox'.", type = 3)
  
#   if (is.null(map_dir)) {
#     map_dir <- paste0(tempdir(), "/moveVis/basemap/")
#     if (!dir.exists(map_dir)) {
#       dir.create(map_dir, recursive = T)
#     }
#   } else {
#     if (!dir.exists(map_dir)) {
#       moveVis:::out("The directory defined with 'map_dir' does not exist.", type = 3)
#     }
#   }
}
  num.args <- c(margin_factor = margin_factor, tail_length = tail_length, 
                tail_size = tail_size, path_size = path_size, path_mitre = path_mitre)
  catch <- sapply(1:length(num.args), function(i) if (!is.numeric(num.args[[i]])) 
    moveVis:::out(paste0("Argument '", names(num.args)[[i]], "' must be of type 'numeric'."), 
                  type = 3))
  char.args <- c(path_end = path_end, path_join = path_join, 
                 path_legend_title = path_legend_title)
  catch <- sapply(1:length(char.args), function(i) if (!is.character(char.args[[i]])) 
    moveVis:::out(paste0("Argument '", names(char.args)[[i]], "' must be of type 'numeric'."), 
                  type = 3))
  if (!is.null(ext)) {
    if (!inherits(ext, "Extent")) 
      moveVis:::out("Argument 'ext' must be of type 'Extent' (see raster::extent), if defined.", 
                    type = 3)
    if (isTRUE(ext < extent(mobj))) 
      moveVis:::out("The frame extent defined using argument 'ext' is smaller than extent(m). Be aware that movements moveVis:::outside of 'ext' will be clipped.", 
                    type = 2)
  }
  if (!is.null(path_arrow)) 
    if (!inherits(path_arrow, "arrow")) 
      moveVis:::out("Argument 'path_arrow' must be of type 'arrrow' (see grid::arrow), if defined.", 
                    type = 3)
  if (is.character(path_colours)) 
    if (length(path_colours) != n.indiv(mobj)) 
      moveVis:::out("Argument 'path_colours' must be of same length as the number of individual tracks of 'm', if defined. Alternatively, use a column 'colour' for individual colouring per coordinate within 'm' (see details of ?frames_spatial).", 
                    type = 3)
  if (!is.logical(path_legend)) 
    moveVis:::out("Argument 'path_legend' must be of type 'logical'.", type = 3)
  if (is.null(equidistant)) 
    if (is.null(ext)) equidistant <- TRUE else equidistant <- FALSE
  if (!is.logical(equidistant)) 
    moveVis:::out("Argument 'equidistant' must be of type 'logical'.", type = 3)
  if (all(as.integer(sf::st_crs(mobj)$epsg) != as.integer(4326), isTRUE(cross_dateline), na.rm = T)) {
    moveVis:::out("Argument 'cross_dateline' is ignored, since the coordinate reference system of 'm' is not geographical (long/lat).", type = 2)
    cross_dateline <- FALSE
  }
  if (all(isTRUE(cross_dateline), !is.null(r_list))) 
    moveVis:::out("Argument 'cross_dateline' only works with default base maps. Arguments 'r_list' and 'r_times' cannot be used, if cross_dateline = TRUE.\nTip: Reproject 'm' to another CRS that better suits the region if you want to use 'r_list' with tracks crossing the dateline.", 
                  type = 3)
  if (isTRUE(cross_dateline)) 
    equidistant <- FALSE
  moveVis:::out("Checking temporal alignment...")
  moveVis:::.time_conform(mobj)
  moveVis:::out("Processing movement data...")
  m.df <- moveVis:::.m2df(mobj, path_colours = path_colours)
  moveVis:::.stats(n.frames = max(m.df$frame))
  gg.ext <- moveVis:::.ext(m.df, m.crs = sf::st_crs(mobj), ext, margin_factor, equidistant, cross_dateline)
  # assign("gg.ext", gg.ext, envir = .GlobalEnv)
  if (isTRUE(cross_dateline)) {
    rg <- c(pos = diff(range(m.df$x[m.df$x >= 0])), neg = diff(range(m.df$x[m.df$x < 
                                                                              0])))
    if (which.max(rg) == 1) {
      m.df$x[m.df$x < 0] <- 180 + m.df$x[m.df$x < 0] + 
        180
    } else {
      m.df$x[m.df$x >= 0] <- -180 + m.df$x[m.df$x >= 0] - 
        180
    }
  }
  if (is.null(r_list)) {
    moveVis:::out("Retrieving and compositing basemap imagery...")
    r_list <- getMap(gg.ext, map_service, map_type, map_res, m.crs = raster::crs(mobj))
    r_type <- "RGB"
  }
  if (isTRUE(cross_dateline)) {
    gg.ext <- .combine_ext(.expand_ext(list(extent(gg.ext$east[[1]], 
                                                   gg.ext$east[[3]], gg.ext$east[[2]], gg.ext$east[[4]]), 
                                            extent(gg.ext$west[[1]], gg.ext$west[[3]], gg.ext$west[[2]], 
                                                   gg.ext$west[[4]])), rg))
    gg.ext <- st_bbox(c(xmin = gg.ext@xmin, xmax = gg.ext@xmax, 
                        ymin = gg.ext@ymin, ymax = gg.ext@ymax), crs = sf::st_crs(mobj))
    m.df$coord <- list(ggplot2::coord_sf(xlim = c(gg.ext$xmin, 
                                                  gg.ext$xmax), ylim = c(gg.ext$ymin, gg.ext$ymax), 
                                         expand = F, clip = "on"))
    m.df$scalex <- list(ggplot2::scale_x_continuous(labels = .x_labels))
    m.df$scaley <- list(ggplot2::scale_y_continuous(labels = .y_labels))
  }
  else {
    m.df$coord <- list(ggplot2::coord_sf(xlim = c(gg.ext$xmin, 
                                                  gg.ext$xmax), ylim = c(gg.ext$ymin, gg.ext$ymax), 
                                         expand = F, crs = sf::st_crs(mobj)$proj4string, datum = sf::st_crs(mobj)$proj4string, 
                                         clip = "on"))
    m.df$scaley <- m.df$scalex <- NULL
  }
  moveVis:::out("Assigning raster maps to frames...")
  r_list <- moveVis:::.rFrames(r_list, r_times = NULL, m.df, gg.ext, fade_raster = FALSE, crop_raster = crop_raster)
  moveVis:::out("Creating frames...")
  frames <- moveVis:::.gg_spatial(r_list = r_list, r_type = r_type, m.df = m.df, 
                                  equidistant = equidistant, path_size = path_size, path_end = path_end, 
                                  path_join = path_join, path_alpha = path_alpha, path_mitre = path_mitre, 
                                  path_arrow = path_arrow, print_plot = F, path_legend = path_legend, 
                                  path_legend_title = path_legend_title, tail_length = tail_length, 
                                  tail_size = tail_size, tail_colour = tail_colour, trace_show = trace_show, 
                                  trace_colour = trace_colour, path_fade = path_fade, ...)
  frames <- mapply(x = frames, y = unique(m.df$time), function(x, y) {
    attr(x, "time") <- y
    return(x)
  }, SIMPLIFY = F)
  return(frames)
}

# PLOTTING ------------------------------------------------------

plot_survival <- function(x = 0:100, a1 = 2, a2 = 0, a3 = 0.01, b1 = 60, b3 = 8, longevity = 69){
  
  plot(x, seq(0,1, length.out = length(x)), ylim = c(0,1), type = 'n', ylab = "p(survival)", xlab = "Age (years)")
  axis(2, at = seq(0,1,0.1))
  abline(v = longevity, lty = 1, col = "firebrick")
  
  lj = exp((-a1/b1) * (1-exp(-b1*x/longevity)))
  lc = exp(-a2*x/longevity)
  ls = exp((a3/b3) * (1-exp(b3*x/longevity)))
  
  lines(x, lj, col = "grey10", lty = 3)
  lines(x, lc, col = "grey20", lty = 4)
  lines(x, ls, col = "grey30", lty = 5)
  
  lines(x, lj*lc*ls, col = "black", lwd = 1.5)
  
}

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

growth_curve <- function(param, obj, cohortID, whaleID, ylabel){

  cohorts <- obj$param$cohorts
  dat <- data.table(obj$sim[[cohorts[id==cohortID, abb]]][whale %in% whaleID, ])
  ind.survived <- dat[day == 365 & alive == 1, whale,]
  dat <- dat[whale %in% ind.survived]
  n.ind <- length(ind.survived)

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
    ggtitle(label = ifelse(grepl(pattern = "calf", param), 
                           cohorts[1, name],
                           ifelse(grepl(pattern = "fetus_l", param), 
                                  "Fetus", cohorts[id == cohortID, name])
                           )) +
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

plot_raster <- function(r, 
                        duke = FALSE, 
                        positive = FALSE, 
                        quantile = FALSE, 
                        breaks = NULL, 
                        title = NULL,
                        xmin = NULL,
                        xmax = NULL,
                        ymin = NULL,
                        ymax = NULL){
  
  world <- targets::tar_read(world)
  
  if(class(r) == "SpatialGridDataFrame") r <- raster::raster(r)
  
  dat <- raster::as.data.frame(r, xy = TRUE)
  names(dat)[3] <- "Nhat"
  dat <- dat[complete.cases(dat),]
  
  if(positive) dat <- dat |> dplyr::filter(Nhat > 0)
  
  if (is.null(breaks)) {
    if (quantile) {
      colour.breaks <- c(0, quantile(dat$Nhat, seq(0, 1, 0.1)))
    } else {
      if (duke) {
        colour.breaks <- colour_breaks(dat)
      } else {
        colour.breaks <- c(0, rgeoda::natural_breaks(10, dat[, "Nhat", drop = FALSE]), max(dat$Nhat))
      }
    }
  } else {
    colour.breaks <- breaks
  }
  
  colour.breaks <- unique(colour.breaks)
  
  all.breaks <- data.frame(Ncol = cut(colour.breaks, breaks = colour.breaks, include.lowest = TRUE))
  all.breaks$mapcol <- pals::viridis(nrow(all.breaks))
  
  dat <- suppressWarnings(dat |> 
    dplyr::mutate(Ncol = cut(Nhat, breaks = colour.breaks, include.lowest = TRUE)) |> 
    dplyr::mutate(Ncol = factor(Ncol)) |> 
    dplyr::left_join(y = all.breaks, by = "Ncol"))
  
  levels(dat$Ncol) <-  gsub(pattern = ",", replacement = " – ", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\[|\\]", replacement = "", levels(dat$Ncol))
  levels(dat$Ncol) <-  gsub(pattern = "\\(|\\)", replacement = "", levels(dat$Ncol))
  
  r.dat <- raster::rasterFromXYZ(dat[, c("x", "y", "Nhat")], res = raster::res(r), crs = narw_crs())
  
  if(all(is.null(c(xmin, xmax, ymin, ymax)))){
    
    x.lim <- raster::extent(r.dat)[1:2]
    y.lim <- raster::extent(r.dat)[3:4]
    
  } else {
    
    x.lim <- c(xmin, xmax)
    y.lim <- c(ymin, ymax)
    
  }
  
  world.sf <- sf::st_as_sf(world)
  
  p <- ggplot2::ggplot(data = dat) +  
    ggplot2::geom_tile(aes(x,y, fill = Ncol)) +
    ggplot2::geom_sf(data = world.sf, fill = "lightgrey", color = "black", size = 0.25) +
    ylab("") + xlab("") +
    theme_narw() + 
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
    ggplot2::coord_sf(xlim = x.lim, ylim = y.lim, expand = FALSE) +
    # 
    # all.breaks$mapcol
    ggplot2::scale_fill_manual(values = pals::viridis(length(levels(dat$Ncol))),
                               guide = guide_legend(reverse = TRUE),
                               drop = FALSE,
                               na.value = "transparent")
  
  if(!is.null(title)) p <- p + ggplot2::ggtitle(title)
  
  return(p)
  
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

# FIGURES -----------------------------------------------------------------

figures <- function(obj = NULL, cex = 1, lwd = 2, scale.width = 1, scale.height = 1){
  
  regions <- targets::tar_read(regions)
  density_narw <- targets::tar_read(density_narw)
  world <- targets::tar_read(world)
  world_sf <- sf::st_as_sf(world)
  
  # ............................................................
  # Distribution of body mass
  # ............................................................
  
  if(!is.null(obj)){
    
    mass <- purrr::map(.x = obj$sim, .f = ~{
      .x[day > 0 & alive == 1, mass]
    }) |> tibble::enframe() |> 
      tidyr::unnest(cols = c(value))
    
    massplot <- ggplot2::ggplot(data = mass, aes(x = value)) + 
     ggplot2::geom_histogram(aes(y =..density..), position = "identity", fill = "grey", color = "black") +
      ggplot2::geom_density(col = "black") + 
      theme_narw() + 
      theme(panel.background = element_rect(fill = "transparent", colour = "black", size = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggplot2::facet_wrap(~name) + 
      xlab("Total mass (kg)")
    
  quiet(ggplot2::ggsave(filename = file.path("out/total_mass.jpeg"), plot = massplot))
    
  }
  
  # ............................................................
  # Estimated density surface within Canadian waters - June ----
  # ............................................................
  
  JunDens <- raster::raster(density_narw$Jun)
  JunDens <- raster::mask(JunDens, mask = regions[regions$country == "Canada",])
  JunDens <- raster::crop(JunDens, regions[regions$country == "Canada",])
  JuneDens.plot <- plot_raster(JunDens, duke = TRUE)
  
  quiet(ggplot2::ggsave(filename = file.path("out/NARW_density_Canada_June.jpeg"), plot = JuneDens.plot))
  
  # ............................................................
  # Feeding / nursing effort ----
  # ............................................................
  
  x <- seq(0, 1, length.out = 100)
  fcolors <- c("#104E8B","#EEB422")
  Lwidth <- 2500 * scale.width
  Swidth <- 1500 * scale.width
  height <- 1500 * scale.height
  
  jpeg(filename = file.path("out/feeding_nursing_effort.jpeg"), res = 300, width = Lwidth, height = height)
  
  plot(x, feeding_effort_vec(10, 0.41718, x), 
       lwd = lwd, 
       col = fcolors[1],
       type = "l", 
       ylab = "Feeding effort",
       lty = 1, 
       xlab = "Body condition (relative blubber mass, %)", 
       cex.axis = cex, 
       cex.lab = cex)
  
  lines(x, feeding_effort_vec(30, 0.6, x), lwd = 2, col = fcolors[2])
  abline(v = 0.41718, lty = 2, col = fcolors[1])
  abline(v = 0.6, lty = 2, col = fcolors[2])
  legend("bottomleft", bty = "n", lty = c(1,1), col = fcolors, legend = c("Other cohorts", "Calves/pregnant females"))
  
  dev.off()
  
  # f <- function(x, A = 1, D = 0, B, C, S) A + (D-A) / (1 + exp(B*(C-x)))^S

 #  plot(x, f(x, B = 29.04823, C = 0.3607178, S = 1.000036),
 #   lwd = 3, col = "#872E1C",
 #   type = "l", ylab = "Feeding effort", lty = 2, xlab = "Body condition (relative blubber mass, %)", cex.axis = 1.5, cex.lab = 1.5
 # )
 # lines(x, f(x, B = 21.02454, C = 0.4408761, S = 1.000066),
 #   lwd = 3, col = "#872E1C", lty = 1,
 #   type = "l", ylab = "Feeding effort", xlab = "Body condition (relative blubber mass, %)"
 # )

  # ............................................................
  # Foraging response to prey concentration ----
  # ............................................................
  
  jpeg(filename = file.path("out/foraging_response.jpeg"), res = 300, width = Lwidth, height = height)
  
  par(mfrow = c(1, 2))

 quiet(draw("tnorm", 1979.970, 1001.192, 0, Inf, 
       xlab = expression(Prey ~ concentration ~ (ind / m^3)), 
       ylab = "Density", 
       col = fcolors[1],
       main = "Minimum prey concentration", lwd = 2))
  
 abline(v = 2500, lty = 2, col = fcolors[1])
 prey <- seq(0, 10000, length.out = 10000)
 plot(prey, feeding_threshold_vec(2500, prey, 1.670151/1000),
      type = "l", 
      col = fcolors[2],
      lwd = lwd,
      ylab = "Feeding response",
      main = "Foraging response",
      xlab = expression(Prey ~ concentration ~ (ind / m^3)))

 par(mfrow = c(1, 1))
 dev.off()
  
  # ............................................................
  # Milk assimilation ----
  # ............................................................
  
  jpeg(filename = file.path("out/milk_assimilation.jpeg"), res = 300, width = Swidth, height = height)
  
  days <- 1:365
  plot(days, milk_assimilation_vec(days, 365, 288, 0.9),
       xlab = "Day", 
       ylab = "Milk assimilation rate", 
       type = "l",
       lwd = lwd, 
       col = fcolors[1],
       main = "Milk assimilation")
  
  dev.off()
  
  # ............................................................
  # Milk supply ----
  # ............................................................
  
  jpeg(filename = file.path("out/milk_supply.jpeg"), res = 300, width = Swidth, height = height)
  
  bc <- seq(0.05, 1, by = 0.001)
  mass <- 30000
  plot(bc, milk_supply_vec(0.05, 0.41718, rep(mass, length(bc)), bc*mass, -2),
       type = 'l', 
       col = fcolors[1],
       lwd = lwd, 
       xlab = "Body condition",
       main = "Milk provisioning", 
       ylab = "")
  
  dev.off()
  
  # ............................................................
  # Mass-at-length ----
  # ............................................................
  
  jpeg(filename = file.path("out/mass_at_length.jpeg"), res = 300, width = Swidth, height = height)

  ages <- seq(0,30,by = 0.1)
  lengths <- age2length_vec(ages)
  plot(ages, length2mass_vec(lengths, lean = 1), 
       type = "l", 
       lwd = lwd, 
       xlab = "Age (years)",
       ylab = "Mass (kg)")
  lines(ages, length2mass_vec(lengths), lty = 2, lwd = 2)
  legend("bottomright", 
         c("Total mass (Fortune et al. 2020)", "Lean mass (narwind)"), 
         lty = c(1,2), bty = "n", y.intersp = 1.25)
  
  dev.off()
  
  # ............................................................
  # Initial body condition ----
  # ............................................................
  
  jpeg(filename = file.path("out/calf_bodycondition.jpeg"), res = 300, width = Lwidth, height = height)
  
  bc_calf <- matrix(nrow = 100000, ncol = 7)
  colnames(bc_calf) <- c("Lg", "Mm", "Mv", "Mbo", "BCi", "Mbu", "BC")
  for(i in 1:10000){
    a <- numeric(4)
    length.b <- runif(1, 4.1, 5.3)
    a[1] <- fetal_tissue_mass(0.282, length.b) # muscle
    a[2] <- fetal_tissue_mass(0.102, length.b) # viscera
    a[3] <- fetal_tissue_mass(0.125, length.b) # bones
    bc <- rtnorm(-0.0061, 0.0049*sqrt(992), -0.4, 0.4)
    a[4] <- fetal_blubber_mass(length.b,bc, a[1], a[2], a[3], 700, 960, 930, 720)
    
    bc_calf[i, 1] <- length.b
    bc_calf[i, 2] <- a[1]
    bc_calf[i, 3] <- a[2]
    bc_calf[i, 4] <- a[3]
    bc_calf[i, 5] <- bc
    bc_calf[i, 6] <- a[4]
    bc_calf[i, 7] <- a[4]/sum(a)
  }
  hist(bc_calf[,"BC"], breaks = 50, xlab = "Relative fat mass (%))", main = "Body condition at birth", freq = FALSE)
  lines(density(bc_calf[,"BC"], adjust = 2, na.rm = TRUE))
  lines(seq(min(bc_calf[,"BC"], na.rm = TRUE), max(bc_calf[,"BC"], na.rm = TRUE), by = 0.01), 
        dnorm(seq(min(bc_calf[,"BC"], na.rm = TRUE), max(bc_calf[,"BC"], na.rm = TRUE), by = 0.01),
              mean = mean(bc_calf[,"BC"], na.rm = TRUE),
              sd = sd(bc_calf[,"BC"], na.rm = TRUE)), col = fcolors[2], lwd = lwd)
  legend("topleft", legend = c("Simulated", "Normal"), lty = 1, lwd = 2, col = c("black", fcolors[2]), bty = "n")
  
  dev.off()
  
  jpeg(filename = file.path("out/initial_bodycondition.jpeg"), res = 300, width = Lwidth, height = height)
  
  pbc <- init_bc()
  pcol <- c("grey20", "firebrick", "goldenrod2", "deepskyblue4")
  
  # Adults (males, resting females, pregnant females)
  plot(density(rtnorm_vec(10000, pbc[["adults"]][1], pbc[["adults"]][2], 0.05, find_maxBC()),            
                adjust = 2,
                from = 0, to = 1),
       main = "",
        xlab = "Body condition (% fat mass)", 
        xlim = c(0,1),
        ylab = "Density", 
        cex.lab = cex, 
        cex.axis = cex,
        cex.main = cex, 
        lwd = lwd, 
        col = pcol[3])
  
  # Calves at birth
  lines(density(rtnorm_vec(10000, pbc[["calves"]][1], pbc[["calves"]][2], 0.05, find_maxBC()), 
               adjust = 2,
               from = 0, to = 1),
       main = "",
       xlab = "Body condition (% fat mass)", 
       xlim = c(0,1),
       ylab = "Density", 
       cex.lab = cex, 
       cex.axis = cex,
       cex.main = cex, 
       lwd = lwd, 
       col = pcol[1])
  

  
  # Adults (males, resting females, pregnant females)
  lines(density(rtnorm_vec(10000, pbc[["juveniles"]][1], pbc[["juveniles"]][2], 0.05, find_maxBC()), 
                adjust = 2,
                from = 0, to = 1),
        xlab = "Body condition (% fat mass)", 
        xlim = c(0,1),
        ylab = "Density", 
        cex.lab = cex, 
        cex.axis = cex,
        cex.main = cex, 
        lwd = lwd, 
        col = pcol[2])
  
  # Adults (males, resting females, pregnant females)
  lines(density(rtnorm_vec(10000, pbc[["lactating"]][1], pbc[["lactating"]][2], 0.05, find_maxBC()), 
                adjust = 2,
                from = 0, to = 1),
        xlab = "Body condition (% fat mass)", 
        xlim = c(0,1),
        ylab = "Density", 
        cex.lab = cex, 
        cex.axis = cex,
        cex.main = cex, 
        lwd = lwd, 
        col = pcol[4])
  
  
  # # Adults (males, resting females, pregnant females)
  # quiet(draw("tnorm", pbc[["adults"]][1], pbc[["adults"]][2], 0.05, find_maxBC(), main = "", 
  #            xlab = "Body condition (% fat mass)", 
  #            xlim = c(0,1),
  #            ylab = "Density", 
  #            cex.lab = cex, 
  #            cex.axis = cex,
  #            cex.main = cex, 
  #            lwd = lwd, 
  #            # add = TRUE,
  #            col = pcol[3],
  #            from = 0, 
  #            to = 1))
  # 
  # # Lactating females
  # lines(seq(0, 1, length.out = 1000),
  #   truncnorm::dtruncnorm(
  #     seq(0, 1, length.out = 1000), 0.05, find_maxBC(),
  #     pbc[["lactating"]][1], pbc[["lactating"]][2]
  #   ),
  #   xlab = "Body condition (% fat mass)",
  #   xlim = c(0, 1),
  #   ylab = "Density",
  #   type = "l",
  #   cex.lab = cex,
  #   cex.axis = cex,
  #   cex.main = cex,
  #   lwd = lwd,
  #   col = pcol[4]
  # )
  # 
  # # Calves
  # quiet(draw("tnorm", pbc[["calves"]][1], pbc[["calves"]][2], 0.05, find_maxBC(),
  #            main = "", 
  #            xlab = "Body condition (% fat mass)", 
  #            xlim = c(0,1), 
  #            ylab = "Density", 
  #            cex.lab = cex, 
  #            cex.axis = cex, 
  #            cex.main = cex, 
  #            lwd = lwd, 
  #            col = pcol[1], 
  #            add = TRUE, 
  #            from = 0, 
  #            to = 1))
  # 
  # # Juveniles
  # quiet(draw("tnorm", pbc[["juveniles"]][1], pbc[["juveniles"]][2], 0.05, find_maxBC(),
  #            main = "", 
  #            xlab = "Body condition (% fat mass)", 
  #            xlim = c(0,1), 
  #            ylab = "Density", 
  #            cex.lab = cex, 
  #            cex.axis = cex, 
  #            cex.main = cex, 
  #            lwd = lwd, 
  #            col = pcol[2],
  #            add = TRUE, 
  #            from = 0, 
  #            to = 1))
  
  legend("topright", 
         bty = "n",
         lty = 1, 
         col = pcol,
         legend = c("Calves", "Juveniles", "Adults", "Late pregnant females"))
  
  dev.off()
  
  # ............................................................
  # Example prey ----
  # ............................................................
  
  preyr <- raster::raster(dummy_prey$Jul)
  ggplot2::ggplot() + 
    ggplot2::geom_raster(
      mapping = ggplot2::aes(x = x, y = y, fill = z),
      data = raster::as.data.frame(preyr, xy = TRUE)
    ) +
    ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", size = 0.25) +
    ggplot2::scale_fill_gradientn(colours = pals::viridis(100), na.value = "transparent")+
    coord_sf(xlim = raster::extent(preyr)[1:2],
             ylim = raster::extent(preyr)[3:4], expand = TRUE) +
    xlab("") + ylab("") +
    theme_narw() +
    ggplot2::labs(fill = expression(paste("Prey conc. (ind/", m^3, ")")))

  ggsave("out/dummy_prey.pdf", height = height, width = Swidth, units = "px")
  
  # ............................................................
  # Entanglement maps ----
  # ............................................................
  
  fishr <- purrr::map(.x = month.abb[7], .f = ~{
      fishing_layer[[.x]] |>
      raster::raster() |>
    raster::as.data.frame(xy = TRUE) |>
    dplyr::mutate(month = .x)}) |>
    do.call(what = "rbind") |> 
    dplyr::mutate(risk = cut(layer, breaks = c(-1, 0,2e-07, 1e-05, 1e-03, 1.5e-3, 1e-2, 0.025, 0.1, 0.3)))

  fishp <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        mapping = ggplot2::aes(x = x, y = y, fill = risk),
        data = fishr) +
    ggplot2::facet_wrap(vars(month), nrow = 2) +
      ggplot2::geom_sf(data = sf::st_as_sf(world), fill = "lightgrey", color = "black", size = 0.25) +
      # ggplot2::scale_fill_gradientn(colours = pals::viridis(100), na.value = "transparent") +
      coord_sf(xlim = raster::extent(raster::as.data.frame(raster::raster(fishing_layer[[1]]), xy = TRUE))[1:2],
               ylim = raster::extent(raster::as.data.frame(raster::raster(fishing_layer[[1]]), xy = TRUE))[3:4], expand = TRUE) +
      xlab("") + ylab("") +
    ggplot2::scale_fill_manual(values = pals::viridis(9), na.value = "transparent", guide = guide_legend(reverse = TRUE)) +
      theme_narw() +
      ggplot2::labs(fill = "Entanglement risk")
  
  ggsave("out/fishing_risk.pdf", height = height, width = Swidth, units = "px")
  
  # ............................................................
  # Vessel strike maps ----
  # ............................................................
  
   vesselr <- map_vessels(obj = scenario_02, 
                          baseline = TRUE, 
                          which.month = 7, 
                          z = "risk", 
                          do.plot = FALSE, 
                          strike_scalar = 1e-7)
  
  z_brks <- raster::as.data.frame(vesselr[[1]][[2]], na.rm = TRUE)["risk"]
  # z_brks <- rbind(z_brks, raster::as.data.frame(vesselr[[1]][[1]], na.rm = TRUE)["risk"])
  z_brks <-  unique(plyr::round_any(
    c(0, rgeoda::natural_breaks(k = 20, z_brks["risk"]
    ), max(z_brks["risk"])),
    accuracy = 1e-05, f = ceiling))
   
  jpeg(filename = file.path("out/strike_risk.jpeg"), res = 300, width = 2400, height = height)

  terra::plot(terra::rast(vesselr[[1]][[2]]),
              pax = list(cex.axis = 0.85),
              legend = "bottomright",
              col = pals::viridis(20),
              breaks = z_brks)
  
  terra::plot(terra::vect(world), col = "grey", add = TRUE)
  dev.off()
  
   jpeg(filename = file.path("out/strike_risk_scenario01.jpeg"), res = 300, width = 2400, height = height)
   terra::plot(terra::rast(vesselr[[1]][[1]]),
               pax = list(cex.axis = 0.85),
               legend = "bottomright",
               col = pals::viridis(20),
               breaks = z_brks)
   
   terra::plot(terra::vect(world), col = "grey", add = TRUE)
   
   dev.off()
  
   jpeg(filename = file.path("out/strike_risk_difference.jpeg"), res = 300, width = 2400, height = height)
   terra::plot(terra::rast(vesselr[[1]][[1]]-vesselr[[1]][[2]]),
               pax = list(cex.axis = 0.85),
               legend = "bottomright",
               col = pals::viridis(20))
   
   terra::plot(terra::vect(world), col = "grey", add = TRUE)
   
   dev.off()
   
  # ............................................................
  # Noise maps ----
  # ............................................................
   
   noisemaps <- map_noise(obj = scenario) 
   jpeg(filename = file.path("out/noise_risk.jpeg"), res = 300, width = 1800, height = 1800)
   plot(raster::projectRaster(raster::raster(noisemaps$`15-01`), crs = narw_crs(TRUE)),
        col = pals::parula(100),
        xlim = c(-72,-69.5), 
        ylim = c(40.2,41.8), xlab = "Longitude", ylab = "Latitude")
   plot(raster::projectRaster(raster::raster(noisemaps$`15-01`), crs = narw_crs(TRUE)),
        col = pals::parula(100), add = TRUE)
   plot(sp::spTransform(world, narw_crs(TRUE)), col = "grey", axes = TRUE, add = TRUE)
   points(scenario_01$locs[scenario_01$locs$windfarm == 1, c("longitude", "latitude")], pch = "+")
   dev.off()
   
   
  # ............................................................
  # Birth probability ----
  # ............................................................
  
  jpeg(filename = file.path("out/prob_birth.jpeg"), res = 300, width = Swidth, height = height)
  
  days <- seq(40, 135, 1)
  plot(days, pbirth_vec(days, 69, 40),
       type = "l", 
       xlab = "Day",
       ylab = "p(birth)", 
       cex.axis = cex, 
       cex.lab = cex,
       lwd = 1.5)
  abline(v = 69, col = fcolors[2], lty = 2, lwd = lwd)
  legend("topleft", legend = "Enter calving\n grounds", lty = 2, lwd = lwd, col = fcolors[2], bty = "n")
  
  dev.off()
  
  # x <- 0:100
  # plot(x, survivorship_vec(x), ylab = "p(survival)", xlab = "Age (years)", type = "l", lwd = 1.5)
  # lines(x, survivorship_vec(x, a1 = 2, a2 = 0, a3 = 0.01, b1 = 60, b3 = 8), lty = 2)
  # # lines(x, survivorship_vec(x, a1 = 0.1, a2 = 0, a3 = 0.01, b1 = 60, b3 = 8), lty = 2)
  # lines(x, survivorship_vec(x, a1 = 10, a2 = 0.5, b1 = 10, b3 = 4), lty = 2)
  # lines(x, survivorship_vec(x, a1 = 14, a2 = 0.1, b1 = 20, b3 = 7), lty = 2)
  # lines(x, survivorship_vec(x, a1 = 1, a2 = 0.12, b1 = 30, b3 = 12), lty = 2)
  # lines(x, survivorship_vec(x, a1 = 1, a2 = 0.12, b1 = 30, b3 = 25), lty = 2)
  # abline(v = 69, lty = 1, lwd = 1.5, col = "firebrick")
  
  # entgl_check <- purrr::map(.x = 1:100000, .f = ~entanglement_event(1)) |> do.call(what = rbind) |> data.table::as.data.table()
  # names(entgl_check) <- c("entgl", "severity", "head", "gape", "duration", "start", "end", "dead")
  # 
  # table(entgl_check$head)
  # entgl_check[head == 1, .(gape = mean(gape), min = min(gape), max = max(gape)), severity]
  
  
  # 
  # estRayleighParams(6.5, 150)
  # # 45 km as target mean
  # # 45/sqrt(pi) = sigma
  # # 
  # target.mean = 5.4
  # rn <- extraDistr::rrayleigh(1000000, sigma = sqrt(2)*target.mean/sqrt(pi))
  # hist(rn, main = "Rayleigh distribution (sigma = 3.807)", xlab = "Step lengths (km)")
  # abline(v = 5, col = "firebrick", lwd = 1.5)
  # mean(rn)
  # max(rn)
  # sum(rn<=5)/length(rn)
  
}

# UTILITIES ------------------------------------------------------

jpeg2gif <- function(input.dir, pattern = NULL, frame.rate = 1, filename){
  imgs <- list.files(input.dir, full.names = TRUE)
  if(!is.null(pattern)){
    imgs <- imgs[grepl(pattern = pattern, x = imgs)]
  }
  
  cat("Reading input files ... \n")
  img_list <- lapply(imgs, magick::image_read)
  
  ## join the images together
  cat("Creating video ... \n")
  img_joined <- magick::image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated <- magick::image_animate(img_joined, fps = frame.rate)
  
  ## save to disk
  cat("Saving file ... \n")
  magick::image_write(image = img_animated, 
                      quality = 100,
                      density = 300,
                      path = file.path(input.dir, paste0(filename, ".gif")))
  cat("Done!")
}

label <- function(obj, lb){
  if(!inherits(obj, "narwsim") & !inherits(obj, "narwproj")) 
    stop("Input object must be of class <narwsim> or <narwproj>")
  obj$param$label <- lb
  return(obj)
}

proj_timeline <- function(schedule){
  
  operator <- " -----> "
  
  # Baseline
  if(any(schedule == 0)){
    if(all(schedule == 0)) base.msg <- paste0("Baseline (yr 1 to ", length(schedule) - 1, ")") else base.msg <- NULL
  } else {
    base.msg <- NULL
  }
  
  # Construction phase
  if(any(schedule == 1)){
    const.msg <- paste0(if(!is.null(base.msg)) operator, "Construction (yr ", min(which(schedule == 1) - 1), ")")
  } else {
    const.msg <- NULL
  }
  
  # O&M phase
  if(any(schedule == 2)){
    if(all(schedule == 2)) om.msg <- paste0("Operation & maintenance (yr 1 to ", length(schedule) - 1, ")") else om.msg <- paste0(operator, "Operation & maintenance (yr ", min(which(schedule == 2) - 1), " to ", max(which(schedule == 2) - 1), ")") 
  } else {
    om.msg <- NULL
  }
  
  cat(base.msg, const.msg, om.msg, sep = "")
  cat("\n\n")
}

add_name <- function(obj, nm){
  obj$param$name <- nm
  return(obj)
}

update_pkg <- function(){
  
  fromPath <- getwd()
  toPath <- "../narwind"
  
  # R files
  R.files <- list.files(path = "R/", recursive = FALSE)
  R.files <- R.files[!R.files %in% c("narwinddev-package.R", "zzz.R", "RcppExports.R", "data_targets.R")]
  R.files <- file.path("R", R.files)
  
  # C++ files
  C.files <- list.files(path = "src/", recursive = FALSE)
  C.files <- C.files[!C.files %in% c("backup", "narwinddev.so", "RcppExports.cpp", "RcppExports.o", "simtools.o")]
  C.files <- file.path("src", C.files)
  
  # Package help files (including datasets)
  man.files <- list.files(path = "man/", recursive = FALSE)
  man.files <- man.files[!man.files %in% c("narwinddev-package.Rd", "narwinddev.Rd")]
  man.files <- file.path("man", man.files)
  
  # Data files
  dat.files <- list.files(path = "data/", recursive = FALSE)
  dat.files <- dat.files[!dat.files %in% c("calanus", "densitymaps", "dirichlet", "elicitation", "gear", "parameters", "regions", "regions", "sightings", "vessels", "windfarms", "basemap", "popviability", "vignettes")]
  dat.files <- file.path("data", dat.files)
  
  # Web files
  pkgdown.files <- list.files(path = "pkgdown/", recursive = FALSE)
  pkgdown.files <- file.path("pkgdown", pkgdown.files)
  
  # Package files
  root.files <- c("README.md", "README.Rmd", "_targets.yaml", "_targets.R", "DESCRIPTION")
  
  # Vignette files
  vig.files <- list.files(path = "vignettes/", recursive = FALSE)
  vig.files <- file.path("vignettes", vig.files)
  
  # All files
  allfiles <- c(R.files, C.files, man.files, dat.files, pkgdown.files, root.files, vig.files)
  
  pb <- progress::progress_bar$new(total = length(allfiles), 
                                   width = 80, 
                                   format = "Processing changes [:bar] :percent",
                                   clear = TRUE)
  
  file.changes <- vector(mode = "list", length = length(allfiles))
  
  for(f in seq_along(allfiles)){
    pb$tick()
    # If file does not exist in target directory
    if(!file.exists(file.path(toPath, allfiles[f]))){
      file.changes[[f]] <- allfiles[f]
    } else {
      isSame <- unname(tools::md5sum(file.path(fromPath, allfiles[f])) == tools::md5sum(file.path(toPath, allfiles[f])))
      if(!isSame){
        file.changes[[f]] <- allfiles[f]
      }
    }
  }
  
  file.changes <- purrr::compact(file.changes)
  
  
  pb <- progress::progress_bar$new(total = length(file.changes), 
                                   width = 80, 
                                   format = "Copying files [:bar] :percent",
                                   clear = TRUE)
  
  for (g in 1:length(file.changes)) {
    pb$tick()
    invisible(file.copy(from = file.path(fromPath, file.changes[[g]]), 
                        to = file.path(toPath, file.changes[[g]]), 
                        overwrite = TRUE))
  }
  
  if(length(file.changes) == 0){
    
    cat("All files are up to date!")
    
  } else {
    
    cat("Updated files:\n")
    cat("--------------------\n")
    purrr::walk(.x = file.changes, .f = ~cat(.x, "\n"))
    
  }
  
  # Update narw.R
  narw.txt <- readLines("R/narw.R")
  stringr::str_replace(string = narw.txt,
                       pattern = ".packages = c\\(\"Rcpp\"\\)",
                       replacement = ".packages = c\\(\"Rcpp\", \"narwind\"\\)") |>
  stringr::str_replace(pattern = "Rcpp::sourceCpp\\(\"src/simtools.cpp\"\\)",
                         replacement = "") |>
    writeLines(file.path(toPath, "R", "narw.R"))

  # Update DESCRIPTION
  desc.txt <- readLines("DESCRIPTION")
  stringr::str_replace(string = desc.txt,
                       pattern = "narwinddev",
                       replacement = "narwind") |>
    writeLines(file.path(toPath, "DESCRIPTION"))
  

}

# outersect <- function(x, y) {
#   sort(c(setdiff(x, y), setdiff(y, x)))
# }

split_byNA <- function(x){
  idx <- 1 + cumsum( is.na( x ) )
  not.na <- ! is.na( x )
  split( x[not.na], idx[not.na] )
}

tickmark <- function(){
  return("\U2714")
}

crossmark <- function(){
  return("\U2717")
}

console <- function(msg, suffix = NULL){
  if(is.null(suffix)){
    cat("\r", paste0(msg, " ..."), sep = "")
  } else {
    if(nchar(suffix) > 1){
      blank <- paste0(rep(" ", max(nchar(suffix)-3,2)), collapse = "")
    } else {
      blank <- "  "
    }
    cat("\r", paste0(msg, " ", suffix, blank), sep = "")
    cat("\n")
  }
}

save_params <- function(){
  rio::convert(in_file = "/Users/philbouchet/OneDrive - University of St Andrews/BOEM/project/WP1_synthesis/BOEM_140M0121C0008_ModelParameters_3.0.xlsx",
               out_file = "/Volumes/GoogleDrive/My Drive/Documents/git/narwind/data/parameters/BOEM_140M0121C0008_ModelParameters.csv")
  targets::tar_delete(params)
  targets::tar_make(params)
  save_object("params", TRUE, TRUE)
}

do.align <- function(var = "preyconc", lyr = dummy_prey){
  
  alive <- which(m$sim[[1]][, alive][-1] == 1)
  fromSim <- unname(unlist(m$sim[[1]][, var, with = FALSE])[-1])[alive]
  fromRast <- numeric(length(fromSim))
  for(j in alive){
    xy <- m$sim[[1]][day==j, list(easting, northing)]
    mo <- m$sim[[1]][day==j, month]
    r <- lyr[[month.abb[mo]]]
    if(!inherits(r, "RasterLayer")) r <- raster::raster(r)
    fromRast[j] <- raster::extract(x = r, y = xy)
  }
  identical(fromSim, fromRast)
}

get_dates <- function(start.month = 10, ndays = 457, strip = TRUE, gantt = FALSE){
  
  thisyear <- lubridate::year(lubridate::now())
  
  if(gantt){
    
    date_seq <- seq(lubridate::ymd(paste0(thisyear, "-01-01")), lubridate::ymd(paste0(thisyear, "-12-31")), by = '1 day')
    if(sum(grepl(pattern = "02-29", x = date_seq))>0) date_seq <- date_seq[-which(grepl(pattern = "02-29", x = date_seq))]
    
  } else {
    
    # Use 2021 as non-leap year, as AIS data were obtained for 2019, which is a non-leap year as well
    start.date <- lubridate::as_date(paste0(2021, "-", stringr::str_pad(start.month, width = 2, pad = "0"), "-01"))
    date_seq <- seq(start.date, by = 'day', length.out = ndays)
    if(strip){
      date_seq <- format(as.Date(date_seq), "%d-%m") # Strip year
      date_seq <- c(paste0("00-", start.month), date_seq)
    } else {
      date_seq <- c(paste0("2021-", start.month, "-00"), as.character(date_seq))
    }
  }
  return(date_seq)
}


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

fix_paths_vignettes <- function(vignette.name = "narwind"){
  vignette.path <- file.path("./docs/articles", paste0(vignette.name, ".html"))
  tx  <- readLines(vignette.path)
  tx_mod  <- gsub(pattern = "../../../../../My%20Drive/Documents/git/narwind/articles/", replace = "", x = tx)
  writeLines(tx_mod, con = vignette.path)
}

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

format_dt <- function(dt, relative = FALSE, N, direction = "col") {
  if (nrow(dt) > 0) {
    if (relative) {
      dt |>
        janitor::adorn_percentages(denominator = direction) |>
        janitor::adorn_pct_formatting() |>
        janitor::adorn_ns()
    } else {
      dt |>
        janitor::adorn_totals("col") |> 
        dplyr::mutate(Total = N - Total) |> 
        janitor::untabyl() |> 
        janitor::adorn_percentages(denominator = "row") |>
        janitor::adorn_pct_formatting() |>
        janitor::adorn_ns() |> 
        dplyr::select(-Total)

    }
  } else {
    dt
  }
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

draw <- function(distrib = "norm", mu, std = NULL, L = -Inf, U = Inf, seed = 215513, title = NULL, 
                 from = NULL, to = NULL,
                 verbose = TRUE, 
                 ...){
  
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
    
    curve(expr = truncnorm::dtruncnorm(x = x, L, U, L, truncN.param), from = ifelse(is.null(from), min(n), from), to = ifelse(is.null(to), max(n), to), ...)
    
  } else if(distrib == "tnorm"){
    
    n <- truncnorm::rtruncnorm(10000, L, U, mu, std)

    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = truncnorm::dtruncnorm(x = x, a = L, b = U, mean = mu, sd = std), from = ifelse(is.null(from), min(n), from), to = ifelse(is.null(to), max(n), to), ...)
    
  } else if(distrib == "norm"){
    
    n <- rnorm(10000, mu, std)
    
    if(verbose) {
    cat("min:", min(n), "\n")
    cat("max:", max(n),  "\n")
    cat("mean:", mean(n),  "\n")
    cat("sd:", sd(n), "\n")
    }
    
    curve(expr = dnorm(x = x, mu, std), from = ifelse(is.null(from), min(n), from), to = ifelse(is.null(to), max(n), to), ...)
    
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
    
    curve(expr = dgamma(x = x, shape = gamma.params$shape, scale = gamma.params$scale), from = ifelse(is.null(from), min(n), from), to = ifelse(is.null(to), max(n), to),
          ...)
    
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
    
    curve(expr = dbeta(x = x, beta.params$alpha, beta.params$beta), from = 0, to = 1, ...)
  }
}

rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  
  if(is.null(x.min)){
    if(inherits(x, "RasterLayer")){
      x.min <- raster::cellStats(x, "min")
    } else if(is.numeric(x)){
      x.min <- min(x, na.rm = TRUE)
    }
  }
  
  if(is.null(x.max)){
    if(inherits(x, "RasterLayer")){
      x.max <- raster::cellStats(x, "max")
    } else if(is.numeric(x)){
      x.max <- max(x, na.rm = TRUE)
    }
  }
  
  out <- new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
  
  return(out)
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
    if(tg) suppressWarnings(targets::tar_load(obj))
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
