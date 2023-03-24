# MODEL ------------------------------------------------------

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

weighted_density <- function(){
  
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

  rp.SEUS <- rp.GSL <- rp.CA <- rp

  cabot.xy <- sp::SpatialPoints(coords = cbind(1400, 1600), proj4string = narw_crs())
  
  circular.d <- raster::distanceFromPoints(object = densr[[1]], xy = cabot.xy) |> 
    raster::mask(mask = densr[[1]])
  
  # Set value to 1 in SEUS
  rp.SEUS[!rp.SEUS$layer == which(sort(reg$region) == "SEUS")] <- 0
  rp.SEUS[rp.SEUS>0] <- 1

  # Set value to 1 in GSL
  rp.GSL[!rp.GSL$layer == which(sort(reg$region) == "GSL")] <- 0
  rp.GSL[rp.GSL>0] <- 1
  
  # Set value to 1 in SEUS
  rp.CA[rp.CA$layer == which(sort(reg$region) == "GSL")] <- 99
  rp.CA[rp.CA$layer == which(sort(reg$region) == "CABOT")] <- 99
  rp.CA[rp.CA$layer == which(sort(reg$region) == "SCOS")] <- 99
  rp.CA[rp.CA<99] <- 0
  rp.CA[rp.CA==99] <- 1
  
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
      
      # OPTION 1
      
      # # Determine mass outside of SEUS
      # mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.SEUS)), na.rm = TRUE)
      # # Redistribute mass across range according to north-south gradient
      # d1 <- densr[[.x]] * rp.SEUS
      # d2 <- mout * r_south/sum(raster::getValues(r_south), na.rm = TRUE) # Sums to mout
      # d.out <- d1 + d2
      
      # OPTION 2
      
      # # Clip density to SEUS
      # d1 <- densr[[.x]] * rp.SEUS
      # # Re-scale gradient outside of SEUS
      # d1[d1==0] <- NA
      # dv <- na.omit(apply(X = as.matrix(d1), MARGIN = 1, FUN = mean, na.rm = TRUE))
      # d2 <- rescale_raster(x = r_south * (1-rp.SEUS), new.max = dv[1])
      # d2[d2==0] <- NA
      # d.out <- raster::merge(d1, d2)
      
      # OPTION 3
      pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs())
      pts <- sp::spsample(pathGSL, n = 1000, type = "regular")
      distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |> 
        raster::mask(mask = densr[[1]])
      dist.r <- (1/distances)^(1/4) # Fourth root
      d.out <- rescale_raster(dist.r, new.min = raster::cellStats(dist.r, "min"), 
                               new.max = raster::cellStats(dist.r*(1-rp.CA), "max")) + r_southbound
      
      as(d.out, "SpatialGridDataFrame")

    } else if(.x == 6) { # Migrating north in June

      # OPTION 1
      
      # Determine mass outside of GSL
      # mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.GSL)), na.rm = TRUE)
      # Redistribute mass across range according to north-south gradient
      # d1 <- densr[[.x]] * rp.GSL
      # d2 <- mout * r_north/sum(raster::getValues(r_north), na.rm = TRUE) # sums to mout
      # d.out <- d1 + d2
      
      # OPTION 2
      
      # # Clip density to GSL
      # d1 <- densr[[.x]] * rp.GSL
      # # Re-scale gradient outside of GSL
      # d1[d1==0] <- NA
      # dv <- na.omit(apply(X = as.matrix(d1), MARGIN = 1, FUN = mean, na.rm = TRUE))
      # d2 <- rescale_raster(x = r_north * (1-rp.GSL), new.max = dv[1])
      # d2[d2==0] <- NA
      # d.out <- raster::merge(d1, d2)
      
      # OPTION 3
      pathGSL <- raster::shapefile("data/pathGSL.shp") |> sp::spTransform(CRSobj = narw_crs()) |> rgeos::gBuffer(width = 50)
      pts <- sp::spsample(pathGSL, n = 10000, type = "regular")
      distances <- raster::distanceFromPoints(object = densr[[1]], xy = pts) |> raster::mask(mask = densr[[1]])
      dist.r <- (1/distances)^(1/4) # Fourth root
      d.out <- rescale_raster(dist.r, new.min = raster::cellStats(dist.r, "min"), 
                              new.max = raster::cellStats(dist.r*rp.CA, "max")) + r_northbound

      as(d.out, "SpatialGridDataFrame")
      
    } else if(.x %in% c(1:2, 12)) { # Breeding season in the SEUS
      
      # Mass outside of SEUS
      mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.SEUS)), na.rm = TRUE)

      # Redistribute mass
      d <- densr[[.x]] * rp.SEUS
      d <- d + d*mout/sum(raster::getValues(d), na.rm = TRUE) + 1e-06
      as(d, "SpatialGridDataFrame")
      
    } else if(.x %in% c(7:10)) { # Foraging season in the GSL
      
      # Mass outside of GSL
      mout <- sum(raster::getValues(densr[[.x]]*abs(1-rp.GSL)), na.rm = TRUE)
      
      # Redistribute mass
      d <- densr[[.x]] * rp.GSL
      d <- d + d*mout/sum(raster::getValues(d), na.rm = TRUE) + 1e-06
      as(d, "SpatialGridDataFrame")

    } else { # March through to May

      dens[[.x]]

    }

  }) |> purrr::set_names(nm = month.abb)

  return(wdens)
  
}

init_xy <- function(maps, maps.weighted, coords, cohort.id, nsim, northSEUS = -12, southGSL = 1400){
  
  if(cohort.id == 5) m <- maps.weighted else m <- maps
  if(cohort.id > 0) m <- maps.weighted else m <- maps
  
  out <- do.call(rbind, lapply(seq_along(m), function(mo) {
    
    p <- as.numeric(m[[mo]])
    p[!is.finite(p)] <- 0
    
    # Foraging season in the GSL (June though to Oct – Crowe et al., 2021)
    if(mo %in% c(7:10)){
      notGSL <- which(coords[,2] < southGSL)
      GSL <- which(coords[,2] >= southGSL)
      N.notGSL <- sum(p[notGSL])
      p[notGSL] <- 0
      p[GSL][p[GSL] > 0] <- p[GSL][p[GSL] > 0] + (N.notGSL / length(p[GSL][p[GSL] > 0]))
      
      # Breeding season (Nov through to February - Krystan et al. 2018)
    } else if (mo %in% c(1:2, 11:12)){
      notSEUS <- which(coords[,2] > northSEUS)
      SEUS <- which(coords[,2] <= northSEUS)
      N.notSEUS <- sum(p[notSEUS])
      p[notSEUS] <- 0
      p[SEUS][p[SEUS] > 0] <- p[SEUS][p[SEUS] > 0] + (N.notSEUS / length(p[SEUS][p[SEUS] > 0]))
      
    } else if (mo == 5){
      
      area <- which(coords[,2] <= 1200)
      notarea <- setdiff(seq_along(p), area)
      N.notarea <- sum(p[notarea])
      p[notarea] <- 0
      p[area][p[area] > 0] <- p[area][p[area] > 0] + (N.notarea / length(p[area][p[area] > 0]))
      
    } else if (mo == 6){
      
      SCOS <- which(coords[,1] >= 1400 & coords[,2] <=1500)
      notSCOS <- setdiff(seq_along(p), SCOS)
      N.notSCOS <- sum(p[notSCOS])
      p[notSCOS] <- 0
      p[SCOS][p[SCOS] > 0] <- p[SCOS][p[SCOS] > 0] + (N.notSCOS / length(p[SCOS][p[SCOS] > 0]))
      
    } else {
      p[p <= 0.01] <- 0
    }
    sample(x = 1:length(p), size = nsim, replace = TRUE, prob = p)
  }))
  
  return(out)
}

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

consolidate <- function(dtl, nsim, cnames, dates){
  purrr::map2(.x = dtl,
              .y = cnames, 
              .f = ~{
                dt <- data.table::data.table(day = rep(0:365, times = nsim), 
                                               date = rep(dates, times = nsim),
                                               whale = rep(1:nsim, each = 366))
                a <- cbind(array2dt(.x), dt)
                a$region <- sort(regions$region)[a$region]
                a$cohort_name <- .y
                data.table::setcolorder(a, c((ncol(a)-3):(ncol(a)-1), 1:(ncol(a)-4))) 
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
  
proxy_prey <- function(seed = 215513){
  
  set.seed(seed)
  
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
  gradient.r$val <- 1/(1+exp(-0.05*(gradient.r$x-550)))
  gradient.r$val[which(is.na(gradient.r$Nhat))] <- NA
  gradient.r <- raster::rasterFromXYZ(xyz = gradient.r[,c("x", "y", "val")], res = raster::res(target.r), crs = narw_crs())
  # plot(gradient.r)
  
  # Create a correlated random field as a dummy prey surface
  grid <- sp::makegrid(support_poly, cellsize = 15) |>
    dplyr::rename(x = x1, y = x2)
  area.ex <- raster::extent(support_poly)
  grid <- expand.grid(seq(area.ex[1], area.ex[2], length.out = 200),
                      seq(area.ex[3]-500, area.ex[4], length.out = 200))
  names(grid) <- c("x", "y")
  
  out <- purrr::map(.x = month.abb, .f = ~{

    g.dummy <- gstat::gstat(formula = 
                            z ~ 1 + y, 
                            locations = ~ x + y, 
                            dummy = TRUE, 
                            beta = 1800, 
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

plot_raster <- function(r, prob = FALSE, zero = FALSE){
  
  if(class(r) == "SpatialGridDataFrame") r <- raster::raster(r)
  dat <- raster::as.data.frame(r, xy = TRUE)
  names(dat)[3] <- "Nhat"
  dat <- dat[complete.cases(dat),]
  
  if(zero){
    dat.c <- dat |> dplyr::filter(Nhat > 0)
    colour.breaks <- c(0, quantile(dat.c$Nhat, seq(0,1, 0.1)))
  }else {
    if(!prob){
      colour.breaks <- colour_breaks(dat)
    } else {
      colour.breaks <- c(0, rgeoda::natural_breaks(10, dat[, "Nhat", drop = FALSE]), max(dat$Nhat))
      # colour.breaks <- c(seq(0, 0.01, length.out = 5), seq(0.02, 0.1, length.out = 5), 0.2, 0.5, 0.8, 1)
    }
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
    coord_sf(xlim = x.lim, ylim = y.lim, expand = FALSE) +
    scale_fill_manual(values = pals::viridis(length(levels(dat$Ncol))),
                      guide = guide_legend(reverse = TRUE))
  
}

# Jason's colour scale
colour_breaks <- function(dat){
  colour.breaks <- 25*c(0,0.016,0.025,0.04,0.063,0.1,0.16,0.25,0.4,0.63,1,1.6,2.5,4,6.3,10)/100
  colour.breaks <- c(colour.breaks, seq(max(colour.breaks), ceiling(max(dat$Nhat, na.rm = TRUE)), length.out = 5))
  colour.breaks <- round(colour.breaks[!duplicated(colour.breaks)],3)
  return(colour.breaks)
}

#' @export
plot.xyinits <- function(obj, month = NULL) {
  if(is.null(month)) month <- 1:12
  for(k in 1:length(obj)){
    if(!"xyinits" %in% class(obj)) stop("Input not recognized")
    par(mfrow = c(3,4))
    for(h in month){
      sp:::plot.SpatialPolygons(world, col = "grey", axes = TRUE, main = month.name[h])
      xpos <- sapply(obj[[k]][h, ], FUN = function(x) raster::xFromCell(raster::raster(density_narw$Feb), x))
      ypos <- sapply(obj[[k]][h, ], FUN = function(x) raster::yFromCell(raster::raster(density_narw$Feb), x))
      points(cbind(xpos, ypos), pch = 16, col = "orange")
      # legend("topleft", legend = month.abb, col = pals::glasbey(12), pch = 16)
    }
  }
  par(mfrow = c(1,1))
}

deg2radians_R <- function(angle){
  angle*pi/180;
}

gape_size_R <- function(L, omega, alpha){
  (deg2radians_R(alpha) * ((omega^2)/4 + (0.2077*L - 1.095)^2))/2
}


# UTILITIES ------------------------------------------------------

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
  return(out$par)
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


save_objects <- function(obj = NULL, redo = FALSE){
  if(is.null(obj)) obj <- c("regions", "regions_m",  "world", "support_poly", "density_narw", "density_support", "density_weighted", "turbines", "dummy_prey", "dummy_noise", "dummy_vessel", "dummy_fishing", "daylight", "params")
  # if(redo & length(obj) == 1) targets::tar_delete(obj); targets::tar_make(obj)
  suppressWarnings(targets::tar_load(obj))
  for(i in obj) {
    file.name <- paste0("data/", i, ".rda")
    if(file.exists(file.name)) file.remove(file.name)
    do.call(save, c(lapply(i, as.name), file = paste0("data/", i, ".rda")))
  }
}

format_table <- function(df, top = TRUE, bottom = TRUE, sign = "-"){
  # dashes <- purrr::map_dbl(.x = names(df), .f = ~nchar(.x))
  # dashes <- purrr::map(.x = dashes, .f = ~paste0(rep(sign, .x), collapse = ""))
  rbind(df[1:nrow(df) - 1,], c("---", rep("",ncol(df)-1)), df[nrow(df),])
}
