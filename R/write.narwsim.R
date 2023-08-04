#' Summary
#'
#' Summary information
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

write.narwsim <- function(obj, ...){
  
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values for optional arguments
  whaleID <- seq_len(obj$param$nsim)
  cohortID <- obj$param$cohortID
  prefix <- "narwsim"
  overwrite <- TRUE

  # Default values
  if(length(args) > 0) {
    if("cohortID" %in% names(args)) cohortID <- args[["cohortID"]]
    if("whaleID" %in% names(args)) whaleID <- args[["whaleID"]]
    if("prefix" %in% names(args)) prefix <- args[["prefix"]]
    if("overwrite" %in% names(args)) overwrite <- args[["overwrite"]]
  }

  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")

  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohortID, abb]
  cohort.names <- cohorts[id %in% cohortID, name]
  sim <- obj$sim

  cat("Saving ...\n")

  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = length(cohort.ab), clear = FALSE, width = 80
  )
  
  for(k in cohort.ab){

    pb$tick() # Update progress bar
    
    const <- sim[[k]][whale %in% whaleID & day > 0, list(day, date, month, whale, cohort)]

    # // ------------------------------------------------------------
    # // LOCATIONS
    # // ------------------------------------------------------------
    
    xy <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(easting, northing, region, d_travel, swimspeed)])
    
    # // ------------------------------------------------------------
    # // ATTRIBUTES
    # // ------------------------------------------------------------
    
    attrib <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(alive, age, starve, died, p_died, bc, length, length_a, length_b, length_c, mass, leanmass, fatmass, mass_a, mass_b, mouth_r, mouth_a, mouth_w, gsl, seus)])
    
    # // ------------------------------------------------------------
    # // ATTRIBUTES (calves)
    # // ------------------------------------------------------------
    
   
    if (cohorts[abb == k, id] == 5) {
      attrib_calf <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(cohort_calf, born, alive_calf, age_calf, starve_calf, bc_calf, length_calf, La_calf, Lb_calf, Lc_calf, mass_calf, leanmass_calf, fatmass_calf, ma_calf, mb_calf, mouth_r_calf, mouth_a_calf, mouth_w_calf)])
    } else {
      attrib_calf <- cbind(const[, !c("cohort")], data.table::data.table(cohort_calf = NA, born = NA, alive_calf = NA, age_calf = NA, starve_calf = NA, bc_calf = NA, length_calf = NA, La_calf = NA, Lb_calf = NA, Lc_calf = NA, mass_calf = NA, leanmass_calf = NA, fatmass_calf = NA, ma_calf = NA, mb_calf = NA, mouth_r_calf = NA, mouth_a_calf = NA, mouth_w_calf = NA))
    }
    
    # // ------------------------------------------------------------
    # // BEHAVIOR
    # // ------------------------------------------------------------
    
    activ <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(glide, glide_feed, glide_echelon, t_travel, t_feed, t_rest_nurse)])

    # // ------------------------------------------------------------
    # // GESTATION
    # // ------------------------------------------------------------
    
    if (cohorts[abb == k, id] == 4) {
      
      gest <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(abort, fetus_l, fetus_m, delta_fetus_m, delta_fetus_l, birth_l, birth_m, muscle_m, viscera_m, bones_m, blubber_m, muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, blubber_lip, blubber_pro,
          prop_mu, prop_visc, prop_bones, dens_blu, dens_mu, dens_visc, dens_bo)])
      
    } else {
      
      gest <- cbind(const, data.table::data.table(abort = NA, fetus_l = NA, fetus_m = NA, delta_fetus_m = NA, delta_fetus_l = NA, birth_l = NA, birth_m = NA, muscle_m = NA, viscera_m = NA, bones_m = NA, blubber_m = NA, muscle_lip = NA, muscle_pro = NA, visc_lip = NA, visc_pro = NA, bone_lip = NA, blubber_lip = NA, blubber_pro = NA))
    }

    # // ------------------------------------------------------------
    # // STRESSORS
    # // ------------------------------------------------------------
    
    stressors <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(gear_risk, entgl_cost, is_entgl, entgl_head, entgl_sev, entgl_d, entgl_start, entgl_end, is_entgl_calf, entgl_head_calf, entgl_sev_calf, entgl_d_calf, entgl_start_calf, entgl_end_calf, strike_risk, strike, strike_calf, noise_resp, noise_lvl, dB_thresh)])

    # // ------------------------------------------------------------
    # // FORAGING
    # // ------------------------------------------------------------
    
    feed <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(feed, preyconc, minprey, gape, feedspeed, captEff, impedance, daylight, feed_effort, eta_lwrBC, eta_upprBC, targetBC, cop_mass, cop_kJ, digestEff, metabEff_juv, metabEff_ad, E_cop)])

    # // ------------------------------------------------------------
    # // SOMATIC GROWTH
    # // ------------------------------------------------------------
    
    growth <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(delta_fat, delta_m, lip_anab, lip_catab)])
    
    if (cohorts[abb == k, id] == 5) {
    
      growth_calf <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(delta_fat_calf, delta_m_calf, lip_anab, lip_catab)])
    
    } else {
      
    growth_calf <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(delta_fat_calf = NA, delta_m_calf = NA, lip_anab = NA, lip_catab = NA)])
      
    }

    # // ------------------------------------------------------------
    # // ENERGY BALANCE
    # // ------------------------------------------------------------
    
    energy <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(E_tot, E_in, E_out, rmr, LC, scalar_LC, stroke, stroke_feed, E_growth, E_gest, fgrowth, placenta, hic, E_lac, perc_muscle, perc_viscera, perc_bones)])
    
    # // ------------------------------------------------------------
    # // ENERGY BALANCE (calves) + LACTATION
    # // ------------------------------------------------------------
    
    if (cohorts[abb == k, id] == 5) {
      
    energy_calf <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(cohort_calf, E_tot_calf, E_in_calf, E_out_calf, rmr_calf, LC_calf, scalar_LC, E_growth_calf)])
    
    lact <- cbind(const[, !c(cohort)], sim[[k]][whale %in% whaleID & day > 0, list(
      t_lac, assim, provision, zeta, milk_drop, eta_milk, mamm_M, mammEff, milk_rate, targetBC_calf, nursing, milk_lip, milk_pro
    )])
    
    } else {
      
      energy_calf <- cbind(const[, !c("cohort")], data.table::data.table(cohort_calf = NA, E_tot_calf = NA, E_in_calf = NA, E_out_calf = NA, rmr_calf = NA, LC_calf = NA, scalar_LC = NA, E_growth_calf = NA))
      
      lact <- cbind(const[, !c(cohort)], sim[[k]][whale %in% whaleID & day > 0, list(
        t_lac = NA, assim = NA, provision = NA, zeta = NA, milk_drop = NA, eta_milk = NA, 
        mamm_M = NA, mammEff = NA, milk_rate = NA, t_suckling = NA, targetBC_calf = NA, nursing = NA, milk_lip = NA, milk_pro = NA
      )])
      
      
    }

    # // ------------------------------------------------------------
    # // EXPORT
    # // ------------------------------------------------------------
    
    file.name <- stringr::str_replace(string = cohorts[abb == k, name], pattern = " \\(", replacement = "_")
    file.name <- stringr::str_replace(string = file.name, pattern = ", ", replacement = "_")
    file.name <-  stringr::str_replace(string = file.name, pattern = "\\)", replacement = "")
    
    sheet.list <- list("Attributes" = attrib,
                       "Attributes_calf" = attrib_calf,
                       "Movements" = xy,
                       "Behavior" = activ,
                       "Stressors" = stressors,
                       "Energy" = energy,
                       "Energy_calf" = energy_calf,
                       "Feeding" = feed,
                       "Growth" = growth,
                       "Growth_calf" = growth_calf,
                       "Gestation" = gest,
                       "Lactation" = lact)
    
    openxlsx::write.xlsx(x = sheet.list,
      file = paste0(prefix, "_", tolower(file.name), ".xlsx"), 
      asTable = TRUE,
      overwrite = overwrite,
      firstRow = TRUE, 
      tableStyle = "TableStyleMedium1",
      bandedRows = TRUE,
      withFilter = FALSE)

  } # End for loop

  cat("Done!")
  
}
