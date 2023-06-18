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

  if(!"narwsim" %in% class(obj)) stop("Input must be of class <narwsim>")

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

    xy <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(easting, northing, region, d_travel, swimspeed)])

    attrib <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(alive, age, starve, bc, length, length_a, length_b, length_c, mass, leanmass, fatmass, mass_a, mass_b, mouth_r, mouth_a, mouth_w, gsl, seus)])
    
    if (cohorts[abb == k, id] == 5) {
      attrib_calf <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(cohort_calf, born, alive_calf, age_calf, starve_calf, bc_calf, length_calf, La_calf, Lb_calf, Lc_calf, mass_calf, leanmass_calf, fatmass_calf, ma_calf, mb_calf, mouth_r_calf, mouth_a_calf, mouth_w_calf)])
    } else {
      attrib_calf <- cbind(const[, !c("cohort")], data.table::data.table(cohort_calf = NA, born = NA, alive_calf = NA, age_calf = NA, starve_calf = NA, bc_calf = NA, length_calf = NA, La_calf = NA, Lb_calf = NA, Lc_calf = NA, mass_calf = NA, leanmass_calf = NA, fatmass_calf = NA, ma_calf = NA, mb_calf = NA, mouth_r_calf = NA, mouth_a_calf = NA, mouth_w_calf = NA))
    }
    
    activ <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(t_travel, t_feed, t_nurse, t_rest, t_sum, t_remain)])

    if (cohorts[abb == k, id] == 4) {
      gest <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(abort, fetus_l, fetus_m, delta_fetus_m, delta_fetus_l, birth_l, birth_m, muscle_m, viscera_m, bones_m, blubber_m, muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, blubber_lip, blubber_pro)])
    } else {
      gest <- cbind(const, data.table::data.table(abort = NA, fetus_l = NA, fetus_m = NA, delta_fetus_m = NA, delta_fetus_l = NA, birth_l = NA, birth_m = NA, muscle_m = NA, viscera_m = NA, bones_m = NA, blubber_m = NA, muscle_lip = NA, muscle_pro = NA, visc_lip = NA, visc_pro = NA, bone_lip = NA, blubber_lip = NA, blubber_pro = NA))
    }

    stressors <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(gear_risk, is_entgl, entgl_head, entgl_sev, entgl_d, entgl_start, entgl_end, is_entgl_calf, entgl_head_calf, entgl_sev_calf, entgl_d_calf, entgl_start_calf, entgl_end_calf, strike_risk, strike, strike_calf, noise_resp, noise_lvl, dB_thresh)])

    feed <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(feed, preyconc, minprey, gape, feedspeed, captEff, impedance, daylight, feed_effort, targetBC, cop_mass, cop_kJ, digestEff, metabEff_juv, metabEff_ad, E_cop)])

    growth <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(delta_fat, ED_lip, lip_anab, lip_catab)])

    energy <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(E_tot, E_in, E_out, rmr, LC, scalar_LC, stroke, E_growth, E_gest, fgrowth, placenta, hic, E_lac, delta_m)])
    
    if (cohorts[abb == k, id] == 5) {
    energy_calf <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(cohort_calf, E_tot_calf, E_in_calf, E_out_calf, rmr_calf, LC_calf, scalar_LC, E_growth_calf, delta_m_calf)])
    } else {
      energy_calf <- cbind(const[, !c("cohort")], data.table::data.table(cohort_calf = NA, E_tot_calf = NA, E_in_calf = NA, E_out_calf = NA, rmr_calf = NA, LC_calf = NA, scalar_LC = NA, E_growth_calf = NA, delta_m_calf = NA))
    }
    
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
                       "Gestation" = gest)
    
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
