#' Export simulation data
#'
#' Writes outputs from the individual-based model to an .xlsx file on disk.
#' @param obj An object of class \code{narwsim}, as returned by \code{\link{narw}}.
#' @param prefix Character string. Output file name. Defaults to \code{"narwsim"}.
#' @note This function 
#' @inheritParams write.narwproj
#' @return One Excel file per population cohort. It is recommended to only use this function for short simulation runs (\code{nsim = 100} or less, or to use the \code{cohort} and \code{whale} arguments to extract data for specific individuals/cohorts. Each output file contains multiple sheets capturing all simulation parameters, as follows:
#' \itemize{
#' \item \code{Attributes}: day, date, month, unique whale ID, cohort ID, alive (yes/no), age, starve (yes/no), mortality from other sources (yes/no), probability of death from other sources, body condition, body length, coefficients of the length-at-age function, total mass, lean mass, fat mass, coefficients of the mass-at-length function, mouth radius, mouth opening angle, mouth width, migration to SEUS (yes/no), migration to GSL (yes/no)
#' \item \code{Attributes_calf} same as \code{Attributes}, for calves
#' \item \code{Movements}: easting,	northing, region, distance traveled, swim speed
#' \item \code{Behavior}: proportion of time spent gliding, proportion of time spent gliding during foraging, proportion of time spent gliding while swimming in echelon position, time spent traveling, time spent feeding, time spent resting/nursing
#' \item \code{Stressors}: Probability of entanglement, daily energetic cost of entanglement, entanglement event (yes/no), site of attachement (head/not head), entanglement serverity, entanglement duration, day of entanglement start, day of entanglement end, vessel strike	risk,	incidence of strike, incidence of behavioral response to noise, received sound level, threshold of response
#' \item \code{Energy}: Daily energy balance, daily energy intake, daily energy expenditure, resting metabolic rate, locomotory costs, scalar applied to locomotory costs, stroke rate, stroke rate during foraging, energetic cost of somatic growth, total energetic cost of gestation, energetic costs of fetus growth and placental maintenance, heat increment of gestation, energetic cost of lactation, relative proportions of bones, viscera and muscles in lean mass
#' \item \code{Energy_calf} same as \code{Energy}, for calves (where relevant)
#' \item \code{Feeding}: incidence of foraging, prey concentration, minimum prey threshold that triggers foraging, size of the mouth gape, swimming speed while foraging, prey capture efficiency, reduction in the mouth gape during a head entanglement, feeding effort, feed_effort, coefficients of the feeding effort curve, target body condition, average mass of copepods, energy content of copepods, digestive efficiency, metabolizing efficiency	for juveniles and adults (heat increment of feeding)
#' \item \code{Growth}: change in fat mass, daily change in lean mass, anabolism efficiency, catabolism efficiency
#' \item \code{Growth_calf} same as \code{Growth}, for calves
#' \item \code{Gestation}: incidence of abortion, fetus length, fetus mass, change in fetus mass, change in fetus length, predicted birth length, predicted birth mass, mass of muscles, viscera, bones, and blubber in fetal tissues, proportions of lipid and protein in fetal muscles, viscera, bones, bones, and blubber.
#' \item \code{Lactation}: duration of lactation, milk assimilation ratem milk provision rate, factor defining the nonlinearity between milk supply and body condition of the mother, age at which milk consumption starts to decrease, factor defining the nonlinearity between milk assimilation and calf age, mass of mammary glands, mammary efficiency, milk production rate, time spent nursing, target body condition (calf), proportions of protein and lipids in milk
#' 
#' }
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' write(m)}

write.narwsim <- function(obj, 
                          prefix = "narwsim",
                          ...){
  
  # Function ellipsis –– optional arguments
  args <- list(...)

  # Default values
  if("cohort" %in% names(args)) cohort <- args[["cohort"]] else cohort <- obj$param$cohort
  if("whale" %in% names(args)) whaleID <- args[["whale"]] else whaleID <- seq_len(obj$param$nsim)

  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")

  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohort, abb]
  cohort.names <- cohorts[id %in% cohort, name]
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
    
    attrib <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(alive, age, p_starve, starve, died, p_surv, bc, length, length_a, length_b, length_c, mass, leanmass, fatmass, mass_a, mass_b, gape, mouth_r, mouth_a, mouth_w, gsl, seus)])
    
    # // ------------------------------------------------------------
    # // ATTRIBUTES (calves)
    # // ------------------------------------------------------------
    
    if (cohorts[abb == k, id] == 5) {
      
      attrib_calf <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(cohort_calf, born, dob, pbirth, alive_calf, age_calf, p_starve_calf, starve_calf, bc_calf, length_calf, La_calf, Lb_calf, Lc_calf, mass_calf, leanmass_calf, fatmass_calf, ma_calf, mb_calf, mouth_r_calf, mouth_a_calf, mouth_w_calf, died_calf, date_died_calf, p_surv_calf)])
      
    } else {
      
      attrib_calf <- cbind(const[, !c("cohort")], data.table::data.table(cohort_calf = NA, born = NA, dob = NA, pbirth = NA, alive_calf = NA, age_calf = NA, p_starve_calf = NA, starve_calf = NA, bc_calf = NA, length_calf = NA, La_calf = NA, Lb_calf = NA, Lc_calf = NA, mass_calf = NA, leanmass_calf = NA, fatmass_calf = NA, ma_calf = NA, mb_calf = NA, mouth_r_calf = NA, mouth_a_calf = NA, mouth_w_calf = NA, died_calf = NA, date_died_calf = NA, p_surv_calf = NA))
      
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
    
    feed <- cbind(const, sim[[k]][whale %in% whaleID & day > 0, list(feed, preyconc, minprey, gape, feedspeed, captEff, impedance,  feed_effort, eta_lwrBC, eta_upprBC, targetBC, cop_mass, cop_kJ, digestEff, metabEff_juv, metabEff_ad, E_cop)])

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
      
      lact <- cbind(const[, !c("cohort")], sim[[k]][whale %in% whaleID & day > 0, list(
        t_lac = NA, assim = NA, provision = NA, zeta = NA, milk_drop = NA, eta_milk = NA, 
        mamm_M = NA, mammEff = NA, milk_rate = NA, targetBC_calf = NA, nursing = NA, milk_lip = NA, milk_pro = NA
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
      file = paste0(tolower(prefix), "_", file.name, ifelse(!is.null(obj$param$label), paste0("_", obj$param$label)), ".xlsx"), 
      asTable = TRUE,
      firstRow = TRUE, 
      tableStyle = "TableStyleMedium1",
      bandedRows = TRUE,
      withFilter = FALSE,
      ...)

  } # End for loop

  cat("Done!")
  
}
