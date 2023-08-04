#' Overview of model results
#'
#' @import data.table
#' @param obj Object returned by run_model
#' @param rows Number of rows to print
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
#'
print.narwsim <- function(obj, ...){

  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default optional arguments
  rowID <- 1:5
  cohortID <- obj$param$cohortID
  whaleID <- 1 
  
  # Default values
  if(length(args) > 0){
    if("rowID" %in% names(args)) rowID <- args[["rowID"]]
    if("cohortID" %in% names(args)) cohortID <- args[["cohortID"]]
    if("whaleID" %in% names(args)) whaleID <- args[["whaleID"]]
  }

  if(!all(rowID %in% 1:365)) stop("<rows> must be between 1 and 365")

  cohorts <- obj$param$cohorts
  cohort.ab <- cohorts[id %in% cohortID, abb]
  # cohort.names <- cohorts[id %in% cohortID, name]

  for (k in cohortID) {
    
    if(k > 1) cat("\n\n\n")
    
    cat("=========================================================================================\n")
    cat(cohorts[id %in% k, name],"\n")
    cat("=========================================================================================\n")
    cat("\n")

    sim_dt <- obj[["sim"]][[cohorts[id==k, abb]]][whale %in% whaleID, ]
    sim_dt <- sim_dt[day %in% rowID, ]

    cat("--------------------------\n")
    cat("Locations \n")
    cat("--------------------------\n")
    cat("\n")
    print(sim_dt[, list(whale, day, date, easting, northing, region)])
    cat("\n")

    cat("--------------------------\n")
    cat("Attributes \n")
    cat("--------------------------\n")

    if (k %in% 4:5) {
      cat("\n +++ Adults +++\n\n")
    }
    print(sim_dt[, list(whale, day,
      cohort, gsl, seus, alive, age, bc, length, length_a, length_b, length_c,
      mass, leanmass, fatmass, mass_a, mass_b, mouth_r, mouth_a, mouth_w, abort, starve, died, date_died, p_died
    )])

    if (k == 5) {
      cat("\n +++ Calves +++\n\n")
      c.att <- sim_dt[, list(whale, day, cohort_calf,
        alive_calf, age_calf, bc_calf, length_calf,
        La_calf, Lb_calf, Lc_calf,
        mass_calf, leanmass_calf, fatmass_calf,
        ma_calf, mb_calf, mouth_r_calf,
        mouth_a_calf, mouth_w_calf, starve_calf, died_calf, date_died_calf 
      )]
      colnames(c.att) <- c("whale", "day", "cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c",
      "mass", "leanmass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", "mouth_w", "starve", "died", "date_died")
      print(c.att)
    }

    if (k == 4) {
      cat("\n +++ Fetus +++\n\n")
      print(sim_dt[, list(whale, day,
        fetus_l, fetus_m, birth_l, birth_m, muscle_m, viscera_m, bones_m,
        blubber_m, muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, blubber_lip, blubber_pro
      )])
    }

    cat("\n")

    cat("--------------------------\n")
    cat("Stressors \n")
    cat("--------------------------\n")
    cat("\n")
    print(sim_dt[, list(whale, day, gear_risk, is_entgl, entgl_head, entgl_sev, entgl_d, entgl_start, entgl_end, 
                        is_entgl_calf, entgl_head_calf, entgl_sev_calf, entgl_d_calf, entgl_start_calf, entgl_end_calf,
                        strike_risk, strike, strike_calf, noise_resp, noise_lvl, dB_thresh)])
    cat("\n")

    cat("--------------------------\n")
    cat("Activity budgets \n")
    cat("--------------------------\n")
    cat("\n")
    print(sim_dt[, list(whale, day, d_travel, swimspeed, glide, glide_feed, glide_echelon, t_travel, t_feed, t_rest_nurse)])
    cat("\n")

    cat("--------------------------\n")
    cat("Growth \n")
    cat("--------------------------\n")

    if (k == 5) {
      cat("\n +++ Adults +++\n\n")
      print(sim_dt[, list(whale, day, delta_fat, EDlip, EDpro, lip_anab, lip_catab, perc_muscle, perc_viscera, perc_bones)])
      cat("\n")
      cat("\n +++ Calves +++\n\n")
      print(sim_dt[, list(whale, day, delta_fat_calf, EDlip, EDpro, lip_anab, lip_catab, perc_muscle, perc_viscera, perc_bones)])
    } else {
      cat("\n")
      print(sim_dt[, list(whale, day, delta_fat, EDlip, EDpro, lip_anab, lip_catab, perc_muscle, perc_viscera, perc_bones)])
    }
    cat("\n")

    cat("--------------------------\n")
    cat("Energy balance \n")
    cat("--------------------------\n")
    
    if (k == 5) {
      cat("\n +++ Adults +++\n\n")
      print(sim_dt[, list(whale, day, E_tot, E_in, E_out)])
      cat("\n")
      cat("\n +++ Calves +++\n\n")
      print(sim_dt[, list(whale, day, E_tot_calf, E_in_calf, E_out_calf)])
    } else {
      cat("\n")
      print(sim_dt[, list(whale, day, E_tot, E_in, E_out)])
    }
    cat("\n")
    
    
    cat("--------------------------\n")
    cat("Energy intake \n")
    cat("--------------------------\n")

    if (k == 5) {
      cat("\n +++ Adults +++\n\n")
    } else{
      cat("\n")
    }
    print(sim_dt[, list(whale, day,
      feed, preyconc, minprey, gape, feedspeed, captEff, impedance, daylight,
      feed_effort, eta_lwrBC, eta_upprBC, targetBC, cop_mass, cop_kJ, digestEff, metabEff_juv, metabEff_ad, E_cop
    )])

    # Lactating females + calves
    if (k == 5) {
      cat("\n")
      cat("\n +++ Calves +++\n\n")
      print(sim_dt[, list(whale, day,
        t_lac, assim, provision, zeta, milk_drop, eta_milk, mamm_M, mammEff, 
        milk_rate, t_suckling, targetBC_calf, nursing, milk_lip, milk_pro
      )])
    }
    cat("\n")

    cat("--------------------------\n")
    cat("Energetic costs \n")
    cat("--------------------------\n")

    if (k == 5) {
      cat("\n +++ Adults +++\n\n")
    } else {
      cat("\n")
    }
    print(sim_dt[, list(whale, day, E_out, rmr, LC, scalar_LC, stroke, stroke_feed, E_growth)])
    
    if (k == 5) {
      cat("\n +++ Calves +++\n\n")
      print(sim_dt[, list(whale, day, E_out_calf, rmr_calf, LC_calf, E_growth_calf, delta_m_calf)])
    }
    
    
    # Pregnant females
    if (k == 4) {
      cat("\n")
      cat("--------------------------\n")
      cat("Gestation\n")
      cat("--------------------------\n")
      cat("\n")
      print(sim_dt[, list(whale, day,
        E_gest, fgrowth, placenta, hic, delta_fetus_m, delta_fetus_l,
        fetus_l, fetus_m, birth_l, birth_m, muscle_m, viscera_m, bones_m, blubber_m,
        muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, bone_pro, blubber_lip, blubber_pro,
        prop_mu, prop_visc, prop_bones, dens_blu, dens_mu, dens_visc, dens_bo
      )])
    }


    
    if (k == 5) {
      cat("\n")
      cat("--------------------------\n")
      cat("Lactation\n")
      cat("--------------------------\n")
      cat("\n")
      print(sim_dt[, list(whale, day, E_lac)])
    }

  }
  cat("\n\n")
}