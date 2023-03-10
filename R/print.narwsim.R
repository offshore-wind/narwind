#' Overview of model results
#'
#' @param obj Object returned by run_model
#' @param n.rows Number of rows to print
#'
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
print.narwsim <- function(obj,
                          rows = NULL,
                          whale.id = 1,
                          n.rows = 4) {

  # obj = m
  # rows = NULL
  # whale.id = 1
  # n.rows = 4

  if (is.null(rows)) rows <- 1:n.rows
  cohort.id <- obj$param$cohort.id
  cohort.ab <- obj$param$cohort.ab
  cohort.names <- obj$param$cohort.names

  for (k in seq_along(cohort.id)) {
    
    if (!cohort.id[k] %in% c(4, 5)) obj[["sim"]][[k]] <- list(obj[["sim"]][[k]])

    cat("\n++++++++++++++++++++++++++++++++++", cohort.names[k], "++++++++++++++++++++++++++++++++++\n\n")

    # purrr::walk(
    # .x = seq_along(obj[["sim"]][[k]]),
    # .f = ~ {
    sim_dt <- obj[["sim"]][[k]][whale == whale.id, ]
    sim_dt <- sim_dt[rows, ]

    # if (.x == 1) {
    cat("+++ Locations +++\n\n")
    print(sim_dt[, list(date, easting, northing, region)])
    cat("\n")

    cat("+++ Attributes +++\n\n")

    if (cohort.id[k] %in% 4:5) {
      cat(" --------------------------\n")
      cat(" Adults\n")
      cat(" --------------------------\n")
    }
    print(sim_dt[, list(
      cohort, alive, age, bc, length, length_a, length_b, length_c,
      mass, leanmass, fatmass, mass_a, mass_b, mouth_r, mouth_a, mouth_w
    )])

    if (cohort.id[k] == 5) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Calves\n")
      cat(" --------------------------\n")
      c.att <- sim_dt[, list(cohort_calf,
        alive_calf, age_calf, bc_calf, length_calf,
        La_calf, Lb_calf, Lc_calf,
        mass_calf, leanmass_calf, fatmass_calf,
        ma_calf, mb_calf, mouth_r_calf,
        mouth_a_calf, mouth_w_calf
      )]
      colnames(c.att) <- c("cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c",
      "mass", "leanmass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", "mouth_w")
      print(c.att)
    }

    if (cohort.id[k] == 4) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Fetus\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(
        fetus_l, fetus_m, birth_l, birth_m, muscle_m, viscera_m, bones_m,
        blubber_m, muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, blubber_lip, blubber_pro
      )])
    }

    cat("\n\n")

    cat("+++ Stressors +++\n\n")
    print(sim_dt[, list(is_entgl, entgl_head, severity, entgl_d, entgl_start, entgl_end, strike, noise, dB_thresh)])
    cat("\n")

    cat("+++ Activity budgets +++\n\n")
    print(sim_dt[, list(d_travel, swimspeed, t_travel, t_feed, t_nurse, t_rest, n_zero, t_sum, t_remain)])
    cat("\n")

    cat("+++ Growth +++\n\n")

    if (cohort.id[k] == 5) {
      cat(" --------------------------\n")
      cat(" Adults\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(E_tot, E_in, E_out, delta_fat, DE_lip, ED_lip, lip_anab, lip_catab)])
      cat("\n")
      cat(" --------------------------\n")
      cat(" Calves\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(E_tot_calf, E_in_calf, E_out_calf, delta_fat_calf, DE_lip, ED_lip, lip_anab, lip_catab)])
    } else {
      print(sim_dt[, list(E_tot, E_in, E_out, delta_fat, DE_lip, ED_lip, lip_anab, lip_catab)])
    }
    cat("\n")

    cat("+++ Energy intake +++\n\n")
    
    if (cohort.id[k] == 5) {
      cat(" --------------------------\n")
      cat(" Adults\n")
      cat(" --------------------------\n")
    }
    print(sim_dt[, list(
      feed, preyconc, minprey, gape, feedspeed, captEff, impedance, daylight,
      feed_effort, targetBC, cop_mass, cop_kJ, digestEff, metabEff, E_cop
    )])

    # Lactating females + calves
    if (cohort.id[k] == 5) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Calves\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(
        assim, provision, mamm_M, milk_rate, Dmilk, t_suckling, targetBC_calf, nursing, milk_lip, milk_pro, EDlip, EDpro
      )])
    }
    cat("\n")

    cat("+++ Energy costs +++\n")

    if (cohort.id[k] == 5) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Adults\n")
      cat(" --------------------------\n")
    }
    print(sim_dt[, list(rmr, LC, stroke, E_growth)])
    
    if (cohort.id[k] == 5) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Calves\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(rmr_calf, LC_calf, E_growth_calf, delta_m_calf)])
    }
    
    
    # Pregnant females
    if (cohort.id[k] == 4) {
      cat(" --------------------------\n")
      cat(" Gestation\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(
        E_gest, fgrowth, placenta, hic, delta_m,
        fetus_l, fetus_m, birth_l, birth_m, muscle_m, viscera_m, bones_m, blubber_m,
        muscle_lip, muscle_pro, visc_lip, visc_pro, bone_lip, blubber_lip, blubber_pro
      )])
    }

    if (cohort.id[k] == 5) {
      cat("\n")
      cat(" --------------------------\n")
      cat(" Lactation\n")
      cat(" --------------------------\n")
      print(sim_dt[, list(E_lac)])
    }

  }
  cat("\n")
}
# ) # End purrr:walk