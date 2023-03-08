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

print.narwsim <- function(obj, rows = NULL, whale.id = 1, n.rows = 4) {
  
  # obj = m
  # rows = NULL
  # whale.id = 1
  # n.rows = 4
  
  if(is.null(rows)) rows = 1:n.rows
  
  for(k in seq_along(obj$param$cohort.id)){
    
    if(!obj$param$cohort.id[k] %in% c(4,5)) obj[["sim"]][[k]] <- list(obj[["sim"]][[k]])
    
    cat("\n++++++++++++++++++++++++++++++++++", obj$param$cohort.names[k], "++++++++++++++++++++++++++++++++++\n\n")
    
    purrr::walk(
      
      .x = seq_along(obj[["sim"]][[k]]),
      .f = ~ {
        
        if(obj$param$cohort.id[k] == 5){
          cat("--------------------------\n")
          cat(c("Adults", "Calves")[.x], "\n")
          cat("--------------------------\n")
        }
        
        if(obj$param$cohort.id[k] == 4){
          cat("--------------------------\n")
          cat(c("Adults", "Fetus")[.x], "\n")
          cat("--------------------------\n")
        }
        
        if(length(rows) > 1){
          
          mlocs <- obj[["locs"]][[k]][rows, , whale.id]
          m_all <- obj[["sim"]][[k]][[.x]][rows, , whale.id]
          
         
        } else {
          
          rn <- rownames(obj[["locs"]][[k]][, , whale.id])[rows]
          cn <- colnames(obj[["locs"]][[k]][, , whale.id])
          mlocs <- matrix(obj[["locs"]][[k]][rows, , whale.id], nrow = length(rows), dimnames = list(rn, cn))
          
          cn <- colnames(obj[["sim"]][[k]][[.x]][, , whale.id])
          m_all <- matrix(obj[["sim"]][[k]][[.x]][rows, ,whale.id], nrow = length(rows), dimnames = list(rn, cn))

        }
        
        if(.x == 1){
        
        m_attrib <- m_all[, c("cohort", "alive", "age", "bc", "length", "length_a", "length_b", "length_c", 
                              "mass", "leanmass", "fatmass", "mass_a", "mass_b", "mouth_r", "mouth_a", 
                              "mouth_w")]
        
        m_stress <- m_all[, c("is_entgl", "entgl_head", "severity", "entgl_d", "entgl_start", "entgl_end", "strike", "noise", "dB_thresh")]
        
        if(obj$param$cohort.id[k] == 5){
          m_energy <- m_all[, c("E_tot_calves", "E_in_calves", "E_out_calves", "delta_fat_calves", "DE_lip", "ED_lip", "lip_anab", "lip_catab")]
        } else {
          m_energy <- m_all[, c("E_tot", "E_in", "E_out", "delta_fat", "DE_lip", "ED_lip", "lip_anab", "lip_catab")]
        }
        
        m_intake <- m_all[, c("feed", "preyconc", "minprey", "gape", "feedspeed", "captEff", "impedance", "daylight",
                              "feed_effort", "targetBC", "cop_mass", "cop_kJ", "digestEff", "metabEff", "E_cop")]

        }
        
        # Lactating females + calves
        if(obj$param$cohort.id[k] == 5 & .x == 2){
          m_lac <-  m_all[, c("assim", "provision", "mamm_M", "milk_rate", "Dmilk", "t_suckling",
                              "targetBC", "nursing", "milk_lip", "milk_pro", "ED_lip", "ED_pro")]
        }
        

        # Pregnant females
        if(obj$param$cohort.id[k] == 4 & .x == 2){
          m_gest <-  m_all[, c("fetus_l", "fetus_m", "birth_l", "birth_m", "muscle_m", "viscera_m", "bones_m", "blubber_m",
                               "muscle_lip", "muscle_pro", "visc_lip", "visc_pro", "bone_lip", "blubber_lip", "blubber_pro")]
          m_costs <- m_all[, c("rmr", "LC", "stroke", "E_growth", "E_gest", "fgrowth", "placenta", "hic", "E_lac", "delta_m")]
        } else {
          m_costs <- m_all[, c("rmr", "LC", "stroke", "E_growth")]
        }
        
       
        m_activ <- m_all[, c("d_travel", "swimspeed", "t_travel", "t_feed", "t_nurse", "t_rest", "n_zero", "t_sum", "t_remain")]
        
        if(.x == 1){
        cat("+++ Locations +++\n\n")
        print(mlocs)
        cat("\n")
        
        cat("+++ Attributes +++\n\n")
        print(m_attrib)
        cat("\n")
        
        cat("+++ Stressors +++\n\n")
        print(m_stress)
        cat("\n")
        
        cat("+++ Energy balance +++\n\n")
        print(m_energy)
        cat("\n")
        
        cat("+++ Energy intake +++\n\n")
        print(m_intake)
        cat("\n")
        
        cat("+++ Energy costs +++\n\n")
        print(m_costs)
        cat("\n")
        
        cat("+++ Activity budgets +++\n\n")
        print(m_activ)
        cat("\n")
        }
        
        if(obj$param$cohort.id[k] == 4 & .x == 2){
          cat("+++ Gestation +++\n\n")
          print(m_gest)
        }
        
        if(obj$param$cohort.id[k] == 5 & .x == 2){
          cat("+++ Lactation +++\n\n")
          print(m_lac)
        }
        
        cat("\n")
      }
    )
  }
}