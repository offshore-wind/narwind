#' Forecast NARW abundance
#'
#' @param obj Object returned by narw()
#' @param n Number of replicate projections
#' @param yrs Number of years for the projection
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
predict.narwsim <- function(obj,
                            n = 1000,
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
  
  p1 <- make_plot(narw.df, narw.conf)
  p2 <- make_plot(tot.df, tot.conf)
  
  print(p1)
  print(p2)
  
  # Find 95% confidence intervals on final population size
  ci <- tot.df[year == current.yr + yrs,  quantile(N, probs = c(0.025, 0.975))]

  return(ci)
  
}