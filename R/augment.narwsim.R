#' Posterior simulation from fitted survival and health GAMs
#'
#' Simulates random draws from the posterior distribution of GAM coefficients using a Metropolis Hastings algorithm.
#' 
#' @param obj An object of class \code{narwsim}, as returned by \code{\link{narw}}.
#' @param n Integer. Defaults to \code{100,000}. Number of posterior samples to draw.
#' @param bc.range Numeric vector of length 2. Range (minimum, maximum) of body conditions used for model predictions.
#' @param ... Additional arguments passed to \code{\link[mgcv]{gam.mh}}.
#' @note In many cases, a Gaussian approximation to the posterior of the model coefficients is largely accurate, and samples generated from it can be treated as samples from the posterior for the coefficients. In other cases, however, this approximation can become poor and it may be useful to simulate from the posterior using a Metropolis Hastings sampler. A simple approach to this alternates fixed proposals, based on the Gaussian approximation to the posterior, with random walk proposals, based on a shrunken version of the approximate posterior covariance matrix. Fixed proposals often promote rapid mixing, while the random walk component ensures that the chain does not become stuck in regions for which the fixed Gaussian proposal density is much lower than the posterior density. The \code{\link[mgcv]{gam.mh}} function implements this, and is harnessed here.
#'
#' The function reports the acceptance rate of the two samplers. As a general rule, the random walk acceptance probability should not be higher than ca. 0.25. If it is, the \code{rw.scale} argument can be used to increase the factor by which the posterior covariance matrix is scaled.
#' @return The original \code{narwsim} object, augmented with posterior samples (list element named \code{post}).
#' @import mgcv
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(10000)
#' m <- augment(m)
#' }

augment.narwsim <- function(obj,
                           n = 100000,
                           bc.range = c(0.05, find_maxBC()),
                           ...){
  
  # Function ellipsis
  args <- list(...)
  
  if("thin" %in% names(args)) thin <- args[["thin"]] else thin <- 10
  if("burn" %in% names(args)) burn <- args[["burn"]] else burn <- 0.5*n
  if("rw.scale" %in% names(args)) rw.scale <- args[["rw.scale"]] else rw.scale <- 0.25 # Default value from mgcv::gam.mh
  
  # Function checks
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")

  # Prediction data
  pred.x <- expand.grid(start_bc = seq(bc.range[1], bc.range[2], length.out = 100), cohort = cohorts$id)
  
  # Total number of iterations
  n.tot <- n/thin
  
  # Extract GAM objects
  modlist <- obj$gam$fit
    
  ##'...............................
  ## Survival / health ----
  ##''...............................
  
  out.survbc <- purrr::map(.x = names(modlist), .f = ~{
    cat(ifelse(.x == "surv", "SURVIVAL MODEL", "\nHEALTH MODEL"),"\n")
    out.mh <- gam.mh(modlist[[.x]], ns = n, burn = burn, thin = thin, rw.scale = rw.scale)
    out.mh$bs
  }) |> purrr::set_names(nm = names(modlist))
  
  # Predictions 
  preds.survbc <- purrr::map(
    .x = names(modlist),
    .f = ~ {
      
      # Inverse link function
      ilink <- family(modlist[[.x]])$linkinv
      
      # Map coefficients to fitted curves
      Xp <- predict(modlist[[.x]], pred.x, type = "lpmatrix")
      
      opt <- matrix(data = NA, nrow = n.tot, ncol = nrow(pred.x))
      
      for (i in seq_len(n.tot)) {
        opt[i, ] <- ilink(Xp %*% out.survbc[[.x]][i, ])
      }
      opt
    }
  ) |> purrr::set_names(nm = names(modlist))
  
  
  #'------------------------------------------------------------
  # Gestation ----
  #'------------------------------------------------------------

  # The below is commented out as we only consider the mean prediction 
  # in the population projections. Uncomment to account for uncertainty in the
  # relationship between total body mass and minimum body condition required
  # for gestation.
  # ilink.mbc <- family(gam_gest)$linkinv
  # xpred <- seq(15000, 45000, length.out = 100)
  # Xp <- predict(gam_gest, data.frame(mass = xpred), type = "lpmatrix")
  # preds.mbc <- matrix(data = NA, nrow = n.tot, ncol = length(xpred))
  # 
  # cat("\nGESTATION MODEL\n")
  # out.minbc <- gam.mh(gam_gest, 
  #                     ns = n, 
  #                     burn = burn, 
  #                     rw.scale = rw.scale[3], 
  #                     thin = thin)
  # 
  # out.minbc <- out.minbc$bs
   
  # for (i in seq_len(n.tot)) {
  #   preds.mbc[i, ] <- ilink.mbc(Xp %*% out.minbc[i, ])
  # }
  
  # Save results
  obj$post <- list(nsamples = n.tot,
                   samples = out.survbc,
                   # samples = append(out.survbc, list("mbc" = out.minbc)),
                   # preds = list(survbc = preds.survbc, mbc = preds.mbc),
                   preds = list(survbc = preds.survbc),
                   bc.range = bc.range,
                   pred.x = pred.x,
                   # mbc.x = xpred,
                   burn = burn,
                   thin = thin)
  
  cat("Done!\n")
  return(obj)
}
