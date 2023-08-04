#' Posterior
#'
#' Posterior sampling
#' @export
#' @import mgcv
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' animals <- run_model(10)
#' plot(animals)
#' }
#' 
expand.narwsim <- function(obj, 
                           nsamples = 10000,
                           xrange = c(0.05,0.6),
                           prop.sd = c(0.00001, 0.01),
                           gauss.approx = FALSE,
                           rw.start = c(1.25e-06, 0.05),
                           max.iter = 500,
                           ntune = 500,
                           plot = TRUE,
                           plot.n = 200,
                           ...){
  
  
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  seed <- NULL
  linewidth <- 0.65
  names(prop.sd) <- names(rw.start) <- c("surv", "bc")
  thin <- 1
  
  # Default values
  if(length(args) > 0) {
    if("seed" %in% names(args)) seed <- args[["seed"]] 
    if("linewidth" %in% names(args)) linewidth <- args[["linewidth"]] 
    if("thin" %in% names(args)) thin <- args[["thin"]] 
  }
  
  set.seed(seed)
  
  pred.x <- expand.grid(start_bc = seq(xrange[1], xrange[2], length.out = 100), cohort = cohorts$id)
  
  modlist <- obj$gam$fit
  
  if(gauss.approx){
    
    out <- purrr::map(.x = modlist,
               .f = ~{
                 
                 # Inverse link function
                 ilink <- family(.x)$linkinv
                 
                 # Map coefficients to fitted curves
                 Xp <- predict(.x, pred.x, type = "lpmatrix") 
                 
                 # Gaussian approximation to the posterior of the model coefficients
                 if(inherits(.x, "gam")){
                   beta <- coef(.x)
                   Vb <- vcov(.x) 
                 } else if(inherits(.x, "scam")){
                   beta <- as.numeric(tsgam:::coef.scam(.x))
                   Vb <- tsgam:::vcov.scam(.x)  
                 }

                 # Simulate nsamples replicate coefficient vectors from posterior
                 mrand <- MASS::mvrnorm(nsamples, beta, Vb) 
                 
                 obj$gam$post <- mrand
                 
                 opt <- matrix(data = NA, nrow = nsamples, ncol = nrow(pred.x))
                 for (i in seq_len(nsamples)) {
                   opt[i, ] <- ilink(Xp %*% mrand[i, ])
                 }
                 opt
               })
    
  } else {
    
    # Bayesian MCMC
    
    # See ?mgcv::gam.mh
    # The function reports the acceptance rate of the two types of step.
    # If the random walk acceptance probability is higher than a quarter 
    # then rw.step should probably be increased. 
    # Similarly if the acceptance rate is too low, it should be decreased. 
    # The random walk steps can be turned off altogether (see above), 
    # but it is important to check the chains for stuck sections if this is done.
    
    #'------------------------------------------------------------
    # Tuning of scale parameter for random walk proposals ----
    #'------------------------------------------------------------
    # Note: The below requires a suitable starting value (rw.start)
    # which is dependent on the data and model used
    
    post.samples <- purrr::map(.x = names(modlist), .f = ~{
      
      if(.x == "surv"){
        cat("\n\n**** SURVIVAL model ****\n")
      } else if(.x == "bc"){
        cat("\n\n**** HEALTH model ****\n")
      }
      
      
      if(inherits(modlist[[.x]], "gam")){
        post.init <- gam_mh(modlist[[.x]], ns = ntune, burn = 0.5*ntune, rw.scale = rw.start[.x], pb = FALSE, verbose = F, thin = thin)
      } else if(inherits(modlist[[.x]], "scam")){
        post.init <- scam_mh(modlist[[.x]], ns = ntune, burn = 0.5*ntune, rw.scale = rw.start[.x], pb = FALSE, verbose = F, thin = thin)
      }

      counter <- 0
      rw <- rw.start[.x]
      df <- matrix(NA, nrow = max.iter, ncol = 3) 
      colnames(df) <- c("iter", "rw", "accept")
      
      cat("\nTuning the Metropolis Hastings sampler ...\n")
      cat("------------------------------------------\n")

      while (post.init$rw.accept > 0.25 | post.init$rw.accept < 0.24 | counter <= max.iter) {
        if(post.init$rw.accept > 0.25){
          rw <- rw + rtnorm(0, prop.sd[.x], 0, Inf)
        } else if(post.init$rw.accept < 0.24){
          rw <- rw - rtnorm(0, prop.sd[.x], 0, Inf)
        } else {
          break
        }
        cat("Iteration ", stringr::str_pad(counter, 3, pad = "0"), " | RW acceptance = ", post.init$rw.accept, "\n")
        
        if (inherits(modlist[[.x]], "gam")) {
          post.init <- gam_mh(modlist[[.x]],
            ns = ntune,
            burn = 0.5 * ntune,
            rw.scale = rw,
            pb = FALSE,
            verbose = FALSE,
            thin = thin
          )
        } else if (inherits(modlist[[.x]], "scam")) {
          post.init <- scam_mh(modlist[[.x]],
            ns = ntune,
            burn = 0.5 * ntune,
            rw.scale = rw,
            pb = FALSE,
            verbose = FALSE,
            thin = thin
          )
        }
        
        df[counter+1,] <- c(counter, rw, post.init$rw.accept)
        counter <- counter + 1
      }
      
      df <- df[complete.cases(df),]
      
      cat("------------------------------------------\n")
      cat("\nExtracting posterior samples ...\n")
      
      if (inherits(modlist[[.x]], "gam")) {
        post.samples <- gam_mh(modlist[[.x]],
          ns = nsamples,
          burn = 0.5 * nsamples,
          rw.scale = df[nrow(df), 2],
          pb = TRUE,
          verbose = FALSE,
          thin = thin
        )
      } else if (inherits(modlist[[.x]], "scam")) {
        post.samples <- scam_mh(modlist[[.x]],
          ns = nsamples,
          burn = 0.5 * nsamples,
          rw.scale = df[nrow(df), 2],
          pb = TRUE,
          verbose = FALSE,
          thin = thin
        )
      }
      
      cat("\nRandom walk proposal acceptance:", post.samples$rw.accept, "\n")
      post.samples$bs
      
    }) |> purrr::set_names(nm = names(modlist))
    
    obj$gam$post <- post.samples
    
    #'------------------------------------------------------------
    # Posterior sampling ----
    #'------------------------------------------------------------
    
    out <- purrr::map(.x = names(modlist),
                      .f = ~{
                        
                        # Inverse link function
                        ilink <- family(modlist[[.x]])$linkinv
                        
                        # Map coefficients to fitted curves
                        Xp <- predict(modlist[[.x]], pred.x, type = "lpmatrix") 
                        
                        opt <- matrix(data = NA, nrow = nsamples, ncol = nrow(pred.x))
                        for (i in seq_len(nsamples)) {
                          opt[i, ] <- ilink(Xp %*% post.samples[[.x]][i, ])
                        }
                        opt
                      }) |> purrr::set_names(nm = names(modlist))
    
  }
  
  #'------------------------------------------------------------
  # Minimum body condition required for gestation ----
  #'------------------------------------------------------------
 
  cat("\n\n**** GESTATION model ****\n")
  
  post.mbc <- gam_mh(gam_gest$gam, ns = nsamples, burn = 0.5*nsamples, rw.scale = 0.625, pb = TRUE, verbose = FALSE, thin = 1)
  
  cat("\nRandom walk proposal acceptance:", post.mbc$rw.accept, "\n")
  ilink.mbc <- family(gam_gest$gam)$linkinv
  xpred <- seq(15000, 45000, length.out = 100)
  Xp <- predict(gam_gest$gam, data.frame(mass = xpred), type = "lpmatrix")
  opt <- matrix(data = NA, nrow = nsamples, ncol = length(xpred))
  for (i in seq_len(nsamples)) {
    opt[i, ] <- ilink.mbc(Xp %*% post.mbc$bs[i, ])
  }
  cat("\n")
  
  obj$gam$post$mbc <- post.mbc$bs
  
  #'------------------------------------------------------------
  # PLOTTING ----
  #'------------------------------------------------------------
  
  if(plot){

    cat("Plotting ...\n")
    pred.df <- expand.grid(start_bc = seq(xrange[1], xrange[2], length.out = 100), cohort = cohorts$id, mcmc = seq_len(nsamples))
    pred.df$pred_surv <- as.numeric(t(out$surv))
    pred.df$pred_bc <- as.numeric(t(out$bc))
    
    p <- purrr::map(
      .x = names(modlist),
      .f = ~ {
        mdf <- cbind(pred.x, pred = predict(modlist[[.x]], pred.x, "response"))
        
        sdf <- pred.df[pred.df$mcmc %in% sample(nsamples, plot.n, replace = FALSE), ] |> 
          dplyr::select(start_bc, cohort, mcmc, paste0("pred_", .x))
        lookup <- c(pred = paste0("pred_", .x))
        sdf <- dplyr::rename(sdf, all_of(lookup))
        
        facet.names <- cohorts$name
        names(facet.names) <- cohorts$id
        
        ggplot2::ggplot(data = sdf, aes(x = start_bc, y = pred)) +
          ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
          ggplot2::geom_line(data = mdf, linewidth = linewidth) +
          ggplot2::facet_wrap(vars(cohort), labeller = ggplot2::labeller(cohort = facet.names)) +
          ggplot2::scale_y_continuous(limits = c(0,1)) +
          theme_narw() +
          xlab("Starting body condition (%)") +
          ylab(dplyr::case_when(
            .x == "surv" ~ "p(survival)",
            .x == "bc" ~ "Final body condition (%)",
            .default = ""
          ))
      }
    ) |> purrr::set_names(nm = names(modlist))
    
    # ggplot2::ggplot(data = pred.df[pred.df$mcmc %in% sample(nsamples, plot.n, replace = FALSE),], aes(x = start_bc, y = pred_bc)) +
    #   ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.5) +
    #   # ggplot2::geom_line(data = cbind(pred.x, pred_bc = predict(modlist$bc, pred.x, "response")), linewidth = 1) +
    #   ggplot2::facet_wrap(vars(cohort)) + 
    #   theme_narw()
    
    purrr::walk(.x = p, .f = ~print(.x))
    
    mbc.df <- expand.grid(mass = xpred, mcmc = seq_len(nsamples))
    mbc.df$min_bc <- as.numeric(t(opt))
    mbc.df <- mbc.df[mbc.df$mcmc %in% sample(nsamples, plot.n, replace = FALSE), ]
    
    pdf <- data.frame(mass = xpred, min_bc = predict(gam_gest$gam, data.frame(mass = xpred), "response"))
    
    gp <- ggplot2::ggplot(data = mbc.df, aes(x = mass, y = min_bc)) +
      ggplot2::geom_line(aes(group = mcmc), colour = "grey", alpha = 0.2) +
      ggplot2::geom_line(data = pdf, linewidth = linewidth) +
      theme_narw() +
      xlab("Total mass (kg)") +
      ylab("Min body condition (%)")
    
    print(gp)
    
  }
  cat("Done!\n")
  return(obj)
}
