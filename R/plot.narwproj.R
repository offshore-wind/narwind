#' Plot population trends
#'
#' Takes one or more \code{narwproj} object(s) and creates plots of the estimated population trajectories.
#'
#' @param ... One or more objects of class \code{narwproj}, as returned by \link{predict.narwsim}.
#' @param diff Logical. If \code{TRUE}, returns plots of the differences between pairs of projections. Defaults to \code{FALSE}.
#' @param interval Logical. If \code{TRUE}, percentile confidence intervals are shown on the plots. Defaults to \code{TRUE}.
#' @param shade Logical. If \code{TRUE}, confidence intervals are displayed as colored ribbons. Otherwise, confidence bounds are shown as dotted lines. Defaults to \code{TRUE}.
#' @param cohort Logical. If \code{TRUE}, separate plots are returned for each cohort. If \code{FALSE}, a single plot for the whole population is shown. Defaults to \code{FALSE}.
#' @param noaa Logical. If \code{TRUE}, the population trajectory predicted as part of NOAA's population viability analysis (Runge et al., 2023) is also plotted. Defaults to \code{FALSE}.
#' @param full Logical. If \code{TRUE}, the trend plot includes historical population size estimates (dating back to 2000).
#' @param timeline Logical. If \code{TRUE}, the plot includes timeline(s) of offshore wind development.
#' @param timeline.y Numeric. Position of the timeline on the y-axis.
#' @param scales Character. Defines whether axis scales should be constant or vary across plots. Can be one of "fixed" or "free" (the default).
#' @param nudge.y Numeric. Vertical adjustment to nudge timeline labels by.
#' @param ncol Integer. Number of columns for the plot layout when \code{cohort = TRUE}. Defaults to \code{3}.
#' @param nx Integer. Desired number of x-axis intervals. Non-integer values are rounded down. Defaults to \code{5}. See \link[base]{pretty} for details. 
#' @param ny Integer. Desired number of y-axis intervals. Non-integer values are rounded down. Defaults to \code{5}. See \link[base]{pretty} for details.
#' @import ggplot2
#' @references Runge MC, Linden DW, Hostetler JA, Borggaard DL, Garrison LP, Knowlton AR, Lesage V, Williams R, Pace III RM (2023). A management-focused population viability analysis for North Atlantic right whales. US Dept Commer Northeast Fish Sci Cent Tech Memo 307, 93 p. Available at \url{https://www.fisheries.noaa.gov/s3/2023-10/TM307-508-1-.pdf}.
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' m <- augment(m)
#' p <- predict(m)
#' plot(p)
#' }
plot.narwproj <- function(...,
                          diff = FALSE,
                          interval = TRUE,
                          shade = TRUE,
                          cohort = FALSE,
                          noaa = FALSE,
                          full = FALSE,
                          timeline = FALSE,
                          timeline.y = 0,
                          scales = "free",
                          nudge.y = 20,
                          ncol = 2,
                          nx = 5,
                          ny = 5){
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  if(cohort) timeline <- FALSE
  if(noaa) full <- TRUE
  start.proj <- 2019
  
  # Identify <narwproj> objects 
  which.obj <- which(purrr::map_lgl(.x = args, .f = ~inherits(.x, "narwproj")))
  if(length(which.obj)==0) stop("No object of class <narwproj> found.")
  obj <- args[which.obj]
  
  if("vignette" %in% names(args)) vignette <- args[["vignette"]] else vignette <- FALSE
  
    # Function to remove some of the axis labels in ggplot faceted plots
    gtable_filter_remove <- function (x, name, trim = TRUE){
      matches <- !(x$layout$name %in% name)
      x$layout <- x$layout[matches, , drop = FALSE]
      x$grobs <- x$grobs[matches]
      if (trim) 
        x <- gtable::gtable_trim(x)
      x
    }
    
    # Color scale for plotting
    if ("colour" %in% names(args)) {
      colour <- args[["colour"]]
    } else {
      colour <- c("#BF0B98", "#E0B200", "#289DE0", "#0098A0", "#FF9B19FA", "#E6124B")
    }
    
    if (length(obj) > length(colour)) colour <- c(colour, 1:(length(obj) - length(colour)))
    
    # Define years for the x-axis
    if(!full){
      start.year <- obj[[1]]$param$start.year
    } else {
      start.year <- 2000
    }
    end.year <- obj[[1]]$param$start.year + obj[[1]]$param$yrs
    
    # Extract the trend data
    proj.labels <- purrr::map_chr(.x = obj, .f = ~.x$param$label)
    if(all(proj.labels == "")) proj.labels <- paste0("Projection ", stringr::str_pad(seq_along(obj), width = 2, pad = "0"))
    
    # Calculate difference between projections, if desired
    if(diff){
      
      n <- obj[[1]]$param$n
      allcombs <- combn(seq_along(obj), m = 2)
      
      # For all combinations
      difftbl <- purrr::map(.x = 1:ncol(allcombs), .f = ~{
        
        # All combinations of individual projections
        ndiff <- unique(expand.grid(1:n, 1:n))
        
        diffproj <- matrix(data = NA, nrow = nrow(ndiff), ncol = end.year-start.year + 1)
        
        pb <- progress::progress_bar$new(total = nrow(ndiff), width = 80, format = "Calculating differences [:bar] :percent eta: :eta", clear = TRUE)
        
        for(i in seq_len(nrow(ndiff))){
          pb$tick()
          diffproj[i, ] <- obj[[allcombs[1,.x]]]$proj$tbl[year >= start.year & cohort == "North Atlantic right whales" & prj == ndiff[i,1], N] -
            obj[[allcombs[2,.x]]]$proj$tbl[year >= start.year &  cohort == "North Atlantic right whales" & prj == ndiff[i,2], N]
        }
        
        data.table::data.table(year = seq(start.year, end.year),
                               mean = colMeans(diffproj),
                               lwr = apply(X = diffproj, MARGIN = 2, FUN = function(x) quantile(x, 0.025)),
                               uppr = apply(X = diffproj, MARGIN = 2, FUN = function(x) quantile(x, 0.975)),
                               label = paste0("[",proj.labels[allcombs[1,.x]], "] vs. [", proj.labels[allcombs[2,.x]], "]"))
        
      }) |> data.table::rbindlist()
      
      p <- ggplot2::ggplot() +
        {if (interval) ggplot2::geom_ribbon(data = difftbl, 
                                            aes(x = year, y = mean, ymin = lwr, ymax = uppr), alpha = 0.15)} +
        ggplot2::geom_line(data = difftbl, aes(x = year, y = mean), linetype = "solid", linewidth = 0.8) +
        ggplot2::scale_y_continuous(limits = c(min(0, min(difftbl$lwr)), max(0, max(difftbl$uppr)))) +
        ggplot2::geom_abline(intercept = 0, slope = 0, linetype = "dotted") +
        xlab("") + ylab ("Difference in predicted NARW abundance") +
        ggplot2::facet_wrap(~label, scales = scales, ncol = ncol) +
        ggplot2::scale_x_continuous(breaks = pretty(c(start.year, end.year), n = 10)) +
        ggplot2::scale_y_continuous(breaks = pretty(c(min(0, min(difftbl$lwr)), max(0, max(difftbl$uppr))), n = 10)) +
        theme_narw()

      suppressWarnings(print(p))

      
    } else {

    df <- purrr::map2(
      .x = obj, 
      .y = proj.labels,
      .f = ~ dplyr::mutate(.data = .x$proj$mean, label = .y)
    ) |> do.call(what = rbind)
    
    if(!length(obj) == length(unique(df$label))) stop("Missing or incorrect labels")
    
    # Range of y-axis
    yrange <- range(pretty(range(0, df$lwr, df$uppr, noaa_pva$low95, noaa_pva$high95), ny))
    
    # Schedule(s) of activities
    schedules <- purrr::map2(
      .x = 1:length(obj),
      .y = proj.labels,
      .f = ~ {
        obj[[.x]]$param$schedule |>
          tibble::enframe() |>
          dplyr::rename(year = name, phase = value) |>
          dplyr::mutate(
            label = .y,
            phase = c("B", "C", "O")[phase + 1],
            ypos = timeline.y + (.x-1)*(max(yrange)/10),
            year = as.numeric(gsub(pattern = "yr ", replacement = "", x = year)) + start.proj
          ) |>
          data.table::as.data.table()
      }
    ) |> data.table::rbindlist() |> 
      dplyr::filter(year >= start.year)
    
    schedules <- schedules[, list(start = min(year), end = max(year), ypos = unique(ypos)), list(label, phase)]
    schedules[, start:=ifelse(start == start.proj, start, start + 0.125), list(label, phase)]
    schedules[, end:=ifelse(end == end.year, end, end + 0.875), list(label, phase)]
    schedules[, xpos := mean(c(start, end)), list(label, phase)]
    
    # Subset the data by cohort if needed
    if(cohort){
      df <- df[cohort != "North Atlantic right whales"]
    } else {
      df <- df[cohort == "North Atlantic right whales"]
    }
    
    p <- ggplot2::ggplot()
    
    if(noaa) p <- p + 
      # 95% credible interval for NOAA trend
      {if (interval) ggplot2::geom_ribbon(data = noaa_pva, aes(x = year, y = median, ymin = low95, ymax = high95, fill = label, group = label), alpha = 0.15)} +
      # 50% credible interval for NOAA trend
      {if (interval) ggplot2::geom_ribbon(data = noaa_pva, aes(x = year, y = median, ymin = low50, ymax = high50, fill = label, group = label), alpha = 0.15)} +
      ggplot2::geom_line(data = noaa_pva, aes(x = year, y = median, group = label, colour = label), linetype = "solid", 
                         linewidth = 0.8, alpha = 0.75) +
      ggplot2::scale_colour_manual(values = "black", name = "") +
      ggplot2::scale_fill_manual(values = "black", name = "")
    
    # Generate plot(s)
    p <- p +
      {if(noaa) ggnewscale::new_scale_fill() } +
      {if(noaa) ggnewscale::new_scale_colour()} + 
      {if (interval & shade) ggplot2::geom_ribbon(data = df, aes(x = year, y = mean, ymin = lwr, ymax = uppr, fill = label, group = label),
                                                  alpha = 0.3)} +
      {if (interval & !shade) ggplot2::geom_line(data = df, aes(x = year, y = lwr, colour = label, group = label), size = 0.5, linetype = "dotdash")} +
      {if (interval & !shade) ggplot2::geom_line(data = df, aes(x = year, y = uppr, colour = label, group = label), size = 0.5, linetype = "dotdash")} +
      ggplot2::geom_path(data = df, aes(x = year, y = mean, colour = label, group = label), linewidth = 0.8, na.rm = TRUE) +
      ggplot2::facet_wrap(~cohort, scales = scales, ncol = ncol) +
      ggplot2::scale_fill_manual(values = colour, name = "") +
      ggplot2::scale_colour_manual(values = colour, name = "") +
      ggplot2::scale_y_continuous(breaks = ~ pretty(.x, n = ny), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
      ggplot2::scale_x_continuous(limits = c(start.year, end.year), breaks = pretty(c(start.year, end.year), n = nx)) +
      xlab("") +
      ylab("Abundance") +
      theme_narw()
    
    # Height of segment ends
    h <- 5
    
    # Add timeline
    if(timeline){
      p <- p +
        # timeline <- ggplot2::ggplot()
        # Main timeline bar
        ggplot2::geom_segment(
          data = schedules,
          aes(x = start, xend = end, y = ypos, yend = ypos, colour = label),
          linetype = 1, linewidth = 0.75
        ) +
        ggplot2::geom_segment(
          data = schedules,
          aes(x = start, xend = start, y = ypos-h, yend = ypos+h, colour = label),
          linetype = 1, linewidth = 0.75
        ) +
        ggplot2::geom_segment(
          data = schedules,
          aes(x = end, xend = end, y = ypos-h, yend = ypos+h, colour = label),
          linetype = 1, linewidth = 0.75
        ) +
        # Phase labels
        ggplot2::geom_label(
          data = schedules,
          aes(x = xpos, y = ypos, label = phase, colour = label),
          fill = "#ECECEC", label.size = NA,
          nudge_y = nudge.y, size = 5, show.legend = FALSE
        ) 
      
    } # End if timeline
    
    # +
    #   ggplot2::scale_fill_manual(values = colour, name = "") +
    #   ggplot2::scale_colour_manual(values = colour, name = "") +
    #   ggplot2::scale_y_continuous(limits = c(-2*5, h*10)) +
    #   ggplot2::scale_x_continuous(limits = c(start.year, end.year), breaks = pretty(c(start.year, end.year), n = nx)) +
    #   xlab("") +
    #   ylab("") +
    #   theme_narw() +
    #   ggplot2::theme(axis.text.x = ggplot2::element_blank(),
    #                  axis.ticks.x = ggplot2::element_blank(),
    #                  axis.text.y = ggplot2::element_blank(),
    #                  axis.ticks.y = ggplot2::element_blank(),
    #                  panel.background = ggplot2::element_blank(),
    #                  panel.grid.major.y = ggplot2::element_blank(),
    #                  panel.grid.minor.y = ggplot2::element_blank(),
    #                  legend.position = "none")
    
    
    # (p / timeline) + patchwork::plot_layout(heights = c(5,1))
    
    # Adapt layout for inclusion in package vignette
    if(vignette){
      p <- p + ggplot2::theme(panel.spacing.y = unit(0.0, "pt"))
      p_tab <- suppressWarnings(ggplot2::ggplotGrob(p))
      if(ncol == 3) a <- 2 else if(ncol == 2) a <- 3
      p.names <- do.call(c, purrr::map(.x = 1:a, .f = ~paste0("axis-b-", 1:6, "-", .x)))
      p_filtered <- gtable_filter_remove(p_tab, name = p.names, trim = FALSE)
      grid::grid.newpage()
      grid::grid.draw(p_filtered)
    } else {
      suppressWarnings(print(p))
    }
    
  } # End diff
  
}
