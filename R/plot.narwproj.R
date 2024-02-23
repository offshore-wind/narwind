#' Plot population trends
#'
#' Takes one or more \code{narwproj} object(s) and creates plots of the estimated population trajectories.
#'
#' @param ... One or more objects of class \code{narwproj}, as returned by \link{predict.narwsim}.
#' @param interval Logical. If \code{TRUE}, percentile confidence intervals are shown on the plots. Defaults to \code{TRUE}.
#' @param cohort Logical. If \code{TRUE}, separate plots are returned for each cohort. If \code{FALSE}, a single plot for the whole population is shown. Defaults to \code{FALSE}.
#' @param noaa Logical. If \code{TRUE}, the population trajectory predicted as part of NOAA's population viability analysis (Runge et al., 2023) is also plotted. Defaults to \code{FALSE}.
#' @param scales Character. Defines whether axis scales should be constant or vary across plots. Can be one of "fixed" or "free" (the default).
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
                          interval = TRUE,
                          cohort = FALSE,
                          noaa = FALSE,
                          scales = "free",
                          ncol = 3,
                          nx = 5,
                          ny = 5){
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
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
    colour <- c("#0098A0", "#E0B200", "#BF0B98", "#FF9B19FA", "#289DE0", "#E6124B")
  }
  
  if (length(obj) > length(colour)) colour <- c(colour, 1:(length(obj) - length(colour)))
  
  # Define years for the x-axis
  start.year <- lubridate::year(lubridate::now())
  end.year <- start.year + obj[[1]]$param$yrs
  
  # Extract the trend data
  df <- purrr::map(.x = obj, .f = ~{
    nm <- .x$param$label
    if(nm == "") nm <- "narwind"
    dplyr::mutate(.data = .x$proj$mean, label = nm)
  }) |> do.call(what = rbind)
  
  if(!length(obj) == length(unique(df$label))) stop("Missing or incorrect labels")

  # Subset the data by cohort if needed
  if(cohort){
    df <- df[cohort != "North Atlantic right whales"]
  } else {
    df <- df[cohort == "North Atlantic right whales"]
  }
  
  p <- ggplot2::ggplot()
  
  if(noaa) p <- p + 
    {if (interval) ggplot2::geom_ribbon(data = noaa_pva, aes(x = year, y = median, ymin = low95, ymax = high95, fill = label, group = label), alpha = 0.15)} +
    {if (interval) ggplot2::geom_ribbon(data = noaa_pva, aes(x = year, y = median, ymin = low50, ymax = high50, fill = label, group = label), alpha = 0.15)} +
    # ggplot2::geom_line(data = noaa_pva, aes(x = year, y = low95, group = label, colour = label), linetype = "dotted") +
    # ggplot2::geom_line(data = noaa_pva, aes(x = year, y = high95, group = label, colour = label), linetype = "dotted") +
    # ggplot2::geom_line(data = noaa_pva, aes(x = year, y = low50, group = label, colour = label), linetype = "dotdash") +
    # ggplot2::geom_line(data = noaa_pva, aes(x = year, y = high50, group = label, colour = label), linetype = "dotdash") +
    ggplot2::geom_line(data = noaa_pva, aes(x = year, y = median, group = label, colour = label), linetype = "solid", 
                       linewidth = 0.8, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = "black", name = "") +
    ggplot2::scale_fill_manual(values = "black", name = "")
  
  # Generate plot(s)
  p <- p +
    {if(noaa) ggnewscale::new_scale_fill() } +
    {if(noaa) ggnewscale::new_scale_colour()} + 
    {if (interval) ggplot2::geom_ribbon(data = df, aes(x = year, y = mean, ymin = lwr, ymax = uppr, fill = label, group = label),
                                        alpha = 0.3)} +
    ggplot2::geom_path(data = df, aes(x = year, y = mean, colour = label, group = label), linewidth = 0.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~cohort, scales = scales, ncol = ncol) +
    ggplot2::scale_fill_manual(values = colour, name = "") +
    ggplot2::scale_colour_manual(values = colour, name = "") +
    ggplot2::scale_y_continuous(breaks = ~ pretty(.x, n = ny), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::scale_x_continuous(limits = c(start.year, end.year), breaks = pretty(c(start.year, end.year), n = nx)) +
    xlab("") +
    ylab("Abundance") +
    theme_narw()
  
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
  
}
