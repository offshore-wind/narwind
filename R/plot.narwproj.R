#' Plot
#'
#' Make plots
#'
#' @param obj Object returned by run_model
#' @import ggplot2
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

plot.narwproj <- function(obj, 
                          cohort = FALSE, 
                          scales = "free", 
                          ncol = 3,
                          nx = 5,
                          n.breaks = 5,
                          vignette = FALSE){
  
  if(!inherits(obj, "narwproj")) stop("Object must be of class <narwproj>")
  
  gtable_filter_remove <- function (x, name, trim = TRUE){
    matches <- !(x$layout$name %in% name)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    if (trim) 
      x <- gtable::gtable_trim(x)
    x
  }
  
  start.year <- lubridate::year(lubridate::now())
  end.year <- start.year + obj$param$yrs
  
  df <- obj$prj$mean
  
  if(cohort){
    df <- df[cohort != "North Atlantic right whales"]
  } else {
    df <- df[cohort == "North Atlantic right whales"]
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(data = df, aes(x = year, y = mean, group = "cohort"), colour = "#1565C0", na.rm = TRUE) +
    ggplot2::geom_ribbon(data = df, aes(x = year, y = mean, ymin = lwr, ymax = uppr), alpha = 0.25, fill = "#1565C0") +
    ggplot2::facet_wrap(~cohort, scales = scales, ncol = ncol) + 
    ggplot2::scale_y_continuous(breaks = ~ pretty(.x, n = n.breaks)) +
    ggplot2::scale_x_continuous(breaks = pretty(c(start.year, end.year), n = nx)) +
    xlab("") + ylab("Abundance") +
    theme_narw()
  
  if(vignette){
    p <- p + ggplot2::theme(panel.spacing.y = unit(0.0, "pt"))
    p_tab <- suppressWarnings(ggplot2::ggplotGrob(p))
    if(ncol == 3) a <- 2 else if(ncol == 2) a <- 3
    p.names <- do.call(c, purrr::map(.x = 1:a, .f = ~paste0("axis-b-", 1:6, "-", .x)))
    p_filtered <- gtable_filter_remove(p_tab, name = p.names, trim = FALSE)
    grid::grid.newpage()
    grid::grid.draw(p_filtered)
  } else {
    print(p)
  }
}
