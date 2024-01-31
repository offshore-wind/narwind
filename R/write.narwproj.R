#' Export population projection data
#'
#' Writes outputs from the stochastic population model to an .xlsx file on disk.
#' 
#' @param obj An object of class \code{narwproj}, as returned by \code{\link{predict.narwsim}}.
#' @param prefix Character string. Output file name. Defaults to \code{"narwproj"}.
#' @param ... Additional arguments passed to \code{\link[openxlsx]{write.xlsx}}.
#' @return An Excel file containing four sheets:
#' \itemize{
#' \item \code{Abundance}: Estimated mean right whale abundance (calculated across replicate projections), summarized by year and cohort. Lower and upper bounds of the associated confidence interval are also reported.
#' \item \code{Projection}: Estimated right whale abundance, by year, cohort, and projection.
#' \item \code{Births}: Numbers of calving events, summarized by year and projection.
#' \item \code{Deaths}: Numbers of mortality events, summarized by year and projection.
#' }
#' @export
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' m <- narw(1000)
#' m <- augment(m)
#' prj <- predict(m)
#' write(m)
#' }

write.narwproj <- function(obj, 
                           prefix = "narwproj",
                           ...){
  
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Function check
  if(!inherits(obj, "narwproj")) stop("Object must be of class <narwproj>")
  
  console("Saving")
  
  # Retrieve column names
  matt.attribs <- colnames(obj$dat$ind[[1]])
  
  # Extract population data
  narw.out <- purrr::map(.x = 1:obj$param$n, .f = ~{
    reshape_array(obj$dat$ind[[.x]], value, yr, attr, whale) |>
      dplyr::mutate(attr = mat.attribs[attr]) |>
      tidyr::pivot_wider(names_from = attr, values_from = value) |>
      dplyr::mutate(prj = .x) |>
      dplyr::relocate(prj, .before = yr)
  }) |> do.call(what = rbind) |>
    data.table::data.table()
  
  # Name sheets in output .xlsx file
  sheet.list <- list(
    "Abundance" = obj$prj$mean,
    "Projection" = obj$prj$proj,
    "Births" = obj$dat$birth$tot,
    "Deaths" = obj$dat$death
  )
  
  # Write to disk
  openxlsx::write.xlsx(x = sheet.list,
                       file = paste0(tolower(prefix), ifelse(!is.null(obj$param$label), paste0("_", obj$param$label)), ".xlsx"), 
                       asTable = TRUE,
                       firstRow = TRUE, 
                       tableStyle = "TableStyleMedium1",
                       bandedRows = TRUE,
                       withFilter = FALSE,
                       ...)
  
  console("Saving", suffix = tickmark())
  
}
