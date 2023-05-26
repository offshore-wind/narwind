#' Summary
#'
#' Summary information
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

write.narwproj <- function(obj, ...){
  
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default values for optional arguments
  filename <- "narwproj"
  overwrite <- TRUE

  # Default values
  if(length(args) > 0) {
    if("filename" %in% names(args)) prefix <- args[["filename"]]
    if("overwrite" %in% names(args)) overwrite <- args[["overwrite"]]
  }

  if(!"narwproj" %in% class(obj)) stop("Input must be of class <narwproj>")

  cat("Saving ...\n")

  matt.attribs <- colnames(obj$dat$ind[[1]])
  
  narw.out <- purrr::map(.x = 1:obj$param$n, .f = ~{
    reshape_array(obj$dat$ind[[.x]], value, yr, attr, whale) |>
      dplyr::mutate(attr = mat.attribs[attr]) |>
      tidyr::pivot_wider(names_from = attr, values_from = value) |>
      dplyr::mutate(prj = .x) |>
      dplyr::relocate(prj, .before = yr)
  }) |> do.call(what = rbind) |>
    data.table::data.table()
  
  births.df <- purrr::map(.x = 1:n, .f = ~{
    m <- matrix(rowSums(obj$dat$ind[[.x]][2:(yrs+1),"birth",], na.rm = TRUE), ncol = 1)
    colnames(m) <- .x
    m
  }) |> do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "birth") |>
    dplyr::select(prj, year, birth) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |> 
    dplyr::mutate(prj = as.numeric(prj))

  deaths.df <- purrr::map(.x = 1:n, .f = ~{
    m <- matrix(apply(X = obj$dat$ind[[.x]][2:(yrs+1),"alive",],
                      MARGIN = 1,
                      FUN = function(x) {
      r <- x[!is.na(x)]
      r <- sum(r == 0)
      r
      }), ncol = 1)
    colnames(m) <- .x
    m
  }) |> do.call(what = cbind) |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "year") |>
    dplyr::mutate(year = as.numeric(year)) |>
    tidyr::pivot_longer(!year, names_to = "prj", values_to = "death") |>
    dplyr::select(prj, year, death) |>
    dplyr::arrange(prj, year) |>
    data.table::data.table() |> 
    dplyr::mutate(prj = as.numeric(prj))
    
    sheet.list <- list(
      "Abundance" = obj$prj$mean,
      "Projection" = obj$prj$proj,
      "Births" = births.df,
      "Deaths" = deaths.df
    )
    
    openxlsx::write.xlsx(x = sheet.list,
      file = paste0(tolower(filename), ".xlsx"), 
      asTable = TRUE,
      overwrite = overwrite,
      firstRow = TRUE, 
      tableStyle = "TableStyleMedium1",
      bandedRows = TRUE,
      withFilter = FALSE)
  
}
