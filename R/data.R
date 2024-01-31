#' North Atlantic right whale density surfaces
#'
#' Monthly rasters giving the estimated mean abundance of right whales (Eubalaena glacialis) per grid cell for the given month, averaged over the period 2010-2020. This era reflects an apparent major shift in right whale distributions around 2010. 
#' 
#' @note Density predictions for Cape Cod Bay for January-May use estimates from Ganley et al. (2019). Predictions for the month of December use all surveys conducted by the Center for Coastal Studies during the month of December from 2003-2020. Additional survey data were added to v12, and the aggregate database of surveys now extends through spring of 2019. The study area was extended farther inshore in certain bays and estuaries, per NOAA’s request.
#'
#' @references 
#'
#' Roberts JJ, Yack TM, Halpin PN (2022) Density Model for North Atlantic right whale (Eubalaena glacialis) for the U.S. East Coast, Version 12, 2022-02-14. Prepared for Naval Facilities Engineering Command, Atlantic by the Duke University Marine Geospatial Ecology Lab, Durham, NC.
#' 
#' Roberts JJ, Best BD, Mannocci L, Fujioka E, Halpin PN, Palka DL, Garrison LP, Mullin KD, Cole TVN, Khan CB, McLellan WM, Pabst DA, Lockhart GG (2016) Habitat-based cetacean density models for the U.S. Atlantic and Gulf of Mexico. Scientific Reports, 6: 22615. DOI: 10.1038/srep22615
#' 
#' Ganley L, Brault S, Mayo CA (2019) What we see is not what there is: Estimating North Atlantic right whale Eubalaena glacialis local abundance. Endangered Species Research, 38: 101–113. DOI: 10.3354/esr00938
#'
#' @format A list of twelve \code{SpatialGridDataFrame} objects.
#' \describe{
#'   \item{Jan}{Density surface for January}
#'   \item{Feb}{Density surface for February}
#'   \item{Mar}{Density surface for March}
#'   \item{Apr}{Density surface for April}
#'   \item{May}{Density surface for May}
#'   \item{Jun}{Density surface for June}
#'   \item{Jul}{Density surface for July}
#'   \item{Aug}{Density surface for August}
#'   \item{Sep}{Density surface for September}
#'   \item{Oct}{Density surface for October}
#'   \item{Nov}{Density surface for November}
#'   \item{Dec}{Density surface for December}
#' }
#'
#' @name density_narw
#' @docType data
#' @source Data provided by Jason Roberts (Marine Geospatial Ecology Lab, Duke University).
#' @keywords datasets
NULL