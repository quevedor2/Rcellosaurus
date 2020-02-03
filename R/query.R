
#' fullpull
#' @description Gets the synonymous cell lines and cell line that they derived from,
#' as well as any indication of MSI or probelmatic cell line (e.g. contamination)
#'
#' @param cvcl Cellosaurus CVCL id (e.g. CVCL_1384)
#' @param melt.cells data.frame made from cellosaurus xml
#'
#' @return Dataframe of all synonymous (SS), similar origins (OI), MSI positive, or problematic
#' lines associated with the CVCL
#' @export
#'
#' @examples data(melt.cells)
#' cvcl <- getCVCL('Hela', melt.cells)
#' fullpull(cvcl, melt.cells)
#'
#' fullpull("HeLa", melt.cells)
fullpull <- function(cvcl, melt.cells){
  if(!substr(cvcl, 1, 5) == "CVCL_"){
    cvcl <- getCVCL(cvcl, melt.cells)
    if(class(cvcl) == 'data.frame' | length(cvcl) == 0) stop("Could not find the CVCL id for given input")
  }
  rbind(.getSynonymous(cvcl, melt.cells),
        .getDerivedFrom(cvcl, melt.cells))
}


#' getCVCL
#' @description get the Cellosaurus CVCL_ id for a given
#' cell line name
#'
#' @param cellid character: cell line name (e.g. MCF-7)
#' @param melt.cells from data(melt.cells)
#'
#' @return CVCL_ style character name
#' @export
#' 
#' @examples data(melt.cells)
#' getCVCL("Hela", melt.cells)
getCVCL <- function(cellid, melt.cells){
  cl.match <- melt.cells[grep(paste0("^", cellid, "$"), melt.cells$ID),]
  if(nrow(cl.match) > 1) {
    warning("Multiple cells found for the given ID (likely contamination), please manully select 1 CVCL for further analysis")
    return(cl.match)
  } else {
    return(cl.match$CVCL)
  }
}
