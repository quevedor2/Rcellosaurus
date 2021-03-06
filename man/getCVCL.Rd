% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{getCVCL}
\alias{getCVCL}
\title{getCVCL}
\usage{
getCVCL(cellid, melt.cells, prioritize.datasets = TRUE)
}
\arguments{
\item{cellid}{character: cell line name (e.g. MCF-7)}

\item{melt.cells}{from data(melt.cells)}
}
\value{
CVCL_ style character name
}
\description{
get the Cellosaurus CVCL_ id for a given
cell line name
}
\examples{
data(melt.cells)
getCVCL("Hela", melt.cells)
}
