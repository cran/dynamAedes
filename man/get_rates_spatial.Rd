% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_rates_spatial.R
\name{get_rates_spatial}
\alias{get_rates_spatial}
\title{Proportion of successful introductions}
\usage{
get_rates_spatial(
  coords = NULL,
  temps.matrix = NULL,
  species = "aegypti",
  rate_fun = NULL,
  spatial = FALSE,
  rate = TRUE,
  n.clusters = 1
)
}
\arguments{
\item{coords}{matrix (optional). A matrix reporting the spatial coordinates of the temperature observations.}

\item{temps.matrix}{matrix, vector. A matrix or a vector of daily average temperature (in Celsius degree /1000) used to fit the life cycle rates. This matrix must be organised with the daily temperature observations as columns and the geographic position of the i-grid cell as rows.}

\item{species}{character. Select what species to model: \code{"aegypti"}, \code{"albopictus"}, \code{"japonicus"}, \code{"koreicus"}. Default \code{species = "aegypti"}.}

\item{rate_fun}{character. Define one of the temperature-dependent functions of the model, e.g. ".e.hatch.rate".}

\item{spatial}{logical. Get a raster as output. Default \code{FALSE}}

\item{rate}{logical. Returns the daily rate. Default \code{TRUE}}

\item{n.clusters}{postive integer. Define the number of parallel processes.}
}
\value{
Returns the estimated value of a given temperature-dependent function. If \code{spatial = TRUE}, it will return a raster.
}
\description{
Compute the proportion of "successful" introductions.
}
\author{
Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
}
