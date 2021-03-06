% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamAedes.R
\name{dynamAedes}
\alias{dynamAedes}
\title{Life cycle simulation of \emph{Aedes} mosquitoes}
\usage{
dynamAedes(
  species = "aegypti",
  intro.eggs = 0,
  intro.adults = 0,
  intro.juveniles = 0,
  scale = "ws",
  intro.cells = NULL,
  ihwv = 1,
  temps.matrix = NULL,
  startd = 1,
  endd = 10,
  cells.coords = NULL,
  lat = 0,
  long = 0,
  road.dist.matrix = NULL,
  avgpdisp = NA,
  intro.year = 2020,
  iter = 1,
  n.clusters = 1,
  cluster.type = "PSOCK",
  sparse.output = FALSE,
  compressed.output = TRUE,
  suffix = NA,
  cellsize = 250,
  maxadisp = 600,
  dispbins = 10,
  verbose = FALSE,
  seeding = FALSE
)
}
\arguments{
\item{species}{character. Select what species to model: \code{"aegypti"}, \code{"albopictus"}, \code{"japonicus"}, \code{"koreicus"}. Default \code{species = "aegypti"}.}

\item{intro.eggs}{positive integer. number of introduced eggs, default \code{intro.eggs = 100}.}

\item{intro.adults}{positive integer. number of introduced adults, default zero.}

\item{intro.juveniles}{positive integer. number of introduced juveniles, default zero.}

\item{scale}{character. Define the model spatial scale: punctual/weather station "ws", local "lc", or regional "rg". Active and passive dispersal is enabled only for \code{scale = "lc"}. Default \code{scale = "ws"}.}

\item{intro.cells}{positive integer. One or more cells (id) where to introduce the population at local ("lc") scale. If intro.cells=NULL, then a random cell is used for introduction; If intro.cells is a vector of cell ids then a cell is drawn at random from the vector (with repetition) for introduction in each model iteration.}

\item{ihwv}{positive integer. Larval-habitat water volume, define the volume (L) of water habitat presents in each spatial unit (parametrised with data retrieved from \doi{10.1111/1365-2664.12620}). Default \code{lhwv = 1}.}

\item{temps.matrix}{matrix. A matrix of daily (average) temperatures (in degrees \strong{Celsius degree x 1000}) used to fit the life cycle rates. This matrix must be organised with the daily temperature observations as columns and the geographic position of the i-grid cell as rows.}

\item{startd}{positive integer. Day of start of the simulations, referring to the column number of temps.matrix.}

\item{endd}{positive integer. Day of end of the simulation, referring to the column number of temps.matrix.}

\item{cells.coords}{matrix. A matrix reporting the spatial coordinates of the temperature observations.}

\item{lat}{numeric. Latitude value of the area of interested used to derive the photoperiod (and thus the diapause eggs allocation function).}

\item{long}{numeric. Longitude value of the area of interested used to derive the photoperiod (and thus the diapause eggs allocation function)}

\item{road.dist.matrix}{matrix. when \code{scale = "lc"}, defines the matrix containing the distances (in meters) between grid cells intersecting the road network for the mosquito passive dispersal process.}

\item{avgpdisp}{optional. when \code{scale = "lc"}, define the average car trip distance for the mosquito passive dispersal process. The value can be set by the users (positive numeric), or the estimates made by \href{https://publications.jrc.ec.europa.eu/repository/handle/JRC77079}{Pasaoglu et al. 2012}) for the following European countries: France "fra", Germany "deu", Italy "ita", Poland "pol", Spain "esp" and the United Kingdom "uk". The average passive dispersal distance must be smaller than the maximum distance in \strong{road.dist.matrix}.}

\item{intro.year}{numeric. Year of the beginning of iterations (for photoperiod calculation).}

\item{iter}{positive integer. Define the number of model iterations.}

\item{n.clusters}{positive integer. Defines the number of parallel processes.}

\item{cluster.type}{character. Defines the type of cluster, default "PSOCK".}

\item{sparse.output}{logical. The output matrix is optimised for sparse-matrix algebra (e.g. zeros are indexed).}

\item{compressed.output}{logical. Default TRUE, if FALSE provide abundance for each model's subcompartiment; if FALSE abundances are summed per compartment.}

\item{suffix}{character. Model output suffix for output RDS.}

\item{cellsize}{(positive integer. When \code{scale = "lc"}, defines the minimal distance for the active dispersal kernel and should match the spatial resolution of temps.matrix to avoid inconsistencies. Default cellsize = 250}

\item{maxadisp}{positive integer. When \code{scale = "lc"}, defines the maximum daily dispersal, default maxadisp = 600.}

\item{dispbins}{positive integer. When scale = "lc", defines the resolution of the dispersal kernel, default dispbins = 10.}

\item{verbose}{logical. if TRUE then an overview of population dynamics is printed in the console.}

\item{seeding}{logical, default \code{FALSE}, if \code{seeding=TRUE} a fixed seed is applied for result reproducibility.}
}
\value{
Matrix or a list of matrices containing, for each iteration, the number of individuals in each life stage per day (and for each grid cell of the study area if scale="lc" or "rg"). If the argument compressed.output=FALSE (default TRUE), the model returns the daily number of individuals in each life stage sub-compartment.
}
\description{
Function to simulate population dynamics of \emph{Aedes} mosquitoes
}
\examples{
## Run dynmAedes at local scale for 5 days
# Make a toy temperature time series
\dontrun{
w <- matrix(seq(20,25,length.out=5),ncol=5)*1000
# Run the model
	dynamAedes(
	species="koreicus", 
	scale="ws",
	intro.eggs=10, 
	ihwv=2, 
	temps.matrix=w, 
	startd=2, 
	endd=5,
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE)
}
}
\seealso{
Beta regression functions were taken from the R package \code{aomisc}, which is available at \url{https://github.com/OnofriAndreaPG/aomisc}.
}
\author{
Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
}
