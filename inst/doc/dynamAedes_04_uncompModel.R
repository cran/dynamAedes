## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10, fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#")
options(rmarkdown.html_vignette.check_title = FALSE)

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
library(dynamAedes)
data(AedeslifeHistoryList)
knitr::kable(AedeslifeHistoryList$speciesheet, align = "ccccc")

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
species_data <- data.frame(
  Species = rep("*Ae. aegypti*", 6),
  `Sub-compartments` = paste("Sub-compartment", 1:6),
  Eggs = c("New layed egg", "2 day egg", "3 day egg", ">= 4 day egg", NA, NA),
  Juveniles = c("1 day juv", "2 day juv", "3 day juv", "4 day juv", "5 day juv", ">= 6 day juv"),
  Adults = c("blood fed", "ovipositing d1", "ovipositing d2", "Host-seeking", "new emerged", NA)
)

knitr::kable(species_data, format = "markdown", align = c('c', 'c', 'c', 'c', 'c'))

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
species_data <- data.frame(
  Species = rep("*Ae. albopictus*", 6),
  `Sub-compartments` = paste("Sub-compartment", 1:6),
  Eggs = c("New layed egg", "2 day egg", "3 day egg", ">= 4 day egg", NA, NA),
  Juveniles = c("1 day juv", "2 day juv", "3 day juv", "4 day juv", "5 day juv", ">= 6 day juv"),
  Adults = c("blood fed", "ovipositing d1", "ovipositing d2", "Host-seeking", "new emerged", NA),
  `Diapausing Eggs` = c("New layed degg", "2 day degg", "3 day degg", ">= 4 day degg", NA, NA)
)
knitr::kable(species_data, format = "markdown", align = c('c', 'c', 'c', 'c', 'c', 'c'))

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
species_data <- data.frame(
  Species = rep("*Ae. japonicus* or *Ae. koreicus*", 12),
  `Sub-compartments` = paste("Sub-compartment", 1:12),
  Eggs = c("New layed egg", "2 day egg", "3 day egg", "4 day egg", "5 day egg", 
           "6 day egg", "7 day egg", ">=8 day egg", NA, NA, NA, NA),
  Juveniles = c("1 day juv", "2 day juv", "3 day juv", "4 day juv", "5 day juv", 
                "6 day juv", "7 day juv", "8 day juv", "9 day juv", "10 day juv", 
                "11 day juv", ">=12 day juv"),
  Adults = c("blood fed", "ovipositing d1", "ovipositing d2", "Host-seeking", 
             "new emerged", NA, NA, NA, NA, NA, NA, NA),
  `Diapausing Eggs` = c("New layed degg", "2 day degg", "3 day degg", "4 day degg", 
                        "5 day degg", "6 day degg", "7 day degg", ">=8 day degg", 
                        NA, NA, NA, NA)
)
knitr::kable(species_data, format = "markdown", align = c('c', 'c', 'c', 'c', 'c', 'c'))

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
#Load packages
# simulate temperatures
library(eesim)
# plotting
library(ggplot2)
Sys.setlocale("LC_TIME", "en_GB.UTF-8")  

## ----warning=FALSE------------------------------------------------------------
ndays <- 365*1 #length of the time series in days
set.seed(123)
sim_temp <- create_sims(n_reps = 1, 
  n = ndays, 
  central = 16, 
  sd = 2, 
  exposure_type = "continuous", 
  exposure_trend = "cos1", exposure_amp = -1.0, 
  average_outcome = 12,
  outcome_trend = "cos1",
  outcome_amp = 0.8, 
  rr = 1.0055)

# Model settings
## Define the day of introduction (July 1st is day 1)
str <- "2000-07-01"
## Define the end-day of life cycle (August 1st is the last day)
endr <- "2000-08-01"
## Define the number of eggs to be introduced
ie <- 1000
## Define the number of model iterations
it <- 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters <- 1
## Define latitude and longitude for the diapause process
myLat <- 42
myLon <- 7
## Define the number of parallel processes (for sequential iterations set nc=1)
cl <- 1
## convert float temperatures to integer
df_temp <- data.frame("Date" = sim_temp[[1]]$date, "temp" = sim_temp[[1]]$x)
w <- t(as.integer(df_temp$temp*1000)[format(as.Date(str)+1,"%j"):format(as.Date(endr)+1,"%j")])

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
simout <- dynamAedes.m(species="albopictus", 
 scale="ws",  
 jhwv=habitat_liters,  
 temps.matrix=w,  
 startd=str, 
 endd=endr,  
 n.clusters=cl, 
 iter=it,  
 intro.eggs=ie,  
 compressed.output=FALSE, 
 lat=myLat, 
 long=myLon,
 verbose=FALSE,
 seeding=TRUE)

## -----------------------------------------------------------------------------
summary(simout)

## ----warning=FALSE, message=FALSE---------------------------------------------
simout@n_iterations

## ----warning=FALSE------------------------------------------------------------
length(simout@simulation[[1]])

## ----warning=FALSE------------------------------------------------------------
class(simout@simulation[[1]][[1]])
simout@simulation[[1]][[1]]
simout@simulation[[1]][[15]]

## ----message=FALSE, warning=FALSE---------------------------------------------
# Retrieve the maximum number of simulated days
dd <- max(simout)
# Compute the inter-quartile of abundances along the iterations
breaks <- c(0.25,0.50,0.75)
ed <- 1:dd
hs <- adci(simout, eval_date=ed, breaks=breaks, 
           stage="Adults",
           sub_stage = "Host-seeking" ) 
head(hs)
tail(hs)

## ----message=FALSE, warning=FALSE---------------------------------------------
ggplot(hs, aes(x=day, y=X50., group=factor(stage), col=factor(stage))) +
ggtitle("Host-seeking Ae. albopictus Interquantile range abundance") +
geom_ribbon(aes(ymin=X25., ymax=X75., fill=factor(stage)), 
  col="white", 
  alpha=0.2, 
  outline.type="full") +
geom_line(linewidth=0.8) +
ylim(0,10)+
labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage") +
theme_classic() +
theme(legend.pos="bottom",  
  text = element_text(size=16), 
  strip.text = element_text(face = "italic"))


## ----results='hide'-----------------------------------------------------------
library(gstat)
library(terra)
gridDim <- 20 # 5000m/250 m = 20 columns and rows
xy <- expand.grid(x=1:gridDim, y=1:gridDim)
varioMod <- vgm(psill=0.5, range=100, model='Exp') # psill = partial sill = (sill-nugget)
# Set up an additional variable from simple kriging
zDummy <- gstat(formula=z~1, 
                locations = ~x+y, 
                dummy=TRUE,
                beta=1, 
                model=varioMod, 
                nmax=1)
# Generate a randomly autocorrelated predictor data field
set.seed(123)
xyz <- predict(zDummy, newdata=xy, nsim=1)
utm32N <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
r <- terra::rast(nrow=gridDim, ncol=gridDim, crs=utm32N, ext=terra::ext(1220000,1225000, 5700000,5705000))
terra::values(r) <- xyz$sim1
# plot(r, main="SAC landscape")

# convert to a data.frame
df <- data.frame("id"=1:nrow(xyz), terra::crds(r))
bbox <- terra::as.polygons(terra::ext(r), crs=utm32N)

# Store Parameters for autocorrelation
autocorr_factor <- terra::values(r)

# "expand onto space" the temperature time series by multiplying it with the autocorrelated surface simulated above. 
mat <- do.call(rbind, lapply(1:ncell(r), function(x) {
	d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
	return(d_t)
}))

# format simulated temperature
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)
w <- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))
# define a two-column matrix of coordinates to identify each cell in the lattice grid.
cc <- df_temp[,c("x","y")]

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
simout <- dynamAedes.m(species="albopictus", 
            scale="rg",  
            jhwv=habitat_liters,  
            temps.matrix=w[,as.numeric(format(as.Date(str),"%j")):as.numeric(format(as.Date(endr),"%j"))],
            coords.proj4=utm32N,
            cells.coords=as.matrix(cc),
            startd=str,
            endd=endr,
            n.clusters=cl,
            iter=it,
            intro.eggs=ie,
            compressed.output=FALSE,
            seeding=TRUE,
            verbose=FALSE)

## ----warning=FALSE------------------------------------------------------------
summary(simout)

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
# Retrieve the maximum number of simulated days
dd <- max(simout)

# Compute a raster with the  median of the iterations
breaks <- c(0.50)
ed <- 1:dd
hs.r <- adci(simout, eval_date=ed, breaks=breaks, 
     stage="Adults",
     sub_stage = "Host-seeking", type="N")

## ----warning=FALSE------------------------------------------------------------
# inspect the raster
hs.r$`Host-seeking_q_0.5`

# plot the raster with the median host-seeking abundace
plot(hs.r$`Host-seeking_q_0.5`$day30)

