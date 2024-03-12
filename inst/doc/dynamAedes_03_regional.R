## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10, fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#") 
options(rmarkdown.html_vignette.check_title = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Libraries
library(gstat)
library(terra)
library(eesim)
library(dynamAedes)
library(ggplot2)
Sys.setlocale("LC_TIME", "en_GB.UTF-8")

## -----------------------------------------------------------------------------
gridDim <- 20 # 5000m/250 m = 20 columns and rows
xy <- expand.grid(x=1:gridDim, y=1:gridDim)

## ----message=FALSE------------------------------------------------------------
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

## -----------------------------------------------------------------------------
utm32N <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
r <- terra::rast(nrow=gridDim, ncol=gridDim, crs=utm32N, ext=terra::ext(1220000,1225000, 5700000,5705000))
terra::values(r) <- xyz$sim1
plot(r, main="SAC landscape")

# convert to a data.frame
df <- data.frame("id"=1:nrow(xyz), terra::crds(r))
bbox <- terra::as.polygons(terra::ext(r), crs=utm32N)

# Store Parameters for autocorrelation
autocorr_factor <- terra::values(r)

## -----------------------------------------------------------------------------
ndays = 365
set.seed(123)
sim_temp <- eesim::create_sims(n_reps = 1, 
                        n = ndays, 
                        central = 16, 
                        sd = 2, 
                        exposure_type = "continuous", 
                        exposure_trend = "cos1", exposure_amp = -1.0, 
                        average_outcome = 12,
                        outcome_trend = "cos1",
                        outcome_amp = 0.8, 
                        rr = 1.0055)

## -----------------------------------------------------------------------------
hist(sim_temp[[1]]$x, 
     xlab="Temperature (°C)", 
     main="Histogram of Simulated Temperatures")

plot(sim_temp[[1]]$date,
     sim_temp[[1]]$x,
     main="Simulated Temperature Seasonal Trend", 
     xlab="Date", ylab="Temperature (°C)"
     )

## -----------------------------------------------------------------------------
mat <- do.call(rbind, lapply(1:ncell(r), function(x) {
	d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
	return(d_t)
}))

## ----message=FALSE, warning=FALSE, hide=TRUE----------------------------------
oldpar <- par(mfrow = c(1,2)) 

## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
hist(mat, xlab="Temperature (°C)", main="Histogram of simulated spatial autocorreled temperature")
hist(sim_temp[[1]]$x, xlab="Temperature (°C)", main="Histogram of simulated temperatures", col="red")
par(mfrow=c(1,1))

## ----message=FALSE, warning=FALSE, hide=TRUE----------------------------------
par(oldpar) 

## -----------------------------------------------------------------------------
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)

## -----------------------------------------------------------------------------
w <- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))

## -----------------------------------------------------------------------------
cc <- df_temp[,c("x","y")]

## ----evaluate=FALSE-----------------------------------------------------------
## Define the day of introduction (May 1st is day 1)
str <- "2000-06-01"
## Define the end-day of life cycle (July 2nd is the last day)
endr <- "2000-07-02"
## Define the number of eggs to be introduced
ie <- 100
## Define the number of model iterations
it <- 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters <- 1
## Define proj4 string for input coordinates
utm32N <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
## Define the number of parallel processes (for sequential iterations set nc=1)
cl <- 1

## ----results='hide', message=FALSE, warning=FALSE, evaluate=FALSE-------------
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
            compressed.output=TRUE,
            seeding=TRUE,
            verbose=FALSE)

## ----warning=FALSE------------------------------------------------------------
summary(simout)

## ----warning=FALSE, message=FALSE, results='hide'-----------------------------
simout@n_iterations

## ----warning=FALSE, message=FALSE, results='hide'-----------------------------
simout@simulation

## ----warning=FALSE------------------------------------------------------------
length(simout@simulation[[1]])

## ----evaluate=FALSE-----------------------------------------------------------
dim(simout@simulation[[1]][[1]])

## ----evaluate=FALSE-----------------------------------------------------------
psi(input_sim = simout, eval_date = 30)

## ----evaluate=FALSE-----------------------------------------------------------
plot(psi_sp(input_sim = simout, eval_date = 30))

## ----message=FALSE, warning=FALSE, evaluate=FALSE-----------------------------
dd <- max(simout) #retrieve the maximum number of simulated days

# Compute the inter-quartile of abundances along the iterations
breaks=c(0.25,0.50,0.75)
ed=1:dd

# type "O" derives a non-spatial time series
outdf <- rbind(
  adci(simout, eval_date=ed, breaks=breaks, stage="Eggs", type="O"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Juvenile", type="O"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Adults", type="O"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Dia", type="O")
  )

## ----message=FALSE, warning=FALSE---------------------------------------------
outdf$stage <- factor(outdf$stage, levels= c('Egg', 'DiapauseEgg', 'Juvenile', 'Adult'))
outdf$Date <- rep(seq.Date(as.Date(str), as.Date(endr) - 2, by="day"), 4)

ggplot(outdf, aes(x=Date, y=X50., group=factor(stage), col=factor(stage))) +
ggtitle("Ae. albopictus Interquantile range abundance") +
geom_ribbon(aes(ymin=X25., ymax=X75., fill=factor(stage)), 
  col="white", 
  alpha=0.2, 
  outline.type="full") +
geom_line(linewidth=0.8) +
labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage") +
facet_wrap(~stage, scales = "free_y") +
theme_light() +
theme(legend.position="bottom",  
  text = element_text(size=16), 
  strip.text = element_text(face = "italic"))

