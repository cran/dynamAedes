## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10, fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#") 

## ----  message=FALSE, warning=FALSE-------------------------------------------
#Load packages
require(spatstat)
require(sp)
require(gstat)
require(parallel)
require(eesim)
require(tidyverse)
require(geosphere)
require(ggplot2)
require(rgeos)
#require(rgdal)
require(dynamAedes)

Sys.setlocale("LC_TIME", "en_GB.UTF-8")  

## -----------------------------------------------------------------------------
gridDim <- 20 # 5000m/250 m = 20 columns and rows
xy <- expand.grid(x=1:gridDim, y=1:gridDim)

## ---- message=FALSE-----------------------------------------------------------
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
r <- raster(nrow=gridDim, ncol=gridDim, crs=utm32N, ext=extent(1220000,1225000, 5700000,5705000))

values(r)=xyz$sim1
plot(r)

df <- data.frame("id"=1:nrow(xyz), raster::coordinates(r))
bbox <- as(extent(r), "SpatialPolygons")

# Store Parameters for autocorrelation
autocorr_factor <- values(r)

## -----------------------------------------------------------------------------
ndays = 365*1 #length of the time series in days
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

## -----------------------------------------------------------------------------
hist(sim_temp[[1]]$x, 
     xlab="Temperature (°C)", 
     main="Histogram of simulated temperatures")

plot(sim_temp[[1]]$date,
     sim_temp[[1]]$x,
     main="Simulated temperatures seasonal trend", 
     xlab="Date", ylab="Temperature (°C)"
     )

## -----------------------------------------------------------------------------
mat <- lapply(1:ncell(r), function(x) {
	d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
	return(d_t)
})

mat <- do.call(rbind,mat)

## ----  message=FALSE, warning=FALSE, hide=TRUE--------------------------------
oldpar <- par(mfrow = c(1,2)) 

## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
hist(mat, xlab="Temperature (°C)", main="Histogram of simulated spatial autocorreled temperatures")
hist(sim_temp[[1]]$x, xlab="Temperature (°C)", main="Histogram of simulated temperatures", col="red")
par(mfrow=c(1,1))

## ----  message=FALSE, warning=FALSE, hide=TRUE--------------------------------
par(oldpar) 

## -----------------------------------------------------------------------------
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)

## -----------------------------------------------------------------------------
w <- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))

## -----------------------------------------------------------------------------
cc <- df_temp[,c("x","y")]

## ---- evaluate=FALSE----------------------------------------------------------
## Define the day of introduction (May 1st is day 1)
str = "2000-06-01"
## Define the end-day of life cycle (July 2nd is the last day)
endr = "2000-07-02"
## Define the number of eggs to be introduced
ie = 100
## Define the number of model iterations
it = 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters=1
## Define proj4 string for input coordinates
utm32N = "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
## Define the number of parallel processes (for sequential iterations set nc=1)
cl = 1

## ----results='hide', message=FALSE, warning=FALSE, evaluate=FALSE-------------
simout=dynamAedes.m(species="albopictus", 
            scale="rg",  
            jhwv=habitat_liters,  
            temps.matrix=w[,as.numeric(format(as.Date(str),"%j")):as.numeric(format(as.Date(endr),"%j"))],
            coords.proj4=utm32N,
            cells.coords=cc,
            startd=str,
            endd=endr,
            n.clusters=cl,
            iter=it,
            intro.eggs=ie,
            compressed.output=TRUE,
            seeding=TRUE,
            verbose=FALSE)

## ---- evaluate=FALSE----------------------------------------------------------
print(it)
print(length(simout))

## ---- evaluate=FALSE----------------------------------------------------------
length(simout[[1]])

## ---- evaluate=FALSE----------------------------------------------------------
dim(simout[[1]][[1]])

## ---- evaluate=FALSE----------------------------------------------------------
psi(input_sim = simout, eval_date = 30)

## ---- evaluate=FALSE----------------------------------------------------------
plot(psi_sp(coords = cc, input_sim = simout, eval_date = 30, n.clusters=cl))

## ----message=FALSE, warning=FALSE, evaluate=FALSE-----------------------------
dd <- max(sapply(simout, function(x) length(x)))#retrieve the maximum number of simulated days
egg <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=1))
juv <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=2))
ad <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=3))
eggd <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=4))

egg$myStage='Egg'
juv$myStage='Juvenile'
ad$myStage='Adult'
eggd$myStage='Diapausing egg'

outdf=bind_rows(egg, juv, ad, eggd) %>% 
  as_tibble()

outdf$Date=rep(seq.Date(as.Date(str),as.Date(str)+dd-1, by="day"),4)

outdf %>% 
  mutate(myStage=factor(myStage, levels= c('Egg', 'Diapausing egg', 'Juvenile', 'Adult'))) %>% 
  ggplot( aes(y=(`50%`),x=Date, group=factor(myStage),col=factor(myStage))) +
  ggtitle("Ae. albopictus Interquantile range abundance")+
  geom_line(linewidth=1.2)+
  geom_ribbon(aes(ymin=`25%`,ymax=(`75%`),fill=factor(myStage)),
              col="white",
              alpha=0.2,
              outline.type="full")+
  labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage")+
  facet_wrap(~myStage, scales = "free")+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14) , strip.text = element_text(face = "italic"))

