## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10,fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#") 
options(mc.cores=2)

## ---- echo=FALSE, results='hide',message=FALSE--------------------------------
Sys.setlocale("LC_TIME", "en_GB.UTF-8")  

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
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
require(dynamAedes)

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
df_temp <- data.frame("Date" = sim_temp[[1]]$date, "temp" = sim_temp[[1]]$x)
w <- t(as.integer(df_temp$temp*1000))

## -----------------------------------------------------------------------------
## Define the day of introduction (August 1st is day 1)
str = 213
## Define the end-day of life cycle (September 1st is the last day)
endr = 213+31
## Define the number of eggs to be introduced
ie = 100
## Define the number of model iterations
it = 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters=1
#Define latitude and longitude for the diapause process
myLat=42
myLon=7
## Define the number of parallel processes (for sequential itarations set nc=1)
cl = 1
## Set output name for the *.RDS output will be saved
#outname= paste0("dynamAedes_albo_ws_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie)

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
simout <- dynamAedes(species="albopictus", 
                     scale="ws",  
                     ihwv=habitat_liters,  
                     temps.matrix=w,  
                     startd=str, 
                     endd=endr,  
                     n.clusters=cl, 
                     iter=it,  
                     intro.eggs=ie,  
                     compressed.output=TRUE, 
                     lat=myLat, 
                     long=myLon,
                     verbose=FALSE,
                     seeding=TRUE)

## -----------------------------------------------------------------------------
print(it)
print(length(simout))

## -----------------------------------------------------------------------------
length(simout[[1]])

## -----------------------------------------------------------------------------
dim(simout[[1]][[1]])
simout[[1]][[1]]
simout[[1]][[31]]

## ----message=FALSE, warning=FALSE---------------------------------------------
psi(input_sim = simout, eval_date = 31)

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
dd <- max(sapply(simout, function(x) length(x)))#retrieve the maximum number of simulated days
egg <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=1))
juv <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=2))
ad <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=3))
eggd <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=4))

egg$myStage="Egg"
egg$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
juv$myStage="Juvenile"
juv$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
ad$myStage="Adult"
ad$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
eggd$myStage="Diapausing egg"
eggd$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")

outdf=bind_rows(egg, juv, ad, eggd) %>% 
  as_tibble()

outdf %>% 
  mutate(myStage=factor(myStage, levels= c("Egg", "Diapausing egg", "Juvenile", "Adult"))) %>% 
  ggplot( aes(y=(`50%`),x=Date, group=factor(myStage),col=factor(myStage))) +
  ggtitle("Ae. albopictus Interquantile range abundance")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=`25%`,ymax=(`75%`),fill=factor(myStage)),
              col="white",
              alpha=0.2,
              outline.type="full")+
  labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage")+
  facet_wrap(~myStage, scales = "free")+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14) , strip.text = element_text(face = "italic"))

## -----------------------------------------------------------------------------
gridDim <- 20 # 5000m/250 m = 20 columns and rows
xy <- expand.grid(x=1:gridDim, y=1:gridDim)

## ---- message=FALSE-----------------------------------------------------------
varioMod <- vgm(psill=0.005, range=100, model='Exp') # psill = partial sill = (sill-nugget)
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
r <- raster(nrow=gridDim, ncol=gridDim, crs=utm32N, ext=extent(0,5000, 0,5000))
values(r)=xyz$sim1
plot(r)

df <- data.frame("id"=1:nrow(xyz), coordinates(r))
bbox <- as(extent(r), "SpatialPolygons")

# Store Parameters for autocorrelation
autocorr_factor <- values(r)

## -----------------------------------------------------------------------------
mat <- lapply(1:ncell(r), function(x) {
	d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
	return(d_t)
})

mat <- do.call(rbind,mat)

## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
hist(mat, xlab="Temperature (°C)", main="Histogram of simulated spatial autocorreled temperature")
hist(sim_temp[[1]]$x, xlab="Temperature (°C)", main="Histogram of simulated temperatures", col="red")
par(mfrow=c(1,1))

# Format temperature data
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)

## -----------------------------------------------------------------------------
set.seed(123)
pts <- spsample(bbox, 5, type="random")
roads <- spLines(pts)

# Check simulated segment
raster::plot(r)
raster::plot(roads, add=T)

## -----------------------------------------------------------------------------
buff <- buffer(roads, width=250)
crs(buff) <- crs(r)
# Check grid, road segment and buffer
raster::plot(r)
raster::plot(buff, add=T)
raster::plot(roads, add=T, col="red")

## ----  message=FALSE----------------------------------------------------------
df_sp <- df
coordinates(df_sp)=~x+y
df_sp <- raster::intersect(df_sp,buff)

# Check selected cells
raster::plot(r)
raster::plot(buff, add=T)
raster::plot(df_sp, add=T, col="red")

## -----------------------------------------------------------------------------
dist_matrix <- as.matrix(dist(coordinates(df_sp)))

## -----------------------------------------------------------------------------
w <- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))

## -----------------------------------------------------------------------------
cc <- df_temp[,c("x","y")]

## -----------------------------------------------------------------------------
colnames(dist_matrix) <- row.names(dist_matrix)

## -----------------------------------------------------------------------------
dist_matrix <- apply(dist_matrix,2,function(x) round(x/1000,1)*1000) 

# An histogram showing the distribution of distances of cells along the road network
hist(dist_matrix, xlab="Distance (meters)")

## -----------------------------------------------------------------------------
set.seed(123)
icellcoords <- df[sample(row.names(dist_matrix),1),c(2:3)]
set.seed(123)
icellid <- df[sample(row.names(dist_matrix),1),1]

raster::plot(r)
raster::plot(buff, add=T)
raster::plot(df_sp, add=T, col="red")
raster::plot(SpatialPoints(icellcoords), add=T, col="blue", pch=21)
raster::plot(SpatialPoints(coords=matrix(coordinates(r)[icellid,],ncol=2)), add=T, col="black", pch=21)

## -----------------------------------------------------------------------------
## Define cells along roads into which introduce propagules on day 1
intro.vector <- icellid
## Define the day of introduction (August 1st is day 1)
str = 213
## Define the end-day of life cycle (October 1st is the last day)
endr = 213+(61)
## Define the number of adult females to be introduced
ia = 5000
## Define the number of model iterations
it = 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters=1
#Define latitude and longitude for the diapause process
myLat=42
myLon=7
##Define average trip distance (in km)
mypDist=2
## Define the number of parallel processes (for sequential iterations set nc=1)
cl = 1

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
simout=dynamAedes(species="albopictus",
            scale="lc",  
            ihwv=habitat_liters,  
            temps.matrix=w,
            cells.coords=cc,
            road.dist.matrix=dist_matrix,
            avgpdisp=mypDist,
            intro.cells=intro.vector,
            startd=str, 
            endd=endr,  
            n.clusters=cl, 
            iter=it,  
            intro.adults=ia,  
            compressed.output=TRUE, 
            lat=myLat, 
            long=myLon, 
            cellsize=250,
            maxadisp=600,
            dispbins=10,
            seeding=TRUE,
            verbose= FALSE
            )

## -----------------------------------------------------------------------------
print(it)
print(length(simout))

## -----------------------------------------------------------------------------
length(simout[[1]])

## -----------------------------------------------------------------------------
dim(simout[[1]][[1]])

## -----------------------------------------------------------------------------
psi(input_sim = simout, eval_date = 61)

## -----------------------------------------------------------------------------
plot(psi_sp(coords = cc, input_sim = simout, eval_date = 61, n.clusters=cl))
raster::plot(buff, add=T)
raster::plot(df_sp, add=T, col="red")
raster::plot(SpatialPoints(icellcoords), add=T, col="blue", pch=21)


## ----message=FALSE, warning=FALSE---------------------------------------------
dd <- max(sapply(simout, function(x) length(x)))#retrieve the maximum number of simulated days
egg <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=1))
juv <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=2))
ad <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=3))
eggd <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=4))

egg$myStage="Egg"
egg$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
juv$myStage="Juvenile"
juv$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
ad$myStage="Adult"
ad$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
eggd$myStage="Diapausing egg"
eggd$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")

outdf=bind_rows(egg, juv, ad, eggd) %>% 
  as_tibble()

outdf %>% 
  mutate(myStage=factor(myStage, levels= c("Egg", "Diapausing egg", "Juvenile", "Adult"))) %>% 
  ggplot( aes(y=`50%`,x=Date, group=factor(myStage),col=factor(myStage))) +
  ggtitle("Ae. albopictus Interquantile range abundance")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=`25%`,ymax=`75%`,fill=factor(myStage)),
              col="white",
              alpha=0.2,
              outline.type="full")+
  labs(x="Date", y="Interquantile range abundance (Log10)", col="Stage", fill="Stage")+
  facet_wrap(~myStage, scales = "free")+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14) , strip.text = element_text(face = "italic"))

## -----------------------------------------------------------------------------
r=adci_sp(simout, coords=cc, eval_date=61, breaks=c(0.5,0.975), stage=3)
plot(r)


## -----------------------------------------------------------------------------
x=icci(simout, eval_date=1:61, breaks=c(0.25,0.50,0.75))
head(x)
tail(x)

## -----------------------------------------------------------------------------
x=dici(simout, coords=cc, eval_date=seq(1,61,length.out=61), breaks=c(0.25,0.50,0.75), space=FALSE)
plot(`75%`~day,x,type="l",ylab="Population dispersal (in meters) from cell of introduction",xlab="days from introduction")
lines(`50%`~day,x,type="l", col="red")
lines(`25%`~day,x,type="l")

## ---- evaluate=FALSE----------------------------------------------------------
w <- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))

## ---- evaluate=FALSE----------------------------------------------------------
cc <- df_temp[,c("x","y")]

## ---- evaluate=FALSE----------------------------------------------------------
## Define the day of introduction (June 1st is day 1)
str = 152
## Define the end-day of life cycle (October 1st is the last day)
endr = 152+(61*2)
## Define the number of eggs to be introduced
ie = 100
## Define the number of model iterations
it = 1 # The higher the number of simulations the better
## Define the number of liters for the larval density-dependent mortality
habitat_liters=1
#Define latitude and longitude for the diapause process
myLat=42
myLon=7
## Define the number of parallel processes (for sequential iterations set nc=1)
cl = 1

## ----results='hide', message=FALSE, warning=FALSE, evaluate=FALSE-------------
simout=dynamAedes(species="albopictus", 
            scale="rg",  
            ihwv=habitat_liters,  
            temps.matrix=w,  
            cells.coords=cc,
            startd=str, 
            endd=endr,  
            n.clusters=cl, 
            iter=it,  
            intro.eggs=ie,  
            compressed.output=TRUE, 
            lat=myLat, 
            long=myLon, 
            verbose = FALSE)

## ---- evaluate=FALSE----------------------------------------------------------
print(it)
print(length(simout))

## ---- evaluate=FALSE----------------------------------------------------------
length(simout[[1]])

## ---- evaluate=FALSE----------------------------------------------------------
dim(simout[[1]][[1]])

## ---- evaluate=FALSE----------------------------------------------------------
psi(input_sim = simout, eval_date = 123)

## ---- evaluate=FALSE----------------------------------------------------------
plot(psi_sp(coords = cc, input_sim = simout, eval_date = 123,n.clusters=cl))

## ----message=FALSE, warning=FALSE, evaluate=FALSE-----------------------------
dd <- max(sapply(simout, function(x) length(x)))#retrieve the maximum number of simulated days
egg <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=1))
juv <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=2))
ad <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=3))
eggd <- as.data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=4))

egg$myStage="Egg"
egg$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
juv$myStage="Juvenile"
juv$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
ad$myStage="Adult"
ad$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")
eggd$myStage="Diapausing egg"
eggd$Date=seq.Date(sim_temp[[1]]$date[str], sim_temp[[1]]$date[endr], by="day")

outdf=bind_rows(egg, juv, ad, eggd) %>% 
  as_tibble()

outdf %>% 
  mutate(myStage=factor(myStage, levels= c("Egg", "Diapausing egg", "Juvenile", "Adult"))) %>% 
  ggplot( aes(y=(`50%`),x=Date, group=factor(myStage),col=factor(myStage))) +
  ggtitle("Ae. albopictus Interquantile range abundance")+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=`25%`,ymax=(`75%`),fill=factor(myStage)),
              col="white",
              alpha=0.2,
              outline.type="full")+
  labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage")+
  facet_wrap(~myStage, scales = "free")+
  theme_light()+
  theme(legend.pos="bottom",  text = element_text(size=14) , strip.text = element_text(face = "italic"))

