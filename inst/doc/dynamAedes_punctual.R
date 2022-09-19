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
## Define the day of introduction (July 1st is day 1)
str = "2000-07-01"
## Define the end-day of life cycle (August 1st is the last day)
endr = "2000-08-01"
## Define the number of eggs to be introduced
ie = 1000
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
#outname= paste0("dynamAedes.m_albo_ws_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie)

## -----------------------------------------------------------------------------
df_temp <- data.frame("Date" = sim_temp[[1]]$date, "temp" = sim_temp[[1]]$x)
w <- t(as.integer(df_temp$temp*1000)[format(as.Date(str)+1,"%j"):format(as.Date(endr)+1,"%j")])

## ----   message=FALSE, warning=FALSE------------------------------------------
simout <- dynamAedes.m(species="albopictus", 
                     scale="ws",  
                     jhwv=habitat_liters,  
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
simout[[1]][[15]]

## ----message=FALSE, warning=FALSE---------------------------------------------
psi(input_sim = simout, eval_date = 30)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Retrieve the maximum number of simulated days
dd <- max(sapply(simout, function(x) length(x)))

# Compute the abundance
outdf <- rbind(data.frame(cbind(
                    data.frame(adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=1)), 
                    myStage = 'Egg')),
               data.frame(cbind(
                    adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=2), 
                    myStage='Juvenile')),
               data.frame(cbind(
                    adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=3), 
                    myStage='Adult')),
               data.frame(cbind(
                    adci(simout, eval_date=1:dd, breaks=c(0.25,0.50,0.75), st=4), 
                    myStage='Diapausing egg'))
               )

# Plot
outdf %>% 
  mutate(myStage=factor(myStage, levels= c('Egg', 'Diapausing egg', 'Juvenile', 'Adult')), 
         X25.=as.numeric(X25.), 
         X50.=as.numeric(X50.),
         X75.=as.numeric(X75.), 
         Date=rep(seq.Date(as.Date(str), as.Date(endr)-2, by="day"),4)) %>% 
ggplot( aes(x=Date, y=X50., group=factor(myStage),col=factor(myStage))) +
  ggtitle("Ae. albopictus Interquantile range abundance")+
  geom_line(size=0.8)+
  geom_ribbon(aes(ymin=X25.,ymax=X75.,fill=factor(myStage)),
              col="white",
              alpha=0.2,
              outline.type="full")+
  labs(x="Date", y="Interquantile range abundance", col="Stage", fill="Stage")+
  facet_wrap(~myStage, scales = "free_y")+
  theme_light()+
  theme(legend.pos="bottom",  
        text = element_text(size=16) , 
        strip.text = element_text(face = "italic"))

