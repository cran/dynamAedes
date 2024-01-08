## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10, fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#")
options(rmarkdown.html_vignette.check_title = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
#Load packages
# simulate temperatures
library(eesim)

# modelling 
library(dynamAedes)

# data manipulation and plotting
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

## ----warning=FALSE------------------------------------------------------------
hist(sim_temp[[1]]$x, 
 xlab="Temperature (°C)", 
 main="Histogram of simulated temperatures")

plot(sim_temp[[1]]$date,
 sim_temp[[1]]$x,
 main="Simulated temperatures seasonal trend", 
 xlab="Date", ylab="Temperature (°C)"
 )

## ----warning=FALSE------------------------------------------------------------
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
#Define latitude and longitude for the diapause process
myLat <- 42
myLon <- 7
## Define the number of parallel processes (for sequential iterations set nc=1)
cl <- 1

## ----warning=FALSE------------------------------------------------------------
df_temp <- data.frame("Date" = sim_temp[[1]]$date, "temp" = sim_temp[[1]]$x)
w <- t(as.integer(df_temp$temp*1000)[format(as.Date(str)+1,"%j"):format(as.Date(endr)+1,"%j")])

## ----message=FALSE, warning=FALSE---------------------------------------------
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

## ----warning=FALSE------------------------------------------------------------
summary(simout)

## ----warning=FALSE, message=FALSE, results='hide'-----------------------------
simout@n_iterations

## ----warning=FALSE, message=FALSE, results='hide'-----------------------------
simout@simulation

## ----warning=FALSE------------------------------------------------------------
length(simout@simulation[[1]])

## ----warning=FALSE------------------------------------------------------------
simout@simulation[[1]][[1]]
simout@simulation[[1]][[15]]

## ----message=FALSE, warning=FALSE---------------------------------------------
psi(input_sim = simout, eval_date = 30)

## ----message=FALSE, warning=FALSE---------------------------------------------

# Retrieve the maximum number of simulated days
dd <- max(simout)

# Compute the inter-quartile of abundances along the iterations
breaks=c(0.25,0.50,0.75)
ed=1:dd

outdf <- rbind(
  adci(simout, eval_date=ed, breaks=breaks, stage="Eggs"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Juvenile"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Adults"),
  adci(simout, eval_date=ed, breaks=breaks, stage="Dia"))

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
theme(legend.pos="bottom",  
  text = element_text(size=16), 
  strip.text = element_text(face = "italic"))


