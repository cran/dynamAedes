## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=10, fig.height=10, fig.asp = 0.618, out.width = "95%", fig.align = "center", fig.dpi = 150, collapse = FALSE, comment = "#") 
options(rmarkdown.html_vignette.check_title = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
# Libraries
library(dplyr)
library(lubridate)
library(stringr)
library(dynamAedes)
library(ggplot2)
Sys.setlocale("LC_TIME", "en_GB.UTF-8")

## ----message=FALSE, warning=FALSE, echo=FALSE---------------------------------
templatedf <- structure(list(year = c(2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014), 
               setting_date = structure(c(16181,16202, 16223, 16251, 16265, 16279, 16293, 16307, 16321, 16335, 16349, 16363, 16377), class = "Date"), 
               sampling_date = structure(c(16202, 16223, 16251, 16265, 16279, 16293, 16307, 16321, 16335, 16349, 16363, 16377, 16391), class = "Date"), 
               value = c(0, 0, 0, 53, 26, 273, 215, 203, 76, 0, 137, 19, 0), 
               lifeStage = c("Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs", "Eggs"), trap_type = c("Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps", "Ovitraps"), 
               species = c("Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", 
                           "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", "Aedes albopictus", 
                           "Aedes albopictus")), 
          row.names = c(NA, -13L), class = c("tbl_df", "tbl", "data.frame"))

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(templatedf, align = "lccrr")

## -----------------------------------------------------------------------------
templatedf <- templatedf %>% 
  mutate(delta_week=lubridate::week(sampling_date) - lubridate::week(setting_date))

ggplot(templatedf, aes(x = delta_week)) + 
  geom_bar(fill="#1A85FF") + 
  scale_x_continuous(breaks = 1:4)+
  labs(x = "Weeks", y="Frequency") +
  ggtitle("Sampling frequency")+
  theme_classic()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

## -----------------------------------------------------------------------------
mySeq <- tibble("date"=seq.Date(as.Date('2014-01-01'), as.Date('2014-12-31'), by="day")) %>% 
          mutate(temporalID = paste0(lubridate::year(date), "_" , lubridate::week(date)), 
                 wday = wday(date) ) %>% 
  filter(wday ==2) %>% 
  select(date, temporalID)  
  
#add temporal ID to the observational dataset and merge the missing dates
tmp <- templatedf %>%
  mutate(temporalID=paste0(year,  "_",  stringr::str_pad(week(sampling_date), 2, pad="0"))) %>% 
  select(temporalID, value, delta_week) %>% 
  full_join(mySeq, by="temporalID" ) %>% 
  arrange(date)

tmp %>% print(n=52)

## -----------------------------------------------------------------------------
# apply spreader function
ex <- spreader(mydf = tmp, 
         date.field = "date", 
         value.field = "value", 
         counter.field = "delta_week",
         seed=123)
ex %>% print(n=52)

## ----warning=FALSE------------------------------------------------------------
cols <- c("Observed" = "#E1BE6A", "Post-processed" = "#40B0A6" )

dplyr::bind_rows(templatedf %>% 
  mutate(week = lubridate::week(sampling_date), 
         field = "Observed") %>% 
  select(week, value, field), 
  ex %>% 
  mutate(week = lubridate::week(date), 
         field = "Post-processed") %>% 
  select(week, value_adj, field) %>%   
  dplyr::rename(  value = value_adj)
  )%>% 
  ggplot(aes(week, value,  col=field, fill=field)) + 
  geom_bar(position="dodge", stat="identity", size=0.5, width = 0.8)+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  ylim(0, 250)+
  facet_wrap(~field,nrow=2)+
  labs(x="Week", y="Egg abundance" )+
  scale_x_continuous( breaks = seq(5, 50, by =5), limits=c(1, 52))+
  theme_classic()+
  theme(legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        text = element_text(size=16), 
        strip.text = element_text(size=16),
        legend.text = element_text(size=16,angle = 0), legend.title = element_text(size=16),
        legend.key.size = unit(1.5, 'cm'))

## -----------------------------------------------------------------------------
# apply spreader function
ex_noSampl <- spreader(mydf = tmp, 
         date.field = "date", 
         value.field = "value", 
         counter.field = NULL,
         seed=123)

## ----warning=FALSE------------------------------------------------------------
plot(ex$value_adj, ex_noSampl$value_adj, 
     xlab ="User-defined counter.field", 
     ylab ="Automatic counter.field")
abline(a=0,b=1,lty=2)

