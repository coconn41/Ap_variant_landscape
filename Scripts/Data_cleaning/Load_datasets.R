#####
# Load Libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
library(moranfast)
#####
# Read in and separate data
#####
remove_private = F
if(remove_private==T){
  Regression_df <- read_csv("Data/Regression_df/Regression_df.csv")[,-1] %>%
    mutate(Site = as.factor(Site),
           UNIT = as.factor(UNIT),
           Year = as.numeric(substring(Date,1,4)),
           competition = ha - v1) %>%
    filter(is.na(latitude)==F) %>%
    group_by(Site,Lifestage,Year) %>%
    summarize(County = unique(County),
              tot_collected = sum(tot_collected),
              ha = sum(ha),
              v1 = sum(v1),
              und = sum(und),
              tot_tested = sum(tot_tested),
              latitude = max(latitude),
              longitude = max(longitude),
              metric = max(metric),
              UNIT = unique(UNIT),
              Landscape_metric = max(Landscape_metric))}
if(remove_private==F){
  Regression_df <- read_csv("Data/Regression_df/Regression_df_w_private.csv")[,-1] %>%
    mutate(Site = as.factor(Site),
           UNIT = as.factor(UNIT),
           Year = as.numeric(substring(Date,1,4)),
           competition = ha - v1) %>%
    filter(is.na(latitude)==F) %>%
    group_by(Site,Lifestage,Year) %>%
    summarize(County = unique(County),
              tot_collected = sum(tot_collected),
              ha = sum(ha),
              v1 = sum(v1),
              und = sum(und),
              tot_tested = sum(tot_tested),
              latitude = max(latitude),
              longitude = max(longitude),
              metric = max(metric),
              UNIT = unique(UNIT),
              Landscape_metric = max(Landscape_metric))
}
Nymph_df = Regression_df %>%
  filter(Lifestage == "Nymph")

Adult_df = Regression_df %>%
  filter(Lifestage == "Adult")

Adult_df_spat = Adult_df %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618)

Nymph_df_spat = Nymph_df %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618)


Regression_df_w_density = read_csv(paste0(getwd(),
                                          "/Data/Regression_df/Regression_df_w_density.csv"))
  
  
  