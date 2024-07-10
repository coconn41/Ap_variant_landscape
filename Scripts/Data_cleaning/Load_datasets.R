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
    filter(is.na(latitude)==F)}
if(remove_private==F){
  Regression_df <- read_csv("Data/Regression_df/Regression_df_w_private.csv")[,-1] %>%
    mutate(Site = as.factor(Site),
           UNIT = as.factor(UNIT),
           Year = as.numeric(substring(Date,1,4)),
           competition = ha - v1) %>%
    filter(is.na(latitude)==F)
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