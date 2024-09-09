#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
library(mgcViz)
#####
# 1.) Adult ha model:
#####
adult_ha_interaction = mgcv::gam(ha ~ s(metric,Landscape_metric) +
                                         s(longitude,latitude) + 
                                         s(Site,bs='re'),
                                       family = ziP,
                                       offset = log(tot_tested),
                                       data = Adult_df)
summary(adult_ha_interaction)
plot(mgcViz::getViz(adult_ha_interaction))

adult_ha_interaction = mgcv::gam(ha ~ s(metric,Closest_UA) +
                                   s(longitude,latitude) + 
                                   s(Site,bs='re'),
                                 family = ziP,
                                 offset = log(tot_tested),
                                 data = Adult_df)
summary(adult_ha_interaction)
plot(mgcViz::getViz(adult_ha_interaction))

#####
# 2.) Adult v1 model:
#####

adult_v1_interaction = mgcv::gam(v1 ~ s(metric,Landscape_metric) +
                                   s(longitude,latitude) + 
                                   s(Site,bs='re'),
                                 family = ziP,
                                 offset = log(tot_tested),
                                 data = Adult_df)
summary(adult_v1_interaction)
plot(mgcViz::getViz(adult_v1_interaction))

#####
# 3.) Nymph ha model:
#####
nymph_ha_interaction = mgcv::gam(ha ~ s(metric,Landscape_metric) +
                                   s(longitude,latitude) + 
                                   s(Site,bs='re'),
                                 family = ziP,
                                 offset = log(tot_tested),
                                 data = Nymph_df)
summary(nymph_ha_interaction)
plot(mgcViz::getViz(nymph_ha_interaction))

#####
# 2.) Nymph v1 model:
#####

nymph_v1_interaction = mgcv::gam(v1 ~ s(metric,Landscape_metric) +
                                   s(longitude,latitude) + 
                                   s(Site,bs='re'),
                                 family = ziP,
                                 offset = log(tot_tested),
                                 data = Nymph_df)
summary(nymph_v1_interaction)
plot(mgcViz::getViz(nymph_v1_interaction))



