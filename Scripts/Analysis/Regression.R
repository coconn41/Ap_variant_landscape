#####
# Load data and libraries
#####
source(paste0(getwd(),'/Scripts/Data_cleaning/Load_datasets.R'))

#####
# Load Moran Function
#####
source(paste0(getwd(),'/Scripts/Analysis/Local_Moran_fn.R'))

#####
# Aspatial patch models
####
source(paste0(getwd(),"/Scripts/Analysis/Aspatial_patch_models.R"))

#####
# Spatial patch models:
#####
source(paste0(getwd(),'/Scripts/Analysis/Spatial_patch_models.R'))
       
#####
# Aspatial landscape models:
#####
source(paste0(getwd(),'/Scripts/Analysis/Aspatial_landscape_models.R'))

#####
# Spatial landscape models;
#####
source(paste0(getwd(),'/Scripts/Analysis/Spatial_landscape_models.R'))

#####
# Combine results:
#####

if((AIC(adult_ha_aspatial)-AIC(adult_ha_spatial))>2){print("Spatial is better")}
if((AIC(adult_v1_aspatial)-AIC(adult_v1_spatial))>2){print("Spatial is better")}
if((AIC(nymph_ha_aspatial)-AIC(nymph_ha_spatial))>2){print("Spatial is better")}
if((AIC(nymph_v1_aspatial)-AIC(nymph_v1_spatial))>2){print("Spatial is better")}
if((AIC(adult_ha_aspatial_landscape)-AIC(adult_ha_spatial_landscape))>2){print("Spatial is better")}
if((AIC(adult_v1_aspatial_landscape)-AIC(adult_v1_spatial_landscape))>2){print("Spatial is better")}
if((AIC(nymph_ha_aspatial_landscape)-AIC(nymph_ha_spatial_landscape))>2){print("Spatial is better")}
if((AIC(nymph_v1_aspatial_landscape)-AIC(nymph_ha_spatial_landscape))>2){print("Spatial is better")}


final_model_results = rbind(final_patch_results,final_landscape_results) %>%
  mutate(exponentiated = exp(Estimate))
