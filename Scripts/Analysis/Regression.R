#####
# Load data and libraries
#####
source(paste0(getwd(),'/Scripts/Data_cleaning/Load_datasets.R'))

#####
# Variant models:
#####
source(paste0(getwd(),'/Scripts/Analysis/Aspatial_models.R'))
print(aspatial_results)
#####
# Density models
#####

source(paste0(getwd(),'/Scripts/Analysis/SCR_to_density_models.R'))
print(final_scr_to_d_models)

source(paste0(getwd(),'/Scripts/Analysis/density_to_pathogen_models.R'))
print(final_d_to_p_models)

#####
# ERI models
#####

source(paste0(getwd(),'/Scripts/Analysis/ERI_models.R'))
print(final_ERI_models)
