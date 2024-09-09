#####
# Load data and libraries
#####
source(paste0(getwd(),'/Scripts/Data_cleaning/Load_datasets.R'))

#####
# Load Moran Function
#####
source(paste0(getwd(),'/Scripts/Analysis/Local_Moran_fn.R'))

#####
# Aspatial models:
#####
# Option for including local moran is within the Aspatial_models script
# As it runs, it tests for global spatial autocorrelation directly in this file
# Spatial autocorrelation is not present, so aspatial models are the remaining best
source(paste0(getwd(),'/Scripts/Analysis/Aspatial_models.R'))
print(aspatial_results)
print(All_models)
#length(which(moran_df$p_adjust<0.05))
#####
# Density models
#####

source(paste0(getwd(),'/Scripts/Analysis/SCR_to_density_models.R'))
print(final_scr_to_d_models)
#length(which(morandf_fin$p_adjust<0.05))/nrow(morandf_fin)


source(paste0(getwd(),'/Scripts/Analysis/density_to_pathogen_models.R'))
print(final_d_to_p_models)
#length(which(moran_df$padjust<0.05))/nrow(moran_df)

#####
# ERI models
#####

source(paste0(getwd(),'/Scripts/Analysis/ERI_models.R'))
print(final_ERI_models)
