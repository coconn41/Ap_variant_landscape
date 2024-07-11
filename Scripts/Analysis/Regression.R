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
source(paste0(getwd(),'/Scripts/Analysis/Aspatial_models.R'))

print(aspatial_results)

print(All_models)
