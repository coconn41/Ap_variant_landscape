library(tidyverse)
library(sf)
library(readxl)

#####
# Load in and clean location table
#####

Location_table <- read_excel("Data/TL_data/Location_table_10_15.xlsx") %>%
  mutate(loc_latitude = ifelse(Location_ID==728,42.903655,loc_latitude)) %>%
  filter(is.na(loc_longitude)==F,
         is.na(loc_latitude)==F,
         loc_public!="Private property",
         !grepl('private', loc_public),
         !grepl('Private', loc_public)) %>%
  mutate(longitude = as.numeric(loc_longitude)*-1,
         latitude = as.numeric(loc_latitude)) %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(4326) %>%
  st_transform(32618)

Location_table_with_private = read_excel("Data/TL_data/Location_table_10_15.xlsx") %>%
  mutate(loc_latitude = ifelse(Location_ID==728,42.903655,loc_latitude)) %>%
  filter(is.na(loc_longitude)==F,
         is.na(loc_latitude)==F) %>%
  mutate(longitude = as.numeric(loc_longitude)*-1,
         latitude = as.numeric(loc_latitude)) %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(4326) %>%
  st_transform(32618)

#####
# Load in metric forest patches for all NYS
#####

# for(i in 1:12){
#   metric2 = read_sf(paste0(getwd(),'/Data/SCR_Patch/Metric_results',i,'.shp'))
#   if(i==1){metric = read_sf(paste0(getwd(),'/Data/SCR_Patch/Metric_results',i,'.shp'))}
#   if(i>1){metric = rbind(metric,metric2)}
# }
# remove(metric2)
# metric = metric %>%
#   st_transform(.,crs=st_crs(Location_table))
metric = read_sf(paste0(getwd(),"/Data/SCR_patch_full/Metric_2024-06-07.shp")) %>%
  st_transform(.,crs=32618)
#####
# Check the distribution of forest patches
#####
mdn = median(metric$area);mdn
mn = mean(metric$area);mn

ggplot(data=metric,
       aes(x=area))+
  geom_histogram(bins=1000)+
  coord_cartesian(xlim=c(0,1000))+
  geom_vline(xintercept=mn,color='blue')+
  geom_vline(xintercept=mdn,color='red')

#####
# Find closest forest patches greater than the mean forest patch size
#####
red_metric = metric %>% filter(area>=mn)
Location_table$layer=rep(1,nrow(Location_table))
for(i in 1:nrow(Location_table)){
Location_table[i,20]=red_metric[st_nearest_feature(Location_table[i,],
                        red_metric),]$layer
}
Location_table = left_join(Location_table,
                            red_metric %>%
                              st_drop_geometry(),
                            by = c('layer')) %>%
  st_drop_geometry() %>%
  mutate(longitude = as.numeric(loc_longitude)*-1,
         latitude = as.numeric(loc_latitude)) %>%
  dplyr::select(Location_ID,loc_county,latitude,longitude,loc_name,layer,metric,area) 

write.csv(Location_table,
          file = paste0(getwd(),'/Data/Metric_locations/Loc_metric_table.csv'))

Location_table_with_private$layer=rep(1,nrow(Location_table_with_private))
for(i in 1:nrow(Location_table_with_private)){
  Location_table_with_private[i,20]=red_metric[st_nearest_feature(Location_table_with_private[i,],
                                                     red_metric),]$layer
}
Location_table_with_private = left_join(Location_table_with_private,
                           red_metric %>%
                             st_drop_geometry(),
                           by = c('layer')) %>%
  st_drop_geometry() %>%
  mutate(longitude = as.numeric(loc_longitude)*-1,
         latitude = as.numeric(loc_latitude)) %>%
  dplyr::select(Location_ID,loc_county,latitude,longitude,loc_name,layer,metric,area) 

write.csv(Location_table_with_private,
          file = paste0(getwd(),'/Data/Metric_locations/Loc_metric_table_w_private.csv'))
