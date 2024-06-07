# source the shapefile loading:
Urban_areas = read_sf(paste0(getwd(),'/Data/Urban_areas/2010UAUC/urban2010.shp')) %>%
  st_transform(.,32618)

LT_spat = Location_table %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618)
for(i in 1:nrow(LT_spat)){
  closest=min(st_distance(LT_spat[i,],Urban_areas)/1000)
  attributes(closest)=NULL
  if(i==1){closest2=closest}
  if(i>1){closest2 = c(closest2,closest)}
}
Location_table$Closest_UA = closest2
