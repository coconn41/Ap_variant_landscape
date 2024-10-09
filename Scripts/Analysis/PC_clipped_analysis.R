library(rlang)
library(sf)
library(sp)
library(terra)
library(FedData)
library(landscapemetrics)
library(leastcostpath)
library(dplyr)
library(tidyr)
library(readr)
library(units)
library(geosphere)
library(foreach)
library(parallel)
library(readxl)
library(igraph)
library(tmap)
library(tmaptools)
#####
# Parameters to change for analysis:
#####
distance = 1675 # Use this to truncate the LCPs
cores = 44
theta = .001788497 # Theta in PC index calculation
tdir=paste0(getwd(),"/LCP_analysis/Input_data/")

#####
# WMU Shapefile Download:
#####

Resource_polygons = read_sf(paste0(getwd(),"/Data/PC_index_clipped/Input_data/New_polygons.shp")) %>%
  dplyr::select(layer,area,geometry)

Loc_table_sf = read_sf(paste0(getwd(),
                              "/Data/PC_index_clipped/Input_data/Collection_site_nodes.shp")) %>%
  select(Lctn_ID,County,Site,geometry)

#####
# Read in NYS
#####

NYS = read_sf(paste0(getwd(),"/Data/PC_index_clipped/Input_data/cb_2018_us_state_500k.shp")) %>%
  filter(NAME=="New York") %>%
  st_transform(.,crs=32618)

#####
# Read in land cover data
#####
LCr = rast(paste0(getwd(),"/Data/PC_index_clipped/Input_data/NYS NLCD_NLCD_Land_Cover_2019.tif"))
LCproj = terra::project(LCr,crs(Resource_polygons))

LCcrop = terra::crop(x = LCproj,
                     y = Resource_polygons |>
                       terra::vect(),
                     mask = T)

#####
# Road Download and reclassification:
#####

minor_roadways = c(7,8,9,17,18,19)
major_roadways = c(1,2,4,6,11,12,14,16)

roadwaydat = read_sf(paste0(getwd(),
                            "/Data/PC_index_clipped/Input_data/NYSDOT Functional Class 07_10_2023/NYSDOT_Functional_Class_07_10_2023.shp")) %>%
  st_zm() %>%
  mutate(Categorization = ifelse(FUNCTIONAL %in% minor_roadways,1,
                                 ifelse(FUNCTIONAL %in% major_roadways,2,0))) %>%
  st_transform(.,crs=st_crs(NYS)) %>% 
  filter(!st_is_empty(.))
road_vect = terra::vect(roadwaydat)
road_pix  = terra::rasterize(x = road_vect,
                             y = LCcrop,
                             field = "Categorization")

#####
# Alpine zone processing:
#####

alpine_regions = read_sf(paste0(getwd(),"/Data/PC_index_clipped/Input_data/ny_eco_14/ny_eco_l4.shp")) %>% 
  filter(US_L4CODE=="58j") %>%
  mutate(US_L4CODE=1000) %>%
  st_transform(crs=32618) %>% 
  filter(!st_is_empty(.))

alpine_regions_vect = terra::vect(alpine_regions)
alpine_pix = terra::rasterize(x = alpine_regions_vect,
                              y = LCcrop,
                              field = "US_L4CODE",
                              background = 0)
#####
# Node processing:
#####
LC_forest_patches = LCcrop
values(LC_forest_patches)[values(LC_forest_patches)==42] = 41
values(LC_forest_patches)[values(LC_forest_patches)==43] = 41
values(LC_forest_patches)[values(LC_forest_patches)!=41] = NA

#####
# Read in habitat
#####

# Habitat = raster::raster(paste0(getwd(),
#                                 "/PCindex/Clipped/Input_data/Habitat.tif"))

y = get_patches(LC_forest_patches,directions=4)
poly1 = as.polygons(terra::rast(y$layer_1$class_41))
poly2 = st_as_sf(poly1)
poly2$area = st_area(poly2)
poly2$area = set_units(poly2$area,ha)
attributes(poly2$area)=NULL

minpoly = min(poly2$area)
poly2 = poly2 %>%
  filter(area >0.09)
minpoly = min(poly2$area)

poly2 = poly2 %>%
  filter(area > 0.18)
nodes = st_centroid(poly2)

#####
# Create conductance grid
#####

LC_forest = LCcrop
forest_values = c(41,42,43,51,52,71)
values(LC_forest)[!(values(LC_forest)%in%forest_values)] = 0
values(LC_forest)[values(LC_forest)%in%forest_values] = 1

LC_cropland = LCcrop
cropland_values = c(81,82)
values(LC_cropland)[!(values(LC_cropland)%in%cropland_values)] = 0
values(LC_cropland)[values(LC_cropland)%in%cropland_values] = 35

LC_wetland = LCcrop
wetland_values = c(90,95)
values(LC_wetland)[!(values(LC_wetland)%in%wetland_values)] = 0
values(LC_wetland)[values(LC_wetland)%in%wetland_values] = 100

LC_water = LCcrop
water_values = c(11)
values(LC_water)[!(values(LC_water)%in%wetland_values)] = 0
values(LC_water)[values(LC_water)%in%wetland_values] = 1000

LC_high_developed = LCcrop
high_developed_values = c(24)
values(LC_high_developed)[!(values(LC_high_developed)%in%high_developed_values)] = 0
values(LC_high_developed)[values(LC_high_developed)%in%high_developed_values] = 1000

LC_med_developed = LCcrop
med_developed_values = c(23)
values(LC_med_developed)[!(values(LC_med_developed)%in%med_developed_values)] = 0
values(LC_med_developed)[values(LC_med_developed)%in%med_developed_values] = 100

LC_low_developed = LCcrop
low_developed_values = c(21,22)
values(LC_low_developed)[!(values(LC_low_developed)%in%low_developed_values)] = 0
values(LC_low_developed)[values(LC_low_developed)%in%low_developed_values] = 27

#Alpine is already classified

Highways = road_pix
highway_values = c(2)
values(Highways)[!(values(Highways)%in%highway_values)] = 0
values(Highways)[values(Highways)%in%highway_values] = 533

minor_roads = road_pix
minor_road_values = c(1)
values(minor_roads)[!(values(minor_roads)%in%minor_road_values)] = 0
values(minor_roads)[values(minor_roads)%in%minor_road_values] = 100

Resistance_grid = max(LC_forest,LC_cropland,LC_wetland,LC_water,LC_high_developed,
                      LC_med_developed,LC_low_developed,Highways,minor_roads,
                      alpine_pix,na.rm=T)
Resistance_grid[Resistance_grid==0]=NA
Resistance_grid=1/Resistance_grid

#####
# Loop through polygons
#####
xind=0
for(x in 1:nrow(Resource_polygons)){
  xind=xind+1
  print(paste0(x," out of ",nrow(Resource_polygons)))
  sing_RP = Resource_polygons[x,]
  sing_RP$area=st_area(sing_RP)
  sing_RP$area = set_units(sing_RP$area,ha)
  attributes(sing_RP$area)=NULL
  
  
  #####
  # Conductance grid crop:
  ####
  
  new_poly = bb_poly(st_bbox(sing_RP))
  R_crop = terra::crop(x = Resistance_grid,
                       y = new_poly |>
                         terra::vect(),
                       mask = T)
  
  #####
  # Node processing:
  #####
  
  wmu_nodes = sf::st_intersection(sing_RP,nodes) %>%
    dplyr::select(layer.1,area.1,geometry) %>%
    mutate(layer = layer.1,
           area = area.1) %>%
    dplyr::select(layer,area,geometry)
  patch_crop = sf::st_intersection(sing_RP,poly2) %>%
    dplyr::select(layer.1,area.1,geometry) %>%
    mutate(layer = layer.1,
           area = st_area(.)) %>%
    dplyr::select(layer,area,geometry)
  patch_crop$area = set_units(patch_crop$area,ha)
  
  #####
  # constant is equal to -.001788497
  
  # Loop through unbuffered WMUs:
  
  
  
  #Create your comparisons  
  comb_mat = matrix(nrow=nrow(wmu_nodes),
                    ncol = nrow(wmu_nodes))  
  combinations <- which(upper.tri(comb_mat,
                                  diag = FALSE),
                        arr.ind = TRUE)
  combinations=as.data.frame(cbind(combinations, comb_mat[combinations]))
  names(combinations)=c('i','j','Nas')
  combinations = combinations %>%
    arrange(i) %>%
    dplyr::select(-Nas)
  
  #####
  # Calculate lcps:
  #####
  comps = combinations  
  
  
  tr1=leastcostpath::create_cs(terra::crop(x=R_crop,
                                           y = terra::vect(sing_RP),
                                           mask=T))
  nodes = nodes %>%
    st_transform(crs=st_crs(sing_RP))
  
  myCluster <- parallel::makeCluster(cores-10)
  doParallel::registerDoParallel(myCluster)
  
  lcp_network <- foreach::foreach(i = 1:nrow(comps),
                                  .verbose = F,
                                  .errorhandling = "remove",
                                  .combine = "rbind",
                                  .packages = c("sf",
                                                "raster",
                                                "gdistance",
                                                "tmaptools",
                                                "dplyr",
                                                "leastcostpath",
                                                "terra")) %dopar% {
                                                  
                                                  
                                                  
                                                  lcp <- leastcostpath::create_lcp(x = tr1,
                                                                                   origin = wmu_nodes[comps[i,1],drop=FALSE],
                                                                                   destination = wmu_nodes[comps[i,2], drop=FALSE],
                                                                                   check_locations = F) %>%
                                                    sf::st_as_sf() %>%
                                                    dplyr::mutate(length = sf::st_length(.),
                                                                  origin_ID = poly2[comps[i,1],]$layer,
                                                                  destination_ID =poly2[comps[i,2],]$layer) %>%
                                                    st_drop_geometry() %>%
                                                    dplyr::select(length,origin_ID,destination_ID)
                                                  
                                                  return(lcp)
                                                }
  parallel::stopCluster(myCluster)    
  attributes(lcp_network$length)=NULL
  print('finished lcp calculation')
  #####    
  # Do the PC calculations
  #####
  
  myCluster <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(myCluster)
  numerator <- foreach::foreach(a = 1:nrow(combinations),
                                .verbose = F,
                                .errorhandling = "remove",
                                .combine = "rbind",
                                .packages = c("sf",
                                              "raster",
                                              "dplyr",
                                              "terra",
                                              "igraph",
                                              "units")) %dopar% {
                                                ipatch = patch_crop[combinations[a,1],drop=F]
                                                jpatch = patch_crop[combinations[a,2],drop=F]
                                                combodf = data.frame(i = combinations[a,1],
                                                                     j = combinations[a,2],
                                                                     i_area = ipatch$area,
                                                                     j_area = jpatch$area,
                                                                     distance = lcp_network[a,1]) #st_distance(fin_poly[combinations[a,1],],fin_poly[combinations[a,2],])
                                                return(combodf)
                                              }
  parallel::stopCluster(myCluster)
  unregister_dopar()
  attributes(numerator$distance)=NULL
  attributes(numerator$i_area)=NULL
  attributes(numerator$j_area)=NULL
  numerator$product_prob=exp(theta*numerator$distance)
  numerator$max_distance = NA
  print('finished numerator calculation')
  g=igraph::graph_from_data_frame(numerator,directed=TRUE,vertices = NULL)
  E(g)$weight = 1/(numerator$product_prob)
  myCluster <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(myCluster)
  full_max_df <- foreach::foreach(a = 1:nrow(wmu_nodes),
                                  .verbose = F,
                                  .errorhandling = "remove",
                                  .combine = "rbind",
                                  .packages = c("sf","raster","dplyr","terra","igraph","units")) %dopar% {
                                    i_ind=0
                                    list_i = shortest_paths(g,from=V(g)[a],to=V(g))
                                    for(b in 1:as.numeric(length(list_i$vpath))){
                                      if(length(list_i$vpath[[b]])==0){next}
                                      plength = length(list_i$vpath[[b]])
                                      if(plength<=1){next}
                                      i_ind=i_ind+1
                                      for(c in 1:plength){
                                        if(c==plength){break}
                                        if(c==1){combos = data.frame(i=as.numeric(list_i$vpath[[b]][c]),
                                                                     j=as.numeric(list_i$vpath[[b]][c+1]))}
                                        if(c>1){combos2 = data.frame(i=as.numeric(list_i$vpath[[b]][c]),
                                                                     j=as.numeric(list_i$vpath[[b]][c+1]))
                                        combos = rbind(combos2,combos)}
                                      }
                                      if(i_ind==1){ndf=left_join(as.data.frame(combos),as.data.frame(numerator),join_by(i,j)) %>%
                                        summarize(max_distance=sum(distance)) %>%#na.rm=T
                                        as.data.frame()
                                      ndf$i = as.numeric(list_i$vpath[[b]][1])
                                      ndf$j = as.numeric(list_i$vpath[[b]][plength])}
                                      if(i_ind>1){ndf2=left_join(as.data.frame(combos),as.data.frame(numerator),join_by(i,j)) %>%
                                        summarize(max_distance=sum(distance)) %>%#na.rm=T
                                        as.data.frame() 
                                      ndf2$i = as.numeric(list_i$vpath[[b]][1])
                                      ndf2$j = as.numeric(list_i$vpath[[b]][plength])
                                      ndf=rbind(ndf2,ndf)}
                                    }
                                    if(exists('ndf')==TRUE){return(ndf)}
                                  }
  parallel::stopCluster(myCluster)
  unregister_dopar()
  if(nrow(as.data.frame(full_max_df))==0){write.csv(x,paste0(getwd(),'/PCindex/Clipped/Output_data/Error_wmus/wmu_',x,'.csv'))
    next}
  numerator = as.data.frame(numerator) %>% dplyr::select(-max_distance)
  numerator = left_join(numerator,
                        as.data.frame(full_max_df) %>%
                          distinct(i,j,max_distance,.keep_all = TRUE),
                        join_by(i,j)) %>%
    distinct()
  numerator$product_prob = exp(-.001788497*numerator$max_distance)
  numerator$act_numerator = numerator$i_area*numerator$j_area*numerator$product_prob
  result_df = data.frame(Lctn_ID = Loc_table_sf[x,1],
                         index = "PC_index",
                         value = sum(numerator$act_numerator*2)/(sing_RP$area^2))
  write.csv(result_df,paste0(getwd(),'/Data/PCindex_clipped/Output_data/Results/Result_df_wmu_',x,".csv"))
}