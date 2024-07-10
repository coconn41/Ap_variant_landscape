#####
# Load libraries
#####
library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)

#####
# Load data
#####
source(paste0(getwd(),'/Scripts/Figures/NYS_shapefile.R'))
Regression_df <- read_csv("Data/Regression_df/Regression_df_w_private.csv")[,-1] %>%
  mutate(Site = as.factor(Site),
         UNIT = as.factor(UNIT),
         Year = as.numeric(substring(Date,1,4))) %>%
  filter(is.na(latitude)==F)

#####
# Spatialize sites and keep only one site:
#####

Regression_df_spat = Regression_df %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618) %>%
  mutate(id = 1:nrow(.)) %>%
  group_by(Site) %>%
  filter(id == max(id))

#####
# Make map of site level SCR
#####
new_bb = c(0,4483094.7,764133.3,4985489.9)
names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
attr(new_bb, "class") = "bbox"


m1 = tm_shape(NYS,bbox = st_bbox(new_bb)) + 
  tm_borders() +
tm_shape(Regression_df_spat %>%
           arrange(metric))+
  tm_dots(col = 'metric',
          title = "Sinuous Connection\nReduction (Patch)",
          shape = 21,
          size = .3,
          legend.hist = T,
          palette = get_brewer_pal("Greens"))


#####
# Make map of WMU level landscape SCR
#####

WMUs = read_sf(paste0(getwd(),'/Data/WMUs/Wildlife_Management_Units.shp')) %>%
  st_transform(.,crs = st_crs(NYS))

WMUs_4_choro = left_join(WMUs,Regression_df %>%
                           mutate(id = 1:nrow(.)) %>%
                           group_by(UNIT) %>%
                           filter(id == max(id)),
                         by = "UNIT") 

m2=tm_shape(WMUs_4_choro %>%
              st_simplify(dTolerance = 200),
            bbox = st_bbox(new_bb))+
  tm_add_legend(type = 'symbol',
                shape = 21,
                col = 'black',
                label = "Sampling sites")+
  tm_polygons(col = 'Landscape_metric',
              title = "Sinuous Connection\nReduction (Landscape)",
              palette = get_brewer_pal("Greens"))+
  tm_shape(Regression_df_spat)+
  tm_dots(col='black')+
  tm_layout(legend.position = c('left','bottom'))

m3 = tmap_arrange(m1,m2,ncol=1)

tmap_save(m3,
          filename = paste0(getwd(),'/Figures/Figure_4.jpeg'))
