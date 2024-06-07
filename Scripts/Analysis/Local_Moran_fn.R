#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
library(tmap)
#####
# Moran's I test
#####
Moran_maps = function(df,model,seed=F,seed_n){
  residuals = abs(residuals(model))
  if(seed==T){set.seed(seed_n)}
  latlong = data.frame(x=df$longitude + runif(n = nrow(df),min = -.00001,max=.00001),
                       y=df$latitude + runif(n = nrow(df),min = -.00001,max=.00001))
  
  xy=latlong[,c(1,2)]
  latlongdf = SpatialPointsDataFrame(coords=xy,data=latlong)
  knn5list = knearneigh(x=latlongdf, k=5, longlat=T)
  knn5 = knn2nb(knn5list, sym = TRUE)
  df = df %>%
    mutate(moran = localmoran(residuals,listw = nb2listw(knn5))[,5]) %>%
    mutate(`P-value` = factor(ifelse(moran<0.05,'sig.','n.s.'),
                              levels = c('n.s.','sig.'))) %>%
    arrange(desc(moran))
tm_shape(df %>%
           st_as_sf(.,coords=c('longitude','latitude')) %>%
           st_set_crs(.,4326) %>%
           st_transform(.,32618))+
    tm_dots(col='P-value',size=1,
            palette=c('white','red'),
            shape=21)+
  tm_layout(main.title = paste0(substring(model$call[2],1,2), " in ",
                                df$Lifestage, "s"))
}
