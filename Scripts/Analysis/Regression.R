#####
# Load Libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)
#####
# Read in and separate data
#####
Regression_df <- read_csv("Data/Regression_df/Regression_df.csv")[,-1] %>%
  mutate(Site = as.factor(Site),
         UNIT = as.factor(UNIT),
         Year = as.numeric(substring(Date,1,4)))
Nymph_df = Regression_df %>%
  filter(Lifestage == "Nymph")

Adult_df = Regression_df %>%
  filter(Lifestage == "Adult")

#####
# A spatial models: (patch_level)
#####

adult_ha_zip = glmmTMB(formula = ha ~ metric + Closest_UA + 
                         (1|Site),
                       zi = ~1,
                       offset = log(tot_tested),
                       family = poisson,
                       data = Adult_df)
summary(adult_ha_zip)
performance::check_overdispersion(adult_ha_zip)
AIC(adult_ha_zip)
# 4857.615 with UA
# 4868.301 without UA

adult_v1_zip = glmmTMB(formula = v1 ~ metric + #Closest_UA + 
                         (1|Site),
                       zi = ~1,
                       offset = log(tot_tested),
                       family = poisson,
                       data = Adult_df)
summary(adult_v1_zip)
performance::check_overdispersion(adult_v1_zip)
AIC(adult_v1_zip)
# 2769.708 with UA
# 2770.345 without UA

nymph_ha_zip = glmmTMB(formula = ha ~ metric + Closest_UA +
                         (1|Site),
                       zi = ~1,
                       offset = log(tot_tested),
                       family = poisson,
                       data = Nymph_df)
summary(nymph_ha_zip)
performance::check_overdispersion(nymph_ha_zip)
AIC(nymph_ha_zip)
#2024.564 with UA
# 2022.879 without UA
# better without

nymph_v1_zip = glmmTMB(formula = v1 ~ metric + Closest_UA +
                         (1|Site),
                       zi = ~1,
                       offset = log(tot_tested),
                       family = poisson,
                       data = Nymph_df)
summary(nymph_v1_zip)
performance::check_overdispersion(nymph_ha_zip)
AIC(nymph_v1_zip)
# 2007.761 with UA
# 2012.735 without UA

a_spatial_parameters = data.frame(Nymph = c(summary(nymph_ha_zip)$coefficients$cond[2],
                                            summary(nymph_v1_zip)$coefficients$cond[2]),
                                  Adult = c(summary(adult_ha_zip)$coefficients$cond[2],
                                            summary(adult_v1_zip)$coefficients$cond[2]))
row.names(a_spatial_parameters)=c('ha','v1')
a_spatial_parameters
a_spatial_p_values = data.frame(Nymph = c(summary(nymph_ha_zip)$coefficients$cond[11],
                                          summary(nymph_v1_zip)$coefficients$cond[11]),
                                Adult = c(summary(adult_ha_zip)$coefficients$cond[11],
                                          summary(adult_v1_zip)$coefficients$cond[11]))
row.names(a_spatial_p_values)=c("ha","v1")
a_spatial_p_values

#####
# Test residuals for spatial autocorrelation:
#####
Adult_df_spat = Adult_df %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618)

Nymph_df_spat = Nymph_df %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,4326) %>%
  st_transform(.,32618)

Adult_df_spat$ha_resid=abs(residuals(adult_ha_zip))
Adult_df_spat$v1_resid=abs(residuals(adult_v1_zip))
Nymph_df_spat$ha_resid=abs(residuals(nymph_ha_zip))
Nymph_df_spat$v1_resid=abs(residuals(nymph_v1_zip))
Adult_df$ha_resid=abs(residuals(adult_ha_zip))
Adult_df$v1_resid=abs(residuals(adult_v1_zip))
Nymph_df$ha_resid=abs(residuals(nymph_ha_zip))
Nymph_df$v1_resid=abs(residuals(nymph_v1_zip))

resid_combined = Adult_df_spat %>%
  pivot_longer(c(ha_resid,
                 v1_resid),
               names_to = "Genotype",
               values_to = "Residuals") %>%
  mutate(Genotype = ifelse(Genotype=="ha_resid","ha","v1")) %>%
  bind_rows(.,Nymph_df_spat %>%
              pivot_longer(c(ha_resid,
                             v1_resid),
                           names_to = "Genotype",
                           values_to = "Residuals") %>%
              mutate(Genotype = ifelse(Genotype=="ha_resid","ha","v1"))) %>%
  arrange(Residuals)

tm_shape(NYS)+
  tm_borders()+
tm_shape(resid_combined)+
  tm_dots(col='Residuals',style='cont',size=1,shape=21)+
  tm_facets(by=c("Genotype","Lifestage"))
set.seed(1)
latlong_adult_ha = data.frame(x=Adult_df$longitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001),
                              y=Adult_df$latitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001))

xy_adult_ha=latlong_adult_ha[,c(1,2)]
latlong_adult_hadf = SpatialPointsDataFrame(coords=xy_adult_ha,data=latlong_adult_ha)
knn5list = knearneigh(x=latlong_adult_hadf, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Adult_df_spat$ha_moran=localmoran(Adult_df$ha_resid,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_adult_v1 = data.frame(x=Adult_df$longitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001),
                              y=Adult_df$latitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001))

xy_adult_v1=latlong_adult_v1[,c(1,2)]
latlong_adult_v1df = SpatialPointsDataFrame(coords=xy_adult_v1,data=latlong_adult_v1)
knn5list = knearneigh(x=latlong_adult_v1df, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Adult_df_spat$v1_moran=localmoran(Adult_df$v1_resid,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_nymph_ha = data.frame(x=Nymph_df$longitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001),
                              y=Nymph_df$latitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001))

xy_nymph_ha=latlong_nymph_ha[,c(1,2)]
latlong_nymph_hadf = SpatialPointsDataFrame(coords=xy_nymph_ha,data=latlong_nymph_ha)
knn5list = knearneigh(x=latlong_nymph_hadf, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Nymph_df_spat$ha_moran=localmoran(Nymph_df$ha_resid,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_nymph_v1 = data.frame(x=Nymph_df$longitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001),
                              y=Nymph_df$latitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001))

xy_nymph_v1=latlong_nymph_v1[,c(1,2)]
latlong_nymph_v1df = SpatialPointsDataFrame(coords=xy_nymph_v1,data=latlong_nymph_v1)
knn5list = knearneigh(x=latlong_nymph_v1df, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Nymph_df_spat$v1_moran=localmoran(Nymph_df$v1_resid,listw = nb2listw(knn5))[,5]

Combine_moran_results = Adult_df_spat %>%
  pivot_longer(c(ha_moran,
                 v1_moran),
               names_to = "Genotype",
               values_to = "P-value") %>%
  mutate(Genotype = ifelse(Genotype=="ha_moran","ha","v1")) %>%
  bind_rows(.,Nymph_df_spat %>%
              pivot_longer(c(ha_moran,
                             v1_moran),
                           names_to = "Genotype",
                           values_to = "P-value") %>%
              mutate(Genotype = ifelse(Genotype=="ha_moran","ha","v1"))) %>%
  mutate(p_val_cat = `P-value`,
          `P-value` = factor(ifelse(`P-value`<0.05,"sig.",'n.s.'),
                            levels=c('n.s.','sig.'))) %>%
  arrange(desc(p_val_cat))

tm_shape(NYS)+
  tm_borders() +
tm_shape(Combine_moran_results)+
  tm_dots(col='P-value',size=1,
          palette=c('white','red'),
          shape=21)+
  tm_facets(by=c("Genotype","Lifestage"))

#####
# Build spatial models: patch-level
#####

adult_ha_gam = gam(formula = ha ~ metric + #Closest_UA +
                     s(longitude,latitude)+
                     s(Site,bs='re'),
                   family = ziP,
                   offset = log(tot_tested),
                   data = Adult_df)
summary(adult_ha_gam)
AIC(adult_ha_gam) 
# 4474.195 with UA
# 4480.064 without UA
adult_v1_gam = gam(formula = v1 ~ metric + #Closest_UA +
                     s(longitude,latitude)+
                     s(Site,bs='re'),
                   family = ziP,
                   offset = log(tot_tested),
                   data = Adult_df)
summary(adult_v1_gam)

AIC(adult_v1_gam) 
# 2706.186 with UA
# 2706.043 without UA

nymph_ha_gam = gam(formula = ha ~ metric + #Closest_UA +
                     s(longitude,latitude)+
                     s(Site,bs='re'),
                   family = ziP,
                   offset = log(tot_tested),
                   data = Nymph_df)
summary(nymph_ha_gam)

AIC(nymph_ha_gam) 
# 1936.963 with UA
# 1935.41 without UA

nymph_v1_gam = gam(formula = v1 ~ metric +# Closest_UA +
                     s(longitude,latitude)+
                     s(Site,bs='re'),
                   family = ziP,
                   offset = log(tot_tested),
                   data = Nymph_df)
summary(nymph_v1_gam)

AIC(nymph_v1_gam)
# 1940.914 with UA
# 1940.998 without UA
spatial_parameters = data.frame(Nymph = c(summary(nymph_ha_gam)$p.coeff[2],
                                          summary(nymph_v1_gam)$p.coeff[2]),
                                Adult = c(summary(adult_ha_gam)$p.coeff[2],
                                          summary(adult_v1_gam)$p.coeff[2]))
row.names(spatial_parameters)=c('ha','v1')
spatial_parameters
spatial_p_values = data.frame(Nymph = c(summary(nymph_ha_gam)$p.pv[2],
                                          summary(nymph_v1_gam)$p.pv[2]),
                                Adult = c(summary(adult_ha_gam)$p.pv[2],
                                          summary(adult_v1_gam)$p.pv[2]))
row.names(spatial_p_values)=c("ha","v1")
spatial_p_values
#####
# Build a-spatial models: Landscape-level
#####

adult_ha_zip_landscape = glmmTMB(formula = ha ~ Landscape_metric +# Closest_UA + 
                                   (1|UNIT) +
                                   (1|UNIT:Site),
                       zi = ~1,
                       offset = log(tot_tested),
                       family = poisson,
                       data = Adult_df)
summary(adult_ha_zip_landscape)
performance::check_overdispersion(adult_ha_zip_landscape)
AIC(adult_ha_zip_landscape)
#4753.169 with UA
#4751.331 without UA


adult_v1_zip_landscape = glmmTMB(formula = v1 ~ Landscape_metric + #Closest_UA + 
                                   (1|UNIT) +
                                   (1|UNIT:Site),
                                 zi = ~1,
                                 offset = log(tot_tested),
                                 family = poisson,
                                 data = Adult_df)
summary(adult_v1_zip_landscape)
performance::check_overdispersion(adult_v1_zip_landscape)
AIC(adult_v1_zip_landscape)
# 2766.925 with UA
# 2766.118 without UA


nymph_ha_zip_landscape = glmmTMB(formula = ha ~ Landscape_metric +# Closest_UA + 
                                   (1|UNIT) +
                                   (1|UNIT:Site),
                                 zi = ~1,
                                 offset = log(tot_tested),
                                 family = poisson,
                                 data = Nymph_df)
summary(nymph_ha_zip_landscape)
performance::check_overdispersion(nymph_ha_zip_landscape)
AIC(nymph_ha_zip_landscape)
# 2000.654 with UA
# 1998.654 without UA
2000.654-1998.654

nymph_v1_zip_landscape = glmmTMB(formula = v1 ~ Landscape_metric + Closest_UA + 
                                   (1|UNIT) +
                                   (1|UNIT:Site),
                                 zi = ~1,
                                 offset = log(tot_tested),
                                 family = poisson,
                                 data = Nymph_df)
summary(nymph_v1_zip_landscape)
performance::check_overdispersion(nymph_v1_zip_landscape)
AIC(nymph_v1_zip_landscape)
# 2008.286 with UA
# 2012.991 without UA

A_spatial_landscape_parameters = data.frame(Nymph = c(summary(nymph_ha_zip_landscape)$coefficients$cond[2],
                                            summary(nymph_v1_zip_landscape)$coefficients$cond[2]),
                                  Adult = c(summary(adult_ha_zip_landscape)$coefficients$cond[2],
                                            summary(adult_v1_zip_landscape)$coefficients$cond[2]))
row.names(A_spatial_landscape_parameters)=c('ha','v1')
A_spatial_landscape_parameters
A_spatial_landscape_p_values = data.frame(Nymph = c(summary(nymph_ha_zip_landscape)$coefficients$cond[11],
                                                    summary(nymph_v1_zip_landscape)$coefficients$cond[11]),
                                          Adult = c(summary(adult_ha_zip_landscape)$coefficients$cond[11],
                                                    summary(adult_v1_zip_landscape)$coefficients$cond[11]))
row.names(A_spatial_landscape_p_values)=c("ha","v1")
A_spatial_landscape_p_values

#####
# Test for spatial autocorrelation: Landscape-level
#####

Adult_df_spat$ha_resid_landscape=abs(residuals(adult_ha_zip_landscape))
Adult_df_spat$v1_resid_landscape=abs(residuals(adult_v1_zip_landscape))
Nymph_df_spat$ha_resid_landscape=abs(residuals(nymph_ha_zip_landscape))
Nymph_df_spat$v1_resid_landscape=abs(residuals(nymph_v1_zip_landscape))
Adult_df$ha_resid_landscape=abs(residuals(adult_ha_zip_landscape))
Adult_df$v1_resid_landscape=abs(residuals(adult_v1_zip_landscape))
Nymph_df$ha_resid_landscape=abs(residuals(nymph_ha_zip_landscape))
Nymph_df$v1_resid_landscape=abs(residuals(nymph_v1_zip_landscape))

resid_combined_landscape = Adult_df_spat %>%
  pivot_longer(c(ha_resid_landscape,
                 v1_resid_landscape),
               names_to = "Genotype",
               values_to = "Residuals") %>%
  mutate(Genotype = ifelse(Genotype=="ha_resid_landscape","ha","v1")) %>%
  bind_rows(.,Nymph_df_spat %>%
              pivot_longer(c(ha_resid_landscape,
                             v1_resid_landscape),
                           names_to = "Genotype",
                           values_to = "Residuals") %>%
              mutate(Genotype = ifelse(Genotype=="ha_resid_landscape","ha","v1"))) %>%
  arrange(Residuals)

tm_shape(NYS)+
  tm_borders()+
  tm_shape(resid_combined_landscape)+
  tm_dots(col='Residuals',style='cont',size=1,shape=21)+
  tm_facets(by=c("Genotype","Lifestage"))

set.seed(1)
latlong_adult_ha = data.frame(x=Adult_df$longitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001),
                              y=Adult_df$latitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001))

xy_adult_ha=latlong_adult_ha[,c(1,2)]
latlong_adult_hadf = SpatialPointsDataFrame(coords=xy_adult_ha,data=latlong_adult_ha)
knn5list = knearneigh(x=latlong_adult_hadf, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Adult_df_spat$ha_moran_landscape=localmoran(Adult_df$ha_resid_landscape,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_adult_v1 = data.frame(x=Adult_df$longitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001),
                              y=Adult_df$latitude + runif(n = nrow(Adult_df),min = -.00001,max=.00001))

xy_adult_v1=latlong_adult_v1[,c(1,2)]
latlong_adult_v1df = SpatialPointsDataFrame(coords=xy_adult_v1,data=latlong_adult_v1)
knn5list = knearneigh(x=latlong_adult_v1df, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Adult_df_spat$v1_moran_landscape=localmoran(Adult_df$v1_resid_landscape,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_nymph_ha = data.frame(x=Nymph_df$longitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001),
                              y=Nymph_df$latitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001))

xy_nymph_ha=latlong_nymph_ha[,c(1,2)]
latlong_nymph_hadf = SpatialPointsDataFrame(coords=xy_nymph_ha,data=latlong_nymph_ha)
knn5list = knearneigh(x=latlong_nymph_hadf, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Nymph_df_spat$ha_moran_landscape=localmoran(Nymph_df$ha_resid_landscape,listw = nb2listw(knn5))[,5]

set.seed(1)
latlong_nymph_v1 = data.frame(x=Nymph_df$longitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001),
                              y=Nymph_df$latitude + runif(n = nrow(Nymph_df),min = -.00001,max=.00001))

xy_nymph_v1=latlong_nymph_v1[,c(1,2)]
latlong_nymph_v1df = SpatialPointsDataFrame(coords=xy_nymph_v1,data=latlong_nymph_v1)
knn5list = knearneigh(x=latlong_nymph_v1df, k=5, longlat=T)
knn5 = knn2nb(knn5list, sym = TRUE)

Nymph_df_spat$v1_moran_landscape=localmoran(Nymph_df$v1_resid_landscape,listw = nb2listw(knn5))[,5]

Combine_moran_results_landscape = Adult_df_spat %>%
  pivot_longer(c(ha_moran_landscape,
                 v1_moran_landscape),
               names_to = "Genotype",
               values_to = "P-value") %>%
  mutate(Genotype = ifelse(Genotype=="ha_moran_landscape","ha","v1")) %>%
  bind_rows(.,Nymph_df_spat %>%
              pivot_longer(c(ha_moran_landscape,
                             v1_moran_landscape),
                           names_to = "Genotype",
                           values_to = "P-value") %>%
              mutate(Genotype = ifelse(Genotype=="ha_moran_landscape","ha","v1"))) %>%
  mutate(p_val_cat = `P-value`,
         `P-value` = factor(ifelse(`P-value`<0.05,"sig.",'n.s.'),
                            levels=c('n.s.','sig.'))) %>%
  arrange(desc(p_val_cat))

tm_shape(NYS)+
  tm_borders() +
  tm_shape(Combine_moran_results_landscape)+
  tm_dots(col='P-value',size=1,
          palette=c('white','red'),
          shape=21)+
  tm_facets(by=c("Genotype","Lifestage"))

#####
# Spatial models landscape:
#####

adult_ha_gam_landscape = gam(formula = ha ~ Landscape_metric + Closest_UA +
                     s(longitude,latitude)+
                     s(Site,bs='re'),
                   family = ziP,
                   offset = log(tot_tested),
                   data = Adult_df)
summary(adult_ha_gam_landscape)
AIC(adult_ha_gam_landscape)
# 4477.074 with UA
# 4480.435 without UA

adult_v1_gam_landscape = gam(formula = v1 ~ Landscape_metric +# Closest_UA +
                               s(longitude,latitude)+
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Adult_df)
summary(adult_v1_gam_landscape)
AIC(adult_v1_gam_landscape)
#2707.11 with UA
# 2706.7 without UA

nymph_ha_gam_landscape = gam(formula = ha ~ Landscape_metric +# Closest_UA +
                               s(longitude,latitude)+
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Nymph_df)
summary(nymph_ha_gam_landscape)
AIC(nymph_ha_gam_landscape)
# 1934.906 with UA
# 1934.29 without UA


nymph_v1_gam_landscape = gam(formula = v1 ~ Landscape_metric +# Closest_UA +
                               s(longitude,latitude)+
                               s(Site,bs='re'),
                             family = ziP,
                             offset = log(tot_tested),
                             data = Nymph_df)
summary(nymph_v1_gam_landscape)
AIC(nymph_v1_gam_landscape)
# 1941.896 with UA
# 1941.896  without UA

spatial_landscape_parameters = data.frame(Nymph = c(summary(nymph_ha_gam_landscape)$p.coeff[2],
                                          summary(nymph_v1_gam_landscape)$p.coeff[2]),
                                Adult = c(summary(adult_ha_gam_landscape)$p.coeff[2],
                                          summary(adult_v1_gam_landscape)$p.coeff[2]))
row.names(spatial_landscape_parameters)=c('ha','v1')
spatial_landscape_parameters
spatial_p_values_landscape = data.frame(Nymph = c(summary(nymph_ha_gam_landscape)$p.pv[2],
                                        summary(nymph_v1_gam_landscape)$p.pv[2]),
                              Adult = c(summary(adult_ha_gam_landscape)$p.pv[2],
                                        summary(adult_v1_gam_landscape)$p.pv[2]))
row.names(spatial_p_values_landscape)=c("ha","v1")
spatial_p_values_landscape
