#####
# Load libraries
#####
library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)

#####
# Load data:
#####
source(paste0(getwd(),'/Scripts/Figures/NYS_shapefile.R'))
Regression_df <- read_csv("Data/Regression_df/Regression_df_w_private.csv")[,-1] %>%
  mutate(Site = as.factor(Site),
         UNIT = as.factor(UNIT),
         Year = as.numeric(substring(Date,1,4))) %>%
  filter(is.na(latitude)==F)

Regression_df_long = Regression_df %>%
  pivot_longer(cols=c(ha,v1)) %>%
  mutate(Weighted_Prevalence = ifelse(tot_tested>=50,(value/tot_tested)*100,
                                      (value/50)*(value/tot_tested)),
         Prevalence = (value/tot_tested)*100)

Adult_df_long = Regression_df_long %>%
  filter(Target_stage == "Adult")

Regression_df_spat = Regression_df %>%
  pivot_longer(cols=c(ha,v1)) %>%
  mutate(Weighted_Prevalence = ifelse(tot_tested>=50,(value/tot_tested)*100,
                                      (value/50)*(value/tot_tested)),
         Prevalence = (value/tot_tested)*100) %>%
  st_as_sf(.,coords=c('longitude','latitude')) %>%
  st_set_crs(.,value=4326) %>%
  st_transform(.,crs=32618)

m1 = tm_shape(NYS)+
  tm_borders()+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Adult",
                    name == "ha") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence) %>%
             filter(mean_w_prevalence==0))+
  tm_dots(col='white',
          style='cont',
          shape = 21,
          size = .3)+
tm_shape(Regression_df_spat %>%
           filter(Target_stage == "Adult",
                  name == "ha") %>%
           group_by(Site) %>%
           summarize(mean_prevalence = mean(Prevalence),
                     mean_w_prevalence = mean(Weighted_Prevalence)) %>%
           arrange(mean_w_prevalence) %>%
           filter(mean_w_prevalence!=0))+
  tm_dots(col='mean_w_prevalence',
          breaks=c(0,5,10,15,20,25),
          shape = 21,
          size = .3,
          palette = tmaptools::get_brewer_pal('Reds',
                                              contrast = c(.3,1)),
          legend.show = F,
          title = "Mean prevalence")+
  tm_add_legend(type = 'fill',
                title = "Mean prevalence",
                labels = c("0",
                           ">0 to 5",
                           "5 to 10",
                           "10 to 15",
                           "15 to 20",
                           "20 to 25"),
                col = c('white',get_brewer_pal("Reds",n=5,
                                               contrast = c(.3,1))))+
  tm_layout(main.title="Adults, Ap-ha")

m2 = tm_shape(NYS)+
  tm_borders()+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Nymph",
                    name == "ha") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence) %>%
             filter(mean_prevalence==0))+
  tm_dots(col='white',
          style='cont',
          shape = 21,
          size = .3)+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Nymph",
                    name == "ha") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence)  %>%
             filter(mean_w_prevalence!=0))+
  tm_dots(col='mean_w_prevalence',
          breaks = c(0,2,4,6,8,9),
          title = "Ap-ha Prevalence",
          shape = 21,
          size = .3,
          palette = tmaptools::get_brewer_pal("Reds",
                                              contrast = c(.3,1)),
          legend.show = F)+
  tm_add_legend(type = 'fill',
                title = "Mean prevalence",
                labels = c("0",
                           ">0 to 2",
                           "2 to 4",
                           "4 to 6",
                           "6 to 8",
                           "8 to 9"),
                col = c('white',get_brewer_pal("Reds",n=5,
                                               contrast = c(.3,1))))+
  tm_layout(main.title="Nymphs, Ap-ha",
            legend.show=T)

m3 = tm_shape(NYS)+
  tm_borders()+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Adult",
                    name == "v1") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence)  %>%
             filter(mean_w_prevalence==0))+
  tm_dots(col='white',
          shape = 21,
          size = .3)+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Adult",
                    name == "v1") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence)  %>%
             filter(mean_w_prevalence!=0))+
  tm_dots(col='mean_w_prevalence',
          breaks = c(0,2,4,6,8,19),
          shape = 21,
          size = .3,
          palette = tmaptools::get_brewer_pal('Blues',
                                              contrast = c(.3,1)),
          legend.show = F,
          title = "Mean prevalence")+
  tm_add_legend(type = 'fill',
                title = "Mean prevalence",
                labels = c("0",
                           ">0 to 2",
                           "2 to 4",
                           "4 to 6",
                           "6 to 8",
                           "8 to 19"),
                col = c('white',get_brewer_pal("Blues",n=5,
                                               contrast = c(.3,1))))+
  tm_layout(main.title="Adults, Ap-v1");m3

m4 = tm_shape(NYS)+
  tm_borders()+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Nymph",
                    name == "v1") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence)  %>%
             filter(mean_w_prevalence==0))+
  tm_dots(col='white',
          shape = 21,
          size = .3)+
  tm_shape(Regression_df_spat %>%
             filter(Target_stage == "Nymph",
                    name == "v1") %>%
             group_by(Site) %>%
             summarize(mean_prevalence = mean(Prevalence),
                       mean_w_prevalence = mean(Weighted_Prevalence)) %>%
             arrange(mean_w_prevalence)  %>%
             filter(mean_w_prevalence!=0))+
  tm_dots(col='mean_w_prevalence',
          breaks = c(0,2,4,6,8,14),
         title = "Ap-v1 Prevalence",
          shape = 21,
          size = .3,
          palette=tmaptools::get_brewer_pal('Blues',
                                            contrast = c(.3,1)),
         legend.show=F)+
  tm_add_legend(type = 'fill',
                title = "Mean prevalence",
                labels = c("0",
                           ">0 to 2",
                           "2 to 4",
                           "4 to 6",
                           "6 to 8",
                           "8 to 14"),
                col = c('white',get_brewer_pal("Blues",n=5,
                                               contrast = c(.3,1))))+
  tm_scale_bar(position=c('right','bottom'))+
  tm_layout(main.title="Nymphs, Ap-v1")

m5 = tmap_arrange(m1,m2,m3,m4,ncol=2)

tmap_save(m5,
          filename = paste0(getwd(),'/Figures/Figure_3.jpeg'),
          dpi = 300,
          width = 8.25,
          height = 8.25)
