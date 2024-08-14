library(tidyverse)

#####
# Load data
#####
Regression_df <- read_csv("Data/Regression_df/Regression_df_w_private.csv")[,-1] %>%
  mutate(Site = as.factor(Site),
         UNIT = as.factor(UNIT),
         Year = as.numeric(substring(Date,1,4))) %>%
  filter(is.na(latitude)==F) %>%
  pivot_longer(.,
               c("ha","v1")) %>%
  mutate(Weighted_Prevalence = ifelse(tot_tested>=50,(value/tot_tested)*100,
                                      (value/50)*(value/tot_tested)),
         Prevalence = (value/tot_tested)*100)

Regression_df_double = Regression_df %>%
  dplyr::select(-metric) %>%
  rename(metric = Landscape_metric) %>%
  mutate(metric_type = "ln(Landscape)",
         metric = log(metric)) %>%
  bind_rows(.,Regression_df %>%
              dplyr::select(-Landscape_metric) %>%
              mutate(metric_type = "Patch")) %>%
  mutate(name = ifelse(name=="ha","Ap-ha","Ap-v1"),
         metric_type = factor(metric_type,
                              levels = c("Patch","ln(Landscape)")))
#####
# Make plot
#####

plot = ggplot(data = Regression_df_double %>%
         filter(Weighted_Prevalence!=0),
       aes(x = metric,
           y = Weighted_Prevalence,
           col = name)) +
  geom_point(alpha=.1)+
  geom_smooth(method='lm',se=T)+
  scale_color_manual(name = "Genotype",
                     values = c("#E73428","#3B8BC2"))+
  facet_grid(Lifestage~metric_type,scales = 'free')+
  xlab("Metric")+
  ylab("Weighted Prevalence")+
  coord_cartesian(ylim=c(0,10))+
  theme_bw()+
  theme(text = element_text(size = 15));plot

ggsave(plot = plot,
       filename = paste0(getwd(),'/Figures/Figure_6.jpeg'),
       dpi = 300)



