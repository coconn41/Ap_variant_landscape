#####
# Load libraries
#####
library(tidyverse)

#####
# Load data
#####
summary_dataset = read.csv(file = paste0(getwd(),'/Data/Regression_df/Regression_df_w_coinf.csv'))

#####
# Pivot longer:
#####
Long_summary = summary_dataset %>%
  mutate(ha = ha+coinf,
         v1 = v1+coinf) %>%
  pivot_longer(c(ha,v1)) %>%
  mutate(Weighted_Prevalence = ifelse(tot_tested>=50,(value/tot_tested)*100,
                                      (value/50)*(value/tot_tested)),
         Prevalence = (value/tot_tested)*100)

#####
# Generate time-series:
#####

Fig2 = ggplot(data = Long_summary,
       aes(x = Year,
           y = Weighted_Prevalence))+
  geom_point(alpha = 0.05,
             size = 5)+
  geom_smooth(data = Long_summary %>%
                mutate(name = ifelse(name=="ha","Ap-ha","Ap-v1")),
              aes(x = Year,
                  y = Weighted_Prevalence,
                  color = name),
              inherit.aes = F,
              method='gam',
              lwd = 3,
              se = F)+
  scale_color_manual(name = "Genotype",
                     values = c("#E73428","#3B8BC2"))+
  facet_wrap(.~Target_stage)+
  ylab("Weighted Prevalence") + 
  theme_bw()+
  theme(text = element_text(size=20));Fig2

ggsave(Fig2,
       filename = paste0(getwd(),'/Figures/Figure_2.jpeg'),
       dpi = 300,
       width = 12)


