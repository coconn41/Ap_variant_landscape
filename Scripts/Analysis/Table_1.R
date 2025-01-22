#####
# Load data
#####
summary_dataset = read.csv(file = paste0(getwd(),'/Data/Regression_df/Regression_df_w_coinf.csv'))

#####
# Calculate summary stats
#####

# Mean and sd of ticks collected:
summary_dataset %>% group_by(Lifestage) %>% summarize(mn = mean(tot_collected,na.rm=T))

# Total unique sites:
length(unique(summary_dataset$Site))

# Total visits:
nrow(summary_dataset)

# Number collected:
sum(summary_dataset$tot_collected)

# Number tested:
sum(summary_dataset$tot_tested)

# Number positive for Anaplasma (overall):
sum(summary_dataset$coinf)+
  sum(summary_dataset$ha)+
  sum(summary_dataset$v1)+
  sum(summary_dataset$und)

# Number positive for Ap-ha:
sum(summary_dataset$coinf)+
  sum(summary_dataset$ha)

# Number positive for Ap-v1:
sum(summary_dataset$coinf)+
  sum(summary_dataset$v1)

# Number undetermined:
sum(summary_dataset$und)

# Full breakdown by Life Stage:

summary_dataset %>%
  group_by(Target_stage) %>%
  summarize(tot_visits = n(),
            total_collected = sum(tot_collected),
            tot_tested = sum(tot_tested),
            test_perc = (tot_tested/total_collected)*100,
            tot_ANA = sum(ha)+sum(v1)+sum(coinf)+sum(und),
            ANA_perc = (tot_ANA/tot_tested)*100,
            tot_ha = sum(ha)+sum(coinf),
            ha_perc = (tot_ha/tot_tested)*100,
            tot_v1 = sum(v1)+sum(coinf),
            v1_perc = (tot_v1/tot_tested)*100,
            tot_coinf = sum(coinf),
            coinf_perc = (tot_coinf/tot_tested)*100,
            tot_und = sum(und),
            und_perc = (tot_und/tot_tested)*100) %>%
  View()



