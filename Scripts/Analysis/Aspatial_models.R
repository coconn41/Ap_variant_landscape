#####
# Load libraries
#####
library(readr)
library(tidyverse)
library(mgcv) 
library(spdep)
library(glmmTMB)
library(performance)

######
# A-spatial adult ha model:
#####
adult_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_df)
adult_ha_aspatial2 = glmmTMB(formula = ha ~ metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
adult_ha_aspatial3 = glmmTMB(formula = ha ~ Landscape_metric + 
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            offset = log(tot_tested),
                            family = poisson,
                            data = Adult_df)
AIC1 = AIC(adult_ha_aspatial)
AIC2 = AIC(adult_ha_aspatial2)
AIC3 = AIC(adult_ha_aspatial3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
adult_ha_best = adult_ha_aspatial3

summary(adult_ha_best)
######
# A-spatial adult v1 model:
#####

adult_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Adult_df)
adult_v1_aspatial2 = glmmTMB(formula = v1 ~ metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
adult_v1_aspatial3 = glmmTMB(formula = v1 ~ Landscape_metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Adult_df)
AIC1 = AIC(adult_v1_aspatial)
AIC2 = AIC(adult_v1_aspatial2)
AIC3 = AIC(adult_v1_aspatial3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
adult_v1_best = adult_v1_aspatial2 

summary(adult_v1_best)

######
# A-spatial nymph ha model:
#####

nymph_ha_aspatial = glmmTMB(formula = ha ~ metric + 
                              Landscape_metric +
                              (1|Site) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_df)
nymph_ha_aspatial2 = glmmTMB(formula = ha ~ metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
nymph_ha_aspatial3 = glmmTMB(formula = ha ~ Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
AIC1 = AIC(nymph_ha_aspatial)
AIC2 = AIC(nymph_ha_aspatial2)
AIC3 = AIC(nymph_ha_aspatial3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
nymph_ha_best = nymph_ha_aspatial3

summary(nymph_ha_best)
######
# A-spatial nymph v1 model:
#####

nymph_v1_aspatial = glmmTMB(formula = v1 ~ metric + 
                              Landscape_metric +
                              (1|UNIT) + (1|UNIT:Site),
                            zi = ~1,
                            family = poisson,
                            offset = log(tot_tested),
                            data = Nymph_df)
nymph_v1_aspatial2 = glmmTMB(formula = v1 ~ metric + 
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
nymph_v1_aspatial3 = glmmTMB(formula = v1 ~ Landscape_metric +
                               (1|UNIT) + (1|UNIT:Site),
                             zi = ~1,
                             offset = log(tot_tested),
                             family = poisson,
                             data = Nymph_df)
AIC1 = AIC(nymph_v1_aspatial)
AIC2 = AIC(nymph_v1_aspatial2)
AIC3 = AIC(nymph_v1_aspatial3)
diff1 = AIC1-AIC2
diff2 = AIC1-AIC3
diff3 = AIC2-AIC3
nymph_v1_best = nymph_v1_aspatial

summary(nymph_v1_best)
#####
# A spatial results
#####
summary(adult_ha_best)$coefficients$cond
summary(adult_v1_best)$coefficients$cond
summary(nymph_ha_best)$coefficients$cond
summary(nymph_v1_best)$coefficients$cond


aspatial_results = data.frame(Lifestage = c(rep("Adult",4),
                                            rep("Nymph",6)),
                              Genotype = c("ha","ha",
                                           "v1","v1",
                                           "ha","ha",
                                           "v1","v1",
                                           "v1","v1"),
                              Coeff_n = c(2,2,
                                          2,2,
                                          2,2,
                                          4,4,4,4),
                              Variable = c("Landscape metric","Landscape metric",
                                           "metric","metric",
                                           "Landscape metric","Landscape metric",
                                           "metric","metric",
                                           "Landscape metric","Landscape metric"),
                              Parameter = rep(c("Estimate","P-value"),5),
                              Value = c(summary(adult_ha_best)$coefficients$cond[2],
                                        summary(adult_ha_best)$coefficients$cond[8],
                                        summary(adult_v1_best)$coefficients$cond[2],
                                        summary(adult_v1_best)$coefficients$cond[8],
                                        summary(nymph_ha_best)$coefficients$cond[2],
                                        summary(nymph_ha_best)$coefficients$cond[8],
                                        summary(nymph_v1_best)$coefficients$cond[2],
                                        summary(nymph_v1_best)$coefficients$cond[11],
                                        summary(nymph_v1_best)$coefficients$cond[3],
                                        summary(nymph_v1_best)$coefficients$cond[12])) %>%
  pivot_wider(names_from = Parameter,
              values_from = Value) %>%
  mutate(Significant = ifelse(`P-value`<0.05,1,0),
         IRR = exp(Estimate))

                        
                        
                        
                        