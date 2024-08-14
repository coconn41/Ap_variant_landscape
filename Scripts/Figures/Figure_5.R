library(readxl)
library(ggplot2)
ANA_inc <- read_excel(paste0(getwd(),"/Data/NYS_ANA_Inc/ANA_inc.xlsx"))

p1=ggplot(data = ANA_inc,
       aes(x = Year,
           y = Incidence))+
  geom_line()+
  ylab("Anaplasmosis incidence per 100,000")+
  scale_x_continuous(breaks=c(min(ANA_inc$Year):max(ANA_inc$Year)))+
  theme_bw()

ggsave(plot = p1,
       filename = paste0(getwd(),"/Figures/Figure_5.jpeg"),
       dpi = 300)
