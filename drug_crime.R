##############################################################
#
#
#
#
##############################################################

rm(list=ls())
library(ggplot2)
library(xtable)
library(foreign)
library(MASS)
library(dummies)
library(stargazer)
library(ggthemes)
library(AER)
library(multiwayvcov)
library('dplyr')
library('tidyr')

# loading the data:
dirCode = "~/Dropbox/crime/recidivism/Rcode"
dirFigures = "~/Dropbox/crime/recidivism/figures"
dirData = "~/Dropbox/crime/recidivism/data"


setwd(dirData)
d = read.csv("drug_arrests_over_time.csv")
head(d)

# Stayalized facts:
cat("Changing in drug arrests 1980 to 1990:",
    d$total_drug_arrests[d$year==1990 & !is.na(d$year)]/d$total_drug_arrests[d$year==1980 & !is.na(d$year)],"\n")

cat("Changing in drug arrests 1990 to 2000:",
    d$total_drug_arrests[d$year==2000 & !is.na(d$year)]/d$total_drug_arrests[d$year==1990 & !is.na(d$year)],"\n")

d$total_drug_arrests[d$year==1980 & !is.na(d$year)]
d$total_drug_arrests[d$year==1990 & !is.na(d$year)]
d$total_drug_arrests[d$year==2000 & !is.na(d$year)]

d$total_arrests[d$year==1980 & !is.na(d$year)]
d$total_arrests[d$year==1990 & !is.na(d$year)]
d$total_arrests[d$year==2000 & !is.na(d$year)]

# total_arrests
d <- d %>% 
  mutate(total_drug_arrests=100*total_drug_arrests/total_arrests,
         total_violent_crime_arrests=100*total_violent_crime_arrests/total_arrests,
         total_property_crime_arrests=100*total_property_crime_arrests/total_arrests)  %>%
  gather(offence_type,arrests,total_drug_arrests,total_violent_crime_arrests,total_property_crime_arrests)

d$offence_type = as.factor(d$offence_type)
levels(d$offence_type) = c("Drug arrests","Property arrests","Violent arrests")

fig <- ggplot(d,aes(x=year,y=arrests,col=offence_type,shape=offence_type))+
  geom_line()+
  geom_point(size=3)+
  labs(x="Year",
       y="Percent (%) from total arrests")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_colour_grey(start=0, end=0.5)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())

fig<- fig+guides(col=guide_legend(title="Cause of arrest \n"),
           shape=guide_legend(title="Cause of arrest \n"))
  
  

setwd(dirFigures)
ggsave(file="drug_arrests_over_time.pdf")































