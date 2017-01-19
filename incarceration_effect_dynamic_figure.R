
###########################################################################################################################
# Dynamic effect of incarceration according to different time measures of recidivism
###########################################################################################################################

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

# Date and time packages:
library(lubridate)

# loading the data:
d=read.dta("~/Dropbox/crime/recidivism/data/GreenWinik_Criminology_2010_LimitedAnonymizedDataset.2.dta")
dirCode = "~/Dropbox/crime/recidivism/Rcode"
dirFigures = "~/Dropbox/crime/recidivism/figures"

cov.names = c("age" ,"agesq" ,"female" ,"nonblack" ,"priorarr" ,"priordrugarr" ,
              "priorfelarr" ,"priorfeldrugarr" ,"priorcon" ,"priordrugcon" ,"priorfelcon" ,
              "priorfeldrugcon" ,"pwid" ,"dist" ,"marijuana" ,"cocaine" ,"crack" ,
              "heroin" ,"pcp" ,"otherdrug" ,"nondrug")
head(d)

####################################################################################
### Define Jackknife IV, leave-one-out estimators
####################################################################################

d$loo.judge = rep(NA,dim(d)[1])
for (i in c(1:dim(d)[1])){
  d$loo.judge[i] = mean(d$incarcerate[-i][d$calendar[-i]==d$calendar[i]])
}

d$loo.judge.serve = rep(NA,dim(d)[1])
for (i in c(1:dim(d)[1])){
  d$loo.judge.serve[i] = mean(d$toserve[-i][d$calendar[-i]==d$calendar[i]])
}

d$loo.judge.conv = rep(NA,dim(d)[1])
for (i in c(1:dim(d)[1])){
  d$loo.judge.conv[i] = mean(d$conviction[-i][d$calendar[-i]==d$calendar[i]])
}

d$loo.judge.prob = rep(NA,dim(d)[1])
for (i in c(1:dim(d)[1])){
  d$loo.judge.prob[i] = mean(d$probat[-i][d$calendar[-i]==d$calendar[i]])
}

### Defining binary IV:

d$z.incar = (d$calendar %in% c(1:9)[order(tapply(d$incarcerate,d$calendar,mean))][6:9])*1
d$z.serve = (d$calendar %in% c(1:9)[order(tapply(d$toserve,d$calendar,mean))][6:9])*1
d$z.conv = (d$calendar %in% c(1:9)[order(tapply(d$conviction,d$calendar,mean))][6:9])*1
d$z.prob = (d$calendar %in% c(1:9)[order(tapply(d$probat,d$calendar,mean))][6:9])*1

# The binary variables are the same except for 1 judge:
c(1:9)[order(tapply(d$incarcerate,d$calendar,mean))][6:9]
c(1:9)[order(tapply(d$toserve,d$calendar,mean))][6:9]
c(1:9)[order(tapply(d$conviction,d$calendar,mean))][6:9]
c(1:9)[order(tapply(d$probat,d$calendar,mean))][6:9]

d$calendar = as.factor(d$calendar)

####################################################################################
# Define date variables
####################################################################################

d$laterarr.date <- as.Date(d$laterarrdate,format="%m/%d/%Y")
d$laterarr.date[is.na(d$laterarr.date)] = "3000-01-01"
cat("Check: ",sum(is.na(d$laterarr.date))==0,"\n")

d$dispdate <- as.Date(d$dispdate,format="%m/%d/%Y")
cat("Check: ",sum(is.na(d$dispdate))==0,"\n")


####################################################################################
# Estimation with incarceration and probation as endogenous sentencing outcomes
####################################################################################

start.month=1
months.vec <- seq(start.month,209,by=1)
tsls.coef.dynamic <- rep(NA,length(months.vec))

for (j in c(start.month:max(months.vec))){
  
  ### Calculating a recidivism measure for a re-arrest in the first "j" months after the disposition date:
  date.disp.plus.j.months <- d[,"dispdate"]
  #month(date.disp.plus.j.months) = month(date.disp.plus.j.months)+j
  week(date.disp.plus.j.months) = week(date.disp.plus.j.months)+j
  
  d$laterarr.dynamic <-  (difftime(d$laterarr.date,date.disp.plus.j.months, units = "weeks")<=0)*1 
  
  formula = as.formula(paste("laterarr.dynamic~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                             "z.incar+z.prob+",
                             paste(cov.names,collapse="+")))
  i = j-start.month+1
  tsls.coef.dynamic[i] <- coef(ivreg(formula,data=d))[2]
  
}


# check that results match the recdivism measure with 4 years:
formula = as.formula(paste("laterarr~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                           "z.incar+z.prob+",
                           paste(cov.names,collapse="+")))

cat("Check: ",tsls.coef.dynamic[length(tsls.coef.dynamic)]==coef(ivreg(formula,data=d))[2])


################################################################################################
# Figure:
################################################################################################

dp = data.frame(months=months.vec, tsls = tsls.coef.dynamic)

p <- ggplot(dp,aes(x=months,y=tsls))+
  geom_point(size=0.8)+geom_line()+
  labs(x="\n Period from disposition date (weeks)",y="Incarceration effect - TSLS estimates \n")+
  theme_bw()+
  theme(legend.position="")+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


#geom_vline(xintercept=min(months.vec[tsls.coef.dynamic>=0])-0.5,lty=2,col="red4",lwd=0.7)
p <- p+
  annotate("text",x = 30, y=0.13-0.02,label="Incapacitation effect \n dominates",col="red4",size=3)+
  annotate("segment",
           x = 6, y = 0.1-0.02, xend = min(months.vec[tsls.coef.dynamic>=0 & months.vec>50])-2, yend = 0.1-0.02,
           arrow=arrow(length=unit(.2, "cm")))+
  annotate("segment",
           x = min(months.vec[tsls.coef.dynamic>=0 & months.vec>50])-2, y = 0.1-0.02, xend=6, yend=0.1-0.02,
           arrow=arrow(length=unit(.2, "cm")))

p <- p+geom_hline(yintercept=coef(ivreg(formula,data=d))[2],col="blue4",lty=1,lwd=0.4)+
  annotate("text",x = 60, y=0.4,label="TSLS estimate using 4 years from disposition",col="blue4",size=3)


setwd(dirFigures)
ggsave("incar_effect_dynamic_probat.pdf",plot=p, width = 6, height = 4.5)

### statistics:
# The week in which the incapacitation effect stops dominating:
print(min(months.vec[tsls.coef.dynamic>=0 & months.vec>50]))



####################################################################################
# Estimation with only incarceration as endogenous sentencing outcomes
####################################################################################

start.month=1
months.vec <- seq(start.month,209,by=1)
tsls.coef.dynamic <- rep(NA,length(months.vec))

for (j in c(start.month:max(months.vec))){
  
  ### Calculating a recidivism measure for a re-arrest in the first "j" months after the disposition date:
  date.disp.plus.j.months <- d[,"dispdate"]
  #month(date.disp.plus.j.months) = month(date.disp.plus.j.months)+j
  week(date.disp.plus.j.months) = week(date.disp.plus.j.months)+j
  
  d$laterarr.dynamic <-  (difftime(d$laterarr.date,date.disp.plus.j.months, units = "weeks")<=0)*1 
  
  formula = as.formula(paste("laterarr.dynamic~","incarcerate+",paste(cov.names,collapse="+"),"|",
                             "z.incar+",
                             paste(cov.names,collapse="+")))
  i = j-start.month+1
  tsls.coef.dynamic[i] <- coef(ivreg(formula,data=d))[2]
  
}


# check that results match the recdivism measure with 4 years:
formula = as.formula(paste("laterarr~","incarcerate+",paste(cov.names,collapse="+"),"|",
                           "z.incar+",
                           paste(cov.names,collapse="+")))

cat("Check: ",tsls.coef.dynamic[length(tsls.coef.dynamic)]==coef(ivreg(formula,data=d))[2])


################################################################################################
# Figure:
################################################################################################

dp = data.frame(months=months.vec, tsls = tsls.coef.dynamic)

p <- ggplot(dp,aes(x=months,y=tsls))+
  geom_point(size=0.8)+geom_line()+
  labs(x="\n Period from disposition date (weeks)",y="Incarceration effect - TSLS estimates \n")+
  theme_bw()+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


#geom_vline(xintercept=min(months.vec[tsls.coef.dynamic>=0])-0.5,lty=2,col="red4",lwd=0.7)
p <- p+
  annotate("text",x = 30, y=0.13-0.02,label="Incapacitation effect \n dominates",col="red4",size=3)+
  annotate("segment",
           x = 6, y = 0.1-0.02, xend = min(months.vec[tsls.coef.dynamic>=0 & months.vec>50])-2, yend = 0.1-0.02,
           arrow=arrow(length=unit(.2, "cm")))+
  annotate("segment",
           x = min(months.vec[tsls.coef.dynamic>=0 & months.vec>50])-2, y = 0.1-0.02, xend=6, yend=0.1-0.02,
           arrow=arrow(length=unit(.2, "cm")))

p <- p+geom_hline(yintercept=coef(ivreg(formula,data=d))[2],col="blue4",lty=1,lwd=0.4)+
  annotate("text",x = 60, y=0.4,label="TSLS estimate using 4 years from disposition",col="blue4",size=3)


setwd(dirFigures)
ggsave("incar_effect_dynamic.pdf",plot=p, width = 6, height = 4.5)

### statistics:
# The week in which the incapacitation effect stops dominating:
print(min(months.vec[tsls.coef.dynamic>=0 & months.vec>50]))


























