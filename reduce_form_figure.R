
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

library(reshape)

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
# Reduce Form  - above median calendar vs. below median calendar
####################################################################################

#start.month=3
#months.vec <- seq(start.month,12*4,by=1)

start.month=4
months.vec <- seq(start.month,209,by=1)
reduce.treat <- reduce.control <- rep(NA,length(months.vec))

for (j in c(start.month:max(months.vec))){
  
  ### Calculating a recidivism measure for a re-arrest in the first "j" months after the disposition date:
  date.disp.plus.j.months <- d[,"dispdate"]
  #month(date.disp.plus.j.months) = month(date.disp.plus.j.months)+j
  week(date.disp.plus.j.months) = week(date.disp.plus.j.months)+j
  
  d$laterarr.dynamic <-  (difftime(d$laterarr.date,date.disp.plus.j.months, units = "weeks")<=0)*1 
  
  i = j-start.month+1
  reduce.treat[i] = mean(d$laterarr.dynamic[d$z.incar==1])
  reduce.control[i] = mean(d$laterarr.dynamic[d$z.incar==0])
  
}



################################################################################################
# Figure:
################################################################################################

dp = data.frame(months=months.vec, reduce.treat, reduce.control)
dp <- melt(dp,id.vars="months")
colnames(dp) <- c("months","group","rearrest")
head(dp)

levels(dp$group) <- c("Above median calendar","Below median calendar")

p <- ggplot(dp,aes(x=months,y=rearrest,col=group,shape=group))+
  geom_line(lwd=0.8)+
  labs(x="\n Period from disposition date (weeks)",y="Ave. re-arrests  \n")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_colour_economist()+
  theme(legend.position="bottom")+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))

p <- p + 
  annotate("segment",
           x = 4, y = 0.2+0.01, xend = 52, yend = 0.2+0.01,
           arrow=arrow(length=unit(.2, "cm")))+
annotate("segment",
         x = 52, y = 0.2+0.01, xend=4, yend=0.2+0.01,
         arrow=arrow(length=unit(.2, "cm")))+
  annotate("text",x = 30, y=0.2+0.04,label="Incapacitation effect \n dominates",col="red4",size=3)
  
p <- p+geom_vline(xintercept=53,col="red4",lty=2)+
  annotate("text",x=70,y=0.6,label="53 weeks",col="red4",size=4)



setwd(dirFigures)
pdf("above_vs_below_median_reduce_form.pdf")
print(p)
dev.off()



####################################################################################
# Incarcerated vs. not-incarcerated
####################################################################################

### Calculating the DFL weights:

ps.model = glm(paste("incarcerate~",paste(cov.names,collapse="+")),data=d,family=binomial(link="logit"))
summary(ps.model)

ps <- ps.model$fit
pd <- mean(d$incarcerate)

# DFL re-weighting
weight0 = ((1-pd)/pd)*(ps/(1-ps))
weight=rep(1,dim(d)[1])
weight=replace(weight,d$incarcerate==0,weight0[d$incarcerate==0])

# Ppropensity score re-weighting
#weight = rep(NA,dim(d)[1])
#weight[d$incarcerate==1]=1/ps[d$incarcerate==1]
#weight[d$incarcerate==0]=1/(1-ps[d$incarcerate==0])

# Validating the DFL re-weighting generates balance on observables between incarcerated and non-incarcerated defendants
balance.model.lm <- lm(paste("incarcerate~",paste(cov.names,collapse="+")),data=d,weights=weight)
summary(balance.model.lm)

stargazer(balance.model.lm)


start.month=4
months.vec <- seq(start.month,209,by=1)
incarcerated <- not.incarcerated <- not.incarcerated.dfl <- rep(NA,length(months.vec))

for (j in c(start.month:max(months.vec))){
  
  ### Calculating a recidivism measure for a re-arrest in the first "j" months after the disposition date:
  date.disp.plus.j.months <- d[,"dispdate"]
  #month(date.disp.plus.j.months) = month(date.disp.plus.j.months)+j
  week(date.disp.plus.j.months) = week(date.disp.plus.j.months)+j
  
  d$laterarr.dynamic <-  (difftime(d$laterarr.date,date.disp.plus.j.months, units = "weeks")<=0)*1 
  
  i = j-start.month+1
  
  incarcerated[i] = mean(d$laterarr.dynamic[d$incarcerate==1])
  not.incarcerated[i] = mean(d$laterarr.dynamic[d$incarcerate==0])
  not.incarcerated.dfl[i] = weighted.mean(d$laterarr.dynamic[d$incarcerate==0],w=weight[d$incarcerate==0])
  
}


dp = data.frame(months=months.vec, incarcerated, not.incarcerated,not.incarcerated.dfl)
dp <- melt(dp,id.vars="months")
dp <- rename(dp,c(variable="group",value="rearrest"))
head(dp)

levels(dp$group) <- c("Incarcerated","Not-Incarcerated","Not-Incarcerated DFL re-weighted ")
#dp$rearrest <- round(100*dp$rearrest,dig=1)

p <- ggplot(dp,aes(x=months,y=rearrest,col=group,shape=group,lty=group))+
  geom_line(lwd=0.8)+
  labs(x="\n Period from disposition date (weeks)",y="Ave. re-arrests  \n")+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_colour_manual(values=c("red4","black","black"))+
  scale_linetype_manual(values=c(1,1,2))+
  theme(legend.position="bottom")+
  guides(col=guide_legend(title=""),
         linetype=guide_legend(title=""))+
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))


setwd(dirFigures)

pdf("incarcerated_vs_non_incarcerated.pdf")
print(p)
dev.off()




