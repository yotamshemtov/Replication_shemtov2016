###########################################################################################
# Reduce form table of multiple instruments using LOO and binary IV (not judge dummies)
###########################################################################################

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

# loading the data:
d=read.dta("~/Dropbox/crime/recidivism/data/GreenWinik_Criminology_2010_LimitedAnonymizedDataset.2.dta")
dirCode = "~/Dropbox/crime/recidivism/Rcode"
dirFigures = "~/Dropbox/crime/recidivism/figures"

cov.names = c("age" ,"agesq" ,"female" ,"nonblack" ,"priorarr" ,"priordrugarr" ,
              "priorfelarr" ,"priorfeldrugarr" ,"priorcon" ,"priordrugcon" ,"priorfelcon" ,
              "priorfeldrugcon" ,"pwid" ,"dist" ,"marijuana" ,"cocaine" ,"crack" ,
              "heroin" ,"pcp" ,"otherdrug" ,"nondrug")
head(d)

data.labels <- cbind(attributes(d)$names,attributes(d)$var.label)


######################################################################################
### Define Jackknife IV, leave-one-out estimators

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

###########################################
### Data arrangings ###
###########################################
d$calendar = as.factor(d$calendar)

###########################################
# LOO mean instrument
###########################################

### Future arresrt ###

formula1 = as.formula(paste("laterarr~loo.judge+",paste(cov.names,collapse="+")))
tsls1 <- lm(formula1,data=d)

formula2 = as.formula(paste("laterarr~loo.judge+loo.judge.prob+",paste(cov.names,collapse="+")))
tsls2 <- lm(formula2,data=d)

formula3 = as.formula(paste("laterarr~loo.judge+loo.judge.prob+loo.judge.conv+",paste(cov.names,collapse="+")))
tsls3 <- lm(formula3,data=d)

formula4 = as.formula(paste("laterarr~loo.judge+loo.judge.serve+loo.judge.conv+loo.judge.prob+",paste(cov.names,collapse="+")))
tsls4 <- lm(formula4,data=d)

### Future conviction ###

formula5 = as.formula(paste("latercon~loo.judge+",paste(cov.names,collapse="+")))
tsls5 <- lm(formula5,data=d)

formula6 = as.formula(paste("latercon~loo.judge+loo.judge.prob+",paste(cov.names,collapse="+")))
tsls6 <- lm(formula6,data=d)

formula7 = as.formula(paste("latercon~loo.judge+loo.judge.prob+loo.judge.conv+",paste(cov.names,collapse="+")))
tsls7 <- lm(formula7,data=d)

formula8 = as.formula(paste("latercon~loo.judge+loo.judge.prob+loo.judge.conv+loo.judge.serve+",paste(cov.names,collapse="+")))
tsls8 <- lm(formula8,data=d)


### Results: ###

# Export table:
stargazer(tsls1,tsls2,tsls3,tsls4,tsls5,tsls6,tsls7,tsls8,
          keep=c("loo.judge","loo.judge.prob","loo.judge.conv","loo.judge.serve"),
          se = list(sqrt(diag(cluster.vcov(tsls1, d$clusterid))),sqrt(diag(cluster.vcov(tsls2,d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls3, d$clusterid))),sqrt(diag(cluster.vcov(tsls4, d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls5, d$clusterid))),sqrt(diag(cluster.vcov(tsls6, d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls7, d$clusterid))),sqrt(diag(cluster.vcov(tsls8, d$clusterid))))
          )


###########################################
# Above median indicator instrument
###########################################


### Future arresrt ###

formula1 = as.formula(paste("laterarr~z.incar+",paste(cov.names,collapse="+")))
tsls1 <- lm(formula1,data=d)

formula2 = as.formula(paste("laterarr~z.incar+z.prob+",paste(cov.names,collapse="+")))
tsls2 <- lm(formula2,data=d)

formula3 = as.formula(paste("laterarr~z.incar+z.prob+z.conv+",paste(cov.names,collapse="+")))
tsls3 <- lm(formula3,data=d)

formula4 = as.formula(paste("laterarr~z.incar+z.prob+z.conv+z.serve+",paste(cov.names,collapse="+")))
tsls4 <- lm(formula4,data=d)

### Future conviction ###

formula5 = as.formula(paste("latercon~z.incar+",paste(cov.names,collapse="+")))
tsls5 <- lm(formula5,data=d)

formula6 = as.formula(paste("latercon~z.incar+z.prob+",paste(cov.names,collapse="+")))
tsls6 <- lm(formula6,data=d)

formula7 = as.formula(paste("latercon~z.incar+z.prob+z.conv+",paste(cov.names,collapse="+")))
tsls7 <- lm(formula7,data=d)

formula8 = as.formula(paste("latercon~z.incar+z.prob+z.conv+z.serve+",paste(cov.names,collapse="+")))
tsls8 <- lm(formula8,data=d)

### Results: ###
# Export table:
stargazer(tsls1,tsls2,tsls3,tsls4,tsls5,tsls6,tsls7,tsls8,
          keep=c("z.incar","z.conv","z.prob","z.serve"),
          se = list(sqrt(diag(cluster.vcov(tsls1, d$clusterid))),sqrt(diag(cluster.vcov(tsls2,d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls3, d$clusterid))),sqrt(diag(cluster.vcov(tsls4, d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls5, d$clusterid))),sqrt(diag(cluster.vcov(tsls6, d$clusterid))),
                    sqrt(diag(cluster.vcov(tsls7, d$clusterid))),sqrt(diag(cluster.vcov(tsls8, d$clusterid))))
)












