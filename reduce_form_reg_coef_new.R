##################################################################################
# Description: 
# Date: January 2017
#
#
##################################################################################







##################################################################################
# Reduce form table of coefficients using LOO and binary IV (not judge dummies)
##################################################################################

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
library(lfe)

# loading the data:
d=read.dta("~/Dropbox/crime/recidivism/data/GreenWinik_Criminology_2010_LimitedAnonymizedDataset.2.dta")
dirCode = "~/Dropbox/crime/recidivism/Rcode"
dirDrafts = "~/Dropbox/crime/recidivism/figures"

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
# Judge dummies - reduce form regresison 
###########################################

formula0 = as.formula(paste("laterarr~","as.factor(calendar)","|0|0|clusterid"))
lm0 <- felm(formula0, data = d)
summary(lm0)

# The specification with covariates was done in STATA

###########################################
# LOO mean - reduce form regresison 
###########################################

###########################################
# Incarcerate or not 
###########################################

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.incarcerate <- dd$loo.judge

### First-stage
formula1 = as.formula(paste("incarcerate~z.incarcerate"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("incarcerate~z.incarcerate+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.incarcerate")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.incarcerate+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.incarcerate")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.incarcerate+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.incarcerate"),
          se = se.coefficient)

###########################################
# Length of time to served 
###########################################

### time to serve as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.toserve <- dd$loo.judge.serve

### First-stage
formula1 = as.formula(paste("toserve~z.toserve"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("toserve~z.toserve+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.toserve")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.toserve+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.toserve")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.toserve+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.toserve"),
          se = se.coefficient)

###########################################
# Convicted or not 
###########################################

### conviction as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.conviction <- dd$loo.judge.conv

### First-stage
formula1 = as.formula(paste("conviction~z.conviction"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("conviction~z.conviction+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.conviction")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.conviction+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.conviction")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.conviction+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.conviction"),
          se = se.coefficient)

###########################################
# Probation length 
###########################################

### conviction as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.probation <- dd$loo.judge.prob

### First-stage
formula1 = as.formula(paste("probat~z.probation"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("probat~z.probation+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.probation")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.probation+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.probation")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.probation+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.probation"),
          se = se.coefficient)


######################################################################################
# Above median indicator - reduce form regresison 
######################################################################################


###########################################
# Incarcerate or not 
###########################################

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.incarcerate <- dd$z.incar

### First-stage
formula1 = as.formula(paste("incarcerate~z.incarcerate"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("incarcerate~z.incarcerate+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.incarcerate")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.incarcerate+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.incarcerate")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.incarcerate+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.incarcerate"),
          se = se.coefficient)

### Table of effects interms of Standard Deviations:

# Effects in terms of SD deviations
coef.incar = vapply(paste("lm",c(1:6),sep=""),
                    function(x){coef(get(x))[2]},FUN.VALUE=numeric(1))
effect.incar <- coef.incar*sd(d$loo.judge)

# The SE of the estimated effects:
se.coef.incar <- vapply(paste("lm",c(1:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))[2]},FUN.VALUE=numeric(1))
se.effect.incar <- se.coef.incar*sd(d$loo.judge)


###########################################
# Length of time to served 
###########################################

### time to serve as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.toserve <- dd$z.serve

### First-stage
formula1 = as.formula(paste("toserve~z.toserve"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("toserve~z.toserve+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.toserve")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.toserve+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.toserve")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.toserve+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.toserve"),
          se = se.coefficient)

###########################################
# Convicted or not 
###########################################

### conviction as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.conviction <- dd$z.conv

### First-stage
formula1 = as.formula(paste("conviction~z.conviction"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("conviction~z.conviction+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.conviction")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.conviction+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.conviction")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.conviction+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.conviction"),
          se = se.coefficient)

### Table of effects interms of Standard Deviations:
# To-ADD


###########################################
# Probation length 
###########################################

### conviction as the treatment:

# This is for the "stargazer export to be one line..." (completely technical)
dd=d

# loo
dd$z.probation <- dd$z.prob

### First-stage
formula1 = as.formula(paste("probat~z.probation"))
lm1 <- lm(formula1,data=dd)

formula2 = as.formula(paste("probat~z.probation+",paste(cov.names,collapse="+")))
lm2 <- lm(formula2,data=dd)

### Reduce-form laterarr
formula3 = as.formula("laterarr~z.probation")
lm3 <- lm(formula3,data=dd)

formula4 = as.formula(paste("laterarr~z.probation+",paste(cov.names,collapse="+")))
lm4 <- lm(formula4,data=dd)

### Reduce-form latercon
formula5 = as.formula("latercon~z.probation")
lm5 <- lm(formula5,data=dd)

formula6 = as.formula(paste("latercon~z.probation+",paste(cov.names,collapse="+")))
lm6 <- lm(formula6,data=dd)

# Export table:
se.coefficient = lapply(paste("lm",c(3:6),sep=""),
                        function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))})
stargazer(lm3,lm4,lm5,lm6,
          keep=c("z.probation"),
          se = se.coefficient)

















