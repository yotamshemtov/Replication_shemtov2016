
#########################################
# First stage of random judge assignment
#########################################

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
dirDrafts = "~/Dropbox/crime/recidivism/figures"

cov.names = c("age" ,"agesq" ,"female" ,"nonblack" ,"priorarr" ,"priordrugarr" ,
              "priorfelarr" ,"priorfeldrugarr" ,"priorcon" ,"priordrugcon" ,"priorfelcon" ,
              "priorfeldrugcon" ,"pwid" ,"dist" ,"marijuana" ,"cocaine" ,"crack" ,
              "heroin" ,"pcp" ,"otherdrug" ,"nondrug")
head(d)


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

d$calendar = as.factor(d$calendar)


# function for generating tables for each endogenous regressor:
f.table <- function(dd,instrument,regressor,regressor.list=regressor){
  # Outcome laterarr
  formula1  <- as.formula(paste("laterarr~",regressor,"+",paste(cov.names,collapse="+")))
  ols1 <- lm(formula1,data=dd)
  
  formula2 = as.formula(paste("laterarr~",regressor,"|",instrument))
  tsls2 <- ivreg(formula2,data=dd)
  
  formula3 = as.formula(paste("laterarr~",regressor,"+",paste(cov.names,collapse="+"),
                              "|", paste(instrument,"+",cov.names,collapse="+")))
  tsls3 <- ivreg(formula3,data=dd)
  
  # Outcome latercon
  formula4  <- as.formula(paste("latercon~",regressor,"+",paste(cov.names,collapse="+")))
  ols4 <- lm(formula4,data=dd)
  
  formula5 = as.formula(paste("latercon~",regressor,"|",instrument))
  tsls5 <- ivreg(formula5,data=dd)
  
  formula6 = as.formula(paste("latercon~",regressor,"+",paste(cov.names,collapse="+"),
                              "|", paste(instrument,"+",cov.names,collapse="+")))
  tsls6 <- ivreg(formula6,data=dd)
  
  se.coefficients <- lapply(list("ols1","tsls2","tsls3",
                                 "ols4","tsls5","tsls6"),
                            function(x){sqrt(diag(cluster.vcov(get(x),dd$clusterid)))}
  )
  
  # Export table:
  stargazer(ols1,tsls2,tsls3,
            ols4,tsls5,tsls6,
            keep=regressor.list,
            se = se.coefficients)
}

###########################################
# LOO mean - TSLS
###########################################

###########################################
# Incarcerate or not 
###########################################
f.table(d=d,instrument="loo.judge",regressor="incarcerate")

###########################################
# Length of time served 
###########################################
f.table(d=d,instrument="loo.judge.serve",regressor="toserve")

###########################################
# Convicted or not 
###########################################
f.table(d=d,instrument="loo.judge.conv",regressor="conviction")

###########################################
# Probation length 
###########################################
f.table(d=d,instrument="loo.judge.prob",regressor="probat")


###########################################
# Above median indicator - TSLS 
###########################################

###########################################
# Incarcerate or not 
###########################################
f.table(d=d,instrument="z.incar",regressor="incarcerate")

###########################################
# Length of time served 
###########################################
f.table(d=d,instrument="z.serve",regressor="toserve")

###########################################
# Convicted or not 
###########################################
f.table(d=d,instrument="z.conv",regressor="conviction")

###########################################
# Probation length 
###########################################
f.table(d=d,instrument="z.prob",regressor="probat")


######################################################################################
# Above median indicator + LOO mean
######################################################################################

###########################################
# Incarcerate or not 
###########################################
f.table(d=d,instrument="loo.judge+z.incar",regressor="incarcerate")

###########################################
# Length of time served 
###########################################
f.table(d=d,instrument="loo.judge.serve+z.serve",regressor="toserve")

###########################################
# Convicted or not 
###########################################
f.table(d=d,instrument="loo.judge.conv+z.conv",regressor="conviction")

###########################################
# Probation length 
###########################################
f.table(d=d,instrument="loo.judge.prob+z.prob",regressor="probat")

######################################################################################
# LOO mean - Multiple endogeneous regressors 
######################################################################################

# Outcome: laterarr
formula1  <- as.formula(paste("laterarr~","incarcerate","+",paste(cov.names,collapse="+")))
ols1 <- lm(formula1,data=d)

formula2  <- as.formula(paste("laterarr~","incarcerate+probat+conviction+toserve","+",paste(cov.names,collapse="+")))
ols2 <- lm(formula2,data=d)

formula3 = as.formula(paste("laterarr~","incarcerate","+",paste(cov.names,collapse="+"),"|",
                            "loo.judge","+",paste(cov.names,collapse="+")))
tsls3 <- ivreg(formula3,data=d)

formula4 = as.formula(paste("laterarr~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                            "loo.judge+loo.judge.prob+",paste(cov.names,collapse="+")))
tsls4 <- ivreg(formula4,data=d)

formula5 = as.formula(paste("laterarr~","incarcerate+probat+conviction+toserve+",paste(cov.names,collapse="+"),"|",
                            "loo.judge+loo.judge.prob+loo.judge.conv+loo.judge.serve+",paste(cov.names,collapse="+")))
tsls5 <- ivreg(formula5,data=d)



# Outcome: latercon
formula6  <- as.formula(paste("latercon~","incarcerate","+",paste(cov.names,collapse="+")))
ols6 <- lm(formula6,data=d)

formula7  <- as.formula(paste("latercon~","incarcerate+probat+conviction+toserve","+",paste(cov.names,collapse="+")))
ols7 <- lm(formula7,data=d)

formula8 = as.formula(paste("latercon~","incarcerate+",paste(cov.names,collapse="+"),"|","loo.judge+",paste(cov.names,collapse="+")))
tsls8 <- ivreg(formula8,data=d)

formula9 = as.formula(paste("latercon~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                            "loo.judge+loo.judge.prob+",paste(cov.names,collapse="+")))
tsls9 <- ivreg(formula9,data=d)

formula10 = as.formula(paste("latercon~","incarcerate+probat+conviction+toserve+",paste(cov.names,collapse="+"),"|",
                            "loo.judge+loo.judge.prob+loo.judge.conv+loo.judge.serve+",paste(cov.names,collapse="+")))
tsls10 <- ivreg(formula10,data=d)

se.coefficients <- lapply(list("ols1","ols2","tsls3","tsls4","tsls5",
                               "ols6","ols7","tsls8","tsls9","tsls10"),
                          function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))}
)

# Export table:
stargazer(ols1,ols2,tsls3,tsls4,tsls5,
          ols6,ols7,tsls8,tsls9,tsls10,
          keep=c("incarcerate","probat","conviction","toserve"),
          se = se.coefficients)


######################################################################################
# Above median indicator - Multiple endogeneous regressors 
######################################################################################

# Outcome: laterarr
formula1  <- as.formula(paste("laterarr~","incarcerate","+",paste(cov.names,collapse="+")))
ols1 <- lm(formula1,data=d)

formula2  <- as.formula(paste("laterarr~","incarcerate+probat+conviction+toserve","+",paste(cov.names,collapse="+")))
ols2 <- lm(formula2,data=d)

formula3 = as.formula(paste("laterarr~","incarcerate","+",paste(cov.names,collapse="+"),"|",
                            "z.incar","+",paste(cov.names,collapse="+")))
tsls3 <- ivreg(formula3,data=d)

formula4 = as.formula(paste("laterarr~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                            "z.incar+z.prob+",paste(cov.names,collapse="+")))
tsls4 <- ivreg(formula4,data=d)

formula5 = as.formula(paste("laterarr~","incarcerate+probat+conviction+toserve+",paste(cov.names,collapse="+"),"|",
                            "z.incar+z.prob+z.conv+z.serve+",paste(cov.names,collapse="+")))
tsls5 <- ivreg(formula5,data=d)


# Outcome: latercon
formula6  <- as.formula(paste("latercon~","incarcerate","+",paste(cov.names,collapse="+")))
ols6 <- lm(formula6,data=d)

formula7  <- as.formula(paste("latercon~","incarcerate+probat+conviction+toserve","+",paste(cov.names,collapse="+")))
ols7 <- lm(formula7,data=d)

formula8 = as.formula(paste("latercon~","incarcerate+",paste(cov.names,collapse="+"),"|","z.incar+",paste(cov.names,collapse="+")))
tsls8 <- ivreg(formula8,data=d)


formula9 = as.formula(paste("latercon~","incarcerate+probat+",paste(cov.names,collapse="+"),"|",
                                                        "z.incar+z.prob+",paste(cov.names,collapse="+")))
tsls9 <- ivreg(formula9,data=d)


formula10 = as.formula(paste("latercon~","incarcerate+probat+conviction+toserve+",paste(cov.names,collapse="+"),"|",
                             "z.incar+z.prob+z.conv+z.serve+",paste(cov.names,collapse="+")))
tsls10 <- ivreg(formula10,data=d)


se.coefficients <- lapply(list("ols1","ols2","tsls3","tsls4","tsls5",
                               "ols6","ols7","tsls8","tsls9","tsls10"),
                          function(x){sqrt(diag(cluster.vcov(get(x),d$clusterid)))}
)

# Export table:
stargazer(ols1,ols2,tsls3,tsls4,tsls5,
          ols6,ols7,tsls8,tsls9,tsls10,
          keep=c("incarcerate","probat","conviction","toserve"),
          se = se.coefficients)















