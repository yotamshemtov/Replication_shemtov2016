##################################################################################
# Characterizing the compliers, alaway-takers and never-takers
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

# loading the data:
d=read.dta("~/Dropbox/Class.test/Results/Judges/GreenWinik_Criminology_2010_LimitedAnonymizedDataset.2.dta")
dirCode = "~/Dropbox/recidivism/Rcode"
dirDrafts = "~/Dropbox/recidivism/figures"

cov.names = c("age" ,"agesq" ,"female" ,"nonblack" ,"priorarr" ,"priordrugarr" ,
              "priorfelarr" ,"priorfeldrugarr" ,"priorcon" ,"priordrugcon" ,"priorfelcon" ,
              "priorfeldrugcon" ,"pwid" ,"dist" ,"marijuana" ,"cocaine" ,"crack" ,
              "heroin" ,"pcp" ,"otherdrug" ,"nondrug")
#head(d)

data.labels <- cbind(attributes(d)$names,attributes(d)$var.labels)
print(d)

### Defining binary IV:

d$z.incar = (d$calendar %in% c(1:9)[order(tapply(d$incarcerate,d$calendar,mean))][6:9])*1

### The characteristics of the compliers: 

# The function below calculates the compliers first moment (mean) of a pre-treatment characteristic "x": 
f.complier.cov <- function(x){
  denominator <- mean(d$incarcerate[which(d$z.incar==1)]) - mean(d$incarcerate[which(d$z.incar==0)])
  numerator <- mean((x*d$incarcerate)[which(d$z.incar==1)]) - mean((x*d$incarcerate)[which(d$z.incar==0)])
  return(numerator/denominator)
}

f.mean.comparison <- function(name){
  x = d[,name]
  
  means <- matrix(NA,ncol=4) 
  colnames(means) = c("All","Incarcerated","Compliers","Not-incarcerated")
  means[,1] <- mean(x)
  means[,2] <- mean(x[which(d$incarcerat==1)])
  means[,3] <- f.complier.cov(x)
  means[,4] <- mean(x[which(d$incarcerat==0)])
  return(means)
}

tab <- t(sapply(cov.names,f.mean.comparison))
colnames(tab) <-  c("All","Incarcerated","Compliers","Not-incarcerated")

tab <- rbind(round(tab[1,],dig=1),
             apply(round(tab[-c(1,2),]*100,dig=1),2,function(x){return(paste(x,"%",sep=""))})
             )
rownames(tab) <- cov.names[-2]

# Exporting to LaTex
stargazer(tab,
          title="The characteristics of the compliers population relative the incarcerated and non-incarcerated populations",
          label="tab: compliers characteristics")


#################################
# Potential outcomes table:
#################################

tab.potential = matrix(NA,ncol=3,nrow=2)
colnames(tab.potential) <- c("Incarcerated","Not-incarcerated","Compliers")
rownames(tab.potential) <- c("mean Y1","MeanY0")
  
tab.potential[2,1]="-"
tab.potential[1,2]="-"

tab.potential[1,1] <- paste(round(100*mean(d$laterarr[d$incarcerate==1]),dig=0),"%",sep="")
tab.potential[2,2] <- paste(round(100*mean(d$laterarr[d$incarcerate==0]),dig=0),"%",sep="")

tab.potential[1,3] <- paste(round(100*f.complier.cov(d$laterarr),dig=0),"%",sep="")

denominator <- mean((1-d$incarcerate)[which(d$z.incar==1)]) - mean((1-d$incarcerate)[which(d$z.incar==0)])
numerator <- mean((d$laterarr*(1-d$incarcerate))[which(d$z.incar==1)]) - mean((d$laterarr*(1-d$incarcerate))[which(d$z.incar==0)])
tab.potential[2,3] <- paste(round(100*(numerator/denominator),dig=0),"%",sep="")

# Exporting to LaTex
stargazer(tab.potential,
          title="The average potential outcomes with/out incarceration",
          label="tab: potential outcomes table")


###################################################################################################
# Characteristics of: (1) compliers, (2) never takers (3) always takers
###################################################################################################

f.complier.always.never <- function(name){
  x = d[,name]
  
  means <- matrix(NA,ncol=4) 
  colnames(means) = c("All","Always-takers","Compliers","Never-takers")
  means[,1] <- mean(x)
  means[,2] <- mean(x[which(d$incarcerat==1 & d$z.incar==0)])
  means[,3] <- f.complier.cov(x)
  means[,4] <- mean(x[which(d$incarcerat==0 & d$z.incar==1)])
  return(means)
}

tab <- t(sapply(cov.names,f.complier.always.never))
colnames(tab) <-  c("All","Always-takers","Compliers","Never-takers")

tab <- rbind(round(tab[1,],dig=1),
             apply(round(tab[-c(1,2),]*100,dig=1),2,function(x){return(paste(x,"%",sep=""))})
)
rownames(tab) <- cov.names[-2]

# Exporting to LaTex
stargazer(tab,
          title="The characteristics of the compliers, always-takers, and never-takers populations",
          label="tab: compliers classes characteristics")






