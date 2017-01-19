###########################################################################################
# Reduce form table of multiple instruments using LOO and binary IV (not judge dummies)
###########################################################################################

rm(list=ls())
set.seed(12345)
library(lfe)
library(ggplot2)
library(xtable)
library(foreign)
library(MASS)
library(dummies)
library(stargazer)
library(ggthemes)
library(AER)
require(nnet)
library(multiwayvcov)

# loading the data:
d=read.dta("~/Dropbox/crime/recidivism/data/GreenWinik_Criminology_2010_LimitedAnonymizedDataset.2.dta")
dirCode = "~/Dropbox/crime/recidivism/Rcode_final"
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

######################################################################################
# Balance test using Multinomial logit with and without permutation based inference
######################################################################################

logit.judge <- multinom(paste("calendar~",paste(cov.names,collapse="+")), data = d,maxit=1000)
logit0 <- multinom(paste("calendar~1"), data = d)

lrtest(logit.judge,logit0)

### Permutation inference P-value using multinomial logit
lr.obs <- lrtest(logit.judge,logit0)$Chisq[2] 

Num.perm=1000
lr.null <- rep(NA,Num.perm)

for (i in c(1:Num.perm)){
  if(i%%20==0){cat("Iteration: ",i ," \n ")}
  d$calendar0 <- sample(d$calendar,dim(d)[1],replace=FALSE)
  logit.judge.null <- multinom(paste("calendar0~",paste(cov.names,collapse="+")), data = d,maxit=1000)
  logit0.null <- multinom(paste("calendar0~1"), data = d)
  
  lr.null[i] <- lrtest(logit.judge.null,logit0.null)$Chisq[2]
}

cat("Permuation based P-value: ",(sum(lr.obs<lr.null)+1)/(Num.perm+1))

################################
### Figure illustration:
################################

teststat.obs <- lr.obs 
teststat.null <- lr.null

#range = c(min(teststat.obs,teststat.null),max(teststat.obs,teststat.null))
range = c(100,220)
fig<-ggplot(data.frame(statistic.null=teststat.null),aes(x=statistic.null))+
  xlim(range)+
geom_histogram(
  fill="lightblue",
  breaks=seq(range[1],range[2],length=20)
)+
  theme_bw()+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  labs(
    title=NULL,
    x="\n Chi-square test statistic under the null",
    y="\n Frequency"
  )

setwd(dirDrafts)

pdf("balance_permutations.pdf")
fig+geom_vline(xintercept=teststat.obs,lty=1,col="red4",size=1)+
  annotate("text", x = 140, y = 150, label = paste("Observed test statistic \n",
                                                   "Permuation based P-value: ",
                                                   round((sum(lr.obs<lr.null)+1)/(Num.perm+1),dig=3)),
                                                   col="red4")
dev.off()


######################################################################################
# F-tests for all the different instruments
######################################################################################

##############
### LOO
##############

formula1 = as.formula(paste("loo.judge~",paste(cov.names,collapse="+")))
lm.incar <- lm(formula1,data=d)
lm0.incar <- lm("loo.judge~1",data=d)

formula2 = as.formula(paste("loo.judge.serve~",paste(cov.names,collapse="+")))
lm.serve <- lm(formula2,data=d)
lm0.serve <- lm("loo.judge.serve~1",data=d)

formula3 = as.formula(paste("loo.judge.conv~",paste(cov.names,collapse="+")))
lm.conv <- lm(formula3,data=d)
lm0.conv <- lm("loo.judge.conv~1",data=d)

formula4 = as.formula(paste("loo.judge.prob~",paste(cov.names,collapse="+")))
lm.prob <- lm(formula4,data=d)
lm0.prob <- lm("loo.judge.prob~1",data=d)


### partial F-stat
f.stat.loo <- c(anova(lm.incar,lm0.incar)$F[2],
                  anova(lm.serve,lm0.serve)$F[2],
                  anova(lm.conv,lm0.conv)$F[2],
                anova(lm.prob,lm0.prob)$F[2])

f.pv.loo <- c(anova(lm.incar,lm0.incar)$P[2],
              anova(lm.serve,lm0.serve)$P[2],
              anova(lm.conv,lm0.conv)$P[2],
              anova(lm.prob,lm0.prob)$P[2])


##############
### Binary
##############

formula1 = as.formula(paste("z.incar~",paste(cov.names,collapse="+")))
lm.incar <- lm(formula1,data=d)
lm0.incar <- lm("z.incar~1",data=d)

formula2 = as.formula(paste("z.serve~",paste(cov.names,collapse="+")))
lm.serve <- lm(formula2,data=d)
lm0.serve <- lm("z.serve~1",data=d)

formula3 = as.formula(paste("z.conv~",paste(cov.names,collapse="+")))
lm.conv <- lm(formula3,data=d)
lm0.conv <- lm("z.conv~1",data=d)

formula4 = as.formula(paste("z.prob~",paste(cov.names,collapse="+")))
lm.prob <- lm(formula4,data=d)
lm0.prob <- lm("z.prob~1",data=d)


### partial F-stat
f.stat.binary <- c(anova(lm.incar,lm0.incar)$F[2],
                anova(lm.serve,lm0.serve)$F[2],
                anova(lm.conv,lm0.conv)$F[2],
                anova(lm.prob,lm0.prob)$F[2])

f.pv.binary <- c(anova(lm.incar,lm0.incar)$P[2],
              anova(lm.serve,lm0.serve)$P[2],
              anova(lm.conv,lm0.conv)$P[2],
              anova(lm.prob,lm0.prob)$P[2])


###########################################
# Results table
###########################################

tab = matrix(NA,ncol=2,nrow=4)
rownames(tab) = c("Incarceration","Sentence length (months)","Conviction","Probation length (months)")
colnames(tab) = c("Leave-one-out mean (Jackknife)","Binary IV")
tab[,1] = f.pv.loo 
tab[,2] = f.pv.binary 

print(tab)

# Export table to LaTex
stargazer(tab)


#########################################################
# Balance table for incarceration binary instrument
#########################################################

f.stat = function(x1,x0){
  #sd.bias = summary(x1-x0)[4]/sd(c(x1,x0),na.rm=T) # need to correct
  ks = ks.test(x1,x0)$p.value
  ttest = t.test(x1,x0)$p.value
  wilcox = wilcox.test(x1,x0)$p.value
  return(c(mean(x1,na.rm=T),mean(x0,na.rm=T),ttest,wilcox,ks))
}

tab1=round(t(mapply(f.stat,as.list(d[d$z.incar==1,cov.names]),as.list(d[d$z.incar==0,cov.names]))),dig=3)
colnames(tab1) = c("Below median judge (mean)","Above median judge  (mean)","T-test (P-value)","Wilcoxon  (P-value)","KS (P-value)")

tab11 <- rbind(round(tab1[1,],dig=1),
             apply(round(tab1[-c(1,2),c(1,2)]*100,dig=1),2,function(x){return(paste(x,"%",sep=""))})
)
tab1=tab1[-2,]
tab1[,c(1,2)]=tab11
rownames(tab1) <- cov.names[-2]

stargazer(tab1)

setwd(dirCode)
save.image(file="balance_test.RData")
































