
#########################################
# Reduce form of random judge assignment
#########################################

rm(list=ls())
library(ggplot2)
library(xtable)
library(foreign)
library(MASS)
library(dummies)
library(stargazer)
library(ggthemes)
library(sandwich)

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


###########################################
### First stage regressions ###
###########################################
d$calendar = as.factor(d$calendar)

###########################################
# First stage
###########################################
setwd(dirDrafts)

### Incarceration

d$loo.judge.res <- lm(as.formula(paste("loo.judge~",paste(cov.names,collapse="+")))
                     ,data=d)$residual
d$incarcerate.res <- lm(as.formula(paste("incarcerate~",paste(cov.names,collapse="+")))
                      ,data=d)$residual
first.stage.reg = lm(as.formula(paste("incarcerate~loo.judge+",paste(cov.names,collapse="+")))
                     ,data=d)

fig.visual.first.stage <- ggplot(d,aes(y=incarcerate.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized incarceration indicator \n"
       ,x="\n Residualized LOO mean incarceration by judge",
       title="Incarceration \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$incarcerate.res),1.3)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme( legend.position = "bottom" )

pdf("visual_first_stage_incar.pdf")
fig.visual.first.stage+annotate("text", x = 0, y = 1.2,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()


### Conviction

d$loo.judge.res <- lm(as.formula(paste("loo.judge.conv~",paste(cov.names,collapse="+")))
                      ,data=d)$residual
d$conviction.res <- lm(as.formula(paste("conviction~",paste(cov.names,collapse="+")))
                        ,data=d)$residual
first.stage.reg = lm(as.formula(paste("conviction.res~loo.judge.conv+",paste(cov.names,collapse="+")))
                     ,data=d)


fig.visual.first.stage <- ggplot(d,aes(y=conviction.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized conviction indicator \n"
       ,x="\n Residualized LOO mean conviction by judge",
       title="Conviction \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$conviction.res),1.3)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme( legend.position = "bottom" )


pdf("visual_first_stage_conv.pdf")
fig.visual.first.stage+annotate("text", x = 0, y = 1.2,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()



### Probation length to serve 

d$loo.judge.res <- lm(as.formula(paste("loo.judge.prob~",paste(cov.names,collapse="+")))
                      ,data=d)$residual
d$probat.res <- lm(as.formula(paste("probat~",paste(cov.names,collapse="+")))
                    ,data=d)$residual
first.stage.reg = lm(as.formula(paste("probat~loo.judge.prob+",paste(cov.names,collapse="+")))
                     ,data=d)


fig.visual.first.stage <- ggplot(d,aes(y=probat.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized probation length to serve \n"
       ,x="\n Residualized LOO mean probation length to serve by judge",
       title="Probation length to serve \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$probat.res),26)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme( legend.position = "bottom" )


pdf("visual_first_stage_probat.pdf")
fig.visual.first.stage+annotate("text", x = 0, y = 26,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()


### Sentence length to serve 

d$loo.judge.res <- lm(as.formula(paste("loo.judge.serve~",paste(cov.names,collapse="+")))
                      ,data=d)$residual
d$toserve.res <- lm(as.formula(paste("toserve~",paste(cov.names,collapse="+")))
                       ,data=d)$residual
first.stage.reg = lm(as.formula(paste("toserve~loo.judge.serve+",paste(cov.names,collapse="+")))
                     ,data=d)


fig.visual.first.stage <- ggplot(d,aes(y=toserve.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized sentence length to serve \n"
       ,x="\n Residualized LOO mean sentence length to serve by judge",
       title="Sentence length to serve \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$toserve.res),26)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  theme( legend.position = "bottom" )


pdf("visual_first_stage_toserve.pdf")
fig.visual.first.stage+annotate("text", x = 2, y = 26,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()


### Sentence length to serve (Conditional on being incarcerated)

d0=d
d = d[d$incarcerate==1,]

d$loo.judge.res <- lm(as.formula(paste("loo.judge.serve~",paste(cov.names,collapse="+")))
                      ,data=d[d$incarcerate==1,])$residual
d$toserve.res <- lm(as.formula(paste("toserve~",paste(cov.names,collapse="+")))
                    ,data=d[d$incarcerate==1,])$residual
first.stage.reg = lm(as.formula(paste("toserve~loo.judge.serve+",paste(cov.names,collapse="+")))
                     ,data=d[d$incarcerate==1,])


fig.visual.first.stage <- ggplot(d[d$incarcerate==1,],aes(y=toserve.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized sentence length to serve \n"
       ,x="\n Residualized LOO mean sentence length to serve by judge",
       title="Sentence length to serve \n among incarcerated defendants \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$toserve.res),26)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))


pdf("visual_first_stage_toserve2.pdf")
fig.visual.first.stage+annotate("text", x = -1, y = 26,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()


### Probation length conditional on incarceration

d$loo.judge.res <- lm(as.formula(paste("loo.judge.prob~",paste(cov.names,collapse="+")))
                      ,data=d[d$incarcerate==1,])$residual
d$probat.res <- lm(as.formula(paste("probat~",paste(cov.names,collapse="+")))
                    ,data=d[d$incarcerate==1,])$residual
first.stage.reg = lm(as.formula(paste("probat~loo.judge.prob+",paste(cov.names,collapse="+")))
                     ,data=d[d$incarcerate==1,])

fig.visual.first.stage <- ggplot(d[d$incarcerate==1,],aes(y=probat.res,x=loo.judge.res))+
  geom_point(alpha=0.7,col="blue4")+
  labs(y="Residualized probation length to serve \n"
       ,x="\n Residualized LOO mean probation length to serve by judge",
       title="Probation length to serve \n among incarcerated defendants \n"
  )+
  theme_bw()+
  geom_smooth(method="lm",col="red4",lwd=1)+
  ylim(min(d$probat.res),26)+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_economist()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(panel.border = element_blank(), 
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))

pdf("visual_first_stage_probat2.pdf")
fig.visual.first.stage+annotate("text", x = -1, y = 26,size=5, 
                                label = paste("coef = ",round(coef(first.stage.reg)[2],dig=3),"\n",
                                              "robust SE = ",round(sqrt(vcovHC(first.stage.reg)[2,2]),dig=3),
                                              sep=""))
dev.off()


d=d0
rm(d0)


