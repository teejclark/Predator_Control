#### Predator Control Meta-analysis - Journal of Applied Ecology
#### T.J. Clark

rm(list=ls(all=T))

library(tidyr)
library(dplyr)
library(ggplot2)

########################################################################
# 1) FORMAT DATA

pred$Type.of.Response[pred$Type.of.Response=="12 week neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="16 week neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="90 day neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="fawn survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="density"] <- "abundance"
pred$Type.of.Response[pred$Type.of.Response=="growth rate"] <- "abundance"
pred$Type.of.Response <- factor(pred$Type.of.Response)

# create variance, used to weight things
pred$var <- ((pred$SE...Experiment^2)/(pred$Mean...Experiment)^2)+((pred$SE...Control^2)/(pred$Mean...Control)^2)

# use lnCVR from Nakagawa et al. (2015)
pred$CVexp <- pred$SE...Experiment/pred$Mean...Experiment
pred$CVcon <- pred$SE...Control/pred$Mean...Control
pred$lnCVR <- log(pred$CVexp/pred$CVcon)

# 1) Experiment Type
rigor <- pred %>%
  mutate(study = recode(pred$Experiment.Type,
                        "BACI"=2,"SEC"=1,"BA"=0)) %>%
  mutate(natural = recode(pred$Natural,
                          "Yes"=0,"No"=1)) %>%
  mutate(response = recode(pred$Type.of.Response,
                           "abundance"=2,
                           "adult survival"=1,
                           "calf-cow ratios"=0,
                           "calf survival"=0,
                           "recruitment"=0)) %>%
  mutate(rep = recode(pred$Replication,
                         "u"=0,"r"=2)) %>%
  mutate(effectiveness = ifelse(pred$X..Change.in..>=50,
                                1,ifelse(is.na(pred$X..Change.in..),
                                         0,0))) %>%
  mutate(length = ifelse(pred$Prey.Gen.Time<=pred$Temporal.Scale..months.,
                         1,0)) %>%
  mutate(size = ifelse(pred$Pred.Max.Range<=pred$Spatial.Scale..km2.,
                       1,0))
rigor$effectiveness[is.na(rigor$effectiveness)] <- 0

# total rigor score up!
rigor$total <- rigor$study+rigor$natural+rigor$response+rigor$rep+rigor$effectiveness+rigor$length+rigor$size
pred$total <- rigor$total

# biodiversity
pred$biodiv <- pred$Predator.Diversity+pred$Prey.Diversity


###############################################################################
#### 2) BASIC SUMMARY STATS

mean(pred$Log.RR.Effect.Size) # 0.25 -> 1.28

# combine dependent variables
# NOTE: "abundance" = all population indices
table(pred$Type.of.Response)
pred$Type.of.Response[pred$Type.of.Response=="12 week neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="16 week neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="90 day neonate survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="fawn survival"] <- "calf survival"
pred$Type.of.Response[pred$Type.of.Response=="density"] <- "abundance"
pred$Type.of.Response[pred$Type.of.Response=="growth rate"] <- "abundance"
pred$Type.of.Response <- factor(pred$Type.of.Response)


windows();
boxplot(pred$Log.RR.Effect.Size~pred$Type.of.Response)
abline(h=0)
title("Means - all values, unweighted")

# prelim plots
blah <- data.frame(mean=c(-0.029,0.044,0.117,0.303,0.143,0.136,0.558),
                   se=c(1.47,0.151,0.36,0.503,0.287,0.2399,0.283),
                   type=c("abund","AS","CC ratios","CS","dens","lambda","recruit"))
blah$lowerCI <- blah$mean-(blah$se*1.96)
blah$highCI <- blah$mean+(blah$se*1.96)

windows();
ggplot(blah, aes(x=type,y=mean))+
  geom_point(size=4)+
  geom_errorbar(aes(ymax=highCI,ymin=lowerCI),size=2)+
  geom_hline(yintercept=0,linetype="dashed")+
  xlab("Type of Response")

############################################################################
# 3) VARIANCE ANALYSIS

# calculate CVs
removed$CVExperiment <- removed$SE...Experiment/removed$Mean...Experiment
removed$CVControl <- removed$SE...Control/removed$Mean...Control

# calculate lnCVR
removed$lnCVR <- log(removed$CVExperiment/removed$CVControl)

# calculate correlations
corexp <- cor(log(removed$Mean...Experiment),log(removed$SE...Experiment))
corcon <- cor(log(removed$Mean...Control),log(removed$SE...Control))
# 0.913 and 0.918, respectively

# calculate var(lnCVR)
removed$varCVR <- (removed$SE...Experiment)^2-(2*corexp*(sqrt((removed$SE...Experiment/removed$Mean...Experiment)^2))+(removed$SE...Experiment)^2-(2*corcon*((sqrt(removed$SE...Control/removed$Mean...Control^2)))))

# do weighted meta-analysis
library(metafor)
library(glmulti)

## 1) using SE's as weights
consv.model.1 <- rma.mv(yi=lnCVR,
                        mods=~Type.of.Response,
                        V=varCVR,
                        random = list(~1|No./Year),
                        data=removed)
summary(consv.model.1)

# backtransform results
predict(consv.model.1,digits=3,transf=exp)

# check profile plots
# this makes sure it's not overparameterized...
# it IS if there's flat bits...
windows();
profile(consv.model.1)

# calculate SE of imputation
pred$SE...Experiment[is.na(pred$SE...Experiment)] <- apply(imp$imp$SE...Experiment,1,mean)
pred$SE...Control[is.na(pred$SE...Control)] <- apply(imp$imp$SE...Control,1,mean)

# calculate CVs
pred$CVExperiment <- pred$SE...Experiment/pred$Mean...Experiment
pred$CVControl <- pred$SE...Control/pred$Mean...Control

# calculate lnCVR
pred$lnCVR <- log(pred$CVExperiment/pred$CVControl)

# calculate correlations
pred <- pred[which(pred$SE...Experiment!=0),]
corexp <- cor(log(pred$Mean...Experiment),log(pred$SE...Experiment))
corcon <- cor(log(pred$Mean...Control),log(pred$SE...Control))

# calculate var(lnCVR)
pred$varCVR <- (pred$SE...Experiment)^2-(2*corexp*(sqrt((pred$SE...Experiment/pred$Mean...Experiment)^2))+(pred$SE...Experiment)^2-(2*corcon*((sqrt(pred$SE...Control/pred$Mean...Control^2)))))

# imputed data
var.impute <- rma.mv(yi=lnCVR,
              #mods=~Prey.type.s.,
              V=varCVR,
              random = list(~1|No./Year),
              data=pred)
summary(var.impute)


# backtransform results
predict(var.impute,digits=3,transf=exp)

########################################################################
# 4) IMPUTATION OF MISSING DATA

library(mice)
library(metafor)

# helper functions
eval(metafor:::.mice)

# set up predictor matrix
imp <- mice(pred,
            maxit=0)
predict <- imp$predictorMatrix # look at pred. matrix

# edit the predictor matrix
# make it all 0s, and just specify exactly the variables you'd like
# add in variables to be imputed (SEs)
predict[,] <- 0
predict["var","Mean...Experiment"] <- 1
predict["var","Mean...Control"] <- 1
predict["var","Log.RR.Effect.Size"] <- 1
predict["var","Study.ID"] <- 1
predict["var","Type.of.Response"] <- 1
predict["var","Natural"] <- 1
predict["var","X..Change.in.."] <- 1
predict["var","Experiment.Type"] <- 1
predict["var","Temporal.Scale..months."] <- 1
predict["var","Spatial.Scale..km2."] <- 1
predict["var","Sample.Size...Experiment.Plots"] <- 1
predict["var","Sample.Size...Control.Plots"] <- 1
predict["var","Sample.Size...experimental.years"] <- 1
predict["var","Sample.size...control.years"] <- 1

# let's impute
imp <- mice(pred, print=F,
            m=25, predictorMatrix = predict,
            seed=1234, maxit = 35)

# list imputations
View(imp$imp$var)

# completed datasets
imp_2 <- complete(imp,1);View(imp_2)

# check convergence/weirdness
densityplot(imp,~var)
plot(imp)

############################################################################
# 5) WEIGHTED GLM ANALYSIS OF EFFECT SIZE

library(metafor)
library(glmulti)

# calculate variances of logRR
removed$var <- ((removed$SE...Experiment^2)/(removed$Mean...Experiment)^2)+((removed$SE...Control^2)/(removed$Mean...Control)^2)
removed2$var <- ((removed2$SE...Experiment^2)/(removed2$Mean...Experiment)^2)+((removed2$SE...Control^2)/(removed2$Mean...Control)^2)
removed3$var <- ((removed3$SE...Experiment^2)/(removed3$Mean...Experiment)^2)+((removed3$SE...Control^2)/(removed3$Mean...Control)^2)
# using the rma.mv function
# nested RE: random=~1|id1/id2

## 1) using SE's as weights
consv.model.1 <- rma.mv(yi=Log.RR.Effect.Size,
                     #mods=~Prey.type.s.,
                     V=var,
                     random = list(~1|No./Year),
                     data=removed3)
summary(consv.model.1)


# backtransform results
predict(consv.model.1,digits=3,transf=exp)

# heterogeneity
W <- diag(1/consv.model.1$vi)
X <- model.matrix(consv.model.1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(consv.model.1$sigma2) / (sum(consv.model.1$sigma2) + (consv.model.1$k-consv.model.1$p)/sum(diag(P)))

# check profile plots
# this makes sure it's not overparameterized...
# it IS if there's flat bits...
windows();
profile(consv.model.1)

## 2) Model Selection - Ecology

# function to fit model
rma.glmulti <- function(formula, data, ...){
  rma.mv(formula, V=var, random = list(~1|No./Year),
         data=data, method="ML",...)
}

# run model
res <- glmulti(Log.RR.Effect.Size ~ Pred.Mass+Type.of.Response+
                 Experiment.Type + Natural + Predator.prey.Ratio+
                 Treat.Score + X..Change.in.. + Temporal.Scale..months.,
               data=removed, level=1, # main effects only
               fitfunction = rma.glmulti,
               crit="aicc",
               confsetsize=100) # models to keep

# summary
print(res)

# look at top models
top <- weightable(res)
top.aic <- top[top$aicc <= min(top$aicc) + 2,]; top.aic

# importance
windows(); plot(res,type="s")

# multimodel inference

# backwards stepwise selection
consv.model.1 <- rma.mv(yi=Log.RR.Effect.Size,
                        mods=~biodiv + Type.of.Response +
                        Pred.Mass + Prey.Mass+
                          Predator.prey.Ratio+Natural+
                          total + K + Near.K. ,
                        V=var,
                        random = list(~1|No./Year),
                        data=removed)

# final model
consv.model.6 <- rma.mv(yi=Log.RR.Effect.Size,
                        mods=~Type.of.Response + Pred.Mass + Natural,
                        V=var,
                        random = list(~1|No./Year),
                        data=removed)

## 3) Model Selection - Experimentation
res2 <- glmulti(Log.RR.Effect.Size ~ Experiment.Type*Natural +
                  Treat.Score + X..Change.in.. + Mult.Pred +
                  Spatial.Scale..km2. + Temporal.Scale..months. +
                  Replication,
               data=removed, level=1, # main effects only
               fitfunction = rma.glmulti,
               crit="aicc",
               confsetsize=100)

windows(); plot(res2,type="s")


# backwards stepwise selection
removed$Treat.Score <- as.factor(removed$Treat.Score)
consv.model.1 <- rma.mv(yi=Log.RR.Effect.Size,
                        mods=~Experiment.Type + Natural +
                          Treat.Score + X..Change.in.. + Mult.Pred +
                          Spatial.Scale..km2. + Temporal.Scale..months. +
                          Replication,
                        V=var,
                        random = list(~1|No./Year),
                        data=removed)

consv.model.6 <- rma.mv(yi=Log.RR.Effect.Size,
                        mods=~Experiment.Type +
                         X..Change.in..+ Treat.Score+
                          Temporal.Scale..months.,
                        V=var,
                        random = list(~1|No./Year),
                        data=removed)

# final model selection - combine things
consv.model.final.1 <- rma.mv(yi=Log.RR.Effect.Size,
                              mods=~Pred.Mass+Type.of.Response+
                                Experiment.Type + Natural +
                                Treat.Score + X..Change.in.. + Temporal.Scale..months.,
                              V=var,
                              random = list(~1|No./Year),
                              data=removed)

consv.model.final.2 <- rma.mv(yi=Log.RR.Effect.Size,
                              mods=~Type.of.Response+
                                Experiment.Type  + Natural +X..Change.in..,
                              V=var,
                              random = list(~1|No./Year),
                              data=removed)

##################################
# Fit the same model above with the imputed data
# assist AIC selection with Wald's test

# Ecological Model

consv.model.orig <- with(imp,rma.mv(yi=Log.RR.Effect.Size,
                                   mods=~biodiv + Type.of.Response +
                                     Pred.Mass + Prey.Mass+
                                     Predator.prey.Ratio+Natural+
                                     total + K,
                                   V=var,
                                   random = list(~1|No./Year)))

consv.model.imp <- with(imp,rma.mv(yi=Log.RR.Effect.Size,
                        mods=~Type.of.Response +
                          Pred.Mass + Natural,
                        V=var,
                        random = list(~1|No./Year)))

consv.model.imp <- with(imp,rma.mv(yi=Log.RR.Effect.Size,
                                   V=var,
                                   random = list(~1|No./Year)))

meanAIC <- mean(unlist(lapply(consv.model.imp$analyses, AIC)))
meanAIC

# pool results
pool <- pool(consv.model.imp)

# summarize results
round(summary(pool),4)

# Experimental Model
consv.model.orig <- with(imp,rma.mv(yi=Log.RR.Effect.Size,
                                    mods=~Experiment.Type+Natural +
                                      Treat.Score + X..Change.in.. + Mult.Pred +
                                      Spatial.Scale..km2. + Temporal.Scale..months. +
                                      Replication,
                                    V=var,
                                    random = list(~1|No./Year)))

consv.model.imp <- with(imp,rma.mv(yi=Log.RR.Effect.Size,
                                    mods=~Experiment.Type+Natural+
                                    Treat.Score + X..Change.in..+
                                    Temporal.Scale..months.,
                                    V=var,
                                    random = list(~1|No./Year)))

meanAIC <- mean(unlist(lapply(consv.model.imp$analyses, AIC)))
meanAIC


#########################################################################
# 6) PUBLICATION BIAS

# funnel plot
windows();funnel(consv.model.1)

# graph of precision vs. effect size
#calc psuedo CI
x <- seq(0.5,20,.01)
ylow <- .1175-(1.96*(1/x))
yupp <- .1175+(1.96*(1/x))


windows();
plot(1/removed$SE.Log.RR.Effect.Size^2,
     removed$Log.RR.Effect.Size,
     ylim=c(-4,4),xlim=c(0.5,20))
abline(h=0,lty=2)
abline(h=.1175,lty=1,col="red")
lines(x,ylow)
lines(x,yupp)

# test
regtest(x=removed$Log.RR.Effect.Size,
        vi=removed$var)
ranktest(x=removed$Log.RR.Effect.Size,
         vi=removed$var)
regtest.rma(consv.model.1)

# trim and fill method
# NOTE: has no random effects stuff...edges estimate up a bit...
blahblah <- rma(yi=Log.RR.Effect.Size,
                        vi=var,
                        data=removed)

taf <- trimfill(blahblah)

windows();funnel(taf, legend=T)

#################################################################################
