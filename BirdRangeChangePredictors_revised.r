## CJA Bradshaw
## British birds - range change ecological correlates (Cagan Sekercioglu's data with contagion stats from Miguel)
## 25/01/2012; Mataelpino, Spain; updated 03/04/2012
## updated December 2012; February 2013; August 2013

## Remove everything
rm(list = ls())

{ ## Load libraries
  library(boot)
  #library(qpcR)
  library(lme4)
  #library(spdep)
  library(sampling)
  library(vegan)
  #library(FactoMineR)
}


## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC

# Set functions
AICc <- function(...) {
	models <- list(...)
	num.mod <- length(models)
	AICcs <- numeric(num.mod)
	ns <- numeric(num.mod)
	ks <- numeric(num.mod)
	AICc.vec <- rep(0,num.mod)
	for (i in 1:num.mod) {
		if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
		if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
		AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
		ns[i] <- n
		ks[i] <- k
		AICc.vec[i] <- AICcs[i]
	}
	return(AICc.vec)
}

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
	fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
	AIC.vec <- c(AICc(fit.full),AICc(fit.null))
	dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
	ER <- wAIC.vec[1]/wAIC.vec[2]
	r.sq.adj <- as.numeric(summary(fit.full)[9])
	return(c(ER,r.sq.adj))
}

## source
source("/Applications/RStudio.app/Contents/Resources/R/new_lmer_AIC_tables3.R") 

## Import data
setwd("~/Dropbox/Bradshaw working papers backup/Bird Range Change/data/")

datwk.orig <- read.table("birdcagan2.csv", header=T, sep=",")
max.cells <- 2861
prc <- 100*(datwk.orig$t2 - datwk.orig$t1)/datwk.orig$t1
lt1prp <- logit(datwk.orig$t1/max.cells)
lt2prp <- logit(datwk.orig$t2/max.cells)
t2prp <- datwk.orig$t2/max.cells
len.dat <- dim(datwk.orig)[1]
lMass <- log(datwk.orig$Avg..Mass)
Clutch <- ((datwk.orig$Min+datwk.orig$Max)/2)
Long <- (datwk.orig$NatL)
cont <- datwk.orig$cont.avg
tFb <- datwk.orig$feb_tmn_Mean
tJl <- datwk.orig$jul_tmn_Mean
pcp <- datwk.orig$an_sum_pre_Mean

threat <- datwk.orig$Threat
UKthreat <- datwk.orig$Ukthreat
drct <- as.factor(ifelse(prc > 1, 1, 0)) # direction of change (contraction = 0; expansion = 1)
natDisp <- datwk.orig$natDisp
lnatDisp <- log10(natDisp)

hab <- datwk.orig$Habitat
datwk <- data.frame(datwk.orig,prc,lt1prp,lt2prp,t2prp,lMass,Clutch,Long,cont,tFb,tJl,pcp,threat,UKthreat,drct,natDisp,lnatDisp,hab)

## remove outlier
datwk <- datwk[-46,]

dat.ext <- subset(datwk, prc < 0)
dat.exp <- subset(datwk, prc >= 0)

datwk.nw <- subset(datwk, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
datwk.nw$Order <- factor(datwk.nw$Order)

dat.ext.nw <- subset(dat.ext, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
dat.ext.nw$Order <- factor(dat.ext.nw$Order)

dat.exp.nw <- subset(dat.exp, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
dat.exp.nw$Order <- factor(dat.exp.nw$Order)

datwk.nw  <- subset(datwk, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
datwk.nw$Order <- factor(datwk.nw$Order)



## correlation matrix for GLM input data
attach(datwk)
clump <- datwk.orig$clump
tFbmin <- datwk.orig$mean_feb6272tmin
tJlmn <- datwk.orig$mean_jul6272tmean
dat.cor.traits <- data.frame(lMass,Clutch,Long,natDisp,tFb,tFbmin,tJl,tJlmn,pcp,cont,clump,lt1prp)
dat.cor.ext <- data.frame(lt1prp,cont,clump,tFb,tFbmin,tJl,tJlmn,pcp)
detach(datwk)

cor.mat.traits <- cor(dat.cor.traits,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(dat.cor.traits)[2]
for (c in 1:lvar) {
	cor.mat.traits[c,c:lvar] <- NA}
cor.mat.traits[,-lvar]

cor.mat.ext <- cor(dat.cor.ext,y=NULL,use="complete.obs",method="spearman")
for (c in 1:8) {
	cor.mat.ext[c,c:8] <- NA}
cor.mat.ext[,-8]

##################################################
## 'contingency' analysis for UK threat
## model set
mod1 <- "t2prp~lt1prp+UKthreat+drct+drct*UKthreat"
mod2 <- "t2prp~lt1prp+UKthreat+drct"
mod3 <- "t2prp~lt1prp+UKthreat"
mod4 <- "t2prp~lt1prp+drct"
mod5 <- "t2prp~lt1prp"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[5]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

n.sample <- fit$df.null+1
n.sample


## repeat to estimate m-a standardised coefficients
mod1 <- "t2prp~lt1prp+drct/(UKthreat)+0"
mod2 <- "t2prp~lt1prp+drct+UKthreat+0"
mod3 <- "t2prp~lt1prp+UKthreat+0"
mod4 <- "t2prp~lt1prp+drct+0"
mod5 <- "t2prp~lt1prp+0"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=12,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","drct0","drct1","UKthreatG","UKthreatA","UKthreatR","drct0:UKthreatG","drct1:UKthreatG","drct0:UKthreatA","drct1:UKthreatA","drct0:UKthreatR","drct1:UKthreatR")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff
  #term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  term.labs[[i]] <- row.names(summ.fit[[i]]$coefficients)
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+1):(2*terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[5]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or # expansion = 1; contraction = 0


####################################################################################
## GLM
## Traits
## all species
## response = t2 (t1 as control variable)
## model set

mod1 <- "t2prp~lt1prp+lMass+Clutch+Long+natDisp"
mod2 <- "t2prp~lt1prp+lMass+Clutch+natDisp"
mod3 <- "t2prp~lt1prp+lMass+Long+natDisp"
mod4 <- "t2prp~lt1prp+Clutch+Long+natDisp"
mod5 <- "t2prp~lt1prp+lMass+Clutch"
mod6 <- "t2prp~lt1prp+lMass+Long"
mod7 <- "t2prp~lt1prp+Clutch+Long"
mod8 <- "t2prp~lt1prp+Clutch+natDisp"
mod9 <- "t2prp~lt1prp+lMass+natDisp"
mod10 <- "t2prp~lt1prp+Long+natDisp"
mod11 <- "t2prp~lt1prp+lMass"
mod12 <- "t2prp~lt1prp+Clutch"
mod13 <- "t2prp~lt1prp+Long"
mod14 <- "t2prp~lt1prp+natDisp"
mod15 <- "t2prp~lt1prp"
mod16 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[15]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="Clutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="Long",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

plot(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),pch=19)
data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]))


## contractors only
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2prp~lt1prp+lMass+Clutch+Long+natDisp"
mod2 <- "t2prp~lt1prp+lMass+Clutch+natDisp"
mod3 <- "t2prp~lt1prp+lMass+Long+natDisp"
mod4 <- "t2prp~lt1prp+Clutch+Long+natDisp"
mod5 <- "t2prp~lt1prp+lMass+Clutch"
mod6 <- "t2prp~lt1prp+lMass+Long"
mod7 <- "t2prp~lt1prp+Clutch+Long"
mod8 <- "t2prp~lt1prp+Clutch+natDisp"
mod9 <- "t2prp~lt1prp+lMass+natDisp"
mod10 <- "t2prp~lt1prp+Long+natDisp"
mod11 <- "t2prp~lt1prp+lMass"
mod12 <- "t2prp~lt1prp+Clutch"
mod13 <- "t2prp~lt1prp+Long"
mod14 <- "t2prp~lt1prp+natDisp"
mod15 <- "t2prp~lt1prp"
mod16 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[15]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="Clutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="Long",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)



## expanders only
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2prp~lt1prp+lMass+Clutch+Long+natDisp"
mod2 <- "t2prp~lt1prp+lMass+Clutch+natDisp"
mod3 <- "t2prp~lt1prp+lMass+Long+natDisp"
mod4 <- "t2prp~lt1prp+Clutch+Long+natDisp"
mod5 <- "t2prp~lt1prp+lMass+Clutch"
mod6 <- "t2prp~lt1prp+lMass+Long"
mod7 <- "t2prp~lt1prp+Clutch+Long"
mod8 <- "t2prp~lt1prp+Clutch+natDisp"
mod9 <- "t2prp~lt1prp+lMass+natDisp"
mod10 <- "t2prp~lt1prp+Long+natDisp"
mod11 <- "t2prp~lt1prp+lMass"
mod12 <- "t2prp~lt1prp+Clutch"
mod13 <- "t2prp~lt1prp+Long"
mod14 <- "t2prp~lt1prp+natDisp"
mod15 <- "t2prp~lt1prp"
mod16 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[15]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="Clutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="Long",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)



#########################################
## climate
## all species
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2prp~lt1prp+tFb+tJl+pcp"
mod2 <- "t2prp~lt1prp+tFb+tJl"
mod3 <- "t2prp~lt1prp+tJl+pcp"
mod4 <- "t2prp~lt1prp+tFb+pcp"
mod5 <- "t2prp~lt1prp+tFb"
mod6 <- "t2prp~lt1prp+tJl"
mod7 <- "t2prp~lt1prp+pcp"
mod8 <- "t2prp~lt1prp"
mod9 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[8]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## contractors only
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2prp~lt1prp+tFb+tJl+pcp"
mod2 <- "t2prp~lt1prp+tFb+tJl"
mod3 <- "t2prp~lt1prp+tJl+pcp"
mod4 <- "t2prp~lt1prp+tFb+pcp"
mod5 <- "t2prp~lt1prp+tFb"
mod6 <- "t2prp~lt1prp+tJl"
mod7 <- "t2prp~lt1prp+pcp"
mod8 <- "t2prp~lt1prp"
mod9 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[8]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## expanders only
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2prp~lt1prp+tFb+tJl+pcp"
mod2 <- "t2prp~lt1prp+tFb+tJl"
mod3 <- "t2prp~lt1prp+tJl+pcp"
mod4 <- "t2prp~lt1prp+tFb+pcp"
mod5 <- "t2prp~lt1prp+tFb"
mod6 <- "t2prp~lt1prp+tJl"
mod7 <- "t2prp~lt1prp+pcp"
mod8 <- "t2prp~lt1prp"
mod9 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[8]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


######################################
## distribution
## all species
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+clump"
mod2 <- "t2prp~clump"
mod3 <- "t2prp~lt1prp"
mod4 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[3]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## contractors only
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+clump"
mod2 <- "t2prp~clump"
mod3 <- "t2prp~lt1prp"
mod4 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[3]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## expanders only
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+clump"
mod2 <- "t2prp~clump"
mod3 <- "t2prp~lt1prp"
mod4 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove outlier 
dat.exp.nw2 <- dat.exp.nw[-20,]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw2)), data=dat.exp.nw2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[3]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw2)), data=dat.exp.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw2)), data=dat.exp.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

data.frame(topfit1$model,residuals(topfit1,type="partial"))


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
## Combined
## all species with drct (direction of change) included; Mass & Feb temp
## response = t2 (with t1 as control)

# update data.frame
datnew <- na.omit(data.frame(datwk.nw$t1,datwk.nw$t2,datwk.nw$t2prp,datwk.nw$lt1prp,datwk.nw$lMass,datwk.nw$tFb,datwk.nw$clump))
colnames(datnew) <- c("t1","t2","t2prp","lt1prp","lMass","tFb","clump")

## test for nonlinear relationship between R1 and clumpiness
plot(datnew$t1,datnew$clump,pch=19,xlab="logit proportional R1",ylab="clumpiness")

model1 <- lm(datwk.nw$clump~datwk.nw$lt1prp,data=datwk.nw)
model2 <- lm(datwk.nw$clump~poly(datwk.nw$lt1prp,2),data=datwk.nw)
model3 <- lm(datwk.nw$clump~1,data=datwk.nw)

AICc.vec <- c(AICc(model1),AICc(model2),AICc(model3))
dAICc.vec <- delta.IC(AICc.vec)
wAICc.vec <- weight.IC(dAICc.vec)
wAICc.vec


## model set
mod1 <- "t2prp~lt1prp+drct+lMass+tFb+clump+drct*lMass+drct*tFb+drct*clump"
mod2 <- "t2prp~lt1prp+drct+lMass+tFb+clump+drct*clump"
mod3 <- "t2prp~lt1prp+drct+lMass+tFb+clump+drct*tFb"
mod4 <- "t2prp~lt1prp+drct+lMass+tFb+clump+drct*lMass"
mod5 <- "t2prp~lt1prp+drct+lMass+tFb+clump"
mod6 <- "t2prp~lt1prp+drct+lMass+tFb"
mod7 <- "t2prp~lt1prp+drct+lMass+clump"
mod8 <- "t2prp~lt1prp+drct+tFb+clump"
mod9 <- "t2prp~lt1prp+drct+lMass+drct*lMass"
mod10 <- "t2prp~lt1prp+drct+tFb+drct*tFb"
mod11 <- "t2prp~lt1prp+drct+clump+drct*clump"
mod12 <- "t2prp~lt1prp+drct+lMass"
mod13 <- "t2prp~lt1prp+drct+tFb"
mod14 <- "t2prp~lt1prp+drct+clump"
mod15 <- "t2prp~lt1prp+drct"
mod16 <- "t2prp~lt1prp"
mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## reparameterise to get m-a standardised betas
## model set
mod1 <- "t2prp~lt1prp+drct/(lMass+tFb+clump)+0"
mod2 <- "t2prp~lt1prp+drct/(clump)+lMass+tFb+0"
mod3 <- "t2prp~lt1prp+drct/(tFb)+lMass+clump+0"
mod4 <- "t2prp~lt1prp+drct/(lMass)+tFb+clump+0"
mod5 <- "t2prp~lt1prp+drct+lMass+tFb+clump+0"
mod6 <- "t2prp~lt1prp+drct+lMass+tFb+0"
mod7 <- "t2prp~lt1prp+drct+lMass+clump+0"
mod8 <- "t2prp~lt1prp+drct+tFb+clump+0"
mod9 <- "t2prp~lt1prp+drct/(lMass)+0"
mod10 <- "t2prp~lt1prp+drct/(tFb)+0"
mod11 <- "t2prp~lt1prp+drct/(clump)+0"
mod12 <- "t2prp~lt1prp+drct+lMass+0"
mod13 <- "t2prp~lt1prp+drct+tFb+0"
mod14 <- "t2prp~lt1prp+drct+clump+0"
mod15 <- "t2prp~lt1prp+drct+0"
mod16 <- "t2prp~lt1prp+0"
#mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=12,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","drct0","drct1","lMass","tFb","clump","drct0:lMass","drct1:lMass","drct0:tFb","drct1:tFb","drct0:clump","drct1:clump")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff
  #term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  term.labs[[i]] <- row.names(summ.fit[[i]]$coefficients)
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+1):(2*terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or # expansion = 1; contraction = 0



## all species with drct (direction of change) included; Mass & precip
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+drct+lMass+pcp+clump+drct*lMass+drct*pcp+drct*clump"
mod2 <- "t2prp~lt1prp+drct+lMass+pcp+clump+drct*clump"
mod3 <- "t2prp~lt1prp+drct+lMass+pcp+clump+drct*pcp"
mod4 <- "t2prp~lt1prp+drct+lMass+pcp+clump+drct*lMass"
mod5 <- "t2prp~lt1prp+drct+lMass+pcp+clump"
mod6 <- "t2prp~lt1prp+drct+lMass+pcp"
mod7 <- "t2prp~lt1prp+drct+lMass+clump"
mod8 <- "t2prp~lt1prp+drct+pcp+clump"
mod9 <- "t2prp~lt1prp+drct+lMass+drct*lMass"
mod10 <- "t2prp~lt1prp+drct+pcp+drct*pcp"
mod11 <- "t2prp~lt1prp+drct+clump+drct*clump"
mod12 <- "t2prp~lt1prp+drct+lMass"
mod13 <- "t2prp~lt1prp+drct+pcp"
mod14 <- "t2prp~lt1prp+drct+clump"
mod15 <- "t2prp~lt1prp+drct"
mod16 <- "t2prp~lt1prp"
mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Precipitation",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## reparameterise to get m-a standardised betas
## model set
mod1 <- "t2prp~lt1prp+drct/(lMass+pcp+clump)+0"
mod2 <- "t2prp~lt1prp+drct/(clump)+lMass+pcp+0"
mod3 <- "t2prp~lt1prp+drct/(pcp)+lMass+clump+0"
mod4 <- "t2prp~lt1prp+drct/(lMass)+pcp+clump+0"
mod5 <- "t2prp~lt1prp+drct+lMass+pcp+clump+0"
mod6 <- "t2prp~lt1prp+drct+lMass+pcp+0"
mod7 <- "t2prp~lt1prp+drct+lMass+clump+0"
mod8 <- "t2prp~lt1prp+drct+pcp+clump+0"
mod9 <- "t2prp~lt1prp+drct/(lMass)+0"
mod10 <- "t2prp~lt1prp+drct/(pcp)+0"
mod11 <- "t2prp~lt1prp+drct/(clump)+0"
mod12 <- "t2prp~lt1prp+drct+lMass+0"
mod13 <- "t2prp~lt1prp+drct+pcp+0"
mod14 <- "t2prp~lt1prp+drct+clump+0"
mod15 <- "t2prp~lt1prp+drct+0"
mod16 <- "t2prp~lt1prp+0"
#mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=12,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","drct0","drct1","lMass","pcp","clump","drct0:lMass","drct1:lMass","drct0:pcp","drct1:pcp","drct0:clump","drct1:clump")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff
  #term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  term.labs[[i]] <- row.names(summ.fit[[i]]$coefficients)
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+1):(2*terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or # expansion = 1; contraction = 0


## all species with drct (direction of change) included; natDisp & Feb temp
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+drct+natDisp+tFb+clump+drct*natDisp+drct*tFb+drct*clump"
mod2 <- "t2prp~lt1prp+drct+natDisp+tFb+clump+drct*clump"
mod3 <- "t2prp~lt1prp+drct+natDisp+tFb+clump+drct*tFb"
mod4 <- "t2prp~lt1prp+drct+natDisp+tFb+clump+drct*natDisp"
mod5 <- "t2prp~lt1prp+drct+natDisp+tFb+clump"
mod6 <- "t2prp~lt1prp+drct+natDisp+tFb"
mod7 <- "t2prp~lt1prp+drct+natDisp+clump"
mod8 <- "t2prp~lt1prp+drct+tFb+clump"
mod9 <- "t2prp~lt1prp+drct+natDisp+drct*natDisp"
mod10 <- "t2prp~lt1prp+drct+tFb+drct*tFb"
mod11 <- "t2prp~lt1prp+drct+clump+drct*clump"
mod12 <- "t2prp~lt1prp+drct+natDisp"
mod13 <- "t2prp~lt1prp+drct+tFb"
mod14 <- "t2prp~lt1prp+drct+clump"
mod15 <- "t2prp~lt1prp+drct"
mod16 <- "t2prp~lt1prp"
mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

## reparameterise to get m-a standardised betas
## model set
mod1 <- "t2prp~lt1prp+drct/(natDisp+tFb+clump)+0"
mod2 <- "t2prp~lt1prp+drct/(clump)+natDisp+tFb+0"
mod3 <- "t2prp~lt1prp+drct/(tFb)+natDisp+clump+0"
mod4 <- "t2prp~lt1prp+drct/(natDisp)+tFb+clump+0"
mod5 <- "t2prp~lt1prp+drct+natDisp+tFb+clump+0"
mod6 <- "t2prp~lt1prp+drct+natDisp+tFb+0"
mod7 <- "t2prp~lt1prp+drct+natDisp+clump+0"
mod8 <- "t2prp~lt1prp+drct+tFb+clump+0"
mod9 <- "t2prp~lt1prp+drct/(natDisp)+0"
mod10 <- "t2prp~lt1prp+drct/(tFb)+0"
mod11 <- "t2prp~lt1prp+drct/(clump)+0"
mod12 <- "t2prp~lt1prp+drct+natDisp+0"
mod13 <- "t2prp~lt1prp+drct+tFb+0"
mod14 <- "t2prp~lt1prp+drct+clump+0"
mod15 <- "t2prp~lt1prp+drct+0"
mod16 <- "t2prp~lt1prp+0"
#mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=12,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","drct0","drct1","natDisp","tFb","clump","drct0:natDisp","drct1:natDisp","drct0:tFb","drct1:tFb","drct0:clump","drct1:clump")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff
  #term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  term.labs[[i]] <- row.names(summ.fit[[i]]$coefficients)
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+1):(2*terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or # expansion = 1; contraction = 0


## all species with drct (direction of change) included; natDisp & precip
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+drct+natDisp+pcp+clump+drct*natDisp+drct*pcp+drct*clump"
mod2 <- "t2prp~lt1prp+drct+natDisp+pcp+clump+drct*clump"
mod3 <- "t2prp~lt1prp+drct+natDisp+pcp+clump+drct*pcp"
mod4 <- "t2prp~lt1prp+drct+natDisp+pcp+clump+drct*natDisp"
mod5 <- "t2prp~lt1prp+drct+natDisp+pcp+clump"
mod6 <- "t2prp~lt1prp+drct+natDisp+pcp"
mod7 <- "t2prp~lt1prp+drct+natDisp+clump"
mod8 <- "t2prp~lt1prp+drct+pcp+clump"
mod9 <- "t2prp~lt1prp+drct+natDisp+drct*natDisp"
mod10 <- "t2prp~lt1prp+drct+pcp+drct*pcp"
mod11 <- "t2prp~lt1prp+drct+clump+drct*clump"
mod12 <- "t2prp~lt1prp+drct+natDisp"
mod13 <- "t2prp~lt1prp+drct+pcp"
mod14 <- "t2prp~lt1prp+drct+clump"
mod15 <- "t2prp~lt1prp+drct"
mod16 <- "t2prp~lt1prp"
mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

## reparameterise to get m-a standardised betas
## model set
mod1 <- "t2prp~lt1prp+drct/(natDisp+pcp+clump)+0"
mod2 <- "t2prp~lt1prp+drct/(clump)+natDisp+pcp+0"
mod3 <- "t2prp~lt1prp+drct/(pcp)+natDisp+clump+0"
mod4 <- "t2prp~lt1prp+drct/(natDisp)+pcp+clump+0"
mod5 <- "t2prp~lt1prp+drct+natDisp+pcp+clump+0"
mod6 <- "t2prp~lt1prp+drct+natDisp+pcp+0"
mod7 <- "t2prp~lt1prp+drct+natDisp+clump+0"
mod8 <- "t2prp~lt1prp+drct+pcp+clump+0"
mod9 <- "t2prp~lt1prp+drct/(natDisp)+0"
mod10 <- "t2prp~lt1prp+drct/(pcp)+0"
mod11 <- "t2prp~lt1prp+drct/(clump)+0"
mod12 <- "t2prp~lt1prp+drct+natDisp+0"
mod13 <- "t2prp~lt1prp+drct+pcp+0"
mod14 <- "t2prp~lt1prp+drct+clump+0"
mod15 <- "t2prp~lt1prp+drct+0"
mod16 <- "t2prp~lt1prp+0"
#mod17 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=12,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","drct0","drct1","natDisp","pcp","clump","drct0:natDisp","drct1:natDisp","drct0:pcp","drct1:pcp","drct0:clump","drct1:clump")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(datwk.nw)), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff
  #term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  term.labs[[i]] <- row.names(summ.fit[[i]]$coefficients)
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+1):(2*terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[16]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or # expansion = 1; contraction = 0





## contractors only
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+Long+natDisp+tFb+pcp+clump"
mod2 <- "t2prp~lt1prp+Long+natDisp"
mod3 <- "t2prp~lt1prp+tFb+pcp"
mod4 <- "t2prp~lt1prp+clump"
mod5 <- "t2prp~lt1prp+Long+natDisp+clump"
mod6 <- "t2prp~lt1prp+tFb+pcp+clump"
mod7 <- "t2prp~lt1prp+Long+natDisp+tFb+pcp"
mod8 <- "t2prp~lt1prp+Long"
mod9 <- "t2prp~lt1prp+natDisp"
mod10 <- "t2prp~lt1prp+tFb"
mod11 <- "t2prp~lt1prp+pcp"
mod12 <- "t2prp~lt1prp"
mod13 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=6,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","Long","natDisp","tFb","pcp","clump")

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
		summ.fit[[i]] <- summary(fit)
		coeffs[[i]] <- fit$coeff[-1]
		term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
		terml <- length(term.labs[[i]])
		coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+3):(terml+2+terml)]
		coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
		sub <- rep(0,terml)
		for (j in 1:terml) {
			sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
		}
		coeff.st.mat[i,sub] <- coeffs.st[[i]]

		print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[12]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or


## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(3,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="precip",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="Long",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample




## expanders only
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+Long+natDisp+tFb+pcp+clump"
mod2 <- "t2prp~lt1prp+Long+natDisp"
mod3 <- "t2prp~lt1prp+tFb+pcp"
mod4 <- "t2prp~lt1prp+clump"
mod5 <- "t2prp~lt1prp+Long+natDisp+clump"
mod6 <- "t2prp~lt1prp+tFb+pcp+clump"
mod7 <- "t2prp~lt1prp+Long+natDisp+tFb+pcp"
mod8 <- "t2prp~lt1prp+Long"
mod9 <- "t2prp~lt1prp+natDisp"
mod10 <- "t2prp~lt1prp+tFb"
mod11 <- "t2prp~lt1prp+pcp"
mod12 <- "t2prp~lt1prp"
mod13 <- "t2prp~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=6,nrow=Modnum)
colnames(coeff.st.mat) <- c("lt1prp","Long","natDisp","tFb","pcp","clump")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"),  weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  summ.fit[[i]] <- summary(fit)
  coeffs[[i]] <- fit$coeff[-1]
  term.labs[[i]] <- attr(summ.fit[[i]]$terms,"term.labels")
  terml <- length(term.labs[[i]])
  coeffs.se[[i]] <- as.numeric(summ.fit[[i]]$coefficients)[(terml+3):(terml+2+terml)]
  coeffs.st[[i]] <- as.numeric(coeffs[[i]]/coeffs.se[[i]])
  sub <- rep(0,terml)
  for (j in 1:terml) {
    sub[j] <- which(term.labs[[i]][j] == colnames(coeff.st.mat))
  }
  coeff.st.mat[i,sub] <- coeffs.st[[i]]
  
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[12]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

coeff.pr.mat <- ifelse(coeff.st.mat == 0, 0, 1)
wt.mat <- as.numeric(sumtable[,5]) * coeff.pr.mat
wt.sum <- apply(wt.mat,2,"sum")
rs.wt.mat <- t(t(wt.mat)/wt.sum)
rs.coeff.st.mat <- coeff.st.mat*rs.wt.mat
rs.coeff.st <- apply(rs.coeff.st.mat,2,"sum")
rs.coeff.st
rs.coeff.st.abs.or <- sort(abs(rs.coeff.st),decreasing=T)
rs.coeff.st.abs.or


## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(3,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.exp.nw)), data=dat.exp.nw, na.action=na.omit)
#plot1 <- termplot(topfit1,terms="lt1prp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="clump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="precip",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="natDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="Long",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


###############################################################################################
###############################################################################################
## Combined GLMM with Order as random effect
## contractors only
## response = t2 (with t1 as control)
## model set
mod1 <- "t2prp~lt1prp+Long+lMass+tFb+pcp+clump+(1|Order)"
mod2 <- "t2prp~lt1prp+Long+lMass+(1|Order)"
mod3 <- "t2prp~lt1prp+tFb+pcp+(1|Order)"
mod4 <- "t2prp~lt1prp+clump+(1|Order)"
mod5 <- "t2prp~lt1prp+Long+lMass+clump+(1|Order)"
mod6 <- "t2prp~lt1prp+tFb+pcp+clump+(1|Order)"
mod7 <- "t2prp~lt1prp+Long+lMass+tFb+pcp+(1|Order)"
mod8 <- "t2prp~lt1prp+Long+(1|Order)"
mod9 <- "t2prp~lt1prp+lMass+(1|Order)"
mod10 <- "t2prp~lt1prp+tFb+(1|Order)"
mod11 <- "t2prp~lt1prp+pcp+(1|Order)"
mod12 <- "t2prp~lt1prp+(1|Order)"
mod13 <- "t2prp~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = mod.list[[12]], order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

best <- lmer(as.formula(mod.vec[1]),family=binomial(link="logit"), weights = rep(max.cells, nrow(dat.ext.nw)), data=dat.ext.nw, na.action=na.omit)
print(best)
hist(resid(best), breaks=300)
qqnorm(resid(best))
qqline(resid(best), col="blue")
coef(best)$Order

