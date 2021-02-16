## CJA Bradshaw
## British birds - range change ecological correlates (Cagan Sekercioglu's data with contagion stats from Miguel)
## 25/01/2012; Mataelpino, Spain; updated 03/04/2012
## updated December 2012; February 2013

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

# { ## Functions
#   
#   scale <- function(x) {
#     x <- (x-min(x,na.rm=T)) / (max(x,na.rm=T)-min(x,na.rm=T))
#     x
#   }
#     
#   meancent = function(x,y){
#     out = {}
#     for (i in unique(y)){
#       out[y==i] = x[y==i]-mean(x[y==i])
#     }
#     out
#   }
# }


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
#source("C:/Program Files/R/new_lmer_AIC_tables3.R") 
source("/Applications/RStudio.app/Contents/Resources/R/new_lmer_AIC_tables3.R") 

## Import data
#setwd("C:\\Users\\Public\\Documents\\Papers\\Extinctions\\Bird Range Change\\")
setwd("~/Documents/Papers/Birds/UK Bird Range Change/data/")

datwk.orig <- read.table("birdcagan.csv", header=T, sep=",")
len.dat <- dim(datwk.orig)[1]
prc <- 100*(datwk.orig$t2 - datwk.orig$t1)/datwk.orig$t1
prc.tr <- prc + 101
lprc <- log10(abs(prc))
lprc.tr <- log10(prc.tr)
t1sq <- datwk.orig$t1^2
diff <- datwk.orig$t1 - datwk.orig$t2
diff.tr <- diff + abs(range(diff)[1]) + 1
frc <- 100*(datwk.orig$PredFutRangeSize3 - datwk.orig$PredPresRangeSize3)/datwk.orig$PredPresRangeSize3
frc.tr <- frc + 45
frdiff <- datwk.orig$PredFutRangeSize3 - datwk.orig$PredPresRangeSize3


lMass <- (as.vector(scale(datwk.orig$Avg..Mass,center=F,scale=T)))
Mass <- (datwk.orig$Avg..Mass)

lClutch <- as.vector(scale((datwk.orig$Min+datwk.orig$Max)/2,center=T,scale=T))
Clutch <- ((datwk.orig$Min+datwk.orig$Max)/2)

lLong <- log10(as.vector(scale(datwk.orig$NatL,center=F,scale=T)))
Long <- (datwk.orig$NatL)

Mat.Avg <- 12*((datwk.orig$MatMin + datwk.orig$MatMax)/2)
MinAge <- datwk.orig$MinAge
lMAB <- log10(ifelse((is.na(Mat.Avg) == T), MinAge, Mat.Avg))
MAB <- (ifelse((is.na(Mat.Avg) == T), MinAge, Mat.Avg))

lcont <- as.vector(scale(datwk.orig$cont.avg,center=F,scale=T))
lclump <- as.vector(scale(datwk.orig$clump,center=F,scale=T))

## compare contagion & clumpiness
plot(lcont,lclump, pch=19)
linreg.ER(lcont,lclump)

tFb <- scale(datwk.orig$feb_tmn_Mean,center=T,scale=T)
tJl <- scale(datwk.orig$jul_tmn_Mean,center=T,scale=T)
pcp <- scale(datwk.orig$an_sum_pre_Mean,center=T,scale=T)

threat <- datwk.orig$Threat
UKthreat <- datwk.orig$Ukthreat
drct <- as.factor(ifelse(prc > 1, 1, 0)) # direction of change (contraction = 0; expansion = 1)
natDisp <- datwk.orig$natDisp
lnatDisp <- log10(natDisp)

# hab <- ifelse((datwk.orig$Habitat == "forest" | datwk.orig$Habitat == "shrub" | datwk.orig$Habitat == "woodland"), "forest", datwk.orig$Habitat)
# hab <- ifelse((hab == "artificial" | hab == "grassland"), "farm", hab)
hab <- datwk.orig$Habitat

datwk <- data.frame(datwk.orig,t1sq,prc,prc.tr,diff,diff.tr,lMass,lClutch,lLong,lMAB,lcont,lclump,tFb,tJl,pcp,frc,frc.tr,frdiff,Clutch,MAB,Long,Mass,lprc,lprc.tr,threat,UKthreat,drct,natDisp,lnatDisp,hab)

dat.ext <- subset(datwk, prc < 0)
dat.ext <- dat.ext[-5, ] ## remove outlier

dat.exp <- subset(datwk, prc >= 0)

datwk.nw <- subset(datwk, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
datwk.nw$Order <- factor(datwk.nw$Order)

dat.ext.nw <- subset(dat.ext, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
dat.ext.nw$Order <- factor(dat.ext.nw$Order)
#write.csv(dat.ext.nw,file="truncators.csv")

dat.exp.nw <- subset(dat.exp, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
dat.exp.nw$Order <- factor(dat.exp.nw$Order)
#write.csv(dat.exp.nw,file="expanders.csv")

datwk.nw  <- subset(datwk, (Order == "Passeriformes") | (Order == "Accipitriformes") | (Order == "Columbiformes") | (Order == "Galliformes") | (Order == "Strigiformes"))
datwk.nw$Order <- factor(datwk.nw$Order)

plot(datwk.nw$frc,datwk.nw$prc,pch=19,xlim=c(-100,100),ylim=c(-100,100))
abline(h=0,lty=2)
abline(v=0,lty=2)

plot(datwk.nw$frc,datwk.nw$prc,pch=19)
abline(h=0,lty=2)
abline(v=0,lty=2)
sub.anom1 <- which(datwk.nw$prc < 0 & datwk.nw$frc > 0)
datwk.nw$Common[sub.anom1]
datwk.nw.pinc.decr <- datwk.nw[sub.anom1,]
bias <- datwk.nw.pinc.decr$frc+abs(datwk.nw.pinc.decr$prc)
plot(datwk.nw.pinc.decr$t1,bias,pch=19,xlab="range t1",ylab="prop bias to agreement")

sub.anom1 <- which(datwk.nw$prc > 0 & datwk.nw$frc < 0)
datwk.nw$Common[sub.anom1]
datwk.nw.pinc.decr <- datwk.nw[sub.anom1,]
bias <- abs(datwk.nw.pinc.decr$frc)+datwk.nw.pinc.decr$prc
plot(datwk.nw.pinc.decr$t1,bias,pch=19,xlab="range t1",ylab="prop bias to agreement")


plot(datwk.nw$frdiff,datwk.nw$diff,pch=19)
abline(h=0,lty=2)
abline(v=0,lty=2)
sub.anom1 <- which(datwk.nw$diff < 0 & datwk.nw$frdiff > 0)
datwk.nw$Common[sub.anom1]
datwk.nw.pinc.decr <- datwk.nw[sub.anom1,]
bias <- abs(datwk.nw.pinc.decr$diff)+datwk.nw.pinc.decr$frdiff
plot(datwk.nw.pinc.decr$t1,bias,pch=19,xlab="range t1",ylab="bias to agreement")

plot(datwk.nw$PredPresRangeSize3,datwk.nw$t1,pch=19)
fit1 <- lm(datwk.nw$t1 ~ datwk.nw$PredPresRangeSize3)
abline(fit1)
summary(fit1)


## 'contingency' analysis for UK threat
## model set
mod1 <- "lprc~UKthreat+drct+UKthreat*drct"
mod2 <- "lprc~UKthreat+drct"
mod3 <- "lprc~UKthreat"
mod4 <- "lprc~drct"
mod5 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)


## Natal dispersal check
plot((dat.ext.nw$lnatDisp),dat.ext.nw$lprc,pch=19,main="contractors",xlab="mean natal dispersal (km)",ylab="delta R")
fit.ext <- lm(dat.ext.nw$lprc ~ dat.ext.nw$lnatDisp)
abline(fit.ext,lty=2)
linreg.ER(dat.ext.nw$lnatDisp,dat.ext.nw$lprc)
plot(fit.ext)
hist(log10(dat.ext.nw$natDisp))
    
plot((dat.exp.nw$lnatDisp),dat.exp.nw$lprc,pch=19,main="expanders",xlab="mean natal dispersal (km)",ylab="delta R")
fit.exp <- lm(dat.exp.nw$lprc ~ dat.exp.nw$lnatDisp)
abline(fit.exp,lty=2)
linreg.ER(dat.exp.nw$lnatDisp,dat.exp.nw$lprc)
plot(fit.ext)
hist(log10(dat.exp.nw$natDisp))


## correlation matrix for GLM input data
attach(datwk)
clump <- datwk.orig$clump
dat.cor.traits <- data.frame(lMass,lClutch,lLong,lnatDisp,tFb,tJl,pcp,lcont,lclump,t1)
dat.cor.ext <- data.frame(t1,lcont,lclump,tFb,tJl,pcp)
detach(datwk)

cor.mat.traits <- cor(dat.cor.traits,y=NULL,use="complete.obs",method="spearman")
lvar <- dim(dat.cor.traits)[2]
for (c in 1:lvar) {
	cor.mat.traits[c,c:lvar] <- NA}
cor.mat.traits[,-lvar]

cor.mat.ext <- cor(dat.cor.ext,y=NULL,use="complete.obs",method="spearman")
for (c in 1:6) {
	cor.mat.ext[c,c:6] <- NA}
cor.mat.ext[,-6]

hist(diff,br=30,main="",xlab="delta cells (t2-t1)",ylab="frequency",col="grey",border="black")
hist.out <- hist(diff,br=30,main="",xlab="delta cells (t2-t1)",ylab="frequency",col="grey",border="black")
hist.out.dat <- data.frame(hist.out$mids,hist.out$counts)
colnames(hist.out.dat) <- c("mids","freq")
write.csv(hist.out.dat,file="histout.csv")


####################################################################################
## GLM
## Traits
## all species
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~lMass+lClutch+lLong"
mod2 <- "lprc~lMass+lClutch"
mod3 <- "lprc~lMass+lLong"
mod4 <- "lprc~lClutch+lLong"
mod5 <- "lprc~lMass"
mod6 <- "lprc~lClutch"
mod7 <- "lprc~lLong"
mod8 <- "lprc~1"

# mod1 <- "prc.tr~lMass+lClutch+lLong"
# mod2 <- "prc.tr~lMass+lClutch"
# mod3 <- "prc.tr~lMass+lLong"
# mod4 <- "prc.tr~lClutch+lLong"
# mod5 <- "prc.tr~lMass"
# mod6 <- "prc.tr~lClutch"
# mod7 <- "prc.tr~lLong"
# mod8 <- "prc.tr~1"
# 
# 
# mod1 <- "lprc.tr~lMass+lClutch+lLong"
# mod2 <- "lprc.tr~lMass+lClutch"
# mod3 <- "lprc.tr~lMass+lLong"
# mod4 <- "lprc.tr~lClutch+lLong"
# mod5 <- "lprc.tr~lMass"
# mod6 <- "lprc.tr~lClutch"
# mod7 <- "lprc.tr~lLong"
# mod8 <- "lprc.tr~1"

## Make model vector
#mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
#plot3 <- termplot(topfit1,terms="lMAB",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Age Primiparity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


###############################################
## response = t2 (t1 as control variable)
## model set

mod1 <- "t2~t1+lMass+lClutch+lLong"
mod2 <- "t2~t1+lMass+lClutch"
mod3 <- "t2~t1+lMass+lLong"
mod4 <- "t2~t1+lClutch+lLong"
mod5 <- "t2~t1+lMass"
mod6 <- "t2~t1+lClutch"
mod7 <- "t2~t1+lLong"
mod8 <- "t2~t1"
mod9 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


####################################################################################
## GLM
## Traits
## decliners only
## response = transformed % range change (prc)
## model set
# mod1 <- "abs(prc)~lMass+lClutch+lMAB+lLong"
# mod2 <- "abs(prc)~lMass+lClutch+lMAB"
# mod3 <- "abs(prc)~lMass+lClutch+lLong"
# mod4 <- "abs(prc)~lMass+lMAB+lLong"
# mod5 <- "abs(prc)~lClutch+lMAB+lLong"
# mod6 <- "abs(prc)~lMass+lClutch"
# mod7 <- "abs(prc)~lMass+lMAB"
# mod8 <- "abs(prc)~lMass+lLong"
# mod9 <- "abs(prc)~lClutch+lMAB"
# mod10 <- "abs(prc)~lClutch+lLong"
# mod11 <- "abs(prc)~lMAB+lLong"
# mod12 <- "abs(prc)~lMass"
# mod13 <- "abs(prc)~lClutch"
# mod14 <- "abs(prc)~lMAB"
# mod15 <- "abs(prc)~lLong"
# mod16 <- "abs(prc)~1"

# mod1 <- "abs(prc)~Mass+Clutch+MAB+Long"
# mod2 <- "abs(prc)~Mass+Clutch+MAB"
# mod3 <- "abs(prc)~Mass+Clutch+Long"
# mod4 <- "abs(prc)~Mass+MAB+Long"
# mod5 <- "abs(prc)~Clutch+MAB+Long"
# mod6 <- "abs(prc)~Mass+Clutch"
# mod7 <- "abs(prc)~Mass+MAB"
# mod8 <- "abs(prc)~Mass+Long"
# mod9 <- "abs(prc)~Clutch+MAB"
# mod10 <- "abs(prc)~Clutch+Long"
# mod11 <- "abs(prc)~MAB+Long"
# mod12 <- "abs(prc)~Mass"
# mod13 <- "abs(prc)~Clutch"
# mod14 <- "abs(prc)~MAB"
# mod15 <- "abs(prc)~Long"
# mod16 <- "abs(prc)~1"

# mod1 <- "abs(prc)~lMass+lClutch+lLong"
# mod2 <- "abs(prc)~lMass+lClutch"
# mod3 <- "abs(prc)~lMass+lLong"
# mod4 <- "abs(prc)~lClutch+lLong"
# mod5 <- "abs(prc)~lMass"
# mod6 <- "abs(prc)~lClutch"
# mod7 <- "abs(prc)~lLong"
# mod8 <- "abs(prc)~1"

mod1 <- "lprc~lMass+lClutch+lLong"
mod2 <- "lprc~lMass+lClutch"
mod3 <- "lprc~lMass+lLong"
mod4 <- "lprc~lClutch+lLong"
mod5 <- "lprc~lMass"
mod6 <- "lprc~lClutch"
mod7 <- "lprc~lLong"
mod8 <- "lprc~1"

## Make model vector
#mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## remove outlier
#dat.ext.nw2 <- dat.ext.nw[-c(32,61),]
dat.ext.nw2 <- dat.ext.nw[-c(32),]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
#plot3 <- termplot(topfit1,terms="lMAB",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Age Primiparity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prlprcM","lClutch","prlprcC","lLong","prlprcL")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


###############################################
## response = t2 (t1 as control variable)
## model set

mod1 <- "t2~t1+lMass+lClutch+lLong"
mod2 <- "t2~t1+lMass+lClutch"
mod3 <- "t2~t1+lMass+lLong"
mod4 <- "t2~t1+lClutch+lLong"
mod5 <- "t2~t1+lMass"
mod6 <- "t2~t1+lClutch"
mod7 <- "t2~t1+lLong"
mod8 <- "t2~t1"
mod9 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

###############################################
## contractors only, including natal dispersal
## model set
mod1 <- "lprc~lMass+lClutch+lLong+lnatDisp"
mod2 <- "lprc~lMass+lClutch"
mod3 <- "lprc~lMass+lLong"
mod4 <- "lprc~lClutch+lLong"
mod5 <- "lprc~lMass"
mod6 <- "lprc~lClutch"
mod7 <- "lprc~lLong"
mod8 <- "lprc~lMass+lClutch+lLong"
mod9 <- "lprc~lMass+lClutch+lnatDisp"
mod10 <- "lprc~lMass+lLong+lnatDisp"
mod11 <- "lprc~lClutch+lLong+lnatDisp"
mod12 <- "lprc~lMass+lnatDisp"
mod13 <- "lprc~lClutch+lnatDisp"
mod14 <- "lprc~lLong+lnatDisp"
mod15 <- "lprc~lnatDisp"
mod16 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## remove outlier
#dat.ext.nw2 <- dat.ext.nw[-c(32,61),]
dat.ext.nw2 <- dat.ext.nw[-c(32),]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,4),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
#plot3 <- termplot(topfit1,terms="lMAB",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Age Primiparity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="lnatDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prlprcM","lClutch","prlprcC","lLong","prlprcL")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


####################################################################################
## GLM
## Traits
## expanders only
## response = transformed % range change (prc)
## model set
# mod1 <- "abs(prc)~lMass+lClutch+lMAB+lLong"
# mod2 <- "abs(prc)~lMass+lClutch+lMAB"
# mod3 <- "abs(prc)~lMass+lClutch+lLong"
# mod4 <- "abs(prc)~lMass+lMAB+lLong"
# mod5 <- "abs(prc)~lClutch+lMAB+lLong"
# mod6 <- "abs(prc)~lMass+lClutch"
# mod7 <- "abs(prc)~lMass+lMAB"
# mod8 <- "abs(prc)~lMass+lLong"
# mod9 <- "abs(prc)~lClutch+lMAB"
# mod10 <- "abs(prc)~lClutch+lLong"
# mod11 <- "abs(prc)~lMAB+lLong"
# mod12 <- "abs(prc)~lMass"
# mod13 <- "abs(prc)~lClutch"
# mod14 <- "abs(prc)~lMAB"
# mod15 <- "abs(prc)~lLong"
# mod16 <- "abs(prc)~1"

# mod1 <- "abs(prc)~Mass+Clutch+MAB+Long"
# mod2 <- "abs(prc)~Mass+Clutch+MAB"
# mod3 <- "abs(prc)~Mass+Clutch+Long"
# mod4 <- "abs(prc)~Mass+MAB+Long"
# mod5 <- "abs(prc)~Clutch+MAB+Long"
# mod6 <- "abs(prc)~Mass+Clutch"
# mod7 <- "abs(prc)~Mass+MAB"
# mod8 <- "abs(prc)~Mass+Long"
# mod9 <- "abs(prc)~Clutch+MAB"
# mod10 <- "abs(prc)~Clutch+Long"
# mod11 <- "abs(prc)~MAB+Long"
# mod12 <- "abs(prc)~Mass"
# mod13 <- "abs(prc)~Clutch"
# mod14 <- "abs(prc)~MAB"
# mod15 <- "abs(prc)~Long"
# mod16 <- "abs(prc)~1"

# mod1 <- "abs(prc)~lMass+lClutch+lLong"
# mod2 <- "abs(prc)~lMass+lClutch"
# mod3 <- "abs(prc)~lMass+lLong"
# mod4 <- "abs(prc)~lClutch+lLong"
# mod5 <- "abs(prc)~lMass"
# mod6 <- "abs(prc)~lClutch"
# mod7 <- "abs(prc)~lLong"
# mod8 <- "abs(prc)~1"

mod1 <- "lprc~lMass+lClutch+lLong"
mod2 <- "lprc~lMass+lClutch"
mod3 <- "lprc~lMass+lLong"
mod4 <- "lprc~lClutch+lLong"
mod5 <- "lprc~lMass"
mod6 <- "lprc~lClutch"
mod7 <- "lprc~lLong"
mod8 <- "lprc~1"

## Make model vector
#mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## remove outlier
dat.exp.nw2 <- dat.exp.nw[-5,]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
#plot3 <- termplot(topfit1,terms="lMAB",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Age Primiparity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


# export partial residulals
out2 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out2) <- c("lMass","prlprcM","lClutch","prlprcC","lLong","prlprcL")
write.table(out2,file="exp.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


##########################################
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lMass+lClutch+lLong"
mod2 <- "t2~t1+lMass+lClutch"
mod3 <- "t2~t1+lMass+lLong"
mod4 <- "t2~t1+lClutch+lLong"
mod5 <- "t2~t1+lMass"
mod6 <- "t2~t1+lClutch"
mod7 <- "t2~t1+lLong"
mod8 <- "t2~t1"
mod9 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residuals
out1 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,3]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out1) <- c("lMass","prR2M","lClutch","prR2C","lLong","prR2L")
write.table(out1,file="ext.lh.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


###############################################
## expanders only, including natal dispersal

## model set
mod1 <- "lprc~lMass+lClutch+lLong+lnatDisp"
mod2 <- "lprc~lMass+lClutch"
mod3 <- "lprc~lMass+lLong"
mod4 <- "lprc~lClutch+lLong"
mod5 <- "lprc~lMass"
mod6 <- "lprc~lClutch"
mod7 <- "lprc~lLong"
mod8 <- "lprc~lMass+lClutch+lLong"
mod9 <- "lprc~lMass+lClutch+lnatDisp"
mod10 <- "lprc~lMass+lLong+lnatDisp"
mod11 <- "lprc~lClutch+lLong+lnatDisp"
mod12 <- "lprc~lMass+lnatDisp"
mod13 <- "lprc~lClutch+lnatDisp"
mod14 <- "lprc~lLong+lnatDisp"
mod15 <- "lprc~lnatDisp"
mod16 <- "lprc~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,4),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch Size",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
#plot3 <- termplot(topfit1,terms="lMAB",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Age Primiparity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="lnatDisp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Natal Dispersal",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## climate
## all species
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~tFb+tJl+pcp"
mod2 <- "lprc~tFb+tJl"
mod3 <- "lprc~tJl+pcp"
mod4 <- "lprc~tFb+pcp"
mod5 <- "lprc~tFb"
mod6 <- "lprc~tJl"
mod7 <- "lprc~pcp"
mod8 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample



## decliners only
## response = transformed % range change (prc)
## model set
# mod1 <- "abs(prc)~tFb+tJl+pcp"
# mod2 <- "abs(prc)~tFb+tJl"
# mod3 <- "abs(prc)~tJl+pcp"
# mod4 <- "abs(prc)~tFb+pcp"
# mod5 <- "abs(prc)~tFb"
# mod6 <- "abs(prc)~tJl"
# mod7 <- "abs(prc)~pcp"
# mod8 <- "abs(prc)~1"

mod1 <- "lprc~tFb+tJl+pcp"
mod2 <- "lprc~tFb+tJl"
mod3 <- "lprc~tJl+pcp"
mod4 <- "lprc~tFb+pcp"
mod5 <- "lprc~tFb"
mod6 <- "lprc~tJl"
mod7 <- "lprc~pcp"
mod8 <- "lprc~1"


## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out3 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out3) <- c("TFb","prlprcFb","TJl","prlprcJl")
write.table(out3,file="ext.clim.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


####################################
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+tFb+tJl+pcp"
mod2 <- "t2~t1+tFb+tJl"
mod3 <- "t2~t1+tJl+pcp"
mod4 <- "t2~t1+tFb+pcp"
mod5 <- "t2~t1+tFb"
mod6 <- "t2~t1+tJl"
mod7 <- "t2~t1+pcp"
mod8 <- "t2~t1"
mod9 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out3 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out3) <- c("TFb","prR2Fb","TJl","prR2Jl","pcp","prR2pcp")
write.table(out3,file="ext.clim.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

## climate
## expanders only
## response = transformed % range change (prc)
## model set
# mod1 <- "abs(prc)~tFb+tJl+pcp"
# mod2 <- "abs(prc)~tFb+tJl"
# mod3 <- "abs(prc)~tJl+pcp"
# mod4 <- "abs(prc)~tFb+pcp"
# mod5 <- "abs(prc)~tFb"
# mod6 <- "abs(prc)~tJl"
# mod7 <- "abs(prc)~pcp"
# mod8 <- "abs(prc)~1"

mod1 <- "lprc~tFb+tJl+pcp"
mod2 <- "lprc~tFb+tJl"
mod3 <- "lprc~tJl+pcp"
mod4 <- "lprc~tFb+pcp"
mod5 <- "lprc~tFb"
mod6 <- "lprc~tJl"
mod7 <- "lprc~pcp"
mod8 <- "lprc~1"


## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8)
mod.vec <- c(mod2,mod3,mod4,mod5,mod6,mod7,mod8)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


# export partial residulals
out4 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out4) <- c("TFb","prlprcFb","TJl","prlprcJl")
write.table(out4,file="exp.clim.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


####################################
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+tFb+tJl+pcp"
mod2 <- "t2~t1+tFb+tJl"
mod3 <- "t2~t1+tJl+pcp"
mod4 <- "t2~t1+tFb+pcp"
mod5 <- "t2~t1+tFb"
mod6 <- "t2~t1+tJl"
mod7 <- "t2~t1+pcp"
mod8 <- "t2~t1"
mod9 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Ann pcp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out3 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]),topfit1$model[,4],as.vector(residuals(topfit1,type="partial")[,3]))
colnames(out3) <- c("TFb","prlR2Fb","TJl","prR2Jl","pcp","prR2pcp")
write.table(out3,file="ext.clim.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


## distribution
## all species
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont"
mod2 <- "lprc~t1"
mod3 <- "lprc~lcont"
mod4 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove outliers
dat.ext.nw2 <- dat.ext.nw[-c(39),]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


#########################################
## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lcont"
mod2 <- "t2~t1"
mod3 <- "t2~lcont"
mod4 <- "t2~1"

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
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="log"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="contagion",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

#########################################
## response = t2 (t1 as control variable)
## clumpiness instead of contagion
## model set
mod1 <- "t2~t1+lclump"
mod2 <- "t2~t1"
mod3 <- "t2~lclump"
mod4 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove negative 'clump' values
datwk.nw2 <- subset(datwk.nw, clump > 0)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## decliners only
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont"
mod2 <- "lprc~t1"
mod3 <- "lprc~lcont"
mod4 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove outliers
dat.ext.nw2 <- dat.ext.nw[-c(39),]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out5 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out5) <- c("t1","prlprct1","lcont","prlprclcont")
write.table(out5,file="ext.rge.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


## response = t2 (t1 as control variable)
## clumpiness instead of contagion
## model set
mod1 <- "t2~t1+lclump"
mod2 <- "t2~t1"
mod3 <- "t2~lclump"
mod4 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove negative 'clump' values
dat.ext.nw2 <- subset(dat.ext.nw, clump > 0)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.ext.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out5 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out5) <- c("t1","prR2t1","clump","prlR2clump")
write.table(out5,file="ext.rge.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


## expanders only
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont"
mod2 <- "lprc~t1"
mod3 <- "lprc~lcont"
mod4 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove outliers
#dat.ext.nw2 <- dat.ext.nw[-c(39),]

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
		assign(paste("fit",i,sep=""), fit)
		mod.list[[i]] <- fit
    	print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out6 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out6) <- c("t1","prlprct1","lcont","prlprclcont")
write.table(out6,file="exp.rge.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)


## response = t2 (t1 as control variable)
## clumpiness instead of contagion
## model set
mod1 <- "t2~t1+lclump"
mod2 <- "t2~t1"
mod3 <- "t2~lclump"
mod4 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

## remove negative 'clump' values
dat.exp.nw2 <- subset(dat.exp.nw, clump > 0)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw2, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(1,2),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="range T1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out6 <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,3],as.vector(residuals(topfit1,type="partial")[,2]))
colnames(out6) <- c("t1","prR2t1","lcont","prR2clump")
write.table(out6,file="exp.rge.out.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
## Combined
## all species
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass"
mod2 <- "lprc~t1+lcont+lLong+lMass"
mod3 <- "lprc~t1+lcont+tFb+tJl"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass"
mod5 <- "lprc~lcont+lLong+lMass"
mod6 <- "lprc~lcont+tFb+tJl"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass"
mod8 <- "lprc~t1+lLong+lMass"
mod9 <- "lprc~t1+tFb+tJl"
mod10 <- "lprc~t1+lcont"
mod11 <- "lprc~t1"
mod12 <- "lprc~lcont"
mod13 <- "lprc~lLong+lMass"
mod14 <- "lprc~tFb+tJl"
mod15 <- "lprc~lMass"
mod16 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=6,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lcont","tFb","tJl","lLong","lMass")

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,3),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp"
mod10 <- "t2~t1+lLong+lMass+lClutch"
mod11 <- "t2~t1+tFb+tJl+pcp"
mod12 <- "t2~t1+lclump"
mod13 <- "t2~t1+tFb"
mod14 <- "t2~t1+tJl"
mod15 <- "t2~t1+pcp"
mod16 <- "t2~t1+lLong"
mod17 <- "t2~t1+lMass"
mod18 <- "t2~t1+lClutch"
mod19 <- "t2~t1"
mod20 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

## remove neg clumpiness
datwk.nw2 <- subset(datwk.nw, clump > 0)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=8,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lclump","tFb","tJl","pcp","lLong","lMass","lClutch")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(3,3))
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=datwk.nw2, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Precip",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot7 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot8 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample



## All species, but include 'direction' for all species to differentiate expanders & contractors
## with interactions
## model set
mod1 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass"
mod2 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass"
mod3 <- "lprc~t1+lcont+lLong+lMass"
mod4 <- "lprc~t1+lcont+tFb+tJl"
mod5 <- "lprc~lcont+tFb+tJl+lLong+lMass"
mod6 <- "lprc~lcont+lLong+lMass"
mod7 <- "lprc~lcont+tFb+tJl"
mod8 <- "lprc~t1+tFb+tJl+lLong+lMass"
mod9 <- "lprc~t1+lLong+lMass"
mod10 <- "lprc~t1+tFb+tJl"
mod11 <- "lprc~t1+lcont"
mod12 <- "lprc~t1"
mod13 <- "lprc~lcont"
mod14 <- "lprc~lLong+lMass"
mod15 <- "lprc~tFb+tJl"
mod16 <- "lprc~lMass"
mod17 <- "lprc~drct+t1+lcont+lLong+lMass"
mod18 <- "lprc~drct+t1+lcont+tFb+tJl"
mod19 <- "lprc~drct+lcont+tFb+tJl+lLong+lMass"
mod20 <- "lprc~drct+lcont+lLong+lMass"
mod21 <- "lprc~drct+lcont+tFb+tJl"
mod22 <- "lprc~drct+t1+tFb+tJl+lLong+lMass"
mod23 <- "lprc~drct+t1+lLong+lMass"
mod24 <- "lprc~drct+t1+tFb+tJl"
mod25 <- "lprc~drct+t1+lcont"
mod26 <- "lprc~drct+t1"
mod27 <- "lprc~drct+lcont"
mod28 <- "lprc~drct+lLong+lMass"
mod29 <- "lprc~drct+tFb+tJl"
mod30 <- "lprc~drct+lMass"
mod31 <- "lprc~drct"
mod32 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass+t1*drct"
mod33 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass+lcont*drct"
mod34 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)


## response = t2 (t1 as control variable)
## including direction & interactions
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp"
mod10 <- "t2~t1+lLong+lMass+lClutch"
mod11 <- "t2~t1+tFb+tJl+pcp"
mod12 <- "t2~t1+lclump"
mod13 <- "t2~t1+tFb"
mod14 <- "t2~t1+tJl"
mod15 <- "t2~t1+pcp"
mod16 <- "t2~t1+lLong"
mod17 <- "t2~t1+lMass"
mod18 <- "t2~t1+lClutch"
mod19 <- "t2~t1"
mod20 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod21 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lMass+lClutch+t1*drct"
mod22 <- "t2~t1+drct+tFb+tJl+pcp+lLong+lMass+lClutch"
mod23 <- "t2~t1+drct+lclump+tJl+pcp+lLong+lMass+lClutch"
mod24 <- "t2~t1+drct+lclump+tFb+tJl+lLong+lMass+lClutch"
mod25 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lMass+lClutch"
mod26 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lClutch"
mod27 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lMass"
mod28 <- "t2~t1+drct+lclump+lLong+lMass+lClutch"
mod29 <- "t2~t1+drct+lclump+tFb+tJl+pcp"
mod30 <- "t2~t1+drct+lLong+lMass+lClutch"
mod31 <- "t2~t1+drct+tFb+tJl+pcp"
mod32 <- "t2~t1+drct+lclump"
mod33 <- "t2~t1+drct+tFb"
mod34 <- "t2~t1+drct+tJl"
mod35 <- "t2~t1+drct+pcp"
mod36 <- "t2~t1+drct+lLong"
mod37 <- "t2~t1+drct+lMass"
mod38 <- "t2~t1+drct+lClutch"
mod39 <- "t2~t1+drct"
mod40 <- "t2~drct"
mod41 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40,mod41)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0,Modnum)
mod.list <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw, na.action=na.omit)
plot(fit)

n.sample <- topfit1$df.null+1
n.sample




## standardised coefficients
## model set
mod1 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass"
mod2 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass"
mod3 <- "lprc~t1+lcont+lLong+lMass"
mod4 <- "lprc~t1+lcont+tFb+tJl"
mod5 <- "lprc~lcont+tFb+tJl+lLong+lMass"
mod6 <- "lprc~lcont+lLong+lMass"
mod7 <- "lprc~lcont+tFb+tJl"
mod8 <- "lprc~t1+tFb+tJl+lLong+lMass"
mod9 <- "lprc~t1+lLong+lMass"
mod10 <- "lprc~t1+tFb+tJl"
mod11 <- "lprc~t1+lcont"
mod12 <- "lprc~t1"
mod13 <- "lprc~lcont"
mod14 <- "lprc~lLong+lMass"
mod15 <- "lprc~tFb+tJl"
mod16 <- "lprc~lMass"
mod17 <- "lprc~drct+t1+lcont+lLong+lMass"
mod18 <- "lprc~drct+t1+lcont+tFb+tJl"
mod19 <- "lprc~drct+lcont+tFb+tJl+lLong+lMass"
mod20 <- "lprc~drct+lcont+lLong+lMass"
mod21 <- "lprc~drct+lcont+tFb+tJl"
mod22 <- "lprc~drct+t1+tFb+tJl+lLong+lMass"
mod23 <- "lprc~drct+t1+lLong+lMass"
mod24 <- "lprc~drct+t1+tFb+tJl"
mod25 <- "lprc~drct+t1+lcont"
mod26 <- "lprc~drct+t1"
mod27 <- "lprc~drct+lcont"
mod28 <- "lprc~drct+lLong+lMass"
mod29 <- "lprc~drct+tFb+tJl"
mod30 <- "lprc~drct+lMass"
mod31 <- "lprc~drct"
#mod32 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass+t1*drct"
#mod33 <- "lprc~drct+t1+lcont+tFb+tJl+lLong+lMass+lcont*drct"
mod32 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=7,nrow=Modnum)
colnames(coeff.st.mat) <- c("drct","t1","lcont","tFb","tJl","lLong","lMass")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
plot(fit)


## response = t2 (t1 as control variable)
## standardised coefficients
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp"
mod10 <- "t2~t1+lLong+lMass+lClutch"
mod11 <- "t2~t1+tFb+tJl+pcp"
mod12 <- "t2~t1+lclump"
mod13 <- "t2~t1+tFb"
mod14 <- "t2~t1+tJl"
mod15 <- "t2~t1+pcp"
mod16 <- "t2~t1+lLong"
mod17 <- "t2~t1+lMass"
mod18 <- "t2~t1+lClutch"
mod19 <- "t2~t1"
mod20 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod21 <- "t2~t1+drct+tFb+tJl+pcp+lLong+lMass+lClutch"
mod22 <- "t2~t1+drct+lclump+tJl+pcp+lLong+lMass+lClutch"
mod23 <- "t2~t1+drct+lclump+tFb+tJl+lLong+lMass+lClutch"
mod24 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lMass+lClutch"
mod25 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lClutch"
mod26 <- "t2~t1+drct+lclump+tFb+tJl+pcp+lLong+lMass"
mod27 <- "t2~t1+drct+lclump+lLong+lMass+lClutch"
mod28 <- "t2~t1+drct+lclump+tFb+tJl+pcp"
mod29 <- "t2~t1+drct+lLong+lMass+lClutch"
mod30 <- "t2~t1+drct+tFb+tJl+pcp"
mod31 <- "t2~t1+drct+lclump"
mod32 <- "t2~t1+drct+tFb"
mod33 <- "t2~t1+drct+tJl"
mod34 <- "t2~t1+drct+pcp"
mod35 <- "t2~t1+drct+lLong"
mod36 <- "t2~t1+drct+lMass"
mod37 <- "t2~t1+drct+lClutch"
mod38 <- "t2~t1+drct"
mod39 <- "t2~drct"
mod40 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21,mod22,mod23,mod24,mod25,mod26,mod27,mod28,mod29,mod30,mod31,mod32,mod33,mod34,mod35,mod36,mod37,mod38,mod39,mod40)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=9,nrow=Modnum)
colnames(coeff.st.mat) <- c("drct","t1","lclump","tFb","tJl","pcp","lLong","lMass","lClutch")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=datwk.nw, na.action=na.omit)
plot(fit)



## decliners only
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass"
mod2 <- "lprc~t1+lcont+lLong+lMass"
mod3 <- "lprc~t1+lcont+tFb+tJl"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass"
mod5 <- "lprc~lcont+lLong+lMass"
mod6 <- "lprc~lcont+tFb+tJl"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass"
mod8 <- "lprc~t1+lLong+lMass"
mod9 <- "lprc~t1+tFb+tJl"
mod10 <- "lprc~t1+lcont"
mod11 <- "lprc~t1"
mod12 <- "lprc~lcont"
mod13 <- "lprc~lLong+lMass"
mod14 <- "lprc~tFb+tJl"
mod15 <- "lprc~lMass"
mod16 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=6,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lcont","tFb","tJl","lLong","lMass")

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,3),pty="s")
topfit1 <- glm(as.formula(mod.vec[8]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample

# export partial residulals
out.ext.comb <- data.frame(topfit1$model[,2],as.vector(residuals(topfit1,type="partial")[,1]),topfit1$model[,6],as.vector(residuals(topfit1,type="partial")[,5]),topfit1$model[,7],as.vector(residuals(topfit1,type="partial")[,6]))
colnames(out.ext.comb) <- c("t1","prlprct1","lLong","prlprcL","lMass","prlprcM")
write.table(out.ext.comb,file="out.ext.comb.csv",sep=",",dec = ".", row.names = F,col.names = TRUE)



## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp"
mod10 <- "t2~t1+lLong+lMass+lClutch"
mod11 <- "t2~t1+tFb+tJl+pcp"
mod12 <- "t2~t1+lclump"
mod13 <- "t2~t1+tFb"
mod14 <- "t2~t1+tJl"
mod15 <- "t2~t1+pcp"
mod16 <- "t2~t1+lLong"
mod17 <- "t2~t1+lMass"
mod18 <- "t2~t1+lClutch"
mod19 <- "t2~t1"
mod20 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=8,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lclump","tFb","tJl","pcp","lLong","lMass","lClutch")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(3,3))
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Precip",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot7 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot8 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample




## expanders only
## response = transformed % range change (prc)
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass"
mod2 <- "lprc~t1+lcont+lLong+lMass"
mod3 <- "lprc~t1+lcont+tFb+tJl"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass"
mod5 <- "lprc~lcont+lLong+lMass"
mod6 <- "lprc~lcont+tFb+tJl"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass"
mod8 <- "lprc~t1+lLong+lMass"
mod9 <- "lprc~t1+tFb+tJl"
mod10 <- "lprc~t1+lcont"
mod11 <- "lprc~t1"
mod12 <- "lprc~lcont"
mod13 <- "lprc~lLong+lMass"
mod14 <- "lprc~tFb+tJl"
mod15 <- "lprc~lMass"
mod16 <- "lprc~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=6,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lcont","tFb","tJl","lLong","lMass")

    for(i in 1:Modnum) {
		fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(2,3),pty="s")
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="lcont",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Contagion",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Range T1",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="prc partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample


## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp"
mod10 <- "t2~t1+lLong+lMass+lClutch"
mod11 <- "t2~t1+tFb+tJl+pcp"
mod12 <- "t2~t1+lclump"
mod13 <- "t2~t1+tFb"
mod14 <- "t2~t1+tJl"
mod15 <- "t2~t1+pcp"
mod16 <- "t2~t1+lLong"
mod17 <- "t2~t1+lMass"
mod18 <- "t2~t1+lClutch"
mod19 <- "t2~t1"
mod20 <- "t2~1"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- terml <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)
coeff.st.mat <- matrix(0,ncol=8,nrow=Modnum)
colnames(coeff.st.mat) <- c("t1","lclump","tFb","tJl","pcp","lLong","lMass","lClutch")

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
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

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
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
fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot(fit)

## Partial residual plot
## select saturated model
par(mfrow=c(3,3))
topfit1 <- glm(as.formula(mod.vec[1]),family=gaussian(link="sqrt"), data=dat.exp.nw, na.action=na.omit)
plot1 <- termplot(topfit1,terms="t1",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="R1",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot2 <- termplot(topfit1,terms="lclump",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clumpiness",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot3 <- termplot(topfit1,terms="tFb",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Feb Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot4 <- termplot(topfit1,terms="tJl",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Jul Temp",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot5 <- termplot(topfit1,terms="pcp",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Precip",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot6 <- termplot(topfit1,terms="lMass",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Mass",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot7 <- termplot(topfit1,terms="lLong",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Longevity",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
plot8 <- termplot(topfit1,terms="lClutch",partial.resid=T,rug=T,se=T,col.term="black",col.res="black",col.se="black",pch=19,xlab="Clutch size",ylab="R2 partial resid",cex.axis=1.5,cex.lab=1.5)
par(mfrow=c(1,1))

n.sample <- topfit1$df.null+1
n.sample



#################################
### GLMMs for final models ######
#################################


## contractors only
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod2 <- "lprc~t1+lcont+lLong+lMass+(1|Order)"
mod3 <- "lprc~t1+lcont+tFb+tJl+(1|Order)"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod5 <- "lprc~lcont+lLong+lMass+(1|Order)"
mod6 <- "lprc~lcont+tFb+tJl+(1|Order)"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass+(1|Order)"
mod8 <- "lprc~t1+lLong+lMass+(1|Order)"
mod9 <- "lprc~t1+tFb+tJl+(1|Order)"
mod10 <- "lprc~t1+lcont+(1|Order)"
mod11 <- "lprc~t1+(1|Order)"
mod12 <- "lprc~lcont+(1|Order)"
mod13 <- "lprc~lLong+lMass+(1|Order)"
mod14 <- "lprc~tFb+tJl+(1|Order)"
mod15 <- "lprc~lMass+(1|Order)"
mod16 <- "lprc~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table


## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod3 <- "t2~t1+lclump+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod4 <- "t2~t1+lclump+tFb+tJl+lLong+lMass+lClutch+(1|Order)"
mod5 <- "t2~t1+lclump+tFb+tJl+pcp+lMass+lClutch+(1|Order)"
mod6 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lClutch+(1|Order)"
mod7 <- "t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+(1|Order)"
mod8 <- "t2~t1+lclump+lLong+lMass+lClutch+(1|Order)"
mod9 <- "t2~t1+lclump+tFb+tJl+pcp+(1|Order)"
mod10 <- "t2~t1+lLong+lMass+lClutch+(1|Order)"
mod11 <- "t2~t1+tFb+tJl+pcp+(1|Order)"
mod12 <- "t2~t1+lclump+(1|Order)"
mod13 <- "t2~t1+tFb+(1|Order)"
mod14 <- "t2~t1+tJl+(1|Order)"
mod15 <- "t2~t1+pcp+(1|Order)"
mod16 <- "t2~t1+lLong+(1|Order)"
mod17 <- "t2~t1+lMass+(1|Order)"
mod18 <- "t2~t1+lClutch+(1|Order)"
mod19 <- "t2~t1+(1|Order)"
mod20 <- "t2~1+(1|Order)"
mod21 <- "t2~1+(1|dummy)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20,mod21)

dat.ext.nw$dummy <- rep(1,nrow(dat.ext.nw))

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table


best <- lmer(as.formula(mod.vec[1]),family=gaussian(link="identity"), data=dat.ext.nw, na.action=na.omit)
print(best)
hist(resid(best), breaks=300)
qqnorm(resid(best))
qqline(resid(best), col="blue")
coef(best)$Order


{ ## 10-fold cross validation
  n <- dim(dat.ext.nw)[1]
  vec.error <- rep(0,100)
  vec.r2 <- rep(0,100)
  plot(range(dat.ext.nw$t2), range(dat.ext.nw$t2), type="n")
  
  for (i in 1:100) {
    random=srswor(round(n/7),n)
    bootstrap=subset(cbind(random,dat.ext.nw), random==0)
    new.data=subset(cbind(random,dat.ext.nw), random==1)
    mod <- lmer(t2~t1+lclump+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order), family=gaussian(link=log), data=bootstrap)
    #mod <- lmer(S ~ MIG_H + MIG_S + MIG_V + MMI_H + MMI_S + MMI_V + (1|Site), family=poisson, data=bootstrap)
    fixed <- fixef(mod)["t1"] * new.data$t1 + fixef(mod)["lclump"] * new.data$lclump + fixef(mod)["tFb"] * new.data$tFb + fixef(mod)["tJl"] * new.data$tJl + fixef(mod)["pcp"] * new.data$pcp + fixef(mod)["lLong"] * new.data$lLong  + fixef(mod)["lMass"] * new.data$lMass + fixef(mod)["lClutch"] * new.data$lClutch
    rand <-  coef(mod)$Order[new.data$Order,"(Intercept)"]
    pred <- exp(fixed + rand)
    points(new.data$t2, pred)
    error <- abs(pred-new.data$t2)*100/new.data$t2
    vec.error[i] <- mean(error)
    vec.r2[i] <- summary(lm(pred~new.data$t2))$r.squared
  }
  mean(vec.error,na.rm=T)
  mean(vec.r2,na.rm=T)
}



## expanders only
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod2 <- "lprc~t1+lcont+lLong+lMass+(1|Order)"
mod3 <- "lprc~t1+lcont+tFb+tJl+(1|Order)"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod5 <- "lprc~lcont+lLong+lMass+(1|Order)"
mod6 <- "lprc~lcont+tFb+tJl+(1|Order)"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass+(1|Order)"
mod8 <- "lprc~t1+lLong+lMass+(1|Order)"
mod9 <- "lprc~t1+tFb+tJl+(1|Order)"
mod10 <- "lprc~t1+lcont+(1|Order)"
mod11 <- "lprc~t1+(1|Order)"
mod12 <- "lprc~lcont+(1|Order)"
mod13 <- "lprc~lLong+lMass+(1|Order)"
mod14 <- "lprc~tFb+tJl+(1|Order)"
mod15 <- "lprc~lMass+(1|Order)"
mod16 <- "lprc~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table


## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod3 <- "t2~t1+clump+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod4 <- "t2~t1+clump+tFb+tJl+lLong+lMass+lClutch+(1|Order)"
mod5 <- "t2~t1+clump+tFb+tJl+pcp+lMass+lClutch+(1|Order)"
mod6 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lClutch+(1|Order)"
mod7 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lMass+(1|Order)"
mod8 <- "t2~t1+clump+lLong+lMass+lClutch+(1|Order)"
mod9 <- "t2~t1+clump+tFb+tJl+pcp+(1|Order)"
mod10 <- "t2~t1+lLong+lMass+lClutch+(1|Order)"
mod11 <- "t2~t1+tFb+tJl+pcp+(1|Order)"
mod12 <- "t2~t1+clump+(1|Order)"
mod13 <- "t2~t1+tFb+(1|Order)"
mod14 <- "t2~t1+tJl+(1|Order)"
mod15 <- "t2~t1+pcp+(1|Order)"
mod16 <- "t2~t1+lLong+(1|Order)"
mod17 <- "t2~t1+lMass+(1|Order)"
mod18 <- "t2~t1+lClutch+(1|Order)"
mod19 <- "t2~t1+(1|Order)"
mod20 <- "t2~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=dat.exp.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table



## all species
## model set
mod1 <- "lprc~t1+lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod2 <- "lprc~t1+lcont+lLong+lMass+(1|Order)"
mod3 <- "lprc~t1+lcont+tFb+tJl+(1|Order)"
mod4 <- "lprc~lcont+tFb+tJl+lLong+lMass+(1|Order)"
mod5 <- "lprc~lcont+lLong+lMass+(1|Order)"
mod6 <- "lprc~lcont+tFb+tJl+(1|Order)"
mod7 <- "lprc~t1+tFb+tJl+lLong+lMass+(1|Order)"
mod8 <- "lprc~t1+lLong+lMass+(1|Order)"
mod9 <- "lprc~t1+tFb+tJl+(1|Order)"
mod10 <- "lprc~t1+lcont+(1|Order)"
mod11 <- "lprc~t1+(1|Order)"
mod12 <- "lprc~lcont+(1|Order)"
mod13 <- "lprc~lLong+lMass+(1|Order)"
mod14 <- "lprc~tFb+tJl+(1|Order)"
mod15 <- "lprc~lMass+(1|Order)"
mod16 <- "lprc~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table


## response = t2 (t1 as control variable)
## model set
mod1 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod2 <- "t2~t1+tFb+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod3 <- "t2~t1+clump+tJl+pcp+lLong+lMass+lClutch+(1|Order)"
mod4 <- "t2~t1+clump+tFb+tJl+lLong+lMass+lClutch+(1|Order)"
mod5 <- "t2~t1+clump+tFb+tJl+pcp+lMass+lClutch+(1|Order)"
mod6 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lClutch+(1|Order)"
mod7 <- "t2~t1+clump+tFb+tJl+pcp+lLong+lMass+(1|Order)"
mod8 <- "t2~t1+clump+lLong+lMass+lClutch+(1|Order)"
mod9 <- "t2~t1+clump+tFb+tJl+pcp+(1|Order)"
mod10 <- "t2~t1+lLong+lMass+lClutch+(1|Order)"
mod11 <- "t2~t1+tFb+tJl+pcp+(1|Order)"
mod12 <- "t2~t1+clump+(1|Order)"
mod13 <- "t2~t1+tFb+(1|Order)"
mod14 <- "t2~t1+tJl+(1|Order)"
mod15 <- "t2~t1+pcp+(1|Order)"
mod16 <- "t2~t1+lLong+(1|Order)"
mod17 <- "t2~t1+lMass+(1|Order)"
mod18 <- "t2~t1+lClutch+(1|Order)"
mod19 <- "t2~t1+(1|Order)"
mod20 <- "t2~1+(1|Order)"

## Make model vector
mod.vec <- c(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8,mod9,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17,mod18,mod19,mod20)

## Define n.mod
n.mod <- length(mod.vec)

# Model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- pc.dev.vec <- k.vec <- rep(0, Modnum)
mod.list <- 0
mod.num <- seq(1, Modnum, 1)

for (i in 1:Modnum) {
  fit <- glmer(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=datwk.nw, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list <- c(mod.list,fit)
  print(i)
}

mod.list <- as.list(mod.list[-1])
sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,4],decreasing=F),1:9]
summary.table
