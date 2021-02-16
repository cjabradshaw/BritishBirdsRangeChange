######################################################################################################
##
## PGLS
##
######################################################################################################

library(caper)
library(MuMIn)

## Data
datwk.nw.phy.sub <- na.omit(datwk.nw.phy[, c('t1', 't2', 'prop.t1', 'spp.phylo', 'prop.t2', 'logit.t1', 'logit.t2', 'clump', 'feb_tmn_Mean', 'jul_tmn_Mean', 'an_sum_pre_Mean', 'lnMass', 'Order', 'NatL', 'Clutch', 'natDisp')])
names(datwk.nw.phy.sub)[c(9, 11)] <- c("tFb", "pcp")

## Phylogeny
bbp.vc.sub <- comparative.data(BritishBirds.tree, datwk.nw.phy.sub, 'spp.phylo', vcv = TRUE, vcv.dim = 3, na.omit = FALSE)
print(bbp.vc.sub)
head(bbp.vc.sub$data)

##----------------------------------------------------------------------------------------------------
## Null model
##----------------------------------------------------------------------------------------------------

pgls.null <- pgls(logit.t2 ~ logit.t1, data = bbp.vc.sub, lambda = 'ML')
summary(pgls.null)

##----------------------------------------------------------------------------------------------------
## 1. Mass, Feb Temp & Clumpiness
##----------------------------------------------------------------------------------------------------

pgls.m1 <- pgls(logit.t2 ~ logit.t1 + ec/(lnMass + tFb + clump) + 0, data = bbp.vc.sub)
summary(pgls.m1)
pgls.m2 <- pgls(logit.t2 ~ logit.t1 + ec/(lnMass) + tFb + clump + 0, data = bbp.vc.sub)
summary(pgls.m2)
pgls.m3 <- pgls(logit.t2 ~ logit.t1 + ec/(clump) + lnMass + tFb + 0, data = bbp.vc.sub)
summary(pgls.m3)
pgls.m4 <- pgls(logit.t2 ~ logit.t1 + ec/(tFb) + lnMass + clump + 0, data = bbp.vc.sub)
summary(pgls.m4)
pgls.m5 <- pgls(logit.t2 ~ logit.t1 + ec + lnMass + tFb + clump + 0, data = bbp.vc.sub)
summary(pgls.m5)
pgls.m6 <- pgls(logit.t2 ~ logit.t1 + ec + lnMass + tFb + 0, data = bbp.vc.sub)
summary(pgls.m6)
pgls.m7 <- pgls(logit.t2 ~ logit.t1 + ec/lnMass + 0, data = bbp.vc.sub)
summary(pgls.m7)
pgls.m8 <- pgls(logit.t2 ~ logit.t1 + ec/tFb + 0, data = bbp.vc.sub)
summary(pgls.m8)
pgls.m9 <- pgls(logit.t2 ~ logit.t1 + ec + tFb + clump + 0, data = bbp.vc.sub)
summary(pgls.m9)
pgls.m10 <- pgls(logit.t2 ~ logit.t1 + ec + tFb + 0, data = bbp.vc.sub)
summary(pgls.m10)
pgls.m11 <- pgls(logit.t2 ~ logit.t1 + ec + lnMass + clump + 0, data = bbp.vc.sub)
summary(pgls.m11)
pgls.m12 <- pgls(logit.t2 ~ logit.t1 + ec + lnMass + 0, data = bbp.vc.sub)
summary(pgls.m12)
pgls.m13 <- pgls(logit.t2 ~ logit.t1 + ec/clump + 0, data = bbp.vc.sub)
summary(pgls.m13)
pgls.m14 <- pgls(logit.t2 ~ logit.t1 + ec + clump + 0, data = bbp.vc.sub)
summary(pgls.m14)
pgls.m15 <- pgls(logit.t2 ~ logit.t1 + ec + 0, data = bbp.vc.sub)
summary(pgls.m15)
pgls.m16 <- pgls(logit.t2 ~ logit.t1, data = bbp.vc.sub)
summary(pgls.m16)
pgls.m17 <- pgls(logit.t2 ~ 1, data = bbp.vc.sub)
summary(pgls.m17)

## Model selection
(zs1 <- model.sel(pgls.m1, pgls.m2, pgls.m3, pgls.m4, pgls.m5, pgls.m6, pgls.m7, pgls.m8, pgls.m9, pgls.m10, pgls.m11, pgls.m12, pgls.m13, pgls.m14, pgls.m15, pgls.m16, pgls.m17))


