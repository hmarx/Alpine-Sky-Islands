#####################################################################################################################
############# Phylogenetic diversity within alpine summits (alpha) ################################################## 
############# Static Null Model ##################################################################################### 
############# Hannah E. Marx, 16 Feb 2016 ########################################################################### 
#####################################################################################################################

source("analysisSkyIsl.R")
source("R/randomizeSourcePoolNull_v1.5.R")
source("R/alps.damocles.chooseClade.R")
source("R/pruneSpeciesPoolsPez.R")

###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool: Regional Ecrins NP ##################### 
###################################################################################### 

## Random resample from phylogeny pool (==Ecrins NP), equal probability random draw from phylogeny pool
ecrins.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
#comm.sesmpd.samplenull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)
write.csv(ecrins.sesmpd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/ecrins.sesmpd.phylonull.csv")

ecrins.sesmntd.phylonull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
write.csv(ecrins.sesmntd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/ecrins.sesmntd.phylonull.csv")


###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool: All Summits ############################
###################################################################################### 

## Source pool = summits, equal probability random draw from phylogeny (pruned to summits)
summits.sesmpd.phylonull <- ses.mpd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
summits.sesmntd.phylonull <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
write.csv(summits.sesmpd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.phylonull.csv")
write.csv(summits.sesmntd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.phylonull.csv")

## Source pool = summits, "abundance" weighted random draw from phylogeny (pruned to summits)
summits.sesmpd.phylonull.abun <- ses.mpd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = TRUE, runs = 999)
summits.sesmntd.phylonull.abun <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = TRUE, runs = 999)
write.csv(summits.sesmpd.phylonull.abun, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.phylonull.abun.csv")
write.csv(summits.sesmntd.phylonull.abun, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.phylonull.abun.csv")

## source pool = summits (persistent + under ice?), weighted for abundance in summits
summits.sesmpd.sourceSummits <- ses.mpd.sourcePool(dist=cophenetic.phylo(pezAlpes.summits$phy), com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)
summits.sesmntd.sourceSummits <- ses.mntd.sourcePool(dist=cophenetic.phylo(pezAlpes.summits$phy), com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)
write.csv(summits.sesmpd.sourceSummits, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmpd.sourceSummits.csv")
write.csv(summits.sesmntd.sourceSummits, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/summits.sesmntd.sourceSummits.csv")


###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool : Persistent through LGM ################
###################################################################################### 

## Source pool = Persistent, equal probability random draw from phylogeny (pruned to persistent species)
persistent.sesmpd.phylonull <- ses.mpd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
persistent.sesmntd.phylonull <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
write.csv(persistent.sesmpd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.phylonull.csv")
write.csv(persistent.sesmntd.phylonull, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.phylonull.csv")

## Source pool = Persistent, "abundance" weighted draw from phylogeny (pruned to persistent species)
persistent.sesmpd.phylonull.abun <- ses.mpd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = TRUE, runs = 999)
persistent.sesmntd.phylonull.abun <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = TRUE, runs = 999)
write.csv(persistent.sesmpd.phylonull.abun, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.phylonull.abun.csv")
write.csv(persistent.sesmntd.phylonull.abun, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.phylonull.abun.csv")

## source pool = persistent above glacier through LGM, weighted for abundance in persistent pool
persistent.sesmpd.sourcePersis <- ses.mpd.sourcePool(dist=cophenetic.phylo(pezAlpes.persistent$phy), com = pezAlpes.persistent$comm, sourcePool = "Persistent", N =999)
persistent.sesmntd.sourcePersis <- ses.mntd.sourcePool(dist=cophenetic.phylo(pezAlpes.persistent$phy), com = pezAlpes.persistent$comm, sourcePool = "Persistent", N =999)
write.csv(persistent.sesmpd.sourcePersis, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.sourcePersis.csv")
write.csv(persistent.sesmntd.sourcePersis, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.sourcePersis.csv")
  
## source pool = summits (persistent + under ice?), weighted for abundance in summits
persistent.sesmpd.sourceSummits <- ses.mpd.sourcePool(dist=cophenetic.phylo(pezAlpes.persistent$phy), com = pezAlpes.persistent$comm, sourcePool = "Summits", N =999)
persistent.sesmntd.sourceSummits <- ses.mntd.sourcePool(dist=cophenetic.phylo(pezAlpes.persistent$phy), com = pezAlpes.persistent$comm, sourcePool = "Summits", N =999)
write.csv(persistent.sesmpd.sourceSummits, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmpd.sourceSummits.csv")
write.csv(persistent.sesmntd.sourceSummits, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/sourcePools/persistent.sesmntd.sourceSummits.csv")


###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: Regional Ecrins NP #####################
###################################################################################### 

# Convert to pres/abs
tmp <- pezAlpes$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa.all <- t(tmp)

# Pruend phylogeny and community matrix for each clade 
ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa.all)
aster.sesmpd.phylonull <- ses.mpd(t(alps.damocles.aster$data), cophenetic.phylo(alps.damocles.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.phylonull <- ses.mntd(t(alps.damocles.aster$data), cophenetic.phylo(alps.damocles.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa.all)
poales.sesmpd.phylonull <- ses.mpd(t(alps.damocles.poales$data), cophenetic.phylo(alps.damocles.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.phylonull <- ses.mntd(t(alps.damocles.poales$data), cophenetic.phylo(alps.damocles.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.rosales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Rosales")
pruned.tree.alps.rosales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.rosales])
alps.damocles.rosales <- treedata(pruned.tree.alps.rosales, alps.sites.pa.all)
rosales.sesmpd.phylonull <- ses.mpd(t(alps.damocles.rosales$data), cophenetic.phylo(alps.damocles.rosales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.phylonull <- ses.mntd(t(alps.damocles.rosales$data), cophenetic.phylo(alps.damocles.rosales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa.all)
lamiales.sesmpd.phylonull <- ses.mpd(t(alps.damocles.lamiales$data), cophenetic.phylo(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.phylonull <- ses.mntd(t(alps.damocles.lamiales$data), cophenetic.phylo(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.alps.caryo <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryo <- treedata(pruned.tree.alps.caryo, alps.sites.pa.all)
caryo.sesmpd.phylonull <- ses.mpd(t(alps.damocles.caryo$data), cophenetic.phylo(alps.damocles.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.phylonull <- ses.mntd(t(alps.damocles.caryo$data), cophenetic.phylo(alps.damocles.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

##### Combine all output of randomization from phylogeny pool
phylogeny.poolSESmntd <- rbind((ecrins.sesmntd.phylonull %>%  mutate(summits = rownames(ecrins.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(ecrins.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | ecrins.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmntd.phylonull %>%  mutate(summits = rownames(aster.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(aster.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | aster.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                               (poales.sesmntd.phylonull %>%  mutate(summits = rownames(poales.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(poales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | poales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (rosales.sesmntd.phylonull %>%  mutate(summits = rownames(rosales.sesmntd.phylonull), metric = "mntd", clade = "Rosales", sig=ifelse(rosales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | rosales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (lamiales.sesmntd.phylonull %>%  mutate(summits = rownames(lamiales.sesmntd.phylonull), metric = "mntd", clade = "Lamiales", sig=ifelse(lamiales.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | lamiales.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (caryo.sesmntd.phylonull %>%  mutate(summits = rownames(caryo.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryo.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | caryo.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.poolSESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.poolSESmpd <- (rbind((ecrins.sesmpd.phylonull %>%  mutate(summits = rownames(ecrins.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(ecrins.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | ecrins.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmpd.phylonull %>%  mutate(summits = rownames(aster.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(aster.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | aster.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                               (poales.sesmpd.phylonull %>%  mutate(summits = rownames(poales.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(poales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | poales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (rosales.sesmpd.phylonull %>%  mutate(summits = rownames(rosales.sesmpd.phylonull), metric = "mpd", clade = "Rosales", sig=ifelse(rosales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | rosales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (lamiales.sesmpd.phylonull %>%  mutate(summits = rownames(lamiales.sesmpd.phylonull), metric = "mpd", clade = "Lamiales", sig=ifelse(lamiales.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | lamiales.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (caryo.sesmpd.phylonull %>%  mutate(summits = rownames(caryo.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryo.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | caryo.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0)))))
names(phylogeny.poolSESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
phylogeny.poolSES <- rbind(phylogeny.poolSESmntd, phylogeny.poolSESmpd)                               
head(phylogeny.poolSES)     
phylogeny.poolSES <- phylogeny.poolSES[!phylogeny.poolSES$summits == "Ecrins NP",]
phylogeny.poolSES$clade <- factor(phylogeny.poolSES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
phylogeny.poolSES <- cbind(phylogeny.poolSES, pool = rep(x = "Ecrins NP", times = nrow(phylogeny.poolSES)))
write.csv(phylogeny.poolSES, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/phylogeny.pool.SES.csv")


###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: All Summits ############################
###################################################################################### 

# Pres/abs for summits pool
tmp2 <- pezAlpes.summits$comm
tmp2[which(tmp2 != 0)] <- 1
alps.summits.pa.all <- t(tmp2)

alps.damocles.aster.summits <- treedata(pruned.tree.alps.aster, alps.summits.pa.all)
aster.sesmpd.summits <- ses.mpd(t(alps.damocles.aster.summits$data), cophenetic.phylo(alps.damocles.aster.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.summits <- ses.mntd(t(alps.damocles.aster.summits$data), cophenetic.phylo(alps.damocles.aster.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.poales.summits <- treedata(pruned.tree.alps.poales, alps.summits.pa.all)
poales.sesmpd.summits <- ses.mpd(t(alps.damocles.poales.summits$data), cophenetic.phylo(alps.damocles.poales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.summits <- ses.mntd(t(alps.damocles.poales.summits$data), cophenetic.phylo(alps.damocles.poales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.rosales.summits <- treedata(pruned.tree.alps.rosales, alps.summits.pa.all)
rosales.sesmpd.summits <- ses.mpd(t(alps.damocles.rosales.summits$data), cophenetic.phylo(alps.damocles.rosales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.summits <- ses.mntd(t(alps.damocles.rosales.summits$data), cophenetic.phylo(alps.damocles.rosales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.lamiales.summits <- treedata(pruned.tree.alps.lamiales, alps.summits.pa.all)
lamiales.sesmpd.summits <- ses.mpd(t(alps.damocles.lamiales.summits$data), cophenetic.phylo(alps.damocles.lamiales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.summits <- ses.mntd(t(alps.damocles.lamiales.summits$data), cophenetic.phylo(alps.damocles.lamiales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.caryo.summits <- treedata(pruned.tree.alps.caryo, alps.summits.pa.all)
caryo.sesmpd.summits <- ses.mpd(t(alps.damocles.caryo.summits$data), cophenetic.phylo(alps.damocles.caryo.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.summits <- ses.mntd(t(alps.damocles.caryo.summits$data), cophenetic.phylo(alps.damocles.caryo.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

##### Combine all output of randomization from phylogeny pool
summits.poolSESmntd <- rbind((summits.sesmntd.phylonull %>%  mutate(summits = rownames(summits.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(summits.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | summits.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmntd.summits %>%  mutate(summits = rownames(aster.sesmntd.summits), metric = "mntd", clade = "Asterales", sig=ifelse(aster.sesmntd.summits[,"mntd.obs.p"] <= 0.05 | aster.sesmntd.summits[,"mntd.obs.p"] > 0.95, 1,0))), 
                               (poales.sesmntd.summits %>%  mutate(summits = rownames(poales.sesmntd.summits), metric = "mntd", clade = "Poales", sig=ifelse(poales.sesmntd.summits[,"mntd.obs.p"] <= 0.05 | poales.sesmntd.summits[,"mntd.obs.p"] > 0.95, 1,0))),
                               (rosales.sesmntd.summits %>%  mutate(summits = rownames(rosales.sesmntd.summits), metric = "mntd", clade = "Rosales", sig=ifelse(rosales.sesmntd.summits[,"mntd.obs.p"] <= 0.05 | rosales.sesmntd.summits[,"mntd.obs.p"] > 0.95, 1,0))),
                               (lamiales.sesmntd.summits %>%  mutate(summits = rownames(lamiales.sesmntd.summits), metric = "mntd", clade = "Lamiales", sig=ifelse(lamiales.sesmntd.summits[,"mntd.obs.p"] <= 0.05 | lamiales.sesmntd.summits[,"mntd.obs.p"] > 0.95, 1,0))),
                               (caryo.sesmntd.summits %>%  mutate(summits = rownames(caryo.sesmntd.summits), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryo.sesmntd.summits[,"mntd.obs.p"] <= 0.05 | caryo.sesmntd.summits[,"mntd.obs.p"] > 0.95, 1,0))))
names(summits.poolSESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

summits.poolSESmpd <- (rbind((summits.sesmpd.phylonull %>%  mutate(summits = rownames(summits.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(summits.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | summits.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                               (aster.sesmpd.summits %>%  mutate(summits = rownames(aster.sesmpd.summits), metric = "mpd", clade = "Asterales", sig=ifelse(aster.sesmpd.summits[,"mpd.obs.p"] <= 0.05 | aster.sesmpd.summits[,"mpd.obs.p"] > 0.95, 1,0))), 
                               (poales.sesmpd.summits %>%  mutate(summits = rownames(poales.sesmpd.summits), metric = "mpd", clade = "Poales", sig=ifelse(poales.sesmpd.summits[,"mpd.obs.p"] <= 0.05 | poales.sesmpd.summits[,"mpd.obs.p"] > 0.95, 1,0))),
                               (rosales.sesmpd.summits %>%  mutate(summits = rownames(rosales.sesmpd.summits), metric = "mpd", clade = "Rosales", sig=ifelse(rosales.sesmpd.summits[,"mpd.obs.p"] <= 0.05 | rosales.sesmpd.summits[,"mpd.obs.p"] > 0.95, 1,0))),
                               (lamiales.sesmpd.summits %>%  mutate(summits = rownames(lamiales.sesmpd.summits), metric = "mpd", clade = "Lamiales", sig=ifelse(lamiales.sesmpd.summits[,"mpd.obs.p"] <= 0.05 | lamiales.sesmpd.summits[,"mpd.obs.p"] > 0.95, 1,0))),
                               (caryo.sesmpd.summits %>%  mutate(summits = rownames(caryo.sesmpd.summits), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryo.sesmpd.summits[,"mpd.obs.p"] <= 0.05 | caryo.sesmpd.summits[,"mpd.obs.p"] > 0.95, 1,0)))))
names(summits.poolSESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
summit.poolSES <- rbind(summits.poolSESmntd, summits.poolSESmpd)                               
head(summit.poolSES)     
summit.poolSES <- summit.poolSES[!summit.poolSES$summits == "Ecrins NP",]
summit.poolSES$clade <- factor(summit.poolSES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
summit.poolSES <- cbind(summit.poolSES, pool = rep(x = "Summits", times = nrow(summit.poolSES)))
write.csv(summit.poolSES, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/summit.pool.SES.csv")


###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: Persistent through LGM #################
###################################################################################### 

# Pres/abs for persistent pool
tmp3 <- pezAlpes.persistent$comm
tmp3[which(tmp3 != 0)] <- 1
alps.persistent.pa.all <- t(tmp3)

alps.damocles.aster.persistent <- treedata(pruned.tree.alps.aster, alps.persistent.pa.all)
aster.sesmpd.persistent <- ses.mpd(t(alps.damocles.aster.persistent$data), cophenetic.phylo(alps.damocles.aster.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.persistent <- ses.mntd(t(alps.damocles.aster.persistent$data), cophenetic.phylo(alps.damocles.aster.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.poales.persistent <- treedata(pruned.tree.alps.poales, alps.persistent.pa.all)
poales.sesmpd.persistent <- ses.mpd(t(alps.damocles.poales.persistent$data), cophenetic.phylo(alps.damocles.poales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.persistent <- ses.mntd(t(alps.damocles.poales.persistent$data), cophenetic.phylo(alps.damocles.poales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.rosales.persistent <- treedata(pruned.tree.alps.rosales, alps.persistent.pa.all)
rosales.sesmpd.persistent <- ses.mpd(t(alps.damocles.rosales.persistent$data), cophenetic.phylo(alps.damocles.rosales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.persistent <- ses.mntd(t(alps.damocles.rosales.persistent$data), cophenetic.phylo(alps.damocles.rosales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.lamiales.persistent <- treedata(pruned.tree.alps.lamiales, alps.persistent.pa.all)
lamiales.sesmpd.persistent <- ses.mpd(t(alps.damocles.lamiales.persistent$data), cophenetic.phylo(alps.damocles.lamiales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.persistent <- ses.mntd(t(alps.damocles.lamiales.persistent$data), cophenetic.phylo(alps.damocles.lamiales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.damocles.caryo.persistent <- treedata(pruned.tree.alps.caryo, alps.persistent.pa.all)
caryo.sesmpd.persistent <- ses.mpd(t(alps.damocles.caryo.persistent$data), cophenetic.phylo(alps.damocles.caryo.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.persistent <- ses.mntd(t(alps.damocles.caryo.persistent$data), cophenetic.phylo(alps.damocles.caryo.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

##### Combine all output of randomization from phylogeny pool
persistent.poolSESmntd <- rbind((persistent.sesmntd.phylonull %>%  mutate(summits = rownames(persistent.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(persistent.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | persistent.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                             (aster.sesmntd.persistent %>%  mutate(summits = rownames(aster.sesmntd.persistent), metric = "mntd", clade = "Asterales", sig=ifelse(aster.sesmntd.persistent[,"mntd.obs.p"] <= 0.05 | aster.sesmntd.persistent[,"mntd.obs.p"] > 0.95, 1,0))), 
                             (poales.sesmntd.persistent %>%  mutate(summits = rownames(poales.sesmntd.persistent), metric = "mntd", clade = "Poales", sig=ifelse(poales.sesmntd.persistent[,"mntd.obs.p"] <= 0.05 | poales.sesmntd.persistent[,"mntd.obs.p"] > 0.95, 1,0))),
                             (rosales.sesmntd.persistent %>%  mutate(summits = rownames(rosales.sesmntd.persistent), metric = "mntd", clade = "Rosales", sig=ifelse(rosales.sesmntd.persistent[,"mntd.obs.p"] <= 0.05 | rosales.sesmntd.persistent[,"mntd.obs.p"] > 0.95, 1,0))),
                             (lamiales.sesmntd.persistent %>%  mutate(summits = rownames(lamiales.sesmntd.persistent), metric = "mntd", clade = "Lamiales", sig=ifelse(lamiales.sesmntd.persistent[,"mntd.obs.p"] <= 0.05 | lamiales.sesmntd.persistent[,"mntd.obs.p"] > 0.95, 1,0))),
                             (caryo.sesmntd.persistent %>%  mutate(summits = rownames(caryo.sesmntd.persistent), metric = "mntd", clade = "Caryophyllales", sig=ifelse(caryo.sesmntd.persistent[,"mntd.obs.p"] <= 0.05 | caryo.sesmntd.persistent[,"mntd.obs.p"] > 0.95, 1,0))))
names(persistent.poolSESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

persistent.poolSESmpd <- (rbind((persistent.sesmpd.phylonull %>%  mutate(summits = rownames(persistent.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(persistent.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | persistent.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                             (aster.sesmpd.persistent %>%  mutate(summits = rownames(aster.sesmpd.persistent), metric = "mpd", clade = "Asterales", sig=ifelse(aster.sesmpd.persistent[,"mpd.obs.p"] <= 0.05 | aster.sesmpd.persistent[,"mpd.obs.p"] > 0.95, 1,0))), 
                             (poales.sesmpd.persistent %>%  mutate(summits = rownames(poales.sesmpd.persistent), metric = "mpd", clade = "Poales", sig=ifelse(poales.sesmpd.persistent[,"mpd.obs.p"] <= 0.05 | poales.sesmpd.persistent[,"mpd.obs.p"] > 0.95, 1,0))),
                             (rosales.sesmpd.persistent %>%  mutate(summits = rownames(rosales.sesmpd.persistent), metric = "mpd", clade = "Rosales", sig=ifelse(rosales.sesmpd.persistent[,"mpd.obs.p"] <= 0.05 | rosales.sesmpd.persistent[,"mpd.obs.p"] > 0.95, 1,0))),
                             (lamiales.sesmpd.persistent %>%  mutate(summits = rownames(lamiales.sesmpd.persistent), metric = "mpd", clade = "Lamiales", sig=ifelse(lamiales.sesmpd.persistent[,"mpd.obs.p"] <= 0.05 | lamiales.sesmpd.persistent[,"mpd.obs.p"] > 0.95, 1,0))),
                             (caryo.sesmpd.persistent %>%  mutate(summits = rownames(caryo.sesmpd.persistent), metric = "mpd", clade = "Caryophyllales", sig=ifelse(caryo.sesmpd.persistent[,"mpd.obs.p"] <= 0.05 | caryo.sesmpd.persistent[,"mpd.obs.p"] > 0.95, 1,0)))))
names(persistent.poolSESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")
persistent.poolSES <- rbind(persistent.poolSESmntd, persistent.poolSESmpd)                               
head(persistent.poolSES)     
persistent.poolSES <- persistent.poolSES[!persistent.poolSES$summits == "Ecrins NP",]
persistent.poolSES$clade <- factor(persistent.poolSES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
persistent.poolSES <- cbind(persistent.poolSES, pool = rep(x = "Persistent LGM", times = nrow(persistent.poolSES)))
write.csv(persistent.poolSES, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/persistent.pool.SES.csv")


master.ses.alpha <- rbind(phylogeny.poolSES, summit.poolSES, persistent.poolSES)
head(master.ses.alpha)
write.csv(master.ses.alpha, file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/master.ses.static.alpha.csv")


