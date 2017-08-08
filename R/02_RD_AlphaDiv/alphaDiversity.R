#####################################################################################################################
############# Phylogenetic diversity within alpine summits (alpha) ################################################## 
############# Static Null Model ##################################################################################### 
############# Hannah E. Marx, 25 April 2017 ########################################################################### 
#####################################################################################################################

source("analysisSkyIsl.R")
source("R/functions/chooseClade.R")

###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool: Regional Ecrins NP ##################### 
###################################################################################### 

## Random resample from phylogeny pool (==Ecrins NP), equal probability random draw from phylogeny pool
ecrins.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
ecrins.sesmntd.phylonull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool: All Summits ############################
###################################################################################### 

## Source pool = all summits, equal probability random draw from phylogeny (pruned to summits)
summits.sesmpd.phylonull <- ses.mpd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
summits.sesmntd.phylonull <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

###################################################################################### 
############################### All Spermatophyta ####################################
############################### Species Pool : Persistent through LGM ################
###################################################################################### 

## Source pool = Persistent, equal probability random draw from phylogeny (pruned to persistent species)
persistent.sesmpd.phylonull <- ses.mpd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
persistent.sesmntd.phylonull <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: Regional Ecrins NP #####################
###################################################################################### 

# Convert to pres/abs
alps.sites.pa.all <- t(pezAlpes$comm)

# Pruend phylogeny and community matrix for each clade 
ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa.all)
aster.sesmpd.phylonull <- ses.mpd(t(alps.aster$data), cophenetic.phylo(alps.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.phylonull <- ses.mntd(t(alps.aster$data), cophenetic.phylo(alps.aster$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa.all)
poales.sesmpd.phylonull <- ses.mpd(t(alps.poales$data), cophenetic.phylo(alps.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.phylonull <- ses.mntd(t(alps.poales$data), cophenetic.phylo(alps.poales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.rosales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Rosales")
pruned.tree.alps.rosales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.rosales])
alps.rosales <- treedata(pruned.tree.alps.rosales, alps.sites.pa.all)
rosales.sesmpd.phylonull <- ses.mpd(t(alps.rosales$data), cophenetic.phylo(alps.rosales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.phylonull <- ses.mntd(t(alps.rosales$data), cophenetic.phylo(alps.rosales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa.all)
lamiales.sesmpd.phylonull <- ses.mpd(t(alps.lamiales$data), cophenetic.phylo(alps.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.phylonull <- ses.mntd(t(alps.lamiales$data), cophenetic.phylo(alps.lamiales$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.alps.caryo <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.caryo <- treedata(pruned.tree.alps.caryo, alps.sites.pa.all)
caryo.sesmpd.phylonull <- ses.mpd(t(alps.caryo$data), cophenetic.phylo(alps.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.phylonull <- ses.mntd(t(alps.caryo$data), cophenetic.phylo(alps.caryo$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

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


###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: All Summits ############################
###################################################################################### 

# Pres/abs for summits pool
alps.summits.pa.all <- t(pezAlpes.summits$comm)

alps.aster.summits <- treedata(pruned.tree.alps.aster, alps.summits.pa.all)
aster.sesmpd.summits <- ses.mpd(t(alps.aster.summits$data), cophenetic.phylo(alps.aster.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.summits <- ses.mntd(t(alps.aster.summits$data), cophenetic.phylo(alps.aster.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.poales.summits <- treedata(pruned.tree.alps.poales, alps.summits.pa.all)
poales.sesmpd.summits <- ses.mpd(t(alps.poales.summits$data), cophenetic.phylo(alps.poales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.summits <- ses.mntd(t(alps.poales.summits$data), cophenetic.phylo(alps.poales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.rosales.summits <- treedata(pruned.tree.alps.rosales, alps.summits.pa.all)
rosales.sesmpd.summits <- ses.mpd(t(alps.rosales.summits$data), cophenetic.phylo(alps.rosales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.summits <- ses.mntd(t(alps.rosales.summits$data), cophenetic.phylo(alps.rosales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.lamiales.summits <- treedata(pruned.tree.alps.lamiales, alps.summits.pa.all)
lamiales.sesmpd.summits <- ses.mpd(t(alps.lamiales.summits$data), cophenetic.phylo(alps.lamiales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.summits <- ses.mntd(t(alps.lamiales.summits$data), cophenetic.phylo(alps.lamiales.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.caryo.summits <- treedata(pruned.tree.alps.caryo, alps.summits.pa.all)
caryo.sesmpd.summits <- ses.mpd(t(alps.caryo.summits$data), cophenetic.phylo(alps.caryo.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.summits <- ses.mntd(t(alps.caryo.summits$data), cophenetic.phylo(alps.caryo.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

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


###################################################################################### 
############################### Five Clades Separately ###############################
############################### Species Pool: Persistent through LGM #################
###################################################################################### 

# Pres/abs for persistent pool
alps.persistent.pa.all <- t(pezAlpes.persistent$comm)

alps.aster.persistent <- treedata(pruned.tree.alps.aster, alps.persistent.pa.all)
aster.sesmpd.persistent <- ses.mpd(t(alps.aster.persistent$data), cophenetic.phylo(alps.aster.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
aster.sesmntd.persistent <- ses.mntd(t(alps.aster.persistent$data), cophenetic.phylo(alps.aster.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.poales.persistent <- treedata(pruned.tree.alps.poales, alps.persistent.pa.all)
poales.sesmpd.persistent <- ses.mpd(t(alps.poales.persistent$data), cophenetic.phylo(alps.poales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
poales.sesmntd.persistent <- ses.mntd(t(alps.poales.persistent$data), cophenetic.phylo(alps.poales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.rosales.persistent <- treedata(pruned.tree.alps.rosales, alps.persistent.pa.all)
rosales.sesmpd.persistent <- ses.mpd(t(alps.rosales.persistent$data), cophenetic.phylo(alps.rosales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
rosales.sesmntd.persistent <- ses.mntd(t(alps.rosales.persistent$data), cophenetic.phylo(alps.rosales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.lamiales.persistent <- treedata(pruned.tree.alps.lamiales, alps.persistent.pa.all)
lamiales.sesmpd.persistent <- ses.mpd(t(alps.lamiales.persistent$data), cophenetic.phylo(alps.lamiales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
lamiales.sesmntd.persistent <- ses.mntd(t(alps.lamiales.persistent$data), cophenetic.phylo(alps.lamiales.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

alps.caryo.persistent <- treedata(pruned.tree.alps.caryo, alps.persistent.pa.all)
caryo.sesmpd.persistent <- ses.mpd(t(alps.caryo.persistent$data), cophenetic.phylo(alps.caryo.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
caryo.sesmntd.persistent <- ses.mntd(t(alps.caryo.persistent$data), cophenetic.phylo(alps.caryo.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

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

master.ses.alpha <- rbind(phylogeny.poolSES, summit.poolSES, persistent.poolSES)
head(master.ses.alpha)
master.ses.alpha$pool <- factor(master.ses.alpha$pool, levels = c( "Ecrins NP",  "Summits", "Persistent LGM"))
master.ses.alpha$pool <- factor(master.ses.alpha$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = "LGM"))
master.ses.alpha$clade <- factor(master.ses.alpha$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
master.ses.alpha$metric <- factor(master.ses.alpha$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))
master.ses.alpha[master.ses.alpha == "Persistent"] <- "LGM"
master.ses.alpha[master.ses.alpha == "Summits"] <-  "All Summits"
write.csv(master.ses.alpha, file="output/8_PhyoDiversity/alpha/static/master.ses.static.alpha.csv")


