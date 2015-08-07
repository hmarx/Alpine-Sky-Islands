############################################## ECRINS ALPINE COMMUNITIIES ################################

source("analysisSkyIsl.R")
source("R/randomizeSourcePoolNull_v1.4.R")

########################################### Basic Phylo Diversity ########################################### 

###################################################################################### 
############################### Contemporary Source Pool : Ecrins NP ##################################### 
###################################################################################### 

### pez will trim data/ phylogeny; no need to use drop tip objects
pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = alps.traits)
pezAlpes.dummy <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = alps.traits.dummy)

pezAlpesNoTrait <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env) ### NEED TO FIGURE THIS OUT: see prepPipeline; tail(comMerge)

pezAlpes
sites(pezAlpes)

comm.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
#write.csv(comm.sesmpd.phylonull, file="output/9_PhyoDiversity/comm.sesmpd.phylonull.csv")

comm.sesmntd <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
#write.csv(comm.sesmntd, file="output/9_PhyoDiversity/comm.sesmntd.phylonull.csv")

ecrins.weighted.summits.sesmntd <- ses.mntd.sourcePool(com = pezAlpes$comm, dist = cophenetic.phylo(pezAlpes$phy), sourcePool = "Summits", N = 999)
ecrins.weighted.summits.sesmpd <- ses.mpd.sourcePool(com = pezAlpes$comm, dist = cophenetic.phylo(pezAlpes$phy), sourcePool = "Summits", N = 999)

ecrins.weighted.persistent.sesmntd <- ses.mntd.sourcePool(com = pezAlpes$comm, dist = cophenetic.phylo(pezAlpes$phy), sourcePool = "persistent", N = 999)
ecrins.weighted.persistent.sesmpd <- ses.mpd.sourcePool(com = pezAlpes$comm, dist = cophenetic.phylo(pezAlpes$phy), sourcePool = "persistent", N = 999)


######### picante manual
## Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate phylogenetic evenness, 
# or a greater phylogenetic distance among co-occurring species than expected. 

########## Helmus et al.:
# Phylogenetic Species Veriability (PSV): the expected variance among spe- cies in a community phylogeny for a trait evolving under Brownian motion
## should be exactly half the mpd value when the phylogeny is ultrametric

# PSR: Phylogenetic Species RIchness: equivalent to multiplying mpd by the number of species in the community 



###################################################################################### 
############################### Contemporary Pool: Summits ###########################
###################################################################################### 

## Contemporary species pool = summits 
summits.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
summits.sites <- filter(summits.sites, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

## Reduce data, without traits
pezAlpes.summits <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(2:3),], env = alps.env)
pezAlpes.summits # 231  taxa

## Source pool = summits, equal probability
summits.sesmpd.phylonull <- ses.mpd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
summits.sesmntd.phylonull <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

## source pool = summits (persistent + under ice?), weighted for abundance
summits.sesmpd.sourceSummits <- ses.mpd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)
summits.sesmntd.sourceSUmmits <- ses.mntd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)

## source pool = persistent abouve glacier through LGM, weighted for abundance...doesn't really make sense because this is a smaller source population
#summits.sesmpd.sourcePersis <- ses.mpd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistent", N =999)
#summits.sesmntd.sourcePersis <- ses.mntd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistent", N =999)


###################################################################################### 
############################### Contemporary Pool : Persistent ############################ 
###################################################################################### 

persistent.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
persistent.sites <- filter(persistent.sites, persistent > 0)
head(persistent.sites)
rownames(persistent.sites) <- persistent.sites$taxa
dim(persistent.sites)
persistent.sites <- t(persistent.sites[-1])
persistent.sites <- data.matrix(persistent.sites)
head(persistent.sites)

pezAlpes.persistent <- comparative.comm(phy = alps.phy, comm = persistent.sites[-c(2:3),], env = alps.env)
pezAlpes.persistent # 175 taxa

## Source pool = Persistent, equal probability
persistent.sesmpd.phylonull <- ses.mpd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
persistent.sesmntd.phylonull <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

## source pool = summits (persistent + under ice?), weighted for abundance
persistent.sesmpd.sourceSummits <- ses.mpd.sourcePool(phy=pezAlpes.persistent$phy, com = pezAlpes.persistent$comm, sourcePool = "Summits", N =999)
persistent.sesmntd.sourceSummits <- ses.mntd.sourcePool(phy=pezAlpes.persistent$phy, com = pezAlpes.persistent$comm, sourcePool = "Summits", N =999)

## source pool = persistent abouve glacier through LGM, weighted for abundance
persistent.sesmpd.sourcePersis <- ses.mpd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistent", N =999)
persistent.sesmntd.sourcePersis <- ses.mntd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistent", N =999)


