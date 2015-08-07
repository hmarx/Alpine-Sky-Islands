############################################## ECRINS ALPINE COMMUNITIIES ################################

source("analysisSkyIsl.R")
source("R/randomizeSourcePoolNull.R")

########################################### Basic Phylo Diversity ########################################### 

###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

### pez will trim data/ phylogeny; no need to use drop tip objects
pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = alps.traits)

pezAlpesNoTrait <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env) ### NEED TO FIGURE THIS OUT: see prepPipeline; tail(comMerge)

pezAlpes
sites(pezAlpes)

## Calculate shape metrics: structure of each community (Phylogenetic Alpha Diveristy)
shape.Alpes <- shape(data = pezAlpes, metric = "all-quick")
coef(shape.Alpes)
#write.csv(coef(shape.Alpes), file="output/9_PhyoDiversity/shape/shape.Alpes.csv")

## Calculate eveness: incorporate species abundances ("abundance" proxy = # of time a spcies was counted in a releve)
evenness.Alpes <- evenness(data = pezAlpes, metric = "all-quick")
coef(evenness.Alpes)
#write.csv(coef(evenness.Alpes), file="output/9_PhyoDiversity/evenness/evenness.Alpes.csv")

## Calculate dissimilarity: compare diversity between communities (Phylogenetic Beta Diveristy)
dissimilarity.Alpes <- dissimilarity(data = pezAlpes, metric = "all") ## done on server
#write.csv(as.matrix(dissimilarity.Alpes$unifrac), file="dissimilarity.Alpes.unifrac.csv")
#write.csv(as.matrix(dissimilarity.Alpes$pcd$PCD), file="dissimilarity.Alpes.pcd.csv")
#write.csv(as.matrix(dissimilarity.Alpes$phylosor), file="dissimilarity.Alpes.phylosor.csv")
#write.csv(as.matrix(dissimilarity.Alpes$comdist), file="dissimilarity.Alpes.comdist.csv")

# Caculate dispersion: does PD differ from random expectation ## isn't working on pez...

######## calculate ses.mpd
comm.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
#write.csv(comm.sesmpd.phylonull, file="output/9_PhyoDiversity/comm.sesmpd.phylonull.csv")

comm.sesmpd.taxanull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)

comm.sesmpd.samplenull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)

comm.sesmpd.swap <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "independentswap", abundance.weighted = FALSE, runs = 999)

######## calculate ses.mntd
comm.sesmntd <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
#write.csv(comm.sesmntd, file="output/9_PhyoDiversity/comm.sesmntd.phylonull.csv")

comm.sesmntd.taxanull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)

comm.sesmntd.samplenull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)

## Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate phylogenetic evenness, 
# or a greater phylogenetic distance among co-occurring species than expected. 

########## Helmus et al.:
# Phylogenetic Species Veriability (PSV): the expected variance among spe- cies in a community phylogeny for a trait evolving under Brownian motion
## should be exactly half the mpd value when the phylogeny is ultrametric

# PSR: Phylogenetic Species RIchness: equivalent to multiplying mpd by the number of species in the community 



###################################################################################### 
############################### Reduced Pool : Summits ##################################### 
###################################################################################### 
summits.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
summits.sites <- filter(summits.sites, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

pezAlpes.summits <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(1:3),], env = alps.env, traits = alps.traits)
pezAlpes.summits # 226 taxa

summits.sesmpd.phylonull <- ses.mpd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

summits.sesmntd.phylonull <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)



###################################################################################### 
############################### Reduced Pool : Persistent ############################ 
###################################################################################### 
persistent.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
persistent.sites <- filter(persistent.sites, persistant > 0)
head(persistent.sites)
rownames(persistent.sites) <- persistent.sites$taxa
dim(persistent.sites)
persistent.sites <- t(persistent.sites[-1])
persistent.sites <- data.matrix(persistent.sites)
head(persistent.sites)

pezAlpes.persistent <- comparative.comm(phy = alps.phy, comm = persistent.sites[-3,], env = alps.env, traits = alps.traits)
pezAlpes.persistent # 172 taxa

persistent.sesmpd.phylonull <- ses.mpd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

persistent.sesmntd.phylonull <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)



###################################################################################### 
############################### Reduced Pool : UnderIce ##############################
###################################################################################### 
UnderIce.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
UnderIce.sites <- filter(UnderIce.sites, underIce > 0)
head(UnderIce.sites)
rownames(UnderIce.sites) <- UnderIce.sites$taxa
dim(UnderIce.sites)
UnderIce.sites <- t(UnderIce.sites[-1])
UnderIce.sites <- data.matrix(UnderIce.sites)
head(UnderIce.sites)

pezAlpes.UnderIce <- comparative.comm(phy = alps.phy, comm = UnderIce.sites[-3,], env = alps.env, traits = alps.traits)
pezAlpes.UnderIce # 109 taxa

UnderIce.sesmpd.phylonull <- ses.mpd(pezAlpes.UnderIce$comm, cophenetic.phylo(pezAlpes.UnderIce$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

UnderIce.sesmntd.phylonull <- ses.mntd(pezAlpes.UnderIce$comm, cophenetic.phylo(pezAlpes.UnderIce$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)


###################################################################################### 
############################### Null : Source Pool ###################################
###################################################################################### 

######################################################################################  
##############  Weight for occurences above glacier / total # occurence in summits

##### Summits Source Pool
comm.sesmpd.sourcePoolNull <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "Summits", N =999)
#write.csv(comm.sesmpd.sourcePoolNull, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.csv")

comm.sesmntd.sourcePoolNull <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "Summits", N =999)
#write.csv(comm.sesmntd.sourcePoolNull, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.csv")

##### Persistant Source Pool
comm.sesmpd.sourcePoolNull.pers <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "persistant", N =999)
#write.csv(comm.sesmpd.sourcePoolNull.pers, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.pers.csv")

comm.sesmntd.sourcePoolNull.pers <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "persistant", N =999)
#write.csv(comm.sesmntd.sourcePoolNull.pers, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.pers.csv")

##### UnderIce Source Pool
comm.sesmpd.sourcePoolNull.under <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "underIce", N =999)
#write.csv(comm.sesmpd.sourcePoolNull.under, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.under.csv")

comm.sesmntd.sourcePoolNull.under <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "underIce", N =999)
#write.csv(comm.sesmntd.sourcePoolNull.under, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.under.csv")


###################################################################################### 
############################### Null : Historic Source Pool = Persistent ###################################
###################################################################################### 

## Contemporary species pool = summits 
## Hisoric source pool = persistent abouve glacier through LGM

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

summits.sesmpd.sourcePersis <- ses.mpd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistant", N =999)
summits.sesmntd.sourcePersis <- ses.mntd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "persistant", N =999)


## Contemporary species pool = summits 
## Hisoric source pool = summits (persistent + under ice?)

summits.sesmpd.sourceSummits <- ses.mpd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)
summits.sesmntd.sourceSUmmits <- ses.mntd.sourcePool(phy=pezAlpes.summits$phy, com = pezAlpes.summits$comm, sourcePool = "Summits", N =999)




