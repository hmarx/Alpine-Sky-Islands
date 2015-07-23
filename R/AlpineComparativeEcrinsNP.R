############################################## ECRINS ALPINE COMMUNITIIES ################################

source("analysisSkyIsl.R")

########################################### Basic Phylo Diversity ########################################### 

###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

### pez will trim data/ phylogeny; no need to use drop tip objects
pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = alps.traits)

pezAlpesNoTrait <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env) ### NEED TO FIGURE THIS OUT: see prepPipeline; tail(comMerge)

pezAlpes
sites(pezAlpes)
summary(pezAlpes$comm)

#Get assemblage phylogenies of all sites
assemblage.phylogenies(data = pezAlpes)
pezAlpes$comm[,pezAlpes$comm["EcrinsNP",] == 0]

## Calculate shape metrics: structure of each community (Phylogenetic Alpha Diveristy)
shape.Alpes <- shape(data = pezAlpes, metric = "all-quick")
coef(shape.Alpes)
#write.csv(coef(shape.Alpes), file="output/9_PhyoDiversity/shape.Alpes.csv")

## Calculate eveness: incorporate species abundances ("abundance" proxy = # of time a spcies was counted in a releve)
evenness.Alpes <- evenness(data = pezAlpes, metric = "all-quick")
coef(evenness.Alpes)
#write.csv(coef(evenness.Alpes), file="output/9_PhyoDiversity/evenness.Alpes.csv")

## Calculate dissimilarity: compare diversity between communities (Phylogenetic Beta Diveristy)
dissimilarity.Alpes <- dissimilarity(data = pezAlpes, metric = "all", permute = 10)
dissimilarity.Alpes$unifrac
dissimilarity.Alpes$pcd$PCD
dissimilarity.Alpes$phylosor
dissimilarity.Alpes$comdist
#write.csv(as.matrix(dissimilarity.Alpes$unifrac), file="dissimilarity.Alpes.unifrac.csv")
#write.csv(as.matrix(dissimilarity.Alpes$pcd$PCD), file="dissimilarity.Alpes.pcd.csv")
#write.csv(as.matrix(dissimilarity.Alpes$phylosor), file="dissimilarity.Alpes.phylosor.csv")
#write.csv(as.matrix(dissimilarity.Alpes$comdist), file="dissimilarity.Alpes.comdist.csv")

# Caculate dispersion: does PD differ from random expectation 
dispersion.Alpes <- dispersion(data = pezAlpes, metric = "all", null.model = "phylogeny.pool", permute = 99)
coef(dispersion.Alpes)
#write.csv(coef(disp.divs), file="disp.phylpool.Alps.csv")

######## calculate ses.mpd
pd(pezAlpes$comm, pezAlpes$phy)
comm.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
comm.sesmpd.phylonull
#write.csv(comm.sesmpd.phylonull, file="output/9_PhyoDiversity/comm.sesmpd.phylonull.csv")

comm.sesmpd.taxanull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
comm.sesmpd.taxanull

comm.sesmpd.samplenull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)
comm.sesmpd.samplenull

comm.sesmpd.swap <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "independentswap", abundance.weighted = FALSE, runs = 999)
comm.sesmpd.swap

######## calculate ses.mntd
comm.sesmntd <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
head(comm.sesmntd)
#write.csv(comm.sesmntd, file="output/9_PhyoDiversity/comm.sesmntd.phylonull.csv")

comm.sesmntd.taxanull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
comm.sesmntd.taxanull
comm.sesmntd.taxanull[comm.sesmntd.taxanull$mntd.obs.p <= 0.05,]

comm.sesmntd.samplenull <- ses.mntd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "sample.pool", abundance.weighted = FALSE, runs = 999)
comm.sesmntd.samplenull

## Positive SES values (mpd.obs.z > 0) and high quantiles (mpd.obs.p > 0.95) indicate phylogenetic evenness, 
# or a greater phylogenetic distance among co-occurring species than expected. 

########## Helmus et al.:
# Phylogenetic Species Veriability (PSV): the expected variance among spe- cies in a community phylogeny for a trait evolving under Brownian motion
## should be exactly half the mpd value when the phylogeny is ultrametric

# PSR: Phylogenetic Species RIchness: equivalent to multiplying mpd by the number of species in the community 


###################################  Plotting  ###################################  

shape.divs.long <- melt(as.matrix(coef(shape.Alpes))) # shape.divs
colnames (shape.divs.long) <- c("summit", "metric", "value")
head(shape.divs.long)

colourCount = length(unique(shape.divs.long$metric))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(shape.divs.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_log10("metric")  +
  scale_x_discrete("Summit (increasing PD)") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Shape Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) #+
  #theme(plot.margin = units(c(0.5,2,0.5,0.5), "cm"))

comm.sesmntd[,"sig"] <- ifelse(comm.sesmntd$mntd.obs.p <= 0.05, "TRUE", "FALSE")
comm.sesmpd.phylonull[,"sig1"] <- ifelse(comm.sesmpd.phylonull$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(comm.sesmntd[c(1,6,9)], comm.sesmpd.phylonull[c(1,6,9)], id=rownames(comm.sesmntd))
dispersion.metrics.long <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (dispersion.metrics.long) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(dispersion.metrics.long)
#11                 mont.pelvoux    17 mntd.obs.z  1.7740514

dispersion.metrics.long <- arrange(dispersion.metrics.long, ntaxa)
dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, levels = unique(dispersion.metrics.long$summit))

# Manual levels
#disp_table <- table(dispersion.metrics.long$summit)
#disp_levels <- c(names(disp_table[c(5, 20, 13, 23)]), (names(disp_table)[order(disp_table)][c(-5, -20, -13, -23)]))
#dispersion.metrics.long$summit <- factor(dispersion.metrics.long$summit, as.character(disp_levels))

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
ggplot(dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nEcrins Species Pool: null.model = phylogeny.pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()

###################################################################################### 
############################### Pool : Summits ##################################### 
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
summits.sesmpd.phylonull

summits.sesmntd.phylonull <- ses.mntd(pezAlpes.summits$comm, cophenetic.phylo(pezAlpes.summits$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
summits.sesmntd.phylonull


summits.dispersion.metrics.long <- melt(as.matrix(cbind(summits.sesmpd.phylonull[6], summits.sesmntd.phylonull[6])))
colnames (summits.dispersion.metrics.long) <- c("summit", "metric", "value")
head(summits.dispersion.metrics.long)

# Manual levels
disp_table <- table(summits.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[18]), (names(disp_table)[order(disp_table)][-18]))
summits.dispersion.metrics.long$summit <- factor(summits.dispersion.metrics.long$summit, levels = disp_levels)

#pdf(file="output/9_PhyoDiversity/DispersionMetrics_summitPool.pdf") 
ggplot(summits.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Dispersion Metrics\nSummit Species Pool: null.model = phylogeny.pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()


###################################################################################### 
############################### Pool : Persistent ##################################### 
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
persistent.sesmpd.phylonull

persistent.sesmntd.phylonull <- ses.mntd(pezAlpes.persistent$comm, cophenetic.phylo(pezAlpes.persistent$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
persistent.sesmntd.phylonull

persistent.dispersion.metrics.long <- melt(as.matrix(cbind(persistent.sesmpd.phylonull[6], persistent.sesmntd.phylonull[6])))
colnames(persistent.dispersion.metrics.long) <- c("summit", "metric", "value")
head(persistent.dispersion.metrics.long)

# Manual levels
disp_table <- table(persistent.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[c(12, 19, 22)]), (names(disp_table)[order(names(disp_table))][c(-12, -19, -22)]))
persistent.dispersion.metrics.long$summit <- factor(persistent.dispersion.metrics.long$summit, levels = disp_levels)

pdf(file="output/9_PhyoDiversity/DispersionMetrics_persistentPool.pdf") 
ggplot(persistent.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Persistent Pool: Dispersion Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()


###################################################################################### 
############################### Pool : UnderIce ##################################### 
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
UnderIce.sesmpd.phylonull

UnderIce.sesmntd.phylonull <- ses.mntd(pezAlpes.UnderIce$comm, cophenetic.phylo(pezAlpes.UnderIce$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
UnderIce.sesmntd.phylonull

UnderIce.dispersion.metrics.long <- melt(as.matrix(cbind(UnderIce.sesmpd.phylonull[6], UnderIce.sesmntd.phylonull[6])))
colnames(UnderIce.dispersion.metrics.long) <- c("summit", "metric", "value")
head(UnderIce.dispersion.metrics.long)

# Manual levels
disp_table <- table(UnderIce.dispersion.metrics.long$summit)
disp_levels <- c(names(disp_table[c(22, 19, 12)]), (names(disp_table)[order(names(disp_table))][c(-22, -12, -19)]))
UnderIce.dispersion.metrics.long$summit <- factor(UnderIce.dispersion.metrics.long$summit, levels = disp_levels)

#pdf(file="output/9_PhyoDiversity/DispersionMetrics_UnderIcePool.pdf") 
ggplot(UnderIce.dispersion.metrics.long, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric)) +
  scale_y_continuous("obs.z : null=phylogeny.pool")  +
  scale_x_discrete("Summit") +
  #facet_wrap(~metric) + 
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("UnderIce Pool: Dispersion Metrics") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()


######## Summarize Source pool effects

mpd.pools <- list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, persistent.sesmpd.phylonull, UnderIce.sesmpd.phylonull)
mpd.pools.melt <- melt(mpd.pools, id=c(colnames(comm.sesmpd.phylonull)))
mpd.pools.melt <- na.omit(mpd.pools.melt)
mpd.pools.melt$L1 <- factor(mpd.pools.melt$L1)

pdf(file="output/9_PhyoDiversity/mpd.sourcepools.pdf") 
p1 <- ggplot(mpd.pools.melt, aes(mpd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mpd for different source pools") 
p1
dev.off()


mntd.pools <- list(comm.sesmntd, summits.sesmntd.phylonull, persistent.sesmntd.phylonull, UnderIce.sesmntd.phylonull)
mntd.pools.melt <- melt(mntd.pools, id=c(colnames(comm.sesmntd)))
mntd.pools.melt <- na.omit(mntd.pools.melt)
mntd.pools.melt$L1 <- factor(mntd.pools.melt$L1)

#pdf(file="output/9_PhyoDiversity/mntd.sourcepools.pdf") 
p1 <- ggplot(mntd.pools.melt, aes(mntd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4"),
                       labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mntd for different source pools") 
p1
#dev.off()

###############################################################################################  
################################  Weight for occurences above glacier / total # occurence in summits
pezAlpes$comm["Summits",]
comm.sesmpd.phylonull <- ses.mpd(pezAlpes$comm, cophenetic.phylo(pezAlpes$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)


comm.summitsPool <- randomizeSourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "Summits", N =5)
dim(comm.summitsPool[[1]])
dim(pezAlpes$comm)

##### Summits Source Pool
comm.sesmpd.sourcePoolNull <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "Summits", N =999)
write.csv(comm.sesmpd.sourcePoolNull, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.csv")

comm.sesmntd.sourcePoolNull <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "Summits", N =999)
write.csv(comm.sesmntd.sourcePoolNull, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.csv")

##### Persistant Source Pool
comm.sesmpd.sourcePoolNull.pers <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "persistant", N =999)
write.csv(comm.sesmpd.sourcePoolNull.pers, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.pers.csv")

comm.sesmntd.sourcePoolNull.pers <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "persistant", N =999)
write.csv(comm.sesmntd.sourcePoolNull.pers, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.pers.csv")

##### UnderIce Source Pool
comm.sesmpd.sourcePoolNull.under <- ses.mpd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "underIce", N =999)
write.csv(comm.sesmpd.sourcePoolNull.under, file="output/9_PhyoDiversity/dispersion/comm.sesmpd.sourcePoolNull.under.csv")

comm.sesmntd.sourcePoolNull.under <- ses.mntd.sourcePool(phy=pezAlpes$phy, com = pezAlpes$comm, sourcePool = "underIce", N =999)
write.csv(comm.sesmntd.sourcePoolNull.under, file="output/9_PhyoDiversity/dispersion/comm.sesmntd.sourcePoolNull.under.csv")


mntd.pools.sourcePool <- list(comm.sesmntd.sourcePoolNull, comm.sesmntd.sourcePoolNull.pers, comm.sesmntd.sourcePoolNull.under)
mntd.pools.sourcePool.melt <- melt(mntd.pools.sourcePool, id=c(colnames(comm.sesmntd.sourcePoolNull)))
mntd.pools.sourcePool.melt <- na.omit(mntd.pools.sourcePool.melt)
mntd.pools.sourcePool.melt$L1 <- factor(mntd.pools.sourcePool.melt$L1)

#pdf(file="output/9_PhyoDiversity/mntd.sourcePoolsNull.pdf") 
p1 <- ggplot(mntd.pools.sourcePool.melt, aes(mntd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3"),
                       labels=c("Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mntd for different source pools\nSource Pool Null") 
p1
#dev.off()



mpd.pools.sourcePool <- list(comm.sesmpd.sourcePoolNull, comm.sesmpd.sourcePoolNull.pers, comm.sesmpd.sourcePoolNull.under)
mpd.pools.sourcePool.melt <- melt(mpd.pools.sourcePool, id=c(colnames(comm.sesmpd.sourcePoolNull)))
mpd.pools.sourcePool.melt <- na.omit(mpd.pools.sourcePool.melt)
mpd.pools.sourcePool.melt$L1 <- factor(mpd.pools.sourcePool.melt$L1)

#pdf(file="output/9_PhyoDiversity/mpd.sourcePoolsNull.pdf") 
p1 <- ggplot(mpd.pools.sourcePool.melt, aes(mpd.obs.z, colour=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3"),
                       labels=c("Summits", "Peristent", "Under Ice")) +
  ggtitle("Distribution of SES mpd for different source pools\nSource Pool Null") 
p1
#dev.off()

colnames(mntd.pools.sourcePool.melt)[1] <- "ntaxa"
levels(mntd.pools.sourcePool.melt$L1) <- c(5,6,7)
all.mntd <- rbind(mntd.pools.melt, mntd.pools.sourcePool.melt)

#pdf(file="output/9_PhyoDiversity/mntd.all.sourcepools.pdf") 
p1 <- ggplot(all.mntd, aes(mntd.obs.z, colour=L1, linetype=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4", "5", "6", "7"),
                       labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  scale_linetype_manual(name="Source Pool", breaks=c("1", "2", "3", "4", "5", "6", "7"), values=c(1,1,1,1, 5,5,5), 
                        labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                 "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  ggtitle("Distribution of SES mntd for different source pools") 
p1
#dev.off()


colnames(mpd.pools.sourcePool.melt)[1] <- "ntaxa"
levels(mpd.pools.sourcePool.melt$L1) <- c(5,6,7)
all.mpd <- rbind(mpd.pools.melt, mpd.pools.sourcePool.melt)

#pdf(file="output/9_PhyoDiversity/mpd.all.sourcepools.pdf") 
p1 <- ggplot(all.mpd, aes(mpd.obs.z, colour=L1, linetype=L1)) + 
  geom_density() +
  scale_color_discrete(name="Source Pool",
                       breaks=c("1", "2", "3", "4", "5", "6", "7"),
                       labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  scale_linetype_manual(name="Source Pool", 
                        breaks=c("1", "2", "3", "4", "5", "6", "7"), values=c(1,1,1,1, 5,5,5), 
                        labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                                 "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull")) +
  ggtitle("Distribution of SES mpd for different source pools") 
p1
#dev.off()


################ spacodiR ################ 
phy.nodetimes(phy.Alpes, time.range = c(0, max(phy.Alpes$edge.length)), proportion = TRUE)

