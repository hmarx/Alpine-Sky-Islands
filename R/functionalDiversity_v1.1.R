############################################## ECRINS ALPINE COMMUNITIIES ################################

source("analysisSkyIsl.R")


################################ Prepare dataset ################################ 

head(pezAlpes$data) #1100
summary(pezAlpes$data)

alps.data.complete <- (pezAlpes$data[complete.cases(pezAlpes$data[c(3,5,7,9,10)]),c(3,5,7,9,10)])
head(alps.data.complete)
str(alps.data.complete)
levels(alps.data.complete$WOODY)
pezAlpes.Traits <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = (alps.data.complete)) #599
summary(pezAlpes.Traits$data)

pezAlpes.Traits$data$PL_H_flora_MAX <- log10(pezAlpes.Traits$data$PL_H_flora_MAX)
pezAlpes.Traits$data$SEEDM <- log10(pezAlpes.Traits$data$SEEDM)
pezAlpes.Traits$data$PLOIDY <- as.factor(pezAlpes.Traits$data$PLOIDY)
pezAlpes.Traits$data$VEG_DISP <- as.factor(pezAlpes.Traits$data$VEG_DISP)
pezAlpes.Traits$data$WOODY <- as.factor(pezAlpes.Traits$data$WOODY)

summary(pezAlpes.Traits$data)


################################ Test for phylogenetic signal of traits in community ################################ 

traits.sig <- pezAlpes.Traits$data[pezAlpes.Traits$phy$tip.label, ]

phylo.signal <- multiPhylosignal(traits.sig, pezAlpes.Traits$phy)

# From picante_Intro:
# Phylogenetic signal is a quantitative measure of the degree to which phylogeny predicts the ecological similarity 
# of species. The K statistic is a measure of phy- logenetic signal that compares the observed signal in a trait to 
# the signal under a Brownian motion model of trait evolution on a phylogeny (Blomberg et al. 2003). 
# K values of 1 correspond to a Brownian motion process, which implies some degree of phylogenetic signal or conservatism. 
# K values closer to zero correspond to a random or convergent pattern of evolution, 
# while K values greater than 1 indicate strong phylogenetic signal and conservatism of traits. 
# The statistical significance of phylo- genetic signal can be evaluated by comparing observed patterns of the 
# variance of independent contrasts of the trait to a null model of shuffling taxa labels across the tips of the phylogeny.

# The higher the K statistic, the more phylogenetic signal in a trait. 
# PIC.variance.P is the quantile of the observed phylogenetically independent contrast variance versus the null 
# distribution, which can be used as a 1-tailed P-value to test for greater phylogenetic signal than expected. 
# Traits with PIC.variance.P < 0.05 have non-random phylogenetic signal.

################################ Analysis of Functional Diversity of Alpine Summits ################################ 


## calculate multivariate distance matric of trait values (combining discrete and continuous measures)
alps.traits.grower <- as.matrix(gowdis(x = pezAlpes.Traits$data))
#daisy(pezAlpes.Traits$data, metric = "gower")

dim(alps.traits.grower) #599

traits.mntd.taxalabels <- ses.mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "taxa.labels", runs = 999)
traits.mpd.taxalabels <- ses.mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "taxa.labels", runs = 999)

traits.mntd.independentswap <- ses.mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "independentswap", runs = 999)
traits.mpd.independentswap <- ses.mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "independentswap", runs = 999)

traits.mntd.frequency <- ses.mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "frequency", runs = 999)
traits.mpd.frequency <- ses.mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, null.model = "frequency", runs = 999)

pdf(file="output/8_TraitDiveristy/FunctionalMNTDEcrinsPool.pdf") 
plotDistributionSES(outputSES = list(traits.mntd.taxalabels, 
                                     traits.mntd.independentswap,
                                     traits.mntd.frequency), 
                    breaks=c("1", "2", "3"),
                    labels=c("Taxa Labels", "Independent Swap", "Frequency"),
                    mainTitle = "Distribution of Functional SES mntd", values=c(1,1,5,5))
dev.off()

pdf(file="output/8_TraitDiveristy/FunctionalMPDEcrinsPool.pdf") 
plotDistributionSES(outputSES = list(traits.mpd.taxalabels, 
                                     traits.mpd.independentswap,
                                     traits.mpd.frequency), 
                    breaks=c("1", "2", "3"),
                    labels=c("Taxa Labels", "Independent Swap", "Frequency"),
                    mainTitle = "Distribution of Functional SES mpd", values=c(1,1,5,5))
dev.off()

#pdf(file="output/8_TraitDiveristy/FunctionalDispersionMetricsEcrinsPool.pdf") 
plotSESdispersion(mpd = traits.mpd.taxalabels, mntd = traits.mntd.taxalabels, 
                  mainTitle = "Functional Dispersion \nnull.model = taxa labels")
#dev.off()

### Convex Hull Volume Metrics: 
alps.Traits.FD <- dbFD(x = pezAlpes.Traits$data, a = as.data.frame(pezAlpes.Traits$comm))

#FRic or Functional Richness
alps.Traits.FD$FRic

#Functional Evenness (FEve): The FEve metric utilizes a minimum spanning tree (MST) 
#to connect all species in functional trait
#space and measures the regularity of the species points along the branches of this
#tree and the regularity of their abundances. Thus, we may expect it to be very similar
#to a nearest neighbor metric.

#The FDiv metric first measures the average distance of all species from the centroid
#of the trait space in a community and then sums the magnitude of the divergences
#from that mean. Thus higher values are supposed to indicate more dispersion
#towards the maximum and minimum of the range of traits.

#The FDis metric calculates the distance of each species from the centroid of the
#community traits. Thus, it is not quite the same calculation as the pairwise distance
#between species we can expect that it might be highly correlated. It is also similar
#to the FDiv metric, though conceptually perhaps easier to understand the biological
#meaning of FDis and it is likely a clearer indicator of trait dispersion in a community.
mpd.traits <- mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower)
mntd.traits <- mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower)
mpd.traits.ab <- mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, abundance.weighted = T)
mntd.traits.ab <- mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, abundance.weighted = T)
outputs <- as.data.frame(cbind(mpd.traits, mpd.traits.ab, mntd.traits, mntd.traits.ab, alps.Traits.FD$FRic, alps.Traits.FD$FEve, alps.Traits.FD$FDiv, alps.Traits.FD$FDis))

plot(outputs, pch=16)
cor(outputs)


########### Constrain Species Pool to Summits
## Contemporary species pool = summits 
summits.sites <- as.data.frame(cbind("taxa" = colnames(alps.sites), as.data.frame(t(alps.sites))))
summits.sites <- filter(summits.sites, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

## Reduce data, with traits (complete cases, log transformed)
pezAlpes.summits <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(2:3),], env = alps.env, traits =  pezAlpes.Traits$data)
pezAlpes.summits # 103  taxa

## calculate multivariate distance matric of trait values (combining discrete and continuous measures)
summits.traits.grower <- as.matrix(gowdis(x = pezAlpes.summits$data))
#daisy(pezAlpes.Traits$data, metric = "gower")

dim(summits.traits.grower) #103

sumits.traits.mntd.taxalabels <- ses.mntd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "taxa.labels", runs = 999)
sumits.traits.mpd.taxalabels <- ses.mpd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "taxa.labels", runs = 999)

sumits.traits.mntd.independentswap <- ses.mntd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "independentswap", runs = 999)
sumits.traits.mpd.independentswap <- ses.mpd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "independentswap", runs = 999)

sumits.traits.mntd.frequency <- ses.mntd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "frequency", runs = 999)
sumits.traits.mpd.frequency <- ses.mpd(samp = pezAlpes.summits$comm, dis = summits.traits.grower, null.model = "frequency", runs = 999)

#pdf(file="output/8_TraitDiveristy/FunctionalMNTDSummitsPool.pdf") 
plotDistributionSES(outputSES = list(sumits.traits.mntd.taxalabels, 
                                     sumits.traits.mntd.independentswap,
                                     sumits.traits.mntd.frequency), 
                    breaks=c("1", "2", "3"),
                    labels=c("Taxa Labels", "Independent Swap", "Frequency"),
                    mainTitle = "Distribution of Functional SES mntd on Summits", values=c(1,1,1))
#dev.off()

#pdf(file="output/8_TraitDiveristy/FunctionalMPDMetricsSummitsPool.pdf") 
plotDistributionSES(outputSES = list(sumits.traits.mpd.taxalabels, 
                                     sumits.traits.mpd.independentswap,
                                     sumits.traits.mpd.frequency), 
                    breaks=c("1", "2", "3"),
                    labels=c("Taxa Labels", "Independent Swap", "Frequency"),
                    mainTitle = "Distribution of Functional SES mpd on Summits", values=c(1,1,1))
#dev.off()

#pdf(file="output/8_TraitDiveristy/FunctionalDispersionMetricsSummitsPool.pdf") 
plotSESdispersion(mpd = sumits.traits.mpd.frequency, mntd = sumits.traits.mntd.frequency, 
                  mainTitle = "Functional Dispersion on Summits \nnull.model = frequency")
#dev.off()


summits.sesmpd.functional.sourceSummits <- ses.mpd.sourcePool(com=pezAlpes.summits$comm, 
                                                              dist = summits.traits.grower, sourcePool = "Summits", N =999)




plotDistributionSES(outputSES = list(sumits.traits.mpd.taxalabels, 
                                     sumits.traits.mpd.independentswap,
                                     sumits.traits.mpd.frequency, summits.sesmpd.functional.sourceSummits), 
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Taxa Labels", "Independent Swap", "Frequency", "Source Summits"),
                    mainTitle = "Distribution of Functional SES mpd on Summits", values=c(1,1,1,5))


######################################## Seed Mass ######################################## 

#### dataset with complete data for this trait 
SeedMass <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(2:3),], env = alps.env, 
                             traits = log10(alps.traits[complete.cases(alps.traits$SEEDM),]["SEEDM"]))

traits.SeedMass <- data.matrix(dist(SeedMass$data))
dim(traits.SeedMass) #143

traits.mntd.SeedMass <- ses.mntd(samp = SeedMass$comm, dis = traits.SeedMass,  null.model = "independentswap", runs = 999)
traits.mpd.SeedMass <- ses.mpd(samp = SeedMass$comm, dis = traits.SeedMass, null.model = "independentswap", runs = 999)

pdf(file="output/8_TraitDiveristy/FunctionalDispersionSeedMass_summits.pdf") 
plotSESdispersion(mpd = traits.mpd.SeedMass, mntd = traits.mntd.SeedMass, 
                  mainTitle = "Functional Dispersion Seed Mass on Summits \nnull.model = independent swap")
#dev.off()
dev.off()




######################################## Max Height ######################################## 

#### dataset with complete data for this trait 
MaxHeight <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(2:3),], env = alps.env, 
                              traits = log10(alps.traits[complete.cases(alps.traits$PL_H_flora_MAX),]["PL_H_flora_MAX"])) 

traits.MaxHeight <- data.matrix(dist((MaxHeight$data)))
dim(traits.MaxHeight) #191

traits.mntd.MaxHeight <- ses.mntd(samp = MaxHeight$comm, dis = traits.MaxHeight, null.model = "independentswap",runs = 999)
traits.mpd.MaxHeight <- ses.mpd(samp = MaxHeight$comm, dis = traits.MaxHeight, null.model = "independentswap", runs = 999)

pdf(file="output/8_TraitDiveristy/FunctionalDispersionMaxHeight_summits.pdf") 
plotSESdispersion(mpd = traits.mpd.MaxHeight, mntd = traits.mntd.MaxHeight, 
                  mainTitle = "Functional Dispersion Seed Max Height on Summits\nnull.model = independent swap")
#dev.off()
dev.off()



######################################## Veg Disp ######################################## 
#### dataset with complete data for this trait 
cushion.data <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                            traits = alps.traits[complete.cases(alps.traits$cushions),]["cushions"]) 

#alps.VegDisp <- data.matrix(gowdis(x = vegDisp$data))
cushion <- vegdist(x = cushion.data$data, method = "bray")
dim(cushion) #1073

traits.mntd.vegDisp <- ses.mntd(samp = vegDisp$comm, dis = alps.VegDisp, runs = 999)
traits.mpd.vegDisp <- ses.mpd(samp = vegDisp$comm, dis = alps.VegDisp, runs = 999)

#pdf(file="output/8_TraitDiveristy/FunctionalDispersionVeg.pdf") 
plotSESdispersion(mpd = traits.mpd.vegDisp, mntd = traits.mntd.vegDisp, 
                  mainTitle = "Functional Dispersion Habit\nnull.model = equal prbability Ecrins")
#dev.off()


######################################## Ploidy ######################################## 
(ploidy$data)

#### dataset with complete data for this trait 
ploidy.data <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                            traits = alps.traits[complete.cases(alps.traits$diploid),]["diploid"]) #599

traits.ploidy <- vegdist(ploidy.data$data, method =  "bray")
dim(traits.ploidy) #721

traits.mntd.ploidy <- ses.mntd(samp = ploidy.data$comm, dis = vegdist(ploidy.data$data, method =  "bray"), runs = 999)
traits.mpd.ploidy <- ses.mpd(samp = ploidy.data$comm, dis = vegdist(ploidy.data$data, method =  "bray"), runs = 999)


pdf(file="output/8_TraitDiveristy/FunctionalDispersionPloidy.pdf") 
plotSESdispersion(mpd = traits.mpd.ploidy, mntd = traits.mntd.ploidy, 
                  mainTitle = "Functional Dispersion Ploidy\nnull.model = equal prbability Ecrins")

dev.off()





######################################## Woodiness ######################################## 

#### dataset with complete data for this trait 
Woody <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                           traits = alps.traits[complete.cases(alps.traits$WOODY),]["WOODY"]) #599

traits.Woody <- data.matrix(daisy(Woody$data))
dim(traits.Woody) #960

traits.mntd.Woody <- ses.mntd(samp = Woody$comm, dis = traits.Woody, runs = 999)
traits.mpd.Woody <- ses.mpd(samp = Woody$comm, dis = traits.Woody, runs = 999)


pdf(file="output/8_TraitDiveristy/FunctionalDispersionWoody.pdf") 

dev.off()




