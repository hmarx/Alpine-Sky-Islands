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
pezAlpes.Traits$data

pezAlpes.Traits$data$PL_H_flora_MAX <- log(pezAlpes.Traits$data$PL_H_flora_MAX)
pezAlpes.Traits$data$SEEDM <- log(pezAlpes.Traits$data$SEEDM)
pezAlpes.Traits$data$PLOIDY <- as.numeric(pezAlpes.Traits$data$PLOIDY)
pezAlpes.Traits$data$VEG_DISP <- as.numeric(pezAlpes.Traits$data$VEG_DISP)
pezAlpes.Traits$data$WOODY <- as.numeric(pezAlpes.Traits$data$WOODY)

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
alps.traits.grower <- daisy(pezAlpes.Traits$data, metric = "gower")

alps.trait.matrix.grower <- data.matrix(alps.traits.grower)
dim(alps.trait.matrix.grower) #599

traits.mntd <- ses.mntd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, runs = 999)
traits.mpd <- ses.mpd(samp = pezAlpes.Traits$comm, dis = alps.trait.matrix.grower, runs = 999)


traits.mntd[,"sig"] <- ifelse(traits.mntd$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd[,"sig1"] <- ifelse(traits.mpd$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd[c(1,6,9)], traits.mpd[c(1,6,9)], id=rownames(traits.mntd))
trait.dispersion <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion)

trait.dispersion <- arrange(trait.dispersion, ntaxa)
trait.dispersion$summit <- factor(trait.dispersion$summit, levels = unique(trait.dispersion$summit))

#pdf(file="output/8_TraitDiveristy/FunctionalDispersionMetricsEcrinsPoolLog.pdf") 
ggplot(trait.dispersion, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion Metrics\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
#dev.off()


######################################## Veg Disp ######################################## 

#### dataset with complete data for this trait 
vegDisp <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                            traits = alps.traits[complete.cases(alps.traits$VEG_DISP),]["VEG_DISP"]) 

alps.VegDisp <- data.matrix(daisy(vegDisp$data))
dim(alps.VegDisp) #1073

traits.mntd.vegDisp <- ses.mntd(samp = vegDisp$comm, dis = alps.VegDisp, runs = 999)
traits.mpd.vegDisp <- ses.mpd(samp = vegDisp$comm, dis = alps.VegDisp, runs = 999)

traits.mntd.vegDisp[,"sig"] <- ifelse(traits.mntd.vegDisp$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd.vegDisp[,"sig1"] <- ifelse(traits.mpd.vegDisp$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd.vegDisp[c(1,6,9)], traits.mpd.vegDisp[c(1,6,9)], id=rownames(traits.mntd.vegDisp))
trait.dispersion.veg <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                            varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                            times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion.veg) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion.veg)

trait.dispersion.veg <- arrange(trait.dispersion.veg, ntaxa)
trait.dispersion.veg$summit <- factor(trait.dispersion.veg$summit, levels = unique(trait.dispersion.veg$summit))

pdf(file="output/8_TraitDiveristy/FunctionalDispersionVeg.pdf") 
ggplot(trait.dispersion.veg, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion of Vegitative Habit\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()


######################################## Ploidy ######################################## 

#### dataset with complete data for this trait 
ploidy <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                            traits = alps.traits[complete.cases(alps.traits$PLOIDY),]["PLOIDY"]) #599

traits.ploidy <- data.matrix(daisy(ploidy$data))
dim(traits.ploidy) #721

traits.mntd.ploidy <- ses.mntd(samp = ploidy$comm, dis = traits.ploidy, runs = 999)
traits.mpd.ploidy <- ses.mpd(samp = ploidy$comm, dis = traits.ploidy, runs = 999)

traits.mntd.ploidy[,"sig"] <- ifelse(traits.mntd.ploidy$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd.ploidy[,"sig1"] <- ifelse(traits.mpd.ploidy$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd.ploidy[c(1,6,9)], traits.mpd.ploidy[c(1,6,9)], id=rownames(traits.mntd.ploidy))
trait.dispersion.ploidy <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion.ploidy) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion.ploidy)

trait.dispersion.ploidy <- arrange(trait.dispersion.ploidy, ntaxa)
trait.dispersion.ploidy$summit <- factor(trait.dispersion.ploidy$summit, levels = unique(trait.dispersion.ploidy$summit))

pdf(file="output/8_TraitDiveristy/FunctionalDispersionPloidy.pdf") 
ggplot(trait.dispersion.ploidy, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion of Ploidy\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()





######################################## Woodiness ######################################## 

#### dataset with complete data for this trait 
Woody <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                           traits = alps.traits[complete.cases(alps.traits$WOODY),]["WOODY"]) #599

traits.Woody <- data.matrix(daisy(Woody$data))
dim(traits.Woody) #960

traits.mntd.Woody <- ses.mntd(samp = Woody$comm, dis = traits.Woody, runs = 999)
traits.mpd.Woody <- ses.mpd(samp = Woody$comm, dis = traits.Woody, runs = 999)

traits.mntd.Woody[,"sig"] <- ifelse(traits.mntd.Woody$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd.Woody[,"sig1"] <- ifelse(traits.mpd.Woody$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd.Woody[c(1,6,9)], traits.mpd.Woody[c(1,6,9)], id=rownames(traits.mntd.Woody))
trait.dispersion.Woody <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                   varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                   times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion.Woody) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion.Woody)

trait.dispersion.Woody <- arrange(trait.dispersion.Woody, ntaxa)
trait.dispersion.Woody$summit <- factor(trait.dispersion.Woody$summit, levels = unique(trait.dispersion.Woody$summit))

pdf(file="output/8_TraitDiveristy/FunctionalDispersionWoody.pdf") 
ggplot(trait.dispersion.Woody, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion of Woody\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()




######################################## Seed Mass ######################################## 

#### dataset with complete data for this trait 
SeedMass <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                          traits = alps.traits[complete.cases(alps.traits$SEEDM),]["SEEDM"]) 

traits.SeedMass <- data.matrix(dist(log10(SeedMass$data)))
dim(traits.SeedMass) #734

traits.mntd.SeedMass <- ses.mntd(samp = SeedMass$comm, dis = traits.SeedMass, null.model = "taxa.labels", runs = 999)
traits.mpd.SeedMass <- ses.mpd(samp = SeedMass$comm, dis = traits.SeedMass, null.model = "taxa.labels", runs = 999)

traits.mntd.SeedMass[,"sig"] <- ifelse(traits.mntd.SeedMass$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd.SeedMass[,"sig1"] <- ifelse(traits.mpd.SeedMass$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd.SeedMass[c(1,6,9)], traits.mpd.SeedMass[c(1,6,9)], id=rownames(traits.mntd.SeedMass))
trait.dispersion.SeedMass <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                  varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                  times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion.SeedMass) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion.SeedMass)

trait.dispersion.SeedMass <- arrange(trait.dispersion.SeedMass, ntaxa)
trait.dispersion.SeedMass$summit <- factor(trait.dispersion.SeedMass$summit, levels = unique(trait.dispersion.SeedMass$summit))

pdf(file="output/8_TraitDiveristy/FunctionalDispersionSeedMass.pdf") 
ggplot(trait.dispersion.SeedMass, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion of SeedMass\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()




######################################## Max Height ######################################## 

#### dataset with complete data for this trait 
MaxHeight <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, 
                             traits = alps.traits[complete.cases(alps.traits$PL_H_flora_MAX),]["PL_H_flora_MAX"]) 

traits.MaxHeight <- data.matrix(dist(log10(MaxHeight$data)))
dim(traits.MaxHeight) #940

traits.mntd.MaxHeight <- ses.mntd(samp = MaxHeight$comm, dis = traits.MaxHeight, null.model = "taxa.labels", runs = 999)
traits.mpd.MaxHeight <- ses.mpd(samp = MaxHeight$comm, dis = traits.MaxHeight, null.model = "taxa.labels", runs = 999)

traits.mntd.MaxHeight[,"sig"] <- ifelse(traits.mntd.MaxHeight$mntd.obs.p <= 0.05, "TRUE", "FALSE")
traits.mpd.MaxHeight[,"sig1"] <- ifelse(traits.mpd.MaxHeight$mpd.obs.p <= 0.05, "TRUE", "FALSE")
tmp <- cbind(traits.mntd.MaxHeight[c(1,6,9)], traits.mpd.MaxHeight[c(1,6,9)], id=rownames(traits.mntd.MaxHeight))
trait.dispersion.MaxHeight <- reshape(tmp, direction = "long", idvar = c("id", "ntaxa"), 
                                     varying = list(c("mntd.obs.z", "mpd.obs.z"), c("sig", "sig1")), 
                                     times = c("mntd.obs.z", "mpd.obs.z"), v.names = c("metric", "value"))
colnames (trait.dispersion.MaxHeight) <- c("ntaxa", "ntaxa.1", "summit", "metric", "value", "sig")
head(trait.dispersion.MaxHeight)

trait.dispersion.MaxHeight <- arrange(trait.dispersion.MaxHeight, ntaxa)
trait.dispersion.MaxHeight$summit <- factor(trait.dispersion.MaxHeight$summit, levels = unique(trait.dispersion.MaxHeight$summit))

pdf(file="output/8_TraitDiveristy/FunctionalDispersionMaxHeight.pdf") 
ggplot(trait.dispersion.MaxHeight, aes(x=summit, y=value)) +
  geom_point(aes(colour = metric, fill= factor(sig), shape=factor(metric)), size = 4) +
  scale_y_continuous("obs.z")  +
  scale_x_discrete("Summit (increasing species richness)") +
  scale_fill_manual(values=c("FALSE"= "white", "TRUE" = "grey"), guide="none") +
  scale_color_manual(values=c("chocolate3", "cyan4")) +
  scale_shape_manual(name=" ", values=c("mntd.obs.z" = 21, "mpd.obs.z" =24),  guide="none") +
  theme_bw() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ggtitle("Functional Dispersion of MaxHeight\nEcrins Species Pool") +
  theme(plot.title=element_text(size=rel(1.5))) 
dev.off()






