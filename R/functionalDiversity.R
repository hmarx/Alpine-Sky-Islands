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

################################ Test for phylogenetic signal of traits in community ################################ 

traits.sig <- pezAlpes.Traits$data[pezAlpes.Traits$phy$tip.label, ]
traits.sig$PLOIDY <- as.numeric(traits.sig$PLOIDY)
traits.sig$VEG_DISP <- as.numeric(traits.sig$VEG_DISP)
traits.sig$WOODY <- as.numeric(traits.sig$WOODY)
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
alps.traits.grower <- daisy(pezAlpes.Traits$data, metric = "grower")

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

#pdf(file="output/8_TraitDiveristy/FunctionalDispersionMetricsEcrinsPool.pdf") 
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








