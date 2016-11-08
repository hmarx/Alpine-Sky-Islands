
##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

alps.phy <- read.tree(file ="phy.Alpes.taxized.tre")

alps.sites <- read.csv(file="alps.sites.csv", row.names=1, header = T)
alps.sites <- data.matrix(alps.sites)

pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites) #traits = alps.traits

## Change community matrix to presence / absence 
tmp <- pezAlpes$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa <- data.matrix(cbind("taxa"=rownames(t(tmp)), t(tmp[4:nrow(tmp),])))
alps.phy <- pezAlpes$phy

####
tax=read.csv(file="fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

pruneCladeTaxonomyLookup <- function(tip.labels, tax, level, taxonomy){
  tips.ecrins=sapply(tip.labels, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips.ecrins, rownames(tax))
  ecrins_tax=tax[ll,]
  rownames(ecrins_tax)=names(tips.ecrins)
  ecrins_tax=as.matrix(ecrins_tax)
  ecrins_tax[is.na(ecrins_tax)]=""
  head(ecrins_tax)
  #length(which(ecrins_tax[,"Angiospermae"] == "Angiospermae")) # 1064 species are in Spermatophyta
  ecrins.clade <- names(which(ecrins_tax[,level] == taxonomy))
  #return(as.data.frame(ecrins_tax))
  return(ecrins.clade)
}

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")

idparsopt = 1:2
initparsopt = c(0.1,0.1)

############ just on asterales
pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa)

damo.bs.aster <- DAMOCLES_bootstrap(
  phy = alps.damocles.aster$phy,
  pa = alps.damocles.aster$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 100,
  estimate_pars = FALSE,
  conf.int = 0.95
)

#write.csv(damo.bs.aster$null_community_data, file="null_community_data.csv")
#write.csv(damo.bs.aster$summary_table, file="summary_table.csv")

damo.output <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/null_community_data.csv")
damo.summary <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/summary_table.csv")

densityplot(damo.output$mntd.RD)
densityplot(damo.output$mntd.DAMOCLES, )

damo.output.melt <- melt(damo.output, measure.vars = c("mntd.RD", "mntd.DAMOCLES", "mpd.RD", "mpd.DAMOCLES"))
head(damo.output.melt)


ggplot(damo.output.melt, aes(x = value, colour = variable, 
                         linetype=variable)) + 
  geom_density() +
  ggtitle(label = "DAMOCLES Asterales")

pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.asterales.pdf")
ggplot(damo.output.melt, aes(x = value, fill = variable)) + 
  geom_histogram(binwidth=1) +
  ggtitle(label = "DAMOCLES Asterales")
dev.off()

#### 100 bootstrap reps

######## ASTERALES
## null communities simulated under DAMOCLES : 
damo.boots.aster <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/04asterales/summaryAsterales.csv")

(mean(as.numeric(as.character(damo.boots.aster[damo.boots.aster$desc == "mntd.obs",3]))) - 
  mean(as.numeric(as.character(damo.boots.aster[damo.boots.aster$desc == "mntd.mean.RD",3])))) / 
  mean(as.numeric(as.character(damo.boots.aster[damo.boots.aster$desc == "mntd.sd.RD",3])))


ses.damo.boots.aster <- rbind(damo.boots.aster[damo.boots.aster$desc == "mntd.obs.z.RD",], damo.boots.aster[damo.boots.aster$desc == "mpd.obs.z.RD",], 
                              damo.boots.aster[damo.boots.aster$desc == "mntd.obs.z.DAMOCLES",], damo.boots.aster[damo.boots.aster$desc == "mpd.obs.z.DAMOCLES",])
ses.damo.boots.aster$value <- as.numeric(as.character(ses.damo.boots.aster$value))

## 
aster.pez <- comparative.comm(comm = t(alps.damocles.aster$data), phy=alps.damocles.aster$phy)
aster.obs.mntd.100 <- lapply(1:100, function(x) ses.mntd(samp =aster.pez$comm[c(2,15),], dis = cophenetic(aster.pez$phy), null.model = "phylogeny.pool", runs = 999))
aster.obs.mpd.100 <- lapply(1:100, function(x) ses.mpd(samp =aster.pez$comm[c(2,15),], dis = cophenetic(aster.pez$phy), null.model = "phylogeny.pool", runs = 999))

aster.obs.mntd.100b <- lapply(1:100, function(x) ses.mntd(samp =aster.pez$comm[c(2,15),], dis = cophenetic(aster.pez$phy), null.model = "phylogeny.pool", runs = 100))
aster.obs.mpd.100b <- lapply(1:100, function(x) ses.mpd(samp =aster.pez$comm[c(2,15),], dis = cophenetic(aster.pez$phy), null.model = "phylogeny.pool", runs = 100))


aster.obs.mntd.100b.summits <- (do.call(rbind, lapply(aster.obs.mntd.100b, "[", 2, )))
hist(aster.obs.mntd.100b.summits$mntd.obs.z)
ast.mean.ses.mntd <- mean(aster.obs.mntd.100b.summits$mntd.obs.z)

aster.obs.mpd.100b.summits <- (do.call(rbind, lapply(aster.obs.mpd.100b, "[", 2, )))
hist(aster.obs.mpd.100b.summits$mpd.obs.z)
ast.mean.ses.mpd <- mean(aster.obs.mpd.100b.summits$mpd.obs.z)

ses.damo.boots.aster <- rbind(cbind(desc = rep("mntd.obs.z", nrow(aster.obs.mntd.100b.summits)), value=-1*(aster.obs.mntd.100b.summits$mntd.obs.z)),
                              damo.boots.aster[damo.boots.aster$desc == "mntd.obs.z.RD", c(2:3)],
                              damo.boots.aster[damo.boots.aster$desc == "mntd.obs.z.DAMOCLES", c(2:3)], 
                              cbind(desc =rep("mpd.obs.z", nrow(aster.obs.mpd.100b.summits)), value=-1*(aster.obs.mpd.100b.summits$mpd.obs.z)), 
                              damo.boots.aster[damo.boots.aster$desc == "mpd.obs.z.RD", c(2:3)], 
                              damo.boots.aster[damo.boots.aster$desc == "mpd.obs.z.DAMOCLES",c(2:3)])
ses.damo.boots.aster$value <- as.numeric(as.character(ses.damo.boots.aster$value))

#pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.asterales2.pdf")
ggplot(ses.damo.boots.aster, aes(x = value, fill = desc)) + 
  geom_histogram(binwidth=.05) +
  geom_vline(aes(xintercept=-1*(ast.mean.ses.mpd))) +
  geom_vline(aes(xintercept=-1*(ast.mean.ses.mntd))) +
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(0,50)) +
  ggtitle(label = "DAMOCLES Asterales SESmetrics")
#dev.off()
t.aster.mpd <- t.test(ses.damo.boots.aster[ses.damo.boots.aster$desc == "mpd.obs.z.RD", "value"], 
                      ses.damo.boots.aster[ses.damo.boots.aster$desc == "mpd.obs.z.DAMOCLES", "value"])
t.aster.mntd <- t.test(ses.damo.boots.aster[ses.damo.boots.aster$desc == "mntd.obs.z.RD", "value"], 
                       ses.damo.boots.aster[ses.damo.boots.aster$desc == "mntd.obs.z.DAMOCLES", "value"])
pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.asterales.boxplot.v2.pdf")
ggplot(ses.damo.boots.aster, aes(x = desc, y=value, fill = desc)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(-2,2.5)) +
  geom_hline(aes(yintercept=0), lty=3) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())  +
  #geom_text(x=1.5, y=2.5, label = "***") +
  geom_hline(aes(yintercept=-1*(ast.mean.ses.mpd)), lty=2, col="violet") +
  geom_hline(aes(yintercept=-1*(ast.mean.ses.mntd)), lty=2, col="coral") +
  ggtitle(label = "DAMOCLES Asterales SESmetrics")
dev.off()

######### POALES
damo.boots.poales <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/04poales/summaryPoalses.csv")
ses.damo.boots.poales <- rbind(damo.boots.poales[damo.boots.poales$desc == "mntd.obs.z.RD",], damo.boots.poales[damo.boots.poales$desc == "mpd.obs.z.RD",], 
                              damo.boots.poales[damo.boots.poales$desc == "mntd.obs.z.DAMOCLES",], damo.boots.poales[damo.boots.poales$desc == "mpd.obs.z.DAMOCLES",])
ses.damo.boots.poales$value <- as.numeric(as.character(ses.damo.boots.poales$value))

poales.pez <- comparative.comm(comm = t(alps.damocles.poales$data), phy=alps.damocles.poales$phy)
poales.obs.mntd.100 <- lapply(1:100, function(x) ses.mntd(samp =poales.pez$comm[c(2,15),], dis = cophenetic(poales.pez$phy), null.model = "phylogeny.pool", runs = 999))
poales.obs.mpd.100 <- lapply(1:100, function(x) ses.mpd(samp =poales.pez$comm[c(2,15),], dis = cophenetic(poales.pez$phy), null.model = "phylogeny.pool", runs = 999))

poales.obs.mntd.100.summits <- (do.call(rbind, lapply(poales.obs.mntd.100, "[", 2, )))
hist(poales.obs.mntd.100.summits$mntd.obs.z)
poales.mean.ses.mntd <- mean(poales.obs.mntd.100.summits$mntd.obs.z)

poales.obs.mpd.100.summits <- (do.call(rbind, lapply(poales.obs.mpd.100, "[", 2, )))
hist(poales.obs.mpd.100.summits$mpd.obs.z)
poales.mean.ses.mpd <- mean(poales.obs.mpd.100.summits$mpd.obs.z)

pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.poales.pdf")
ggplot(ses.damo.boots.poales, aes(x = value, fill = desc)) + 
  geom_histogram(binwidth=.05) +
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(0,50)) +
  ggtitle(label = "DAMOCLES poales SESmetrics")
dev.off()
t.poales.mpd <- t.test(ses.damo.boots.poales[ses.damo.boots.poales$desc == "mpd.obs.z.RD", "value"], 
                      ses.damo.boots.poales[ses.damo.boots.poales$desc == "mpd.obs.z.DAMOCLES", "value"])
t.poales.mntd <- t.test(ses.damo.boots.poales[ses.damo.boots.poales$desc == "mntd.obs.z.RD", "value"], 
                       ses.damo.boots.poales[ses.damo.boots.poales$desc == "mntd.obs.z.DAMOCLES", "value"])
pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.poales.boxplot.pdf")
ggplot(ses.damo.boots.poales, aes(x = desc, y=value, fill = desc)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(-2,2.5)) +
  geom_hline(aes(yintercept=0), lty=3) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())  +
  geom_hline(aes(yintercept=-1*(poales.mean.ses.mpd)), lty=2, col="violet") +
  geom_hline(aes(yintercept=-1*(poales.mean.ses.mntd)), lty=2, col="coral") +
  #geom_text(x=1.5, y=2.5, label = "***") +
  ggtitle(label = "DAMOCLES Poales SESmetrics")
dev.off()

####### LAMIALES
damo.boots.lamiales <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/06lamiales/summaryLamiales.csv")
ses.damo.boots.lamiales <- rbind(damo.boots.lamiales[damo.boots.lamiales$desc == "mntd.obs.z.RD",], damo.boots.lamiales[damo.boots.lamiales$desc == "mpd.obs.z.RD",], 
                               damo.boots.lamiales[damo.boots.lamiales$desc == "mntd.obs.z.DAMOCLES",], damo.boots.lamiales[damo.boots.lamiales$desc == "mpd.obs.z.DAMOCLES",])
ses.damo.boots.lamiales$value <- as.numeric(as.character(ses.damo.boots.lamiales$value))

lamiales.pez <- comparative.comm(comm = t(alps.damocles.lamiales$data), phy=alps.damocles.lamiales$phy)
lamiales.obs.mntd.100 <- lapply(1:100, function(x) ses.mntd(samp =lamiales.pez$comm[c(2,15),], dis = cophenetic(lamiales.pez$phy), null.model = "phylogeny.pool", runs = 999))
lamiales.obs.mpd.100 <- lapply(1:100, function(x) ses.mpd(samp =lamiales.pez$comm[c(2,15),], dis = cophenetic(lamiales.pez$phy), null.model = "phylogeny.pool", runs = 999))

lamiales.obs.mntd.100.summits <- (do.call(rbind, lapply(lamiales.obs.mntd.100, "[", 2, )))
hist(lamiales.obs.mntd.100.summits$mntd.obs.z)
lamiales.mean.ses.mntd <- mean(lamiales.obs.mntd.100.summits$mntd.obs.z)

lamiales.obs.mpd.100.summits <- (do.call(rbind, lapply(lamiales.obs.mpd.100, "[", 2, )))
hist(lamiales.obs.mpd.100.summits$mpd.obs.z)
lamiales.mean.ses.mpd <- mean(lamiales.obs.mpd.100.summits$mpd.obs.z)

pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.lamiales.pdf")
ggplot(ses.damo.boots.lamiales, aes(x = value, fill = desc)) + 
  geom_histogram(binwidth=.05) +
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(0,50)) +
  ggtitle(label = "DAMOCLES lamiales SESmetrics")
dev.off()
t.lamiales.mpd <- t.test(ses.damo.boots.lamiales[ses.damo.boots.lamiales$desc == "mpd.obs.z.RD", "value"], 
                       ses.damo.boots.lamiales[ses.damo.boots.lamiales$desc == "mpd.obs.z.DAMOCLES", "value"])
t.lamiales.mntd <- t.test(ses.damo.boots.lamiales[ses.damo.boots.lamiales$desc == "mntd.obs.z.RD", "value"], 
                        ses.damo.boots.lamiales[ses.damo.boots.lamiales$desc == "mntd.obs.z.DAMOCLES", "value"])
pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.lamiales.boxplot.pdf")
ggplot(ses.damo.boots.lamiales, aes(x = desc, y=value, fill = desc)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(-2,2.5)) +
  geom_hline(aes(yintercept=0), lty=3) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())  +
  geom_hline(aes(yintercept=-1*(lamiales.mean.ses.mpd)), lty=2, col="violet") +
  geom_hline(aes(yintercept=-1*(lamiales.mean.ses.mntd)), lty=2, col="coral") +
  geom_text(x=1.5, y=2.5, label = "***") +
  ggtitle(label = "DAMOCLES Lamiales SESmetrics")
dev.off()

##### CARYOPHYLLALALES
damo.boots.caryophyllales <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/clades/06cary/summaryCary.csv")
ses.damo.boots.caryophyllales <- rbind(damo.boots.caryophyllales[damo.boots.caryophyllales$desc == "mntd.obs.z.RD",], damo.boots.caryophyllales[damo.boots.caryophyllales$desc == "mpd.obs.z.RD",], 
                                 damo.boots.caryophyllales[damo.boots.caryophyllales$desc == "mntd.obs.z.DAMOCLES",], damo.boots.caryophyllales[damo.boots.caryophyllales$desc == "mpd.obs.z.DAMOCLES",])
ses.damo.boots.caryophyllales$value <- as.numeric(as.character(ses.damo.boots.caryophyllales$value))

caryoph.pez <- comparative.comm(comm = t(alps.damocles.caryophyllales$data), phy=alps.damocles.caryophyllales$phy)
caryoph.obs.mntd.100 <- lapply(1:100, function(x) ses.mntd(samp =caryoph.pez$comm[c(2,15),], dis = cophenetic(caryoph.pez$phy), null.model = "phylogeny.pool", runs = 999))
caryoph.obs.mpd.100 <- lapply(1:100, function(x) ses.mpd(samp =caryoph.pez$comm[c(2,15),], dis = cophenetic(caryoph.pez$phy), null.model = "phylogeny.pool", runs = 999))

caryoph.obs.mntd.100.summits <- (do.call(rbind, lapply(caryoph.obs.mntd.100, "[", 2, )))
hist(caryoph.obs.mntd.100.summits$mntd.obs.z)
caryoph.mean.ses.mntd <- mean(caryoph.obs.mntd.100.summits$mntd.obs.z)

caryoph.obs.mpd.100.summits <- (do.call(rbind, lapply(caryoph.obs.mpd.100, "[", 2, )))
hist(caryoph.obs.mpd.100.summits$mpd.obs.z)
caryoph.mean.ses.mpd <- mean(caryoph.obs.mpd.100.summits$mpd.obs.z)


pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.caryophyllales.pdf")
ggplot(ses.damo.boots.caryophyllales, aes(x = value, fill = desc)) + 
  geom_histogram(binwidth=.05) +
  scale_x_continuous(limits=c(-3,3)) +
  scale_y_continuous(limits=c(0,50)) +
  ggtitle(label = "DAMOCLES caryophyllales SESmetrics")
dev.off()
t.caryophyllales.mpd <- t.test(ses.damo.boots.caryophyllales[ses.damo.boots.caryophyllales$desc == "mpd.obs.z.RD", "value"], 
                         ses.damo.boots.caryophyllales[ses.damo.boots.caryophyllales$desc == "mpd.obs.z.DAMOCLES", "value"])
t.caryophyllales.mntd <- t.test(ses.damo.boots.caryophyllales[ses.damo.boots.caryophyllales$desc == "mntd.obs.z.RD", "value"], 
                          ses.damo.boots.caryophyllales[ses.damo.boots.caryophyllales$desc == "mntd.obs.z.DAMOCLES", "value"])
pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/DAMOCLES.SES.caryophyllales.boxplot.pdf")
ggplot(ses.damo.boots.caryophyllales, aes(x = desc, y=value, fill = desc)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(-2,2.5)) +
  geom_hline(aes(yintercept=0), lty=3) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())  +
  geom_text(x=1.5, y=0.25, label = "***") +
  geom_text(x=3.5, y=0.25, label = "***") +
  ggtitle(label = "DAMOCLES caryophyllales SESmetrics")
dev.off()


