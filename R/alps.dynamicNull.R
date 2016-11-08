
## Dynamic null model of communtiy assembly in the Ecrin NP, France

install.packages("DAMOCLES")
require(DAMOCLES)
data(NWPrimates_data)

phy = rcoal(10)
pa = matrix(c(phy$tip.label,sample(c(0,1),Ntip(phy),replace = T)),
            nrow = Ntip(phy),ncol = 2)

demo.bs <- DAMOCLES_bootstrap(
  phy = phy,
  pa = pa,
  initparsopt = c(0.1,0.1),
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 2,
  estimate_pars = TRUE,
  conf.int = 0.95
)

demo.bs <- DAMOCLES_bootstrap(
  phy = rcoal(10),
  pa = matrix(c(phy$tip.label,sample(c(0,1),Ntip(phy),replace = T)),nrow = Ntip(phy),ncol = 2),
  initparsopt = c(0.1,0.1),
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 9,
  estimate_pars = TRUE,
  conf.int = 0.95
)

plot(demo.bs$null_community_data$gamma_0.DAMOCLES, demo.bs$null_community_data$mu.DAMOCLES)
lines(3.76491777424178, 2.26989532387427, col = "red")



###################################################### Test on Alps dataset
tmp <- pezAlpes$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa <- data.matrix(cbind("taxa"=rownames(t(tmp)), t(tmp[4:nrow(tmp),])))
dim(alps.sites.pa)
head(alps.sites.pa)
alps.phy <- pezAlpes$phy

alps.damocles <- treedata(pezAlpes$phy, alps.sites.pa)

DAMOCLES_ML(phy = alps.damocles$phy, pa = alps.damocles$data[,c(1,16)])

##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)

#computes the maximum likelihood estimates of colonisation and local extinction
#rate for a given phylogeny and presence-absence data under the DAMOCLES model
damo.bs <- DAMOCLES_bootstrap(
  phy = alps.damocles$phy,
  pa = alps.damocles$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 2,
  estimate_pars = TRUE,
  conf.int = 0.95
)

write.csv(damo.bs$null_community_data, file="null_community_data.csv")
write.table(damo.bs$summary_table, file="summary_table")


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

write.csv(damo.bs.aster$null_community_data, file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/null_community_data.csv")
write.csv(damo.bs.aster$summary_table, file="output/9_PhyoDiversity/Spermatophyta/dynamic/clades/summary_table.csv")

