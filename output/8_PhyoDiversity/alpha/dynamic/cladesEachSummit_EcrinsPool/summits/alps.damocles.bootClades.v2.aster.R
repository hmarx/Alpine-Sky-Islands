##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)

source("alps.damocles.chooseClade.R")

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############ just on asterales
pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa)

obs.mntd.aster <- ses.mntd(samp =t(alps.damocles.aster$data)[c(3,16),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.aster <- ses.mpd(samp =t(alps.damocles.aster$data)[c(3,16),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.aster, file="obs.mntd.aster.csv")
write.csv(obs.mpd.aster, file="obs.mpd.aster.csv")

damo.bs.aster <- DAMOCLES_bootstrap(
  phy = alps.damocles.aster$phy,
  pa = alps.damocles.aster$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.aster$null_community_data, file="asterales.null_community_data.csv")
write.csv(damo.bs.aster$summary_table, file="asterales.summary_table.csv")


