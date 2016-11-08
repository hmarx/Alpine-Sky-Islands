##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)

source("alps.damocles.chooseClade.R")

ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############ just on poales
pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa)

obs.mntd.poales <- ses.mntd(samp =t(alps.damocles.poales$data)[c(3,10),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.poales <- ses.mpd(samp =t(alps.damocles.poales$data)[c(3,10),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.poales, file="obs.mntd.persis.poales.csv")
write.csv(obs.mpd.poales, file="obs.mpd.persis.poales.csv")

damo.bs.poales <- DAMOCLES_bootstrap(
  phy = alps.damocles.poales$phy,
  pa = alps.damocles.poales$data[,c(1, 10)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.poales$null_community_data, file="poales.persis.null_community_data.csv")
write.csv(damo.bs.poales$summary_table, file="poales.persis.summary_table.csv")

