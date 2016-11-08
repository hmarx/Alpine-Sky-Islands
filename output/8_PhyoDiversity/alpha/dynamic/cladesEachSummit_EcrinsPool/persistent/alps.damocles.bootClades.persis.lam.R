##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)

source("alps.damocles.chooseClade.R")

ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############ just on lamiales
pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa)


obs.mntd.lamiales <- ses.mntd(samp =t(alps.damocles.lamiales$data)[c(3,10),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.lamiales <- ses.mpd(samp =t(alps.damocles.lamiales$data)[c(3,10),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.lamiales, file="obs.mntd.persis.lamiales.csv")
write.csv(obs.mpd.lamiales, file="obs.mpd.persis.lamiales.csv")

damo.bs.lamiales <- DAMOCLES_bootstrap(
  phy = alps.damocles.lamiales$phy,
  pa = alps.damocles.lamiales$data[,c(1, 10)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.lamiales$null_community_data, file="lamiales.persis.null_community_data.csv")
write.csv(damo.bs.lamiales$summary_table, file="lamiales.persis.summary_table.csv")

