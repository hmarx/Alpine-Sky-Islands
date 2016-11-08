##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)

source("alps.damocles.chooseClade.R")

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
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
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.aster$null_community_data, file="asterales.null_community_data.csv")
write.csv(damo.bs.aster$summary_table, file="asterales.summary_table.csv")

obs.mntd.aster <- ses.mntd(samp =t(alps.damocles.aster$data)[c(3,16),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.aster <- ses.mpd(samp =t(alps.damocles.aster$data)[c(3,16),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.aster, file="obs.mntd.aster.csv")
write.csv(obs.mpd.aster, file="obs.mpd.aster.csv")


############ just on poales
pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa)

damo.bs.poales <- DAMOCLES_bootstrap(
  phy = alps.damocles.poales$phy,
  pa = alps.damocles.poales$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.poales$null_community_data, file="poales.null_community_data.csv")
write.csv(damo.bs.poales$summary_table, file="poales.summary_table.csv")

obs.mntd.poales <- ses.mntd(samp =t(alps.damocles.poales$data)[c(3,16),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.poales <- ses.mpd(samp =t(alps.damocles.poales$data)[c(3,16),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.poales, file="obs.mntd.poales.csv")
write.csv(obs.mpd.poales, file="obs.mpd.poales.csv")


############ just on lamiales
pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa)

damo.bs.lamiales <- DAMOCLES_bootstrap(
  phy = alps.damocles.lamiales$phy,
  pa = alps.damocles.lamiales$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.lamiales$null_community_data, file="lamiales.null_community_data.csv")
write.csv(damo.bs.lamiales$summary_table, file="lamiales.summary_table.csv")

obs.mntd.lamiales <- ses.mntd(samp =t(alps.damocles.lamiales$data)[c(3,16),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.lamiales <- ses.mpd(samp =t(alps.damocles.lamiales$data)[c(3,16),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.lamiales, file="obs.mntd.lamiales.csv")
write.csv(obs.mpd.lamiales, file="obs.mpd.lamiales.csv")


############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.sites.pa)

damo.bs.caryophyllales <- DAMOCLES_bootstrap(
  phy = alps.damocles.caryophyllales$phy,
  pa = alps.damocles.caryophyllales$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.caryophyllales$null_community_data, file="caryophyllales.null_community_data.csv")
write.csv(damo.bs.caryophyllales$summary_table, file="caryophyllales.summary_table.csv")

obs.mntd.caryophyllales <- ses.mntd(samp =t(alps.damocles.caryophyllales$data)[c(3,16),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.caryophyllales <- ses.mpd(samp =t(alps.damocles.caryophyllales$data)[c(3,16),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.caryophyllales, file="obs.mntd.caryophyllales.csv")
write.csv(obs.mpd.caryophyllales, file="obs.mpd.caryophyllales.csv")

          


