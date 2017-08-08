#####################################################################################################################
############# Dynamic null model of communtiy assembly ############################################################## 
############# Each order within the LGM community, All Summits Species Pool  ########################################
############# Hannah E. Marx, 11 Feb 2017 ###########################################################################
#####################################################################################################################

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)
require(dplyr)

source("R/chooseClade.R")  #generates alps.summits.pa.all

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
ecrins.spermatophyta <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "Spermatophyta", taxonomy = "Spermatophyta")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############# Spermatophyta order within the LGM species pool ###########################################################
pruned.tree.alps.sperma <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.spermatophyta])
alps.damocles.sperma <- treedata(pruned.tree.alps.sperma, alps.summits.pa.all)

# Random Draw Null 
obs.mntd.sperma <- ses.mntd(samp =t(alps.damocles.sperma$data)[c(12,18),], dis = cophenetic(alps.damocles.sperma$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.sperma <- ses.mpd(samp =t(alps.damocles.sperma$data)[c(12,18),], dis = cophenetic(alps.damocles.sperma$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.sperma, file="output/obs.mntd.persis.sperma.csv")
write.csv(obs.mpd.sperma, file="output/obs.mpd.persis.sperma.csv")

# DAMOCLES null 
damo.bs.sperma <- DAMOCLES_bootstrap(
  phy = alps.damocles.sperma$phy,
  pa = alps.damocles.sperma$data[,c(1, 12)],  # LGM (persistent) species pool only
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.sperma$null_community_data, file="output/Spermatophyta.persis.null_community_data.csv")
write.csv(damo.bs.sperma$summary_table, file="output/Spermatophyta.persis.summary_table.csv")



############# Asterales order within the LGM species pool ###########################################################
pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.summits.pa.all)

# Random Draw Null 
obs.mntd.aster <- ses.mntd(samp =t(alps.damocles.aster$data)[c(12,18),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.aster <- ses.mpd(samp =t(alps.damocles.aster$data)[c(12,18),], dis = cophenetic(alps.damocles.aster$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.aster, file="output/obs.mntd.persis.aster.csv")
write.csv(obs.mpd.aster, file="output/obs.mpd.persis.aster.csv")

# DAMOCLES null 
damo.bs.aster <- DAMOCLES_bootstrap(
  phy = alps.damocles.aster$phy,
  pa = alps.damocles.aster$data[,c(1, 12)],  # LGM (persistent) species pool only
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.aster$null_community_data, file="output/asterales.persis.null_community_data.csv")
write.csv(damo.bs.aster$summary_table, file="output/asterales.persis.summary_table.csv")



############# Poales order within the LGM species pool ######################################################
pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.summits.pa.all)

# Random Draw Null 
obs.mntd.poales <- ses.mntd(samp =t(alps.damocles.poales$data)[c(12,18),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.poales <- ses.mpd(samp =t(alps.damocles.poales$data)[c(12,18),], dis = cophenetic(alps.damocles.poales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.poales, file="output/obs.mntd.persis.poales.csv")
write.csv(obs.mpd.poales, file="output/obs.mpd.persis.poales.csv")

# DAMOCLES null 
damo.bs.poales <- DAMOCLES_bootstrap(
  phy = alps.damocles.poales$phy,
  pa = alps.damocles.poales$data[,c(1, 12)],  # LGM (persistent) species pool only
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.poales$null_community_data, file="output/poales.persis.null_community_data.csv")
write.csv(damo.bs.poales$summary_table, file="output/poales.persis.summary_table.csv")



############# Lamiales order within the LGM species pool ######################################################
pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.summits.pa.all)

# Random Draw Null 
obs.mntd.lamiales <- ses.mntd(samp =t(alps.damocles.lamiales$data)[c(12,18),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.lamiales <- ses.mpd(samp =t(alps.damocles.lamiales$data)[c(12,18),], dis = cophenetic(alps.damocles.lamiales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.lamiales, file="output/obs.mntd.persis.lamiales.csv")
write.csv(obs.mpd.lamiales, file="output/obs.mpd.persis.lamiales.csv")

# DAMOCLES null 
damo.bs.lamiales <- DAMOCLES_bootstrap(
  phy = alps.damocles.lamiales$phy,
  pa = alps.damocles.lamiales$data[,c(1, 12)],   # LGM (persistent) species pool only
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.lamiales$null_community_data, file="output/lamiales.persis.null_community_data.csv")
write.csv(damo.bs.lamiales$summary_table, file="output/lamiales.persis.summary_table.csv")


############# Caryophyllales order within the LGM species pool ######################################################
pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.summits.pa.all)

# Random Draw Null 
obs.mntd.caryophyllales <- ses.mntd(samp =t(alps.damocles.caryophyllales$data)[c(12,18),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.caryophyllales <- ses.mpd(samp =t(alps.damocles.caryophyllales$data)[c(12,18),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.caryophyllales, file="output/obs.mntd.persis.caryophyllales.csv")
write.csv(obs.mpd.caryophyllales, file="output/obs.mpd.persis.caryophyllales.csv")

# DAMOCLES null 
damo.bs.caryophyllales <- DAMOCLES_bootstrap(
  phy = alps.damocles.caryophyllales$phy,
  pa = alps.damocles.caryophyllales$data[,c(1, 12)],   # LGM (persistent) species pool only
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.caryophyllales$null_community_data, file="output/caryophyllales.persis.null_community_data.csv")
write.csv(damo.bs.caryophyllales$summary_table, file="output/caryophyllales.persis.summary_table.csv")




