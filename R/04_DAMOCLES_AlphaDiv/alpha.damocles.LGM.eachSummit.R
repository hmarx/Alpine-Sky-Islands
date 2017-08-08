#####################################################################################################################
############# Dynamic null model of communtiy assembly ############################################################## 
############# Each order within each summit community, 'LGM' species pool ###########################################
############# Hannah E. Marx, 26 Feb 2016 ###########################################################################
#####################################################################################################################

################################ Loop DAMOCLES over all summits
args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

source("R/functions/chooseClade.R")

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



############ just on asterales in summit source pool
### NOTE: i = 8 never converges to optimized parameters...skiped to i = 9 (each clade run as separate code)

pruned.tree.alps.aster <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.persistent.pa)

#for (i in 2:ncol(alps.damocles.aster$data))

damo.bs.aster <- DAMOCLES_bootstrap(
  phy = alps.damocles.aster$phy,
  pa = alps.damocles.aster$data[,c(1, i)],  #LGM P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.aster$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Asterales/asterales.persistent.null_community_data.", i, ".csv", sep=""))
write.csv(damo.bs.aster$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Asterales/asterales.persistent.summary_table.", i, ".csv", sep=""))

############ just on poales
pruned.tree.alps.poales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.persistent.pa)

#for (i in 2:ncol(alps.damocles.poales$data))
damo.bs.poales <- DAMOCLES_bootstrap(
  phy = alps.damocles.poales$phy,
  pa = alps.damocles.poales$data[,c(1, i)],  #LGM P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95)

write.csv(damo.bs.poales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Poales/poales.persistent.null_community_data.", i, ".csv", sep=""))
write.csv(damo.bs.poales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Poales/poales.persistent.summary_table.", i, ".csv", sep=""))


############ just on lamiales
pruned.tree.alps.lamiales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.persistent.pa)

#for (i in 2:ncol(alps.damocles.lamiales$data))
damo.bs.lamiales  <- DAMOCLES_bootstrap(
  phy = alps.damocles.lamiales$phy,
  pa = alps.damocles.lamiales$data[,c(1, i)],  #LGM P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95)

write.csv(damo.bs.lamiales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Poales/lamiales.persistent.null_community_data.", i, ".csv", sep=""))
write.csv(damo.bs.lamiales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Poales/lamiales.persistent.summary_table.", i, ".csv", sep=""))

############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.persistent.pa)


#for (i in 2:ncol(alps.damocles.caryophyllales$data))

damo.bs.caryophyllales <- DAMOCLES_bootstrap(
  phy = alps.damocles.caryophyllales$phy,
  pa = alps.damocles.caryophyllales$data[,c(1, i)],  #LGM P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95)

write.csv(damo.bs.caryophyllales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Caryo/caryophyllales.persistent.null_community_data.", i, ".csv", sep=""))
write.csv(damo.bs.caryophyllales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Caryo/caryophyllales.persistent.summary_table.", i, ".csv", sep=""))

############ just on spermatophyta in persistent source pool
pruned.tree.alps.sperma <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.spermatophyta])
alps.damocles.sperma <- treedata(pruned.tree.alps.sperma, alps.persistent.pa)

#for (i in 2:ncol(alps.damocles.sperma$data)) #2:16

damo.bs.sperma <- DAMOCLES_bootstrap(phy = alps.damocles.sperma$phy, pa = alps.damocles.sperma$data[,c(1, i)], idparsopt = 1:length(initparsopt), parsfix = 0, idparsfix = (1:3)[-idparsopt], pars2 = c(1E-3,1E-4,1E-5,1000), pchoice = 0, runs = 1000, estimate_pars = FALSE, conf.int = 0.95)

write.csv(damo.bs.sperma$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Spermatophyta/spermatophyta.persistent.null_community_data.", i, ".csv", sep=""))
write.csv(damo.bs.sperma$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/3_LGM/Spermatophyta/spermatophyta.persistent.summary_table.", i, ".csv", sep=""))


