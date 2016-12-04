#####################################################################################################################
############# Dynamic null model of communtiy assembly ############################################################## 
############# Each order within each summit community, 'All Summits' species pool ###################################
############# Hannah E. Marx, 16 Feb 2016 ###########################################################################
#####################################################################################################################

############################## Loop DAMOCLES across all summits
require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

source("R/chooseClade.R")

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)



############ just on asterales in summit source pool
### NOTE: i = 8 never converges to optimized parameters...skiped to i = 9

pruned.tree.alps.aster <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.summits.pa)


for (i in 9:ncol(alps.damocles.aster$data)){
  damo.bs.aster <- DAMOCLES_bootstrap(
    phy = alps.damocles.aster$phy,
    pa = alps.damocles.aster$data[,c(1, i)],  #summits P/A,
    idparsopt = 1:length(initparsopt),
    parsfix = 0,
    idparsfix = (1:3)[-idparsopt],
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    runs = 1000,
    estimate_pars = FALSE,
    conf.int = 0.95
  )
  
  write.csv(damo.bs.aster$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Asterales/asterales.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.aster$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Asterales/asterales.summary_table.", i, ".csv", sep=""))

}


############ just on poales
pruned.tree.alps.poales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.summits.pa)

for (i in 2:ncol(alps.damocles.poales$data)){
  damo.bs.poales <- DAMOCLES_bootstrap(
    phy = alps.damocles.poales$phy,
    pa = alps.damocles.poales$data[,c(1, i)],  #summits P/A,
    idparsopt = 1:length(initparsopt),
    parsfix = 0,
    idparsfix = (1:3)[-idparsopt],
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    runs = 1000,
    estimate_pars = FALSE,
    conf.int = 0.95)
  
  write.csv(damo.bs.poales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Poales/poales.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.poales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Poales/poales.summary_table.", i, ".csv", sep=""))
  
}

############ just on lamiales
pruned.tree.alps.lamiales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.summits.pa)


for (i in 2:ncol(alps.damocles.lamiales$data)){
  damo.bs.lamiales <- DAMOCLES_bootstrap(
    phy = alps.damocles.lamiales$phy,
    pa = alps.damocles.lamiales$data[,c(1, i)],  #summits P/A,
    idparsopt = 1:length(initparsopt),
    parsfix = 0,
    idparsfix = (1:3)[-idparsopt],
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    runs = 1000,
    estimate_pars = FALSE,
    conf.int = 0.95)
  
  write.csv(damo.bs.lamiales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Lamiales/lamiales.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.lamiales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Lamiales/lamiales.summary_table.", i, ".csv", sep=""))
  
}    


############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.summits.pa)


for (i in 2:ncol(alps.damocles.caryophyllales$data)){
  
  damo.bs.caryophyllales <- DAMOCLES_bootstrap(
    phy = alps.damocles.caryophyllales$phy,
    pa = alps.damocles.caryophyllales$data[,c(1, i)],  #summits P/A,
    idparsopt = 1:length(initparsopt),
    parsfix = 0,
    idparsfix = (1:3)[-idparsopt],
    pars2 = c(1E-3,1E-4,1E-5,1000),
    pchoice = 0,
    runs = 1000,
    estimate_pars = FALSE,
    conf.int = 0.95
  )
  
  write.csv(damo.bs.caryophyllales$null_community_data, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Caryo/caryophyllales.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.caryophyllales$summary_table, file=paste("output/8_PhyoDiversity/alpha/dynamic/speciesPools/2_allSummits/Caryo/caryophyllales.summary_table.", i, ".csv", sep=""))
  
}

