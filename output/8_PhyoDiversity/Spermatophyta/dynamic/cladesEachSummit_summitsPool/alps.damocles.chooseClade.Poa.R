
## Loop DAMOCLES over all summits

##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)

source("alps.damocles.chooseClade.R")

ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

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
  
  write.csv(damo.bs.poales$null_community_data, file=paste("poales.null_community_data.", i, ".csv", sep=""))
  write.csv(damo.bs.poales$summary_table, file=paste("poales.summary_table.", i, ".csv", sep=""))
  
}

