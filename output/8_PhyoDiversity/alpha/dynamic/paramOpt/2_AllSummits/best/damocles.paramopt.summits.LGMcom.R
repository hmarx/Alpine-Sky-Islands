#####################################################################################################################
############# Explore parameter space for rates of colonization (gamma) and local extinction (mu) for dynamic DAMOCLES null model ###################################################### 
############# Community of persistent (LGM) species within the ‘Regional’ Écrins NP species pool
############# Hannah E. Marx, 07 Dec 2015 ###########################################################################
#####################################################################################################################

## Dynamic null model of communtiy assembly in the Ecrin NP, France
#args <- commandArgs(TRUE)

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(dplyr)

source("R/chooseClade.R") #generates alps.summits.pa

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.rosales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Rosales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")


## Repeat optimization of parameters 10 times to see if they converge 
# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs

for (i in 1:30){
  initparsopt <- sample(c(1:100), size=2, replace = T)
  
  ## Spermatophyta
  alps.damocles <- treedata(pezAlpes.summits$phy, alps.sites.pa)
  
  opt.params.spermatophyta <- DAMOCLES_ML(phy = alps.damocles$phy,
                                          pa = alps.damocles$data[,c(1, 10)], ## Just LGM (Persistent) community
                                          initparsopt = initparsopt,
                                          idparsopt = 1:2,
                                          parsfix = 0,
                                          idparsfix = 3,
                                          pars2 = c(1E-3,1E-4,1E-5,1000),
                                          pchoice = 0)
  
  write.csv(opt.params.spermatophyta, file=paste("output/opt.params.spermatophyta.", i, ".csv", sep=""))
  
  
  ############ just on asterales
  pruned.tree.alps.aster <-drop.tip(pezAlpes.summits$phy, pezAlpes.summits$phy$tip.label[!pezAlpes.summits$phy$tip.label %in% ecrins.asterales])
  alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa)
  
  opt.params.asterales <- DAMOCLES_ML(phy = alps.damocles.aster$phy,
                                      pa = alps.damocles.aster$data[,c(1, 10)],
                                      initparsopt = initparsopt,
                                      idparsopt = 1:2,
                                      parsfix = 0,
                                      idparsfix = 3,
                                      pars2 = c(1E-3,1E-4,1E-5,1000),
                                      pchoice = 0)
  
  write.csv(opt.params.asterales, file=paste("output/opt.params.asterales.", i, ".csv", sep=""))
  
  ############ just on poales
  pruned.tree.alps.poales <-drop.tip(pezAlpes.summits$phy, pezAlpes.summits$phy$tip.label[!pezAlpes.summits$phy$tip.label %in% ecrins.poales])
  alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa)
  
  opt.params.poales <- DAMOCLES_ML(phy = alps.damocles.poales$phy,
                                   pa = alps.damocles.poales$data[,c(1, 10)],
                                   initparsopt = initparsopt,
                                   idparsopt = 1:2,
                                   parsfix = 0,
                                   idparsfix = 3,
                                   pars2 = c(1E-3,1E-4,1E-5,1000),
                                   pchoice = 0)
  
  write.csv(opt.params.poales, file=paste("output/opt.params.poales.", i, ".csv", sep=""))
  
  
  ############ just on rosales
  pruned.tree.alps.rosales <-drop.tip(pezAlpes.summits$phy, pezAlpes.summits$phy$tip.label[!pezAlpes.summits$phy$tip.label %in% ecrins.rosales])
  alps.damocles.rosales <- treedata(pruned.tree.alps.rosales, alps.sites.pa)
  
  opt.params.rosales <- DAMOCLES_ML(phy = alps.damocles.rosales$phy,
                                    pa = alps.damocles.rosales$data[,c(1, 10)],
                                    initparsopt = initparsopt,
                                    idparsopt = 1:2,
                                    parsfix = 0,
                                    idparsfix = 3,
                                    pars2 = c(1E-3,1E-4,1E-5,1000),
                                    pchoice = 0)
  
  write.csv(opt.params.rosales, file=paste("output/opt.params.rosales.", i, ".csv", sep=""))
  
  
  ############ just on lamiales
  pruned.tree.alps.lamiales <-drop.tip(pezAlpes.summits$phy, pezAlpes.summits$phy$tip.label[!pezAlpes.summits$phy$tip.label %in% ecrins.lamiales])
  alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa)
  
  opt.params.lamiales <- DAMOCLES_ML(phy = alps.damocles.lamiales$phy,
                                     pa = alps.damocles.lamiales$data[,c(1, 10)],
                                     initparsopt = initparsopt,
                                     idparsopt = 1:2,
                                     parsfix = 0,
                                     idparsfix = 3,
                                     pars2 = c(1E-3,1E-4,1E-5,1000),
                                     pchoice = 0)
  
  write.csv(opt.params.lamiales, file=paste("output/opt.params.lamiales.", i, ".csv", sep=""))
  
  
  ############ just on caryophyllales
  pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes.summits$phy, pezAlpes.summits$phy$tip.label[!pezAlpes.summits$phy$tip.label %in% ecrins.caryophyllales])
  alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.sites.pa)
  
  opt.params.caryophyllales <- DAMOCLES_ML(phy = alps.damocles.caryophyllales$phy,
                                           pa = alps.damocles.caryophyllales$data[,c(1, 10)],
                                           initparsopt = initparsopt,
                                           idparsopt = 1:2,
                                           parsfix = 0,
                                           idparsfix = 3,
                                           pars2 = c(1E-3,1E-4,1E-5,1000),
                                           pchoice = 0)
  
  write.csv(opt.params.caryophyllales, file=paste("output/opt.params.caryophyllales.", i, ".csv", sep=""))
  
}  





