#####################################################################################################################
############# Explore parameter space for rates of colonization (gamma) and local extinction (mu) for dynamic DAMOCLES null model ###################################################### 
############# Community of persistent (LGM) species within the ‘Regional’ Écrins NP species pool
############# Hannah E. Marx, 07 Dec 2015 ###########################################################################
#####################################################################################################################

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(dplyr)

## Taxonomy lookup table
tax=read.csv(file="data/AnalysesDatasets/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

pruneCladeTaxonomyLookup <- function(tip.labels, tax, level, taxonomy){
  tips.ecrins=sapply(tip.labels, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips.ecrins, rownames(tax))
  ecrins_tax=tax[ll,]
  rownames(ecrins_tax)=names(tips.ecrins)
  ecrins_tax=as.matrix(ecrins_tax)
  ecrins_tax[is.na(ecrins_tax)]=""
  head(ecrins_tax)
  #length(which(ecrins_tax[,"Angiospermae"] == "Angiospermae")) # 1064 species are in Spermatophyta
  ecrins.clade <- names(which(ecrins_tax[,level] == taxonomy))
  #return(as.data.frame(ecrins_tax))
  return(ecrins.clade)
}



############### Loop parameter optimization over all 100 bootstrap trees

for (i in 1:100){
  
  ### read tree
  alps.phy <- read.tree(file = paste("data/AnalysesDatasets/boots/scale", i, "dated.rename.tre", sep="."))
  roottaxa <- c("Pinus_sylvestris", "Pinus_nigra", "Pinus_mugo", "Picea_abies", "Larix_decidua", "Abies_alba", "Cedrus_atlantica", "Juniperus_communis", "Juniperus_sabina", "Taxus_baccata")
  alps.phy.root <- root(alps.phy, roottaxa, resolve.root = TRUE)
  
  alps.sites <- read.csv(file="data/AnalysesDatasets/alps.sites.csv", row.names=1, header = T)
  alps.sites <- data.matrix(alps.sites)
  
  pezAlpes <- comparative.comm(phy = alps.phy.root, comm = alps.sites) #traits = alps.traits
  
  ###
  
  ## REGIONAL community matrix to presence / absence 
  tmp <- pezAlpes$comm
  tmp[which(tmp != 0)] <- 1
  alps.sites.pa <- data.matrix(cbind("taxa"=rownames(t(tmp)), t(tmp[4:nrow(tmp),])))
  alps.phy <- pezAlpes$phy
  
  ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
  ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
  ecrins.rosales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Rosales")
  ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
  ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
  
  
  ## Repeat optimization of parameters 10 times to see if they converge 
  # inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
  # == 227.2309 * 1/10 = 22 / mil yrs
  
  for (j in 1:30){
    initparsopt <- sample(c(1:100), size=2, replace = T)
    
    ## Spermatophyta
    alps.damocles <- treedata(pezAlpes$phy, alps.sites.pa)
    
    opt.params.spermatophyta <- DAMOCLES_ML(phy = alps.damocles$phy,
                                            pa = alps.damocles$data[,c(1, 16)], ## Just summits community
                                            initparsopt = initparsopt,
                                            idparsopt = 1:2,
                                            parsfix = 0,
                                            idparsfix = 3,
                                            pars2 = c(1E-3,1E-4,1E-5,1000),
                                            pchoice = 0)
    
    write.csv(opt.params.spermatophyta, file=paste("output/scale", i, "opt.params.spermatophyta", j, "csv", sep="."))
    
    
    ############ just on asterales
    pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
    alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa)
    
    opt.params.asterales <- DAMOCLES_ML(phy = alps.damocles.aster$phy,
                                        pa = alps.damocles.aster$data[,c(1, 16)],
                                        initparsopt = initparsopt,
                                        idparsopt = 1:2,
                                        parsfix = 0,
                                        idparsfix = 3,
                                        pars2 = c(1E-3,1E-4,1E-5,1000),
                                        pchoice = 0)
    
    write.csv(opt.params.asterales, file=paste("output/scale", i, "opt.params.asterales", j, "csv", sep="."))
    
    ############ just on poales
    pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
    alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa)
    
    opt.params.poales <- DAMOCLES_ML(phy = alps.damocles.poales$phy,
                                     pa = alps.damocles.poales$data[,c(1, 16)],
                                     initparsopt = initparsopt,
                                     idparsopt = 1:2,
                                     parsfix = 0,
                                     idparsfix = 3,
                                     pars2 = c(1E-3,1E-4,1E-5,1000),
                                     pchoice = 0)
    
    write.csv(opt.params.poales, file=paste("output/scale", i, "opt.params.poales", j, "csv", sep="."))
    
    
    ############ just on rosales
    pruned.tree.alps.rosales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.rosales])
    alps.damocles.rosales <- treedata(pruned.tree.alps.rosales, alps.sites.pa)
    
    opt.params.rosales <- DAMOCLES_ML(phy = alps.damocles.rosales$phy,
                                      pa = alps.damocles.rosales$data[,c(1, 16)],
                                      initparsopt = initparsopt,
                                      idparsopt = 1:2,
                                      parsfix = 0,
                                      idparsfix = 3,
                                      pars2 = c(1E-3,1E-4,1E-5,1000),
                                      pchoice = 0)
    
    write.csv(opt.params.rosales, file=paste("output/scale", i, "opt.params.rosales", j, "csv", sep="."))
    
    
    ############ just on lamiales
    pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
    alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa)
    
    opt.params.lamiales <- DAMOCLES_ML(phy = alps.damocles.lamiales$phy,
                                       pa = alps.damocles.lamiales$data[,c(1, 16)],
                                       initparsopt = initparsopt,
                                       idparsopt = 1:2,
                                       parsfix = 0,
                                       idparsfix = 3,
                                       pars2 = c(1E-3,1E-4,1E-5,1000),
                                       pchoice = 0)
    
    write.csv(opt.params.lamiales, file=paste("output/scale", i, "opt.params.lamiales", j, "csv", sep="."))
    
    
    ############ just on caryophyllales
    pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
    alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.sites.pa)
    
    opt.params.caryophyllales <- DAMOCLES_ML(phy = alps.damocles.caryophyllales$phy,
                                             pa = alps.damocles.caryophyllales$data[,c(1, 16)],
                                             initparsopt = initparsopt,
                                             idparsopt = 1:2,
                                             parsfix = 0,
                                             idparsfix = 3,
                                             pars2 = c(1E-3,1E-4,1E-5,1000),
                                             pchoice = 0)
    
    write.csv(opt.params.caryophyllales, file=paste("output/scale", i, "opt.params.caryophyllales", j, "csv", sep="."))
    
  }  
  
  
}

