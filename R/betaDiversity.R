#####################################################################################################################
############# Phylogenetic diversity between alpine summits (beta) ################################################## 
############# Decomposed beta diveristy (PhyloSor) to get 'true' turnover (independent of Species Richness) #########
############# Species pool = Ecrins NP ##############################################################################
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################

source("R/beta_pd_decomp.r")

alps.phy <- read.tree(file ="data/AnalysesDatasets/phy.Alpes.taxized.tre")

alps.sites.df <- read.csv(file="data/AnalysesDatasets/alps.sites.csv", row.names=1, header=T)

###### Spermatophyta
decompo_beta <- beta.pd.decompo(com = alps.sites.df, tree=alps.phy, type="both", output.dist=T, random=999)

#write.csv(decompo_beta_df$betadiv, file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/Dryad_decompo_beta_betadiv.csv")
#write.csv(decompo_beta_df$util.pd, file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/decompo_beta_utilpd.csv")

#### Distance Matricies take forever to run

#Observed 
write.csv(data.matrix(decompo_beta$betadiv$PhyloSor), file= "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_PhloSor.csv")
write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_turn), file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_PhyloSor_turn.csv")
write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_PD), file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_PhyloSor_PD.csv")
write.csv(data.matrix(decompo_beta$betadiv$UniFrac), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_UniFrac.csv")
write.csv(data.matrix(decompo_beta$betadiv$UniFrac_turn), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta$betadiv$UniFrac_PD), file= "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_observed/dist_UniFrac_PD.csv")

#SES
write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_PhyloSor.csv")
write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_turn), file= "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_PhyloSor_turn.csv")
write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_PD), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_PhyloSor_PD.csv")
write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_UniFrac.csv")
write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_turn), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_UniFrac_turn.csv")
write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_PD), file = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_UniFrac_PD.csv")

########### Specific clades ########### 

## Clade = Asterales
betaDecompClade(phy= alps.phy, sites = alps.sites.pa, clade = "Asterales")

## Clade = Caryphyllales
betaDecompClade(phy= alps.phy, sites = alps.sites.pa, clade = "Caryophyllales")

## Clade = Lamiales
betaDecompClade(phy= alps.phy, sites = alps.sites.pa, clade = "Lamiales")

## Clade = Poales
betaDecompClade(phy= alps.phy, sites = alps.sites.pa, clade = "Poales")

## Clade = Rosales
betaDecompClade(phy= alps.phy, sites = alps.sites.pa, clade = "Rosales")

#### Calculate Beta decomposed for clade:
betaDecompClade <- function(phy, sites, clade){
  ecrins.clade <- pruneCladeTaxonomyLookup(tip.labels = phy$tip.label, tax, level = "order", taxonomy = clade)
  pruned.tree.alps.clade <-drop.tip(phy, phy$tip.label[!phy$tip.label %in% ecrins.clade])
  alps.damocles.clade <- treedata(pruned.tree.alps.clade, sites)
  pruned.sites.alps.clade <- as.data.frame(t(alps.damocles.clade$data[,-1]), strings.as.factors=F)
  pruned.sites.alps.clade[] <- apply(pruned.sites.alps.clade, 2, function(x) as.numeric(as.character(x)))
  
  decompo_beta <- beta.pd.decompo(com = pruned.sites.alps.clade, tree=pruned.tree.alps.clade, type="both", output.dist=T, random=999)
  
  write.csv(data.matrix(decompo_beta$betadiv$PhyloSor), file= paste(clade, "dist_PhloSor.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_turn), file=paste(clade, "dist_PhyloSor_turn.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_PD), file=paste(clade, "dist_PhyloSor_PD.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$UniFrac), file = paste(clade, "dist_UniFrac.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$UniFrac_turn), file = paste(clade, "dist_UniFrac_turn.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$UniFrac_PD), file= paste(clade, "dist_UniFrac_PD.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor), file = paste(clade, "SES_PhyloSor.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_turn), file= paste(clade, "SES_PhyloSor_turn.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_PD), file = paste(clade, "SES_PhyloSor_PD.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac), file = paste(clade, "SES_UniFrac.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_turn), file = paste(clade, "SES_UniFrac_turn.csv", sep="_"))
  write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_PD), file = paste(clade, "SES_UniFrac_PD.csv", sep="_"))
  
  #decompo_beta_df <- beta.pd.decompo(com = pruned.sites.alps.clade, tree=pruned.tree.alps.clade, type="both", output.dist=F, random=999)
  #write.csv(decompo_beta_df$betadiv, file = paste(clade, "decompo_beta_betadiv.csv", sep="_"))
  #write.csv(decompo_beta_df$util.pd, file = paste(clade, "decompo_beta_utilpd.csv", sep="_"))
  
}




