
###### Beta diveristy of Ecrins National Park Summits 
## Source pool = Ecrins NP
## Clade = Spermatophyta

source("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/beta_pd_decomp.r")
source("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/alps.damocles.chooseClade.ecrins.R")


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




