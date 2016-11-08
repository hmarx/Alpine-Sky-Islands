
###### Beta diveristy of Ecrins National Park Summits 
## Source pool = Ecrins NP
## Clade = Spermatophyta

source("beta_pd_decomp.r")
## Spermatophyta
alps.phy <- read.tree(file ="phy.Alpes.taxized.tre")

alps.sites.df <- read.csv(file="alps.sites.csv", row.names=1, header=T)

decompo_beta <- beta.pd.decompo(com = alps.sites.df, tree=alps.phy, type="both", output.dist=T, random=999)

#write.csv(decompo_beta_df$betadiv, file = "decompo_beta_betadiv3.csv")
          
#write.csv(decompo_beta_df$util.pd, file = "decompo_beta_utilpd3.csv")

#### Distance Matricies take forever to run
write.csv(data.matrix(decompo_beta$betadiv$PhyloSor), file= "dist_PhloSor.csv")

write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_turn), file="dist_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$PhyloSor_PD), file="dist_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac), file = "dist_UniFrac.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac_turn), file = "dist_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$UniFrac_PD), file= "dist_UniFrac_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor), file = "SES_PhyloSor.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_turn), file= "SES_PhyloSor_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_PhyloSor_PD), file = " SES_PhyloSor_PD.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac), file = "SES_UniFrac.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_turn), file = "SES_UniFrac_turn.csv")

write.csv(data.matrix(decompo_beta$betadiv$SES_UniFrac_PD), file = " SES_UniFrac_PD.csv")


