########################### Beta diveristy (turnover) of Ecrins Alpine Summits for Clades ########################### 
#### 15 March, 2016 H. Marx
##########################################################################################################
source("R/chooseClade.ecrinsPool.R")
source("R/spacodiWrap.R")
source("R/plotBetaPair.R")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Asterales", 
            output = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/clades")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Caryophyllales", 
            output = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/clades")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Lamiales", 
            output = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/clades")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Poales", 
            output = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/clades")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Rosales", 
            output = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/clades")

##################################
beta.Sperma <- plotBetaPair(clade = "Spermatophyta", 
                           outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/")

beta.Aster <- plotBetaPair(clade = "Asterales", 
             outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Asterales/")

beta.Caryo <-plotBetaPair(clade = "Caryophyllales", 
             outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Caryophyllales/")

beta.Lam <-plotBetaPair(clade = "Lamiales", 
             outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Lamiales/")

beta.Poa <-plotBetaPair(clade = "Poales", 
             outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Poales/")

beta.Ros <-plotBetaPair(clade = "Rosales", 
             outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Rosales/")

write.csv(rbind(beta.Sperma, beta.Aster, beta.Caryo, beta.Lam, beta.Poa, beta.Ros), file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/decomposedBetaClades.csv")



