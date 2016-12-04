#####################################################################################################################
############# Environmental drivers of phylogenetic beta diversity patterns #########################################
############# multiple regression (MRM) on SES beta diversity (decomposed)  #########################################
############# Hannah E. Marx, 15 Mar 2016 ###########################################################################
#####################################################################################################################

source("R/envCategory.R")
source("R/betaDivEnvSummary.R") 


############################################### Geographic : Spatial Distance ############################################### 
spatial.dist <- vegdist(na.omit(pezAlpes$env[c("X_WGS84", "Y_WGS84")]), method = "euclid")

MRM_geog_dist_sperm <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_PhyloSor.csv",
                                      geog.dist <- spatial.dist,
                                      nperm = 1000,
                                      clade = "Spermatophyta",
                                      independent = "geog_dist",
                                      dependent = "PhyloSorTurn")

MRM_geog_dist_aster <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Asterales/distance_SES/Asterales_SES_PhyloSor.csv",
                                      geog.dist <- spatial.dist,
                                      nperm = 1000,
                                      clade = "Asterales",
                                      independent = "geog_dist",
                                      dependent = "PhyloSorTurn")

MRM_geog_dist_cary <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Caryophyllales/distance_SES/Caryophyllales_SES_PhyloSor.csv",
                                     geog.dist <- spatial.dist,
                                     nperm = 1000,
                                     clade = "Caryophyllales",
                                     independent = "geog_dist",
                                     dependent = "PhyloSorTurn")

MRM_geog_dist_lam <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Lamiales/distance_SES/Lamiales_SES_PhyloSor.csv",
                                    geog.dist <- spatial.dist,
                                    nperm = 1000,
                                    clade = "Lamiales",
                                    independent = "geog_dist",
                                    dependent = "PhyloSorTurn")

MRM_geog_dist_pao <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Poales/distance_SES/Poales_SES_PhyloSor.csv",
                                    geog.dist <- spatial.dist,
                                    nperm = 1000,
                                    clade = "Poales",
                                    independent = "geog_dist",
                                    dependent = "PhyloSorTurn")

MRM_geog_dist_ros <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Rosales/distance_SES/Rosales_SES_PhyloSor.csv",
                                    geog.dist <- spatial.dist,
                                    nperm = 1000,
                                    clade = "Rosales",
                                    independent = "geog_dist",
                                    dependent = "PhyloSorTurn")

MRM_geog_dist <- rbind(MRM_geog_dist_sperm, MRM_geog_dist_aster, MRM_geog_dist_cary, MRM_geog_dist_lam, MRM_geog_dist_pao, MRM_geog_dist_ros)

############################################### Geographic : Historical Connectivity ############################################### 


############################################# Climatic: Available Energy ############################################# 
edispc.availener.score <- vegdist(pc.availener.score, method = "euclid")

MRM_env_energy_sperm <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_PhyloSor.csv",
                                       geog.dist <- edispc.availener.score,
                                       nperm = 1000,
                                       clade = "Spermatophyta",
                                       independent = "env_energy",
                                       dependent = "PhyloSorTurn")

MRM_env_energy_aster <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Asterales/distance_SES/Asterales_SES_PhyloSor.csv",
                                       geog.dist <- edispc.availener.score,
                                       nperm = 1000,
                                       clade = "Asterales",
                                       independent = "env_energy",
                                       dependent = "PhyloSorTurn")

MRM_env_energy_cary <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Caryophyllales/distance_SES/Caryophyllales_SES_PhyloSor.csv",
                                      geog.dist <- edispc.availener.score,
                                      nperm = 1000,
                                      clade = "Caryophyllales",
                                      independent = "env_energy",
                                      dependent = "PhyloSorTurn")

MRM_env_energy_lam <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Lamiales/distance_SES/Lamiales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.availener.score,
                                     nperm = 1000,
                                     clade = "Lamiales",
                                     independent = "env_energy",
                                     dependent = "PhyloSorTurn")

MRM_env_energy_pao <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Poales/distance_SES/Poales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.availener.score,
                                     nperm = 1000,
                                     clade = "Poales",
                                     independent = "env_energy",
                                     dependent = "PhyloSorTurn")

MRM_env_energy_ros <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Rosales/distance_SES/Rosales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.availener.score,
                                     nperm = 1000,
                                     clade = "Rosales",
                                     independent = "env_energy",
                                     dependent = "PhyloSorTurn")

MRM_env_energy <- rbind(MRM_env_energy_sperm, MRM_env_energy_aster, MRM_env_energy_cary, MRM_env_energy_lam, MRM_env_energy_pao, MRM_env_energy_ros)

############################################# Climatic: Environmental stress ############################################# 
edispc.pc.stress.score <- vegdist(pc.stress.score, method = "euclid")
MRM_env_stress_sperm <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_PhyloSor.csv",
                                       geog.dist <- edispc.pc.stress.score,
                                       nperm = 1000,
                                       clade = "Spermatophyta",
                                       independent = "env_stress",
                                       dependent = "PhyloSorTurn")

MRM_env_stress_aster <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Asterales/distance_SES/Asterales_SES_PhyloSor.csv",
                                       geog.dist <- edispc.pc.stress.score,
                                       nperm = 1000,
                                       clade = "Asterales",
                                       independent = "env_stress",
                                       dependent = "PhyloSorTurn")

MRM_env_stress_cary <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Caryophyllales/distance_SES/Caryophyllales_SES_PhyloSor.csv",
                                      geog.dist <- edispc.pc.stress.score,
                                      nperm = 1000,
                                      clade = "Caryophyllales",
                                      independent = "env_stress",
                                      dependent = "PhyloSorTurn")

MRM_env_stress_lam <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Lamiales/distance_SES/Lamiales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stress.score,
                                     nperm = 1000,
                                     clade = "Lamiales",
                                     independent = "env_stress",
                                     dependent = "PhyloSorTurn")

MRM_env_stress_pao <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Poales/distance_SES/Poales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stress.score,
                                     nperm = 1000,
                                     clade = "Poales",
                                     independent = "env_stress",
                                     dependent = "PhyloSorTurn")

MRM_env_stress_ros <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Rosales/distance_SES/Rosales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stress.score,
                                     nperm = 1000,
                                     clade = "Rosales",
                                     independent = "env_stress",
                                     dependent = "PhyloSorTurn")

MRM_env_stress <- rbind(MRM_env_stress_sperm, MRM_env_stress_aster, MRM_env_stress_cary, MRM_env_stress_lam, MRM_env_stress_pao, MRM_env_stress_ros)


############################################# Climatic: Environmental stability ############################################# 
edispc.pc.stabil.score <- vegdist(pc.stabil.score, method = "euclid")

MRM_env_stabil_sperm <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/SES_PhyloSor.csv",
                                       geog.dist <- edispc.pc.stabil.score,
                                       nperm = 1000,
                                       clade = "Spermatophyta",
                                       independent = "env_stabil",
                                       dependent = "PhyloSorTurn")

MRM_env_stabil_aster <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Asterales/distance_SES/Asterales_SES_PhyloSor.csv",
                                       geog.dist <- edispc.pc.stabil.score,
                                       nperm = 1000,
                                       clade = "Asterales",
                                       independent = "env_stabil",
                                       dependent = "PhyloSorTurn")

MRM_env_stabil_cary <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Caryophyllales/distance_SES/Caryophyllales_SES_PhyloSor.csv",
                                      geog.dist <- edispc.pc.stabil.score,
                                      nperm = 1000,
                                      clade = "Caryophyllales",
                                      independent = "env_stabil",
                                      dependent = "PhyloSorTurn")

MRM_env_stabil_lam <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Lamiales/distance_SES/Lamiales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stabil.score,
                                     nperm = 1000,
                                     clade = "Lamiales",
                                     independent = "env_stabil",
                                     dependent = "PhyloSorTurn")

MRM_env_stabil_pao <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Poales/distance_SES/Poales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stabil.score,
                                     nperm = 1000,
                                     clade = "Poales",
                                     independent = "env_stabil",
                                     dependent = "PhyloSorTurn")

MRM_env_stabil_ros <- betaDivSummary(betaOutput = "output/8_PhyoDiversity/beta/decomposedBeta/clades/Rosales/distance_SES/Rosales_SES_PhyloSor.csv",
                                     geog.dist <- edispc.pc.stabil.score,
                                     nperm = 1000,
                                     clade = "Rosales",
                                     independent = "env_stabil",
                                     dependent = "PhyloSorTurn")

MRM_env_stabil <- rbind(MRM_env_stabil_sperm, MRM_env_stabil_aster, MRM_env_stabil_cary, MRM_env_stabil_lam, MRM_env_stabil_pao, MRM_env_stabil_ros)


########################## Combine output:

write.csv(rbind(MRM_geog_dist, MRM_env_energy, MRM_env_stress, MRM_env_stabil), file="output/9_Environment/beta/PhyloSor_env.csv")  #MRM_env_hetero = MRM_env_hetero$r.squared



