
source("R/sesPlotFunctions.R")

###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
## Contemporary species pool = Summits 
## Hisoric source pool = Ecrins NP
plotSESdispersion(mntd = comm.sesmntd, mpd = comm.sesmpd.phylonull, 
                  mainTitle = "Dispersion Metrics Ecrins NP\nnull.model = equal probability draw from Ecrins Pool")
#dev.off()

###################################################################################### 
############################### Reduced Pool : Summits ############################### 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/DispersionMetricsSummitsPool.pdf") 
plotSESdispersion(mpd = summits.sesmpd.phylonull, mntd = summits.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = equal probability draw from Summit Pool")
#dev.off()
 
###################################################################################### 
############################### Reduced Pool : Persistent ############################ 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/DispersionMetrics_persistentPool.pdf") 
plotSESdispersion(mpd = persistent.sesmpd.phylonull, mntd = persistent.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Persistent\nnull.model = equal probability draw from Persistent Pool")
#dev.off()


###################################################################################### 
########################   Summarize Reduced Source pool effects ########################   
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/mpd.sourcepools.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull, UnderIce.sesmpd.phylonull), 
                    breaks=c("1", "2", "3", "4"), labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice"),
                    mainTitle = "Distribution of SES mpd for different source pools", values=c(1,1,1,1))
#dev.off()

#pdf(file="output/9_PhyoDiversity/mntd.sourcepools.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull, UnderIce.sesmntd.phylonull), 
                    breaks=c("1", "2", "3", "4"), labels=c("Ecrins NP", "Summits", "Peristent", "Under Ice"),
                    mainTitle = "Distribution of SES mntd for different source pools", values=c(1,1,1,1))
#dev.off()




########################################################################################################### 
################# Distribution of SES under reduced source pools and  sourcePool Null ##################### 
########################################################################################################### 

#pdf(file="output/9_PhyoDiversity/mpd.all.sourcepools.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull, UnderIce.sesmpd.phylonull, 
                                     comm.sesmpd.sourcePoolNull, comm.sesmpd.sourcePoolNull.pers, comm.sesmpd.sourcePoolNull.under), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7"),
                    labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                             "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull"),
                    mainTitle = "Distribution of SES mpd for different source pools", values=c(1,1,1,1,5,5,5))
#def.off()

#pdf(file="output/9_PhyoDiversity/mntd.all.sourcepools.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull, UnderIce.sesmntd.phylonull,
                                     comm.sesmntd.sourcePoolNull, comm.sesmntd.sourcePoolNull.pers, 
                                     comm.sesmntd.sourcePoolNull.under), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7"),
                    labels=c("Ecrins NP_phyNull", "Summits_phyNull", "Peristent_phyNull", "Under Ice_phyNull", 
                             "Summits_sourceNull", "Peristent_sourceNull", "Under Ice_sourceNull"),
                    mainTitle = "Distribution of SES mntd for different source pools", values=c(1,1,1,1,5,5,5))
#dev.off()


###################################################################################### 
############################### Contemporary Pool: Sumits  ###########################
###################################################################################### 

## Contemporary species pool = summits 
## Hisoric source pool = persistent abouve glacier through LGM
plotSESdispersion(mpd = summits.sesmpd.sourcePersis,mntd = summits.sesmntd.sourcePersis, 
                  mainTitle = "Dispersion Metrics Persistent above LGM\nnull.model = weighted draw from summits source pool")


## Contemporary species pool = summits 
## Hisoric source pool = summits (persistent + under ice?)
plotSESdispersion(mpd = summits.sesmpd.sourceSummits, mntd = summits.sesmntd.sourceSUmmits, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = weighted draw from summits source pool")

#pdf(file="output/9_PhyoDiversity/summitPoolnullSourceDistribMPD.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, summits.sesmpd.sourceSummits, summits.sesmpd.sourcePersis), 
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Summits, weighted for abundance", "Persistent, weighted for abundance"),
                    mainTitle = "Distribution of SES mpd for Summit Species Pool", values=c(1,1,5,5))
#dev.off()

#pdf(file="output/9_PhyoDiversity/summitPoolnullSourceDistribMNTD.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, summits.sesmntd.sourceSUmmits, summits.sesmntd.sourcePersis), 
                    breaks=c("1", "2", "3", "4"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Summits, weighted for abundance", "Persistent, weighted for abundance"),
                    mainTitle = "Distribution of SES mntd for Summit Species Pool", values=c(1,1,5,5))
#dev.off()


###################################################################################### 
############################### Contemporary Pool: Persistent ########################
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/persistentPoolnullSourceDistribMPD.pdf") 
plotDistributionSES(outputSES = list(persistent.sesmpd.phylonull, persistent.sesmpd.sourceSummits, persistent.sesmpd.sourcePersis), 
                    breaks=c("1", "2", "3"),
                    labels=c("Persistent, equal probability", "Summits, weighted for abundance", "Persistent, weighted for abundance"),
                    mainTitle = "Distribution of SES mpd for Persistent Species Pool", values=c(1,5,5))
#dev.off()

#pdf(file="output/9_PhyoDiversity/persistentPoolnullSourceDistribMNTD.pdf") 
plotDistributionSES(outputSES = list(persistent.sesmntd.phylonull, persistent.sesmntd.sourceSummits, persistent.sesmntd.sourcePersis), 
                    breaks=c("1", "2", "3"),
                    labels=c("Persistent, equal probability", "Summits, weighted for abundance", "Persistent, weighted for abundance"),
                    mainTitle = "Distribution of SES mntd for Persistent Species Pool", values=c(1,5,5))
#dev.off()




