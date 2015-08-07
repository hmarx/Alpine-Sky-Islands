
source("R/sesPlotFunctions.R")

###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/DispersionMetricsEcrinsPool.pdf") 
plotSESdispersion(mntd = comm.sesmntd, mpd = comm.sesmpd.phylonull, 
                  mainTitle = "Dispersion Metrics Ecrins NP\nnull.model = equal probability draw from Ecrins Pool")
#dev.off()

###################################################################################### 
############################### Reduced Pool : Summits ############################### 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/summits_sourceEcrins.pdf") 
plotSESdispersion(mpd = summits.sesmpd.phylonull, mntd = summits.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = equal probability draw from Summit Pool")
#dev.off()
 
###################################################################################### 
############################### Reduced Pool : Persistent ############################ 
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/persistent_sourceEcrins.pdf") 
plotSESdispersion(mpd = persistent.sesmpd.phylonull, mntd = persistent.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Persistent\nnull.model = equal probability draw from Persistent Pool")
#dev.off()








###################################################################################### 
############################### Contemporary Pool: Sumits  ###########################
###################################################################################### 

## Contemporary species pool = summits 
## Hisoric source pool = persistent abouve glacier through LGM
#pdf(file="output/9_PhyoDiversity/summits_sourcePersistent.pdf") 
plotSESdispersion(mpd = summits.sesmpd.sourcePersis,mntd = summits.sesmntd.sourcePersis, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = weighted draw from persistent source pool")
#dev.off()

## Contemporary species pool = summits 
## Hisoric source pool = summits (persistent + under ice?)
#pdf(file="output/9_PhyoDiversity/summits_sourceSummits.pdf") 
plotSESdispersion(mpd = summits.sesmpd.sourceSummits, mntd = summits.sesmntd.sourceSUmmits, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = weighted draw from summits source pool")
#dev.off()

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



###################################################################################### 
############################### ALL ########################
###################################################################################### 

pdf(file="output/9_PhyoDiversity/mpd_distribution.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull,  persistent.sesmpd.phylonull, 
                                     summits.sesmpd.sourceSummits, summits.sesmpd.sourcePersis, 
                                     persistent.sesmpd.sourceSummits, persistent.sesmpd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Persistent, equal probability", 
                             "Summits, weighted for summits", "Summits, weighted for persistent", 
                              "Persistent, weighted for summits", "Persistent, weighted for persistent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,5,5,5))
dev.off()


pdf(file="output/9_PhyoDiversity/mntd_distribution.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, summits.sesmntd.sourceSUmmits, persistent.sesmntd.phylonull, 
                                     summits.sesmntd.sourcePersis, persistent.sesmntd.sourceSummits, persistent.sesmntd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Persistent, equal probability",
                             "Summits, weighted for summits", "Summits, weighted for persistent", 
                              "Persistent, weighted for summits", "Persistent, weighted for persistent"),
                    mainTitle = "Distribution of SES mntd", values=c(1,1,1,5,5,5,5))
dev.off()



###################################################################################### 
########################   Summarize Reduced Source pool effects ########################   
###################################################################################### 

pdf(file="output/9_PhyoDiversity/mpd.equalprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull), 
                    breaks=c("1", "2", "3" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1), 
                    colors=c("blue", "magenta3", "green4"))
dev.off()

pdf(file="output/9_PhyoDiversity/mntd.equalprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull), 
                    breaks=c("1", "2", "3" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1), 
                    colors=c("blue", "magenta3", "green4"))
dev.off()


###################################################################################### 
############################### Just of Interest ########################
###################################################################################### 


pdf(file="output/9_PhyoDiversity/mpd.weightedprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull, summits.sesmpd.sourceSummits, persistent.sesmpd.sourceSummits, persistent.sesmpd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5","6" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent","Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,5,5), 
                    colors=c("blue", "magenta3", "green4", "blue", "magenta3", "green4"))
dev.off()


pdf(file="output/9_PhyoDiversity/mntd.weightedprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull,summits.sesmntd.sourceSUmmits, persistent.sesmntd.sourceSummits, persistent.sesmntd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5","6" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent","Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,5,5), 
                    colors=c("blue", "magenta3", "green4", "blue", "magenta3", "green4"))
dev.off()




