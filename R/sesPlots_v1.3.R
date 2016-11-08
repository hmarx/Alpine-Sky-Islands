
source("R/sesPlotFunctions.R")

###################################################################################### 
############################### Pool : Ecrins NP ##################################### 
###################################################################################### 

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/EcrinsPool.phylonull.pdf") 
plotSESdispersion(mntd = ecrins.sesmntd.phylonull, mpd = ecrins.sesmpd.phylonull, 
                  mainTitle = "Dispersion Metrics Ecrins NP\nnull.model = equal probability draw from Ecrins Pool")
dev.off()

###################################################################################### 
############################### Reduced Pool : Summits ############################### 
###################################################################################### 

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/Summits.phylonull.pdf") 
plotSESdispersion(mpd = summits.sesmpd.phylonull, mntd = summits.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = equal probability draw from Summit Pool")
dev.off()
 
###################################################################################### 
############################### Reduced Pool : Persistent ############################ 
###################################################################################### 

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/Persistent.phylonull.pdf") 
plotSESdispersion(mpd = persistent.sesmpd.phylonull, mntd = persistent.sesmntd.phylonull, 
                  mainTitle = "Dispersion Metrics Persistent\nnull.model = equal probability draw from Persistent Pool")
dev.off()








###################################################################################### 
############################### Contemporary Pool: Sumits  ###########################
###################################################################################### 


## Contemporary species pool = summits 
## Hisoric source pool = summits (persistent + under ice?)
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/Summits.sourceSummits.pdf") 
plotSESdispersion(mpd = summits.sesmpd.sourceSummits, mntd = summits.sesmntd.sourceSummits, 
                  mainTitle = "Dispersion Metrics Summits\nnull.model = weighted draw from summits source pool")
dev.off()


###################################################################################### 
############################### Contemporary Pool: Persistent ########################
###################################################################################### 

## Contemporary species pool = persistent 
## Hisoric source pool = summits (persistent + under ice?)
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/Persistent.sourceSummits.pdf") 
plotSESdispersion(mpd = persistent.sesmpd.sourcePersis, mntd = persistent.sesmntd.sourcePersis, 
                  mainTitle = "Dispersion Metrics Persistent\nnull.model = weighted draw from summits source pool")
dev.off()

## Contemporary species pool = persistent 
## Hisoric source pool = summits (persistent + under ice?)
pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/Persistent.sourcePersistent.pdf") 
plotSESdispersion(mpd = persistent.sesmpd.sourceSummits, mntd = persistent.sesmntd.sourceSummits, 
                  mainTitle = "Dispersion Metrics Persistent\nnull.model = weighted draw from persistent source pool")
dev.off()


###################################################################################### 
############################### ALL ########################
###################################################################################### 

pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/mpd_distribution.pdf") 
plotDistributionSES(outputSES = list(ecrins.sesmpd.phylonull, summits.sesmpd.phylonull,  persistent.sesmpd.phylonull, 
                                     summits.sesmpd.phylonull.abun, summits.sesmpd.sourceSummits, 
                                     persistent.sesmpd.phylonull.abun, persistent.sesmpd.sourcePersis, persistent.sesmpd.sourceSummits), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7", "8"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Persistent, equal probability", 
                             "Summits, abundance weighted", "Summits, weighted for persistent", 
                             "Persistent, abundance weighted", "Persistent, weighted for persistent", "Persistent, weighted for summits"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,6,5,6,4), 
                    colors=c("blue","magenta3","green4","magenta3","magenta3","green4", "green4","green4"),
                    max.y=1)
dev.off()


pdf(file="output/9_PhyoDiversity/Spermatophyta/static/MaxLikelihood/figs/mntd_distribution.pdf") 
plotDistributionSES(outputSES = list(ecrins.sesmntd.phylonull, summits.sesmntd.phylonull, persistent.sesmntd.phylonull, 
                                     summits.sesmntd.phylonull.abun, summits.sesmntd.sourceSummits,
                                     persistent.sesmntd.phylonull.abun, persistent.sesmntd.sourcePersis, persistent.sesmntd.sourceSummits), 
                    breaks=c("1", "2", "3", "4", "5", "6", "7","8"),
                    labels=c("Ecrins, equal probability", "Summits, equal probability", "Persistent, equal probability", 
                             "Summits, abundance weighted", "Summits, weighted for persistent", 
                             "Persistent, abundance weighted", "Persistent, weighted for persistent", "Persistent, weighted for summits"),
                    mainTitle = "Distribution of SES mntd", values=c(1,1,1,5,6,5,6,4), 
                    colors=c("blue","magenta3","green4","magenta3","magenta3","green4", "green4","green4"),
                    max.y=1)
dev.off()









###################################################################################### 
########################   Summarize Reduced Source pool effects ########################   
###################################################################################### 

#pdf(file="output/9_PhyoDiversity/mpd.equalprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull), 
                    breaks=c("1", "2", "3" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1), 
                    colors=c("blue", "magenta3", "green4"))
#dev.off()

#pdf(file="output/9_PhyoDiversity/mntd.equalprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull), 
                    breaks=c("1", "2", "3" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1), 
                    colors=c("blue", "magenta3", "green4"))
#dev.off()


###################################################################################### 
############################### Just of Interest ########################
###################################################################################### 


#pdf(file="output/9_PhyoDiversity/mpd.weightedprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmpd.phylonull, summits.sesmpd.phylonull, 
                                     persistent.sesmpd.phylonull, summits.sesmpd.sourceSummits, persistent.sesmpd.sourceSummits, persistent.sesmpd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5","6" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent","Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,5,5), 
                    colors=c("blue", "magenta3", "green4", "blue", "magenta3", "green4"))
#dev.off()


#pdf(file="output/9_PhyoDiversity/mntd.weightedprob.pdf") 
plotDistributionSES(outputSES = list(comm.sesmntd, summits.sesmntd.phylonull, 
                                     persistent.sesmntd.phylonull,summits.sesmntd.sourceSUmmits, persistent.sesmntd.sourceSummits, persistent.sesmntd.sourcePersis), 
                    breaks=c("1", "2", "3", "4", "5","6" ),
                    labels=c( "Ecrins NP", "Summits", "Peristent","Ecrins NP", "Summits", "Peristent"),
                    mainTitle = "Distribution of SES mpd", values=c(1,1,1,5,5,5), 
                    colors=c("blue", "magenta3", "green4", "blue", "magenta3", "green4"))
#dev.off()




