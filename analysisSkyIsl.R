#####################################################################################################################
############# Code for: Riders in the sky (islands): using a mega-phylogenetic approach to understand plant species #
###################### distribution and coexistence at the altitudinal limits of angiosperm plant life ##############
############# Core datasets and required R packages #################################################################
############# Hannah E. Marx, 25 April 2017 ###########################################################################
#####################################################################################################################

## Install and load the following required packages
#devtools::install_github("richfitz/storr")
#devtools::install_github("wcornwell/TaxonLookup", force = TRUE)
#devtools::install_github("taxize", "ropensci")
#devtools::install_github("oschwery/MonoPhy")
#devtools::install_github("GuangchuangYu/ggtree")

library(FD)
library(geiger)
library(picante)
library(spacodiR)
library(pez)
library(ape)
library(ggplot2)
library(gplots)
library(pheatmap)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(taxize)
library(plyr)
library(dplyr)
library(Biostrings)
library(MonoPhy)
library(taxonlookup)
library(cluster)
library(diversitree)
library(grid)
require(DAMOCLES)
library(ade4)
library(sp)
library(rgeos)
library(ecodist)
library(broom)
library(tidyr)
library(stargazer)
library(maps)
library(mapdata)
library(sp)
library(rgdal)
library(rgeos)
library(ggmap)
library(mapplots)
library(maptools)
library(raster)
library(paleotree)
library(vegan)
library(xtable)
library(BSDA)
library(ggtree)
library(data.table)
library(stringr)
require(phytools)

#rm(list=ls())

# Deal with French properly
options(encoding="latin1")

source("R/functions/TaxonomyHarmony.R")

############################################################################################################# 
############# See prepPipelineSkyIsl.R, dataWrangling.R for origin of the following datasets ################
############################################################################################################# 

################################################# Phylogeny ################################################# 

alps.phy <- read.tree(file ="data/AnalysesDatasets/phy.Alpes.taxized.taxref.tre")
alps.phy #1070 taxa

################################################# Community ################################################# 

alps.sites.df <- read.csv(file="data/AnalysesDatasets/v5tax_edit.csv", row.names=1, header=T)
colnames(alps.sites.df) <- (gsub(colnames(alps.sites.df), pattern = "[.]", replacement = " ")) # replace periods with spaces in summit names
head(alps.sites.df)
alps.sites <- data.matrix(alps.sites.df[-c(1,2)]) # convert to data matrix and remove unnecessary columns
rownames(alps.sites) #19 (4 species pools + 15 summits)
alps.sites <- t(alps.sites) # transform data frame (rows = summits; columns = species)
dim(alps.sites) # 1345 species

head(t(alps.sites))
## The number of species collected in the All Summits pool
length(which((t(alps.sites)[,3]) == 1)) #306
## The number of species collected in the LGM pool
length(which((t(alps.sites)[,4]) == 1)) #235

################################################# Metadata ################################################# 

alps.env <- read.csv(file="data/AnalysesDatasets/alps.env.csv", row.names=1, header=T)
alps.env
names(alps.env) #climate and geological variables

alps.env.sprich <- merge(cbind(ntax = (rowSums((alps.sites)))), alps.env,  by=0)

alps.env.sprich <- alps.env.sprich[order(as.numeric(alps.env.sprich$ntax)),]
rownames(alps.env.sprich) <- alps.env.sprich$Row.names
head(alps.env.sprich)

alps.env.sprich.summits <- cbind(na.omit(alps.env.sprich), elevation = c("Brevoort" = 3765,
                                                                "Pelvoux" = 3946,
                                                                "Choisy" = 3671,
                                                                "Burlan" = 3207,
                                                                "Plat de la Selle" = 3597,
                                                                "Occidentale Ailefroide" = 3954,
                                                                "La Meije" = 3983,
                                                                "Barre des Ecrins" = 4102,
                                                                "Rouies" =3589,
                                                                "Muraillette" = 3019,
                                                                "Olan" = 3564,
                                                                "Sirac" = 3440,
                                                                "Lauvitel" = 2904,
                                                                "Rateau" = 3809,
                                                                "Rocher de la Selle" = 2791))

#write.csv(alps.env.sprich.summits, "data/OutputDatasets/AppendixS2_alps.env.csv")


######################################################################################################################## 
############################ Combined Data Objects for each species pool ###############################################
######################################################################################################################## 

###################################################################################### 
############################### Contemporary Pool: Regional Ecrins NP ################
###################################################################################### 

### pez will trim data/phylogeny so that all datasets are the same size 
pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env.sprich) 
pezAlpes$dropped
pezAlpes #1065 taxa
1065/1345 #0.7918216 taxa in GenBank
sites(pezAlpes)

###################################################################################### 
############################### Contemporary Pool: All Summits #######################
###################################################################################### 

alps.sites.df <- as.data.frame(alps.sites)
summits.sites <- as.data.frame(cbind("taxa" = names(alps.sites.df), as.data.frame(t(alps.sites.df))))
summits.sites <- filter(summits.sites, Summits > 0) # filter for only speices occuring on alpine summits
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

# Remove: Ecrins NP, Persistent
rownames(summits.sites)
pezAlpes.summits <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(2),], env = alps.env) 
pezAlpes.summits # 215 taxa
215/306 #0.7026144% summit taxa in GenBank

sites(pezAlpes.summits)
# Pres/abs for all summits species pool
alps.summits.pa.all <- t( pezAlpes.summits$comm)

###################################################################################### 
############################### Contemporary Pool: LGM Persistent ####################
###################################################################################### 

persistent.sites <- as.data.frame(cbind("taxa" = names(alps.sites.df), as.data.frame(t(alps.sites.df))))
persistent.sites <- filter(persistent.sites, Persistent > 0) # filter for only speices occuring in areas that persisted through the LGM
head(persistent.sites)
rownames(persistent.sites) <- persistent.sites$taxa
dim(persistent.sites)
persistent.sites <- t(persistent.sites[-1])
persistent.sites <- data.matrix(persistent.sites)
head(persistent.sites)

#Remove: Ecrins NP
rownames(persistent.sites)
pezAlpes.persistent <- comparative.comm(phy = alps.phy, comm = persistent.sites[-c(2,3),], env = alps.env) 
pezAlpes.persistent # 167 taxa
167/235 #0.7106383% LGM species in GenBank

# Pres/abs for persistent pool
alps.persistent.pa.all <- t(pezAlpes.persistent$comm)


###################################################################################### 
############################### Contemporary Pool: Regional Ecrins NP ################
############################### Without Endemic Species               ################
###################################################################################### 

alps.sites.t <- as.data.frame(t(alps.sites))

endemics <- alps.sites.t[alps.sites.t$Endemic == 1, ]
head(endemics)
dim(endemics) # 101 endemic species 

dim(endemics[which(rownames(endemics) %in% pezAlpes$phy$tip.label),]) #58 endemics in phylogeny 

treedata.endemics <- treedata(phy  = alps.phy, data = as.matrix(endemics)) 

length(which(treedata.endemics$phy$tip.label %in% pezAlpes.summits$phy$tip.label)) #24

phy.NOendemics <- drop.tip(phy = alps.phy, tip = treedata.endemics$phy$tip.label) 

pezAlpes.NOendmics <- comparative.comm(phy = phy.NOendemics, comm = alps.sites, env = alps.env.sprich) 
rownames(pezAlpes.NOendmics$comm) #1005

pezAlpes.NOendmics.summits <- comparative.comm(phy = phy.NOendemics, comm = summits.sites, env = alps.env.sprich)
rownames(pezAlpes.NOendmics.summits$comm) #190

pezAlpes.persistent.NOendemics <- comparative.comm(phy = phy.NOendemics, comm = persistent.sites, env = alps.env.sprich) #Remove: Ecrins NP, Under Ice  
pezAlpes.persistent.NOendemics # 145 taxa

