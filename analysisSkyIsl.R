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
#devtools::install_github("taxize", "ropensci")
library(plyr)
library(dplyr)
library(Biostrings)
#devtools::install_github("oschwery/MonoPhy")
library(MonoPhy)
#devtools::install_github("richfitz/storr")
#devtools::install_github("wcornwell/TaxonLookup")
library(TaxonLookup)
#install.packages("stringi")
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
#devtools::install_github("GuangchuangYu/ggtree")
library(ggtree)
#rm(list=ls())

# Deal with French properly
options(encoding="latin1")

source("R/TaxonomyHarmony.R")

################################################# 
################################################# 
############# See prepPipeline, dataWrangling for origin of the following datasets ###########
################################################# 
################################################# Phylogeny ################################################# 

## Spermatophyta
alps.phy <- read.tree(file ="data/AnalysesDatasets/phy.Alpes.taxized.tre")
alps.phy #1084

################################################# Community ################################################# 

alps.sites.df <- read.csv(file="data/AnalysesDatasets/alps.sites.csv", row.names=1, header=T)
#colnames(alps.sites.df)
alps.sites <- data.matrix(alps.sites.df)
dim(alps.sites) 
head(alps.sites) 
rownames(alps.sites) #19 (15 summits + 4 pools)

tmp <- alps.sites
tmp[which(tmp != 0)] <- 1
alps.sites.pa.all <- t(tmp)

################################################# Metadata ################################################# 

alps.env <- read.csv(file="data/AnalysesDatasets/alps.env.csv", row.names=1, header=T)
alps.env
names(alps.env)

alps.env.sprich <- merge(cbind(ntax = (rowSums(t(alps.sites.pa.all)))), alps.env,  by=0)

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

#df.t <- matrix(nrow=4,ncol=ncol(alps.env.sprich))

#rownames(df.t) <- c("Persistent", "Summits", "Ecrins NP", "Under Ice")
#colnames(df.t) <- colnames(alps.env.sprich)
  
#alps.env.sprich <- rbind(alps.env.sprich, df.t)

#write.csv(alps.env.sprich.summits, "data/OutputDatasets/AppendixS2_alps.env.csv")


################################################# Endemics ################################################# 

endemics <- read.csv("data/OutputDatasets/endemics.taxize.edit.csv", header =T, as.is=T, row.names=1)
head(endemics)
endemics.df <- cbind(endemics[1],"Endemic"= rep("1", nrow(endemics)))
com.endemic.tmp <- merge(endemics.df[2], alps.sites.pa.all, by=0, all.x=F, all.y=T)
com.endemic.tmp$Endemic <- as.character(as.numeric(com.endemic.tmp$Endemic))
com.endemic.tmp$Endemic[!is.na(com.endemic.tmp$Endemic)] <- 1
com.endemic.tmp$Endemic[is.na(com.endemic.tmp$Endemic)] <- 0
head(com.endemic.tmp)
#write.csv(com.endemic.tmp, "data/OutputDatasets/ecrins.releve.tnrs.community.csv")


dim(endemics[which(rownames(endemics) %in% pezAlpes$phy$tip.label),]) #59

dim(endemics[which(rownames(endemics) %in% ecrins.include.sperma),]) #98

treedata.endemics <- treedata(phy  = alps.phy, data = as.matrix(endemics)) 

length(which(treedata.endemics$phy$tip.label %in% pezAlpes.summits$phy$tip.label)) #25

phy.NOendemics <- drop.tip(phy = alps.phy, tip = treedata.endemics$phy$tip.label) 

pezAlpes.NOendmics <- comparative.comm(phy = phy.NOendemics, comm = alps.sites, env = alps.env.sprich) #traits = alps.traits
rownames(pezAlpes.NOendmics$comm) #1022

################################################# Combined Data Object ################################################# 
### pez will trim data/ phylogeny; no need to use drop tip objects
pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env.sprich) #traits = alps.traits

pezAlpes #1081
sites(pezAlpes)

