library(FD)
library(geiger)
library(picante)
library(spacodiR)
library(pez)
library(ape)
library(ggplot2)
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

# Deal with French properly
options(encoding="latin1")

source("R/TaxonomyHarmony.R")

################################################# 
################################################# 
############# See prepPipeline, dataWrangling for origin of the following datasets ###########
################################################# 
################################################# Phylogeny ################################################# 

alps.phy <- read.tree(file ="data/AnalysesDatasets/phy.Alpes.taxized.tre")
alps.phy #1154

################################################# Community ################################################# 

alps.sites <- read.csv(file="data/AnalysesDatasets/alps.sites.csv", row.names=1)
alps.sites <- data.matrix(alps.sites)
                          
################################################# Traits ################################################# 

alps.traits <- read.csv(file="data/AnalysesDatasets/alps.traits.csv", row.names=1)
alps.traits

alps.traits.dummy <- read.csv(file="data/AnalysesDatasets/alps.traits.dummy.csv", row.names=1)
head(alps.traits.dummy)

################################################# Metadata ################################################# 

alps.env <- read.csv(file="data/AnalysesDatasets/alps.env.csv", row.names=1)
alps.env

