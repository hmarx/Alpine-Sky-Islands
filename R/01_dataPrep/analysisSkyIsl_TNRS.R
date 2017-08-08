#####################################################################################################################
############# Core to load datasets with TNRS Taxonomy ##############################################################
############# Taxonomy changed to TAXREF v.9 for final version ######################################################
############# Some of the data prep. depends on these version, so this is kept for reference ########################
############# Hannah E. Marx, 25 April 2017 #########################################################################
#####################################################################################################################

source("R/functions/TaxonomyHarmony.R")

################################################# 
################################################# 
############# See prepPipelineSkyIsl, dataWrangling for origin of the following datasets ###########
################################################# 
################################################# Phylogeny ################################################# 

## Spermatophyta
alps.TNRS.phy <- read.tree(file ="data/TNRS/phy.Alpes.taxized.tre")
alps.TNRS.phy #1084

################################################# Community ################################################# 

alps.TNRS.sites.df <- read.csv(file="data/TNRS/alps.sites.csv", row.names=1, header=T)
#colnames(alps.TNRS.sites.df)
alps.TNRS.sites <- data.matrix(alps.TNRS.sites.df)
dim(alps.TNRS.sites) 
head(alps.TNRS.sites)
rownames(alps.TNRS.sites) #19 (15 summits + 4 pools)

tmp <- alps.TNRS.sites
tmp[which(tmp != 0)] <- 1
alps.TNRS.sites.pa.all <- t(tmp)
dim(alps.TNRS.sites.pa.all) #1414

################################################# Metadata ################################################# 

alps.TNRS.env <- read.csv(file="data/AnalysesDatasets/alps.env.csv", row.names=1, header=T)
alps.TNRS.env
names(alps.TNRS.env)

alps.TNRS.env.sprich <- merge(cbind(ntax = (rowSums(t(alps.TNRS.sites.pa.all)))), alps.TNRS.env,  by=0)

alps.TNRS.env.sprich <- alps.TNRS.env.sprich[order(as.numeric(alps.TNRS.env.sprich$ntax)),]
rownames(alps.TNRS.env.sprich) <- alps.TNRS.env.sprich$Row.names
head(alps.TNRS.env.sprich)

alps.TNRS.env.sprich.summits <- cbind(na.omit(alps.TNRS.env.sprich), elevation = c("Brevoort" = 3765,
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

#df.t <- matrix(nrow=4,ncol=ncol(alps.TNRS.env.sprich))

#rownames(df.t) <- c("persistent.TNRS", "Summits", "Ecrins NP", "Under Ice")
#colnames(df.t) <- colnames(alps.TNRS.env.sprich)
  
#alps.TNRS.env.sprich <- rbind(alps.TNRS.env.sprich, df.t)

#write.csv(alps.TNRS.env.sprich.summits, "data/OutputDatasets/AppendixS2_alps.TNRS.env.csv")


######################################################################################################################## 
################################################# Combined Data Object #################################################
####################################################################################################################### 

###################################################################################### 
############################### Contemporary Pool: Regional Ecrins NP ################
###################################################################################### 
### pez will trim data/ phylogeny; no need to use drop tip objects
pezalps.TNRS <- comparative.comm(phy = alps.TNRS.phy, comm = alps.TNRS.sites, env = alps.TNRS.env.sprich) #traits = alps.TNRS.traits
pezalps.TNRS$dropped
pezalps.TNRS #1081
sites(pezalps.TNRS)
pezalps.TNRS$phy

###################################################################################### 
############################### Contemporary Pool: All Summits #######################
###################################################################################### 
alps.TNRS.sites.df <- as.data.frame(alps.TNRS.sites)
## Contemporary species pool = summits 
summits.sites <- as.data.frame(cbind("taxa" = names(alps.TNRS.sites.df), as.data.frame(t(alps.TNRS.sites.df))))
summits.sites <- filter(summits.sites, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

## Reduce data, without traits
pezalps.TNRS.summits <- comparative.comm(phy = alps.TNRS.phy, comm = summits.sites[-c(1, 3:4),], env = alps.TNRS.env) # Remove: Ecrins NP, persistent.TNRS, Under Ice
pezalps.TNRS.summits # 215  taxa

tmp2 <- pezalps.TNRS.summits$comm
tmp2[which(tmp2 != 0)] <- 1
alps.TNRS.summits.pa.all <- t(tmp2)

###################################################################################### 
############################### Contemporary Pool: LGM persistent.TNRS ####################
###################################################################################### 

persistent.TNRS.sites <- as.data.frame(cbind("taxa" = names(alps.TNRS.sites.df), as.data.frame(t(alps.TNRS.sites.df))))
persistent.TNRS.sites <- filter(persistent.TNRS.sites, Persistent > 0)
head(persistent.TNRS.sites)
rownames(persistent.TNRS.sites) <- persistent.TNRS.sites$taxa
dim(persistent.TNRS.sites)
persistent.TNRS.sites <- t(persistent.TNRS.sites[-1])
persistent.TNRS.sites <- data.matrix(persistent.TNRS.sites)
head(persistent.TNRS.sites)

pezalps.TNRS.persistent.TNRS <- comparative.comm(phy = alps.TNRS.phy, comm = persistent.TNRS.sites[-c(1, 4),], env = alps.TNRS.env) #Remove: Ecrins NP, Under Ice  
pezalps.TNRS.persistent.TNRS # 164 taxa

# Pres/abs for persistent.TNRS pool
tmp3 <- pezalps.TNRS.persistent.TNRS$comm
tmp3[which(tmp3 != 0)] <- 1
alps.TNRS.persistent.TNRS.pa.all <- t(tmp3)


################################################# endemic.TNRSs.TNRS ################################################# 

endemic.TNRSs.TNRS <- read.csv("data/OutputDatasets/endemics.taxize.edit.csv", header =T, as.is=T, row.names=1)
head(endemic.TNRSs.TNRS)
endemic.TNRSs.TNRS.df <- cbind(endemic.TNRSs.TNRS[1],"endemic.TNRS"= rep("1", nrow(endemic.TNRSs.TNRS)))
com.endemic.TNRS.tmp <- merge(endemic.TNRSs.TNRS.df[2], alps.TNRS.sites.pa.all, by=0, all.x=F, all.y=T)
com.endemic.TNRS.tmp$endemic.TNRS <- as.character(as.numeric(com.endemic.TNRS.tmp$endemic.TNRS))
com.endemic.TNRS.tmp$endemic.TNRS[!is.na(com.endemic.TNRS.tmp$endemic.TNRS)] <- 1
com.endemic.TNRS.tmp$endemic.TNRS[is.na(com.endemic.TNRS.tmp$endemic.TNRS)] <- 0
head(com.endemic.TNRS.tmp)
dim(com.endemic.TNRS.tmp)
#write.csv(com.endemic.TNRS.tmp, "data/OutputDatasets/ecrins.releve.tnrs.community.csv")

dim(endemic.TNRSs.TNRS[which(rownames(endemic.TNRSs.TNRS) %in% pezalps.TNRS$phy$tip.label),]) #59

#dim(endemic.TNRSs.TNRS[which(rownames(endemic.TNRSs.TNRS) %in% ecrins.include.sperma),]) #98

treedata.endemic.TNRSs.TNRS <- treedata(phy  = alps.TNRS.phy, data = as.matrix(endemic.TNRSs.TNRS)) 

length(which(treedata.endemic.TNRSs.TNRS$phy$tip.label %in% pezalps.TNRS.summits$phy$tip.label)) #25

phy.NOendemic.TNRSs.TNRS <- drop.tip(phy = alps.TNRS.phy, tip = treedata.endemic.TNRSs.TNRS$phy$tip.label) 

pezalps.TNRS.NOendmics <- comparative.comm(phy = phy.NOendemic.TNRSs.TNRS, comm = alps.TNRS.sites, env = alps.TNRS.env.sprich) #traits = alps.TNRS.traits
rownames(pezalps.TNRS.NOendmics$comm) #1022

pezalps.TNRS.NOendmics.summits <- comparative.comm(phy = phy.NOendemic.TNRSs.TNRS, comm = summits.sites[-c(1, 3:4),], env = alps.TNRS.env.sprich) #traits = alps.TNRS.traits
rownames(pezalps.TNRS.NOendmics.summits$comm) #190

pezalps.TNRS.persistent.TNRS.NOendemic.TNRSs.TNRS <- comparative.comm(phy = phy.NOendemic.TNRSs.TNRS, comm = persistent.TNRS.sites[-c(1, 4),], env = alps.TNRS.env.sprich) #Remove: Ecrins NP, Under Ice  
pezalps.TNRS.persistent.TNRS.NOendemic.TNRSs.TNRS # 144 taxa

