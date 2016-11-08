###################################################################################### 
############################### Contemporary Pool: Summits ###########################
###################################################################################### 
alps.sites.df <- as.data.frame(alps.sites)
## Contemporary species pool = summits 
summits.sites <- as.data.frame(cbind("taxa" = names(alps.sites.df), as.data.frame(t(alps.sites.df))))
summits.sites <- filter(summits.sites, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)

## Reduce data, without traits
pezAlpes.summits <- comparative.comm(phy = alps.phy, comm = summits.sites[-c(1, 3:4),], env = alps.env) # Remove: Ecrins NP, Persistent, Under Ice
pezAlpes.summits # 215  taxa

pezAlpes.NOendmics.summits <- comparative.comm(phy = phy.NOendemics, comm = summits.sites[-c(1, 3:4),], env = alps.env.sprich) #traits = alps.traits
rownames(pezAlpes.NOendmics.summits$comm) #190

###################################################################################### 
############################### Contemporary Pool : Persistent ############################ 
###################################################################################### 

persistent.sites <- as.data.frame(cbind("taxa" = names(alps.sites.df), as.data.frame(t(alps.sites.df))))
persistent.sites <- filter(persistent.sites, Persistent > 0)
head(persistent.sites)
rownames(persistent.sites) <- persistent.sites$taxa
dim(persistent.sites)
persistent.sites <- t(persistent.sites[-1])
persistent.sites <- data.matrix(persistent.sites)
head(persistent.sites)

pezAlpes.persistent <- comparative.comm(phy = alps.phy, comm = persistent.sites[-c(1, 4),], env = alps.env) #Remove: Ecrins NP, Under Ice  
pezAlpes.persistent # 164 taxa

pezAlpes.persistent.NOendemics <- comparative.comm(phy = phy.NOendemics, comm = persistent.sites[-c(1, 4),], env = alps.env.sprich) #Remove: Ecrins NP, Under Ice  
pezAlpes.persistent.NOendemics # 144 taxa

