
##################################################################################################  
##### Pre processesing of RAW datasets, and post-processing of prepPipeline output
################################################################################################## 

################################################# Community ################################################# 

com.Alpes <- read.csv("data/OutputDatasets/comMerge.csv", as.is=T, row.names=2)
## NOTE: this is pre-removal of species not resolved to genus
head(com.Alpes)
tail(com.Alpes)
#com.Alpes <- data.matrix(com.Alpes)
dim(com.Alpes) #persistent + underIce + NP + 19 summits + NP = 23, 1413 species

## Add a coulum for the summit counts
com.Alpes.tmp <- com.Alpes[order(rownames(com.Alpes)),]
com.Alpes.tmp <- as.data.frame(t(com.Alpes[-c(1)]), stringsAsFactors = F)
head(com.Alpes.tmp)
com.Alpes.tmp <- data.matrix(com.Alpes.tmp)
com.Alpes2 <- as.data.frame(rbind(com.Alpes.tmp, colSums(com.Alpes.tmp[5:nrow(com.Alpes.tmp),])))
rownames(com.Alpes2)[23] <- "Summits"
com.Alpes2 <- data.matrix(com.Alpes2)
head(com.Alpes2)
tail(com.Alpes2)

colnames(com.Alpes2)
colnames(com.Alpes2) <- sub("-", "_", colnames(com.Alpes2),)
colnames(com.Alpes2) <- sub(" ", "_", colnames(com.Alpes2),)
colSums(t(com.Alpes2) > 0) #summary of community matrix
dim(com.Alpes2)
CommunityMatrixFinal <- com.Alpes2
head(t(CommunityMatrixFinal))

#data.Alpes <- as.data.frame(cbind( "phyloTips" = colnames(CommunityMatrixFinal), t(CommunityMatrixFinal)))
#dim(data.Alpes) #1413

alps.sites <- com.Alpes2
head(alps.sites)
row.names(alps.sites)[3] <- "EcrinsNP"
#write.csv(alps.sites, file="data/AnalysesDatasets/alps.sites.csv")


##### Don't actually need to drop b/c pez does...moved to /Extra
## Need to drop the taxa in tree that are not in the dataset (STILL DON'T KNOW WHY THIS IS>>>>ASK CODY)
drops <- phy.Alpes$tip.label[!phy.Alpes$tip.label %in% rownames(data.Alpes)]
phy.Alpes.drop <- drop.tip(phy = phy.Alpes, tip = drops)
phy.Alpes.drop

tmp <- treedata(phy=phy.Alpes.drop, data=data.Alpes)
tmp$phy
data.Alpes.drop <- as.data.frame(tmp$data)
head(data.Alpes.drop)
#write.csv(data.Alpes.drop, file="data/Extra/data.Alpes.drop.csv")

phy.Alpes.drop$tip.label[!phy.Alpes.drop$tip.label %in% rownames(data.Alpes.drop)]
rownames(data.Alpes.drop)[!rownames(data.Alpes.drop) %in% rownames(data.Alpes.drop)]

## Convert community matrix to pres/absence 
data.Alpes.drop.pa <- data.Alpes.drop[2:ncol(data.Alpes.drop)]
data.Alpes.drop.pa[data.Alpes.drop.pa != 0 ] <- 1
data.Alpes.drop.pa <- cbind("phyloTips" = rownames(data.Alpes.drop.pa), data.Alpes.drop.pa)
head(data.Alpes.drop.pa)
#write.csv(data.Alpes.drop.pa, file="data/Extra/data.Alpes.drop.preAbs.csv")


################################################# Traits ################################################# 
traitnames <- read.csv(file="data/RAW/Traits/sp_list.csv", sep=";")
head(traitnames)
str(traitnames)

traits <- read.csv(file="data/RAW/Traits/traits.csv")
head(traits)
colnames(traits)
unique(traits$code)
str(traits)

traits.data <- as.data.frame(traits[1:3], strings.as.factors=F) #get just the species code, trait name, and value
head(traits.data) 
traits.data[traits.data == "NULL"] <- "NA"
dim(traits.data) #168,942      
traits.data$valeur <- as.numeric(as.character(traits.data$valeur))
str(traits.data)

# Take mean of traits with multiple values for a species:
traits.data.avg <- traits.data %>% group_by(code_cbna, code) %>% select(valeur) %>% 
  summarise(avgVal = mean(valeur, na.rm=T))
traits.data.avg

# Reshape dataset to have species in rows, traits$code in columns 
traits.data.melt <- dcast(traits.data.avg, code_cbna ~ code)
head(traits.data.melt)
dim(traits.data.melt) #1421 species, 130 traits

# Match cbna code to species name
all.traits.alps <- right_join(traitnames, traits.data.melt, by="code_cbna")
head(all.traits.alps)

# taxize to match species occurence data, phylogeny
all.traits.alps.taxize <- add.taxized.column(df = as.data.frame(all.traits.alps), colnum = 2, spliton = " ", sepas = " ", source="iPlant_TNRS")
head(all.traits.alps.taxize) 
head(all.traits.alps.taxize[duplicated(all.traits.alps.taxize$taxized),])

# reduce down to traits of interest
colnames(all.traits.alps.taxize)

cont.traits <- as.data.frame(all.traits.alps.taxize[c("taxized", "L_AREA", "L_DRYM", "PL_H_flora_MAX", "PLOIDY (VALUE)", "SEEDM", "SLA")], strings.as.factors=F)
head(cont.traits)
cols.num <- c(colnames(cont.traits))
cont.traits[2:7]  <- sapply(cont.traits[2:7], as.numeric)
sapply(cont.traits, class)

# Take mean of traits with multiple values for a species again:
all.traits.alps.taxize.avg <- cont.traits %>% group_by(taxized) %>% summarise_each(funs(mean))
head(all.traits.alps.taxize.avg)

cont.traits <- as.data.frame(all.traits.alps.taxize.avg[2:length(all.traits.alps.taxize.avg)])
rownames(cont.traits) <- all.traits.alps.taxize.avg$taxized
head(cont.traits)
dim(cont.traits) # 1400 species 

######################### get discrete traits 
disc.traits.tmp <- traits[c(1:2,4)]

# take first of duplicated values
traits.data.first <- disc.traits.tmp %>% group_by(code_cbna, code) %>% 
  summarise(first(nom))
traits.data.first

# reshape the data to have single species in row
traits.data.melt.disc <- dcast(traits.data.first, code_cbna ~ code)
head(traits.data.melt.disc)
dim(traits.data.melt.disc) #1421 species, 130 traits

# Match cbna code to species name
traits.data.melt.disc <- right_join(traitnames, traits.data.melt.disc, by="code_cbna")
head(traits.data.melt.disc)

# taxize to match species occurence data, phylogeny
disc.traits.alps.taxize <- add.taxized.column(df = as.data.frame(traits.data.melt.disc), colnum = 2, spliton = " ", sepas = " ", source="iPlant_TNRS")
head(disc.traits.alps.taxize) 

# reduce down to traits of interest
disc.traits <- disc.traits.alps.taxize[c("taxized", "PLOIDY", "TAPROOT", "VEG_DISP", "WOODY", "LIFESPAN", "GROWTHF")]
dim(disc.traits)
duplicated(disc.traits$taxized)

# Take first of traits with multiple values for a species again:
disc.traits.alps.taxize.avg <- disc.traits %>% group_by(taxized) %>% summarise_each(funs(first))
head(disc.traits.alps.taxize.avg)

disc.traits <- as.data.frame(disc.traits.alps.taxize.avg[2:ncol(disc.traits.alps.taxize.avg)])
rownames(disc.traits) <- disc.traits.alps.taxize.avg$taxized

## Merge the two 
trait.Alps.tmp <- merge(cont.traits, disc.traits, by=0)
trait.Alps<- trait.Alps.tmp[2:ncol(trait.Alps.tmp)]
rownames(trait.Alps) <- sub(trait.Alps.tmp$Row.names, pattern =" ", replacement="_")
head(trait.Alps)

alps.traits <- trait.Alps
write.csv(alps.traits, file="data/AnalysesDatasets/alps.traits.csv")


################################################# Metadata ################################################# 
#### Metadata: characteristics of summits
########## topographic indicies
########## geological composition (geo-alp: http://www.geol-alp.com)
########## mean climate gradients 
climateVariables <- load("data/RAW/climatic_variables")
list(climateVariables)
head(summits)
dim(summits) #20 summits 85 variables
names(summits)

## rename to match community matrix
summits$NOM[4] <- "brèche.du.râteau"
summits$NOM[6] <- "brèche.de.valsenestre"
summits$NOM[8] <- "aiguille.du.plat.de.la.selle"
summits$NOM[10] <- "rocher.de.la.selle"
summits$NOM[12] <- "pic.du.clapier.du.peyron"
summits$NOM[14] <- "tête.de.la.muraillette"
summits$NOM[18] <- "pic.ouest..le.rateau."
summits$NOM[20] <- "le.râteau"

summits$NOM[13] # not in the community matrix


### Just checking coverage
summits$NOM %in% rownames(com.Alpes2)
rownames(com.Alpes2) %in% summits$NOM 

envAlps <- summits[-3]
rownames(envAlps) <- summits$NOM
envAlps

### Add fill rows for analyses:
envAlps2 <- rbind("persistant" = rep("NA", times = ncol(envAlps)), "underIce" = rep("NA", times = ncol(envAlps)), 
                  "EcrinsNP" = rep("NA", times = ncol(envAlps)), envAlps, "Summits" = rep("NA", times = ncol(envAlps)))

rownames(alps.sites) %in% rownames(envAlps2)

alps.env <- envAlps2
write.csv(alps.env, file="data/AnalysesDatasets/alps.env.csv")

