#####################################################################################################################
############# Pre processesing of RAW datasets, and post-processing of prepPipelineSkyIsl output ####################
############# Both Static & Dynamic Null Models #####################################################################
############# Hannah E. Marx, 1 June 2016 ###########################################################################
#####################################################################################################################

################################################# Community ################################################# 

com.Alpes <- read.csv("data/OutputDatasets/comMerge.csv", as.is=T, row.names=2)
## NOTE: this is pre-removal of species not resolved to genus
head(com.Alpes)
tail(com.Alpes)
#com.Alpes <- data.matrix(com.Alpes)
dim(com.Alpes) #Persistent + underIce + NP + 19 summits + NP = 23, 1414 species

## Add a coulum for the summit counts
com.Alpes.tmp <- com.Alpes[order(rownames(com.Alpes)),]
com.Alpes.tmp <- as.data.frame(t(com.Alpes[-c(1)]), stringsAsFactors = F)
head(com.Alpes.tmp)
com.Alpes.tmp <- data.matrix(com.Alpes.tmp)
com.Alpes2 <- as.data.frame(rbind(com.Alpes.tmp, colSums(com.Alpes.tmp[5:nrow(com.Alpes.tmp),])))
rownames(com.Alpes2)[23] <- "Summits"
rownames(com.Alpes2)
rownames(com.Alpes2) <- sub("-", "_", rownames(com.Alpes2),)
rownames(com.Alpes2) <- sub(" ", "_", rownames(com.Alpes2),)

com.Alpes2 <- data.matrix(com.Alpes2)
head(com.Alpes2)
tail(com.Alpes2)


colSums(t(com.Alpes2) > 0) #summary of community matrix
dim(com.Alpes2)
CommunityMatrixFinal <- com.Alpes2
head(t(CommunityMatrixFinal))

#data.Alpes <- as.data.frame(cbind( "phyloTips" = colnames(CommunityMatrixFinal), t(CommunityMatrixFinal)))
#dim(data.Alpes) #1413

alps.sites <- com.Alpes2
head(alps.sites)
row.names(alps.sites)[3] <- "Ecrins NP"
 
## Combine certain area to one peak: 
alps.sites <- rbind(alps.sites, "Râteau" = colSums(alps.sites[c("brèche.du.râteau", "le.râteau", "pic.ouest..le.rateau."),]))
alps.sites <- rbind(alps.sites, "Lauvitel" = colSums(alps.sites[c("brèche.de.valsenestre", "signal.du.lauvitel", "pic.du.clapier.du.peyron"),]))
tail(alps.sites)

rownames(alps.sites)

alps.sites <- alps.sites[-c(6,7,9,15,16,20),]

rownames(alps.sites) <- c("Persistent", "Under Ice", "Ecrins NP", "Plat de la Selle", "Barre des Ecrins", "La Meije",
                          "Sirac", "Rouies", "Olan", "Pelvoux", "Occidentale Ailefroide", "Brevoort", "Burlan", 
                          "Rocher de la Selle", "Choisy", "Muraillette", "Summits", "Rateau", "Lauvitel")

#write.csv(alps.sites, file="data/OutputDatasets/alps.sites.v1.csv") 

## Edited in excel to sum summit counts that == 0 for Ecrins NP
which(alps.sites["Ecrins NP",] == 0)

# Saved as alps.sites_editEcrins.csv -> data/AnalysesDatasets/alps.sites.csv
# Edited "Persistent", "Under Ice", "Ecrins NP" to capitals


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

### Just checking coverage
summits$NOM %in% rownames(com.Alpes2)
rownames(com.Alpes2) %in% summits$NOM 

envAlps <- summits

## rename to match community matrix
#summits$NOM[4] <- "brèche.du.râteau"
#summits$NOM[6] <- "brèche.de.valsenestre"
#summits$NOM[8] <- "aiguille.du.plat.de.la.selle"
#summits$NOM[10] <- "rocher.de.la.selle"
#summits$NOM[12] <- "pic.du.clapier.du.peyron"
#summits$NOM[14] <- "tête.de.la.muraillette"
#summits$NOM[18] <- "pic.ouest..le.rateau."
#summits$NOM[20] <- "le.râteau"
#summits$NOM[13] # not in the community matrix
levels(envAlps$NOM)
 
levels(envAlps$NOM)[match("occidentale l'ailefroide",levels(envAlps$NOM))] <- "Occidentale Ailefroide"
levels(envAlps$NOM)[match("tour choisy",levels(envAlps$NOM))] <- "Choisy"
levels(envAlps$NOM)[match("l'olan",levels(envAlps$NOM))] <- "Olan"
levels(envAlps$NOM)[match("br\xe8che du r\xe2teau",levels(envAlps$NOM))] <- "Rateau"  ## "brèche.du.râteau" Combined into Rateau
levels(envAlps$NOM)[match("signal du lauvitel",levels(envAlps$NOM))] <- "Lauvitel"
levels(envAlps$NOM)[match("br\xe8che de valsenestre",levels(envAlps$NOM))] <- "Lauvitel" ##"brèche.de.valsenestre" ## Combined into lauvitel
levels(envAlps$NOM)[match("pointe brevoort",levels(envAlps$NOM))] <- "Brevoort"
levels(envAlps$NOM)[match("aiguille du plat de la selle",levels(envAlps$NOM))] <- "Plat de la Selle"
levels(envAlps$NOM)[match("les rouies",levels(envAlps$NOM))] <- "Rouies"
levels(envAlps$NOM)[match("rocher de la selle",levels(envAlps$NOM))] <- "Rocher de la Selle"
levels(envAlps$NOM)[match("pointes de burlan",levels(envAlps$NOM))] <- "Burlan"
levels(envAlps$NOM)[match("pic du clapier du peyron",levels(envAlps$NOM))] <- "Lauvitel" ##"Peyron" ## Combined into lauvitel
levels(envAlps$NOM)[match("aiguille dibona",levels(envAlps$NOM))] <- "Dibona"
levels(envAlps$NOM)[match("t\xeate de la muraillette",levels(envAlps$NOM))] <- "Muraillette"
levels(envAlps$NOM)[match("le sirac",levels(envAlps$NOM))] <- "Sirac"
levels(envAlps$NOM)[match("mont pelvoux",levels(envAlps$NOM))] <- "Pelvoux"
levels(envAlps$NOM)[match("la meije",levels(envAlps$NOM))] <- "La Meije"
levels(envAlps$NOM)[match("pic ouest (le rateau)",levels(envAlps$NOM))] <- "Rateau" ## "pic ouest (le rateau)" ## Combined into Rateau
levels(envAlps$NOM)[match("barre des \xe9crins",levels(envAlps$NOM))] <- "Barre des Ecrins"
levels(envAlps$NOM)[match("le r\xe2teau",levels(envAlps$NOM))] <- "Rateau"

## Summarize environ data for grouped summits
envAlps.tmp <- as.data.frame(envAlps %>% group_by(NOM) %>% summarise_each(funs(mean)))
head(envAlps.tmp)

envAlps.tmp <- envAlps.tmp[-c(2:5)]

rownames(envAlps.tmp) <- envAlps.tmp$NOM

### Add fill rows for analyses:
envAlps2 <- rbind("Persistent" = rep("NA", times = ncol(envAlps.tmp)), "Under Ice" = rep("NA", times = ncol(envAlps.tmp)), 
                  "Ecrins NP" = rep("NA", times = ncol(envAlps.tmp)), envAlps.tmp, "Summits" = rep("NA", times = ncol(envAlps.tmp)))

rownames(alps.sites) %in% rownames(envAlps2) #ALL TRUE

alps.env <- envAlps2

#### Read in elevation abouve LGM, and slope 
lgm.elevation <- read.csv(file="data/RAW/stat_elevation.csv")
colnames(lgm.elevation) <- c("obg", "Summit", "zone", "ct", "area", "min_elev", "max_elev", "range_elev",
                             "mean_elev", "std_elev","median_elev")

slope <- read.csv(file="data/RAW/stat_slope.csv")
colnames(slope) <- c("obg", "Summit", "zone", "ct", "min_slope", "max_slope", "range_slope",
                             "mean_slope", "std_slope")
env.new <- full_join(lgm.elevation[c(2,5:ncol(lgm.elevation))], slope[c(2,5:ncol(slope))])
head(env.new)
env.new$Summit

levels(env.new$Summit)[match("Valsenestre",levels(env.new$Summit))] <- "Lauvitel" ##"brèche.de.valsenestre" ## Combined into lauvitel
levels(env.new$Summit)[match("Peyron",levels(env.new$Summit))] <- "Lauvitel" ##"Peyron" ## Combined into lauvitel
str(env.new)

env.new.tmp <- as.data.frame(env.new %>% group_by(Summit) %>% summarise_each(funs(mean)))

area <- read.csv("data/RAW/sud_slope_sup40.csv")
head(area)
area.df <- as.data.frame(area %>% group_by(summit = sommet.C.50) %>% summarise(S.area = mean(SArea.N.19.11)))
env.new.tmp <- full_join(env.new.tmp, area.df, by= c("Summit" = "summit"))
env.new.tmp$Summit <- as.factor(env.new.tmp$Summit)

rownames(env.new.tmp) <- env.new.tmp[,1] 
head(env.new.tmp)

envAlps3 <- merge(env.new.tmp[,-1], envAlps2[,-1], by=0, all.x=T, all.y=T)
head(envAlps3)
rownames(envAlps3) <- envAlps3[,1]
alps.env <- envAlps3[-1]
names(alps.env)
head(alps.env)

lith.tmp <- read.csv("data/RAW/Geol_stats_edit.csv", row.names=1)
head(lith.tmp)
str(lith.tmp)
levels(lith.tmp$sommet.C.50)[match("Valsenestre",levels(lith.tmp$sommet.C.50))] <- "Lauvitel" ##"brèche.de.valsenestre" ## Combined into lauvitel
levels(lith.tmp$sommet.C.50)[match("Peyron",levels(lith.tmp$sommet.C.50))] <- "Lauvitel" ##"Peyron" ## Combined into lauvitel
lith.tmp$LITHO_SIMP.C.50 = factor(lith.tmp$LITHO_SIMP.C.50, levels=c(levels(lith.tmp$LITHO_SIMP.C.50), "other"))
lith.tmp[is.na(lith.tmp$LITHO_SIMP.C.50), "LITHO_SIMP.C.50"] <- "other"
lith.tmp[lith.tmp$LITHO_SIMP.C.50 == "", "LITHO_SIMP.C.50"] <- "other"

tdf <- as.data.frame(lith.tmp %>% group_by(sommet.C.50, LITHO_SIMP.C.50) %>% summarise(area.lith = sum(AREA.N.19.11))) 

head(tdf)
lith <- as.data.frame(tdf %>% group_by(sommet.C.50) %>% mutate(freq.lith = area.lith  / sum(area.lith )))
colnames(lith) <- c("summit", "lithology", "area.lith", "frequency")

lith <- merge(lith, alps.env.sprich.summits[2], by.x=1, by.y=0, all.x=T)

lith_plot <- ggplot(lith[order(lith$lithology), ], aes(x=reorder(factor(summit), as.numeric(as.character(ntax))), y= frequency, fill=lithology)) +
  geom_bar(stat="identity") +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = base_size *1.5, angle = -45, hjust=0, colour = "black"))
pdf(file="output/10_Environment/env_lith.pdf")
lith_plot
dev.off()

lith.sum <- lith %>% group_by(summit) %>% spread(lithology, frequency)
lith.sum[is.na(lith.sum)] <- 0
lith.sum.df <- as.data.frame(lith.sum %>% summarise_each(funs(max)))[,-c(2:3)]
rownames(lith.sum.df) <- lith.sum.df[,1]

lith.sum.df.com <- cbind(lith.sum.df[,-1], simpson.d.geol = vegan::diversity((lith.sum.df)[-1], index = "simpson"), sum.geol.type = rowSums(lith.sum.df[,-1] != 0))

alps.env.4 <- merge(alps.env, lith.sum.df.com, by=0, all.x=T)
rownames(alps.env.4) <- alps.env.4[,1]
alps.env.4 <- alps.env.4[,-1]
#write.csv(alps.env.4, file="data/AnalysesDatasets/alps.env.csv")

################################################# Endemics ################################################# 

endemics <- read.csv("data/RAW/AlpineEndemics_FloraAlpina_edit.csv", header =T, as.is=T)
head(endemics)

#endemics.taxize <- add.taxized.column(df = endemics, colnum = 1, spliton = "_", sepas = " ", source = "iPlant_TNRS")
head(endemics.taxize)
endemics.taxize.df <- endemics.taxize %>% group_by(taxized) %>% summarise_each(funs(first))
#write.csv(endemics.taxize.df, file="data/OutputDatasets/endemics.taxize.csv")


