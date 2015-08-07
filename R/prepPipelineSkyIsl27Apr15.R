## FRENCH ALPES ## 

source("R/treeFunctions.R")
source("R/TaxonomyHarmony.R")
source("R/dataWrangling.R")

####################################### Read in Raw datasets, homogonize taxonomy ####################################### 
############# ALL RELEVES : from Ecrin National Park, 2009-2014 == Species Pool for French Alpes project (allows for null models of entire flora)
rel09 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2009.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel10 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2010.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel11 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2011.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel12 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2012.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel13 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2013.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel14 <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/French Alpes/Species Lists/Cedric/ecrins_sky islands/fs_2014-11-14_2014.csv", as.is=T, header=T, stringsAsFactors=FALSE)
head(rel09)

## Combine all releves into one dataframe
releves <- rbind(rel09, rel10, rel11, rel12, rel13, rel14)
head(releves)
#write.csv(releves, file="data/RAW/CombinedReleveEcrin.csv")

## add a column that matche iPlant TNRS taxonomy
combinedname.releve.tnrs <- add.taxized.column(df = as.data.frame(releves), colnum = 4, spliton = " ", sepas = " ", source="iPlant_TNRS")
head(combinedname.releve.tnrs)
combinedname.releve.tnrs[16,] ## take a look to make sure it works with a known synonym 

#combinedname.releve.ncbi <- add.taxized.column(df = as.data.frame(releves), colnum = 4, spliton = " ", sepas = " ", source="NCBI")


## Count the number of times each  species is in a releve ~aboundance
releve.count <- combinedname.releve.tnrs %>% group_by(Id_station) %>% count(taxized)
dim(releve.count) #1442...1413 with taxonomy fixed ... 1390 after infraspecific taxa combined

## Species list for entire Ecrins releves: counts = number of times each species was found in a releve from 2009-2014
ecrins.releve.count <- as.data.frame(releve.count)
dim(ecrins.releve.count) #1390

#write.csv(combinedname.releve.tnrs, file="data/OutputDatasets/ecrins.releve.tnrs.csv")


############## Summit Matrix (from Julien) 
#### Summit island "communities": GPS points from each releve assigned to summit
############# filtered above tree-line (from forestry national inventory layers) and elevation (>2500)
summitMatrix <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/AlpinePD/Alpine-Sky-Islands/data/RAW/final_matrix.csv", sep=";", row.names=1, as.is=T, header = T, stringsAsFactors=FALSE)
head(summitMatrix)
dim(summitMatrix) #19 summits 363 species

summitMatrix <- t(summitMatrix)
summitMatrix[is.na(summitMatrix)] <- 0
head(summitMatrix)
summitMatrix <- as.data.frame(cbind(rownames(summitMatrix), summitMatrix)) 

combinedname.summitMatrix <- add.taxized.column(df = as.data.frame(summitMatrix), colnum = 1, spliton = "_", sepas = " ", source="iPlant_TNRS")
head(combinedname.summitMatrix)
combinedname.summitMatrix$taxized

#combinedname.summitMatrix[2:ncol(combinedname.summitMatrix)] <- sapply(combinedname.summitMatrix[2:ncol(combinedname.summitMatrix)], as.character)
summitMatrix.count <- combinedname.summitMatrix %>% group_by(taxized) %>% summarise_each(funs(max))
head(summitMatrix.count)
dim(summitMatrix.count) #322

#write.csv(summitMatrix.count, file="data/OutputDatasets/summit.count.tnrs.csv")


############## LGM count lists (from Julien) 
underIce.tmp <- read.csv('data/RAW/sp_pool_under_ice.csv')
dim(underIce.tmp) #165
persistant.tmp <- read.csv('data/RAW/sp_pool_persistant.csv')
dim(persistant.tmp) #261

LGM.tmp <- merge(persistant.tmp, underIce.tmp,  by=1, all=T) # number of releves with each species above or below LGM layer limit
LGM <- as.data.frame(LGM.tmp, , strings.as.factors=F)
rownames(LGM) <- LGM.tmp[,1]
LGM[is.na(LGM)] <- 0
colnames(LGM) <- c("species", "persistant", "underIce")
LGM <- LGM[sort(rownames(LGM)),]
dim(LGM) #298
head(LGM)
length(which(LGM$persistant != 0)) #261
length(which(LGM$underIce != 0)) #165

LGM.taxize <- add.taxized.column(df = LGM, colnum = 1, spliton = "_", sepas = " ", source = "iPlant_TNRS")
head(LGM.taxize)

#write.csv(LGM.taxize, file="data/OutputDatasets/LGM.count.tnrs.csv")

################################################# Trait Species List -> Julien ################################################# 
## List of all of the species in each dataset, before names are corrected --> email to Julien to get traits from database
head(releves)
releves.gs <- get.genus.species(releves, colnum = 3, spliton = " ", sepas = "_")

head(summitMatrix)
summitMatrix.gs <- rownames(summitMatrix)

head(LGM.tmp)
LGM.tmp.gs <- as.character(LGM.tmp[,1])

tmp <- sort(c(releves.gs, summitMatrix.gs, LGM.tmp.gs))
length(tmp)

spList.traits <- (tmp[!duplicated(tmp)]) #1584
#write.csv(spList.traits, file="data/RAW/Traits/spListEcrins.csv")

################################################# Prepare Clean Summit Matrix  ################################################# 

#######  merge a column of all species from Ecrins releves counts to pres/abs matrix from Julien
comMerge <- as.data.frame(merge(as.data.frame(ecrins.releve.count), summitMatrix.count[c(1, 4:ncol(summitMatrix.count))], by.x=1, by.y=1, all=T))
comMerge[is.na(comMerge)] <- 0
head(comMerge)
dim(comMerge) #1412   21
comMerge$taxized
tail(comMerge) ##### WHY ARE THERE NOT OBSERVATIONS FROM THE ECRINS HERE???? taxize issues?? 
comMerge$n
comMerge[1] ### SPECIES LIST

# Combine with LGM species lists
comMerge2 <- as.data.frame(merge(LGM.taxize[c(1, 4:ncol(LGM.taxize))], comMerge, by=1, all=T))
comMerge2[is.na(comMerge2)] <- 0
head(comMerge2)
tail(comMerge2)
dim(comMerge2) #1437
comMerge2$taxized

## Colapse counts into each species
comMerge2.count <- comMerge2 %>% group_by(taxized) %>% summarise_each(funs(max))
head(comMerge2.count)
dim(comMerge2.count) #1413   23
which(duplicated(comMerge2.count[1]) == T) ## Should be no duplicated names (= 0)
#write.csv(comMerge2.count, file = "data/OutputDatasets/comMerge.csv")


################################################# 0_IncludeFile: Write Species list -> includefile  ################################################# 
comMerge2.count$taxized

# Remove the genera not resolved to species
comMatrixGenus <- comMerge2.count[grepl(pattern = " NA$",comMerge2.count$taxized),]
dim(comMatrixGenus) #37 23

## The resulting species list: use for includefile, to search trait database, etc.
comMatrix <- comMerge2.count[!grepl(pattern = " NA$",comMerge2.count$taxized),]
which(duplicated(comMatrix[1]) == T) ## Should be no duplicated names (= 0)
tail(comMatrix)
comMatrix$taxized
length(comMatrix$taxized) #1376

#write.csv(comMatrix[1], file="output/0_IncludeFile/EcrinTotalSpeciesList.csv", fileEncoding = "UTF-8") ### SPECIES LIST
## saved in textwrangelr as "SpeciesPoolEcrinInclude"
#read.csv(file = "output/0_IncludeFile/EcrinTotalSpeciesList.csv")

includefile <- comMatrix$taxized #1376, the list that was used as the 'includefile'

#################################################  1_PHLAWD 
## for each alignment output for each gene region from PHLAWD, collapse intra-spacific taxa to species level, then keep longest sequence --> concatenate

source("R/ParsePHLAWD.R")
parsePHLAWD("output/1_PHLAWD/ITS/ITS.FINAL.aln.full") #1293...937
parsePHLAWD("output/1_PHLAWD/atpB/atpB.FINAL.aln.full") #144...141
parsePHLAWD("output/1_PHLAWD/matK/matK.FINAL.aln.full") #3009..819
parsePHLAWD("output/1_PHLAWD/rbcL/rbcL.FINAL.aln.full") #1043...865
parsePHLAWD("output/1_PHLAWD/trnTLF/trnTLF.FINAL.aln.full") #1115...912


#################################################  4_Trees
######## Use Congruifier (geiger) to ASSESS Trees 
### Modified from Jon's Congruify_HannahIsWonderful.R

## Concatenated Trees:
#Concat/RAxML_bipartitions.concat.EcrinSpPool.12052015.1000.nex ## This tree was not very good -- try with more stringent gap removal
genetree=read.nexus("output/4_Trees/clean75/RAxML_bipartitions.concat.EcrinSpPool.clean75.16072015.1000.nex") # saved as rooted in FigTree

tax=read.csv(file="output/5_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)
tips=sapply(genetree$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""
atol=read.tree("output/5_Scaling/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
ftax=tax[match(atol$tip.label, rownames(tax)),]
ftax[,2]="Spermatophyta"
fatol=subset(atol, ftax, "family")
phy=genetree
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
out=res$phy
congruif=out$node.label%in%res$calibrations$MRCA
out$node.label=NULL
out=nodelabel.phylo(out, tax, strict=FALSE)
out$node.label=ifelse(congruif, paste(out$node.label, "*", sep=""), out$node.label)
out$tip.label=genetree$tip.label[match(swaptips, res$phy$tip.label)] ### CHANGE tip labels back
out$tip.label=paste(genetree$tip.label[match(swaptips, res$phy$tip.label)], tax[,"family"], sep="=") ### ADD Family to tip labels

#pdf("output/5_Scaling/Figures/EcrinSpPool.12052015.1000.clean75.congruifyTaxonomy.fan.pdf", width=10, height=10) 
tree <- ladderize(out, right=F)
plot.phylo(tree, type="fan", cex=0.05, label.offset = .05) 
nodelabels(out$node.label, frame="n", col="red", cex=0.2)
#dev.off()

###help you evaluate where the tree is inconsistent with the taxonomy 
#see how the best node for lineage differs from the clade definition:
#out$FUN("Poaceae")
#missing from the clade in your tree OR unexpected within clade (but is found here in the tree)
out$FUN("Orobanchaceae")[[1]][[1]] #unexpected
[1] "739_Phelipanche"

out$FUN("Orobanchaceae")

# * = Congruified Nodes
# Noded labels with ""  implies some inconsistency between the tips expected as defined in fleshed_genera.csv and the subtended tips in the tree at the best matching node

######### Assess monophyletic clades (accuracy of PHLAWD)

## Get a taxonomic lookup table: https://github.com/wcornwell/TaxonLookup
devtools::session_info()
head(plant_lookup())
head(plant_lookup(include_counts = TRUE))
tx <- lookup_table(genetree$tip.label, by_species=TRUE)

#tx <- read.csv("~/Desktop/TaxonLookup/plant_lookup.csv")
head(tx)
dim(tx)

## Assess monophyly: MonoPhy
dat <- as.data.frame(treedata(genetree, tx)$dat)
dat <- cbind(rownames(dat), dat[1:3])
head(dat)
tr <- (treedata(genetree, tx)$phy)

mono <- AssessMonophyly(tree = tr, taxonomy = dat)
head(mono[1:5])
mono$taxonomy1$result
mono$taxonomy2$result



### Label tips that are on summmits : 
################### Label bootstrap support on unscaled ML tree
#Concat/RAxML_bipartitions.concat.EcrinSpPool.12052015.1000
treeBS <- read.tree(file="output/4_Trees/clean75/RAxML_bipartitions.concat.EcrinSpPool.clean75.16072015.1000") 
treeBS
## Find the node to root on
getMRCA(treeBS, tip=c("Equisetum_hyemale", "Huperzia_selago")) #2283

treeBSroot <- root(treeBS, node=1519, resolve.root=T)  #2283
is.rooted(treeBSroot)
treeBSroot
#write.tree(treeBSroot, file="output/4_Trees/clean75/EcrinSpPoolclean75.16072015.1000.unscaled.root.tre")
#Concat/EcrinSpPool12052015.1000.unscaled.root.tre

p1 <- character(length(treeBSroot$node.label)) # create a vector of colors that correspond to different support ranges
p1[treeBSroot$node.label >= 95] <- "turquoise3"
p1[treeBSroot$node.label < 95 & treeBSroot$node.label >= 75] <- "slategray"
p1[treeBSroot$node.label < 75] <- "red1"
p1[treeBSroot$node.label == 100] <- "green"
p1 
co = c("green", "turquoise3", "slategray", "red1")

#pdf(file="output/4_Trees/EcrinSpPool.clean75.unscaled.BS.pdf") 
#EcrinSpPool.unscaled.BS.pdf
par(mar = rep(2, 4))
plot.phylo(ladderize(treeBSroot), type="fan", cex=0.2, label.offset=.01)
nodelabels(pch = 21, cex = 0.3, col=p1, bg = p1) #p
#title(main="Unmasked")
# Printing a legend:
legend("bottomleft", legend = expression(BS == 100, BS >= 95, 75 <= BS * " < 95", BS < 75), pch = 21, pt.bg = co, pt.cex = 1, cex=.5, bty = "n")
#dev.off()


########################################################## 5_Scaling  ############################################################

##### WRITE treePL for TreeScaling: Re-run congruifier, changing scale = "NA", and writing out .treePL files

## Saved rooted in FigTree: rooted on all ferns & alies
genetree 

res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 

## WRITE OUT TREEPL FILE -- you'll need to do more than just run the exported file, but this gives you a start
# 4600 
nsites=3509 #SHOULD be number of nucleotides in alignment 
#write.treePL(res1$target, res1$calibrations, nsites, base="output/5_Scaling/clean75/EcrinSpPool.EcrinSpPool.clean75.16072015.1000.treePL", opts=list(prime=TRUE))
##had to modify .infile to adjust location of .intree

##Run on fortytwo: treePL EcrinSpPool.120515.1000.treePL.infile

#### After prime, changed end of block to estimated parameters
#PLACE THE LINES BELOW IN THE CONFIGURATION FILE
#opt = 5
#optad = 5
#optcvad = 0
#moredetailcvad

### For the clean 75
#opt = 1
#optad = 1
#moredetailad
#optcvad = 3
#moredetailcvad

### Read in treePL output to change tip labels back
treePL <- read.tree("output/5_Scaling/clean75/EcrinSpPool.EcrinSpPool.clean75.16072015.1000.treePL.dated.tre")
treePL$tip.label=genetree$tip.label[match(swaptips, res1$target$tip.label)] ### CHANGE tip labels back

#Write .tre file
#write.tree(treePL, file="output/5_Scaling/clean75/EcrinSpPool.EcrinSpPool.clean75.16072015.1000.treePL.dated.rename.tre")


########### Match bootstrap node labels from ML tree to treePL scaled tree: Save as final Tree ###############
#Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
#library(ape)

treeBSroot #<- read.tree(file="output/4_Trees/Concat/EcrinSpPool12052015.1000.unscaled.root.tre") 
treePL #<- read.tree("output/5_Scaling/EcrinSpPool.120515.1000.treePL.dated.rename.tre")
treePL #1154

##### Match bootstrap node labels from ML tree to treePL scaled tree ####
## Takes a while for large trees...
plotTreePLBoot(treePL=treePL, bootTree=treeBSroot, file="output/5_Scaling/EcrinSpPool.bootstrap.tre") 
## Deleted "Root" label from nexus in TextWrangler

########################################################### 6_Visualize Tree ######################################################### 

AlpineFinalTree <- read.tree(file="output/5_Scaling/EcrinSpPool.bootstrap.tre")
AlpineFinalTree$tip.label #note: no hyphens
#1154 tips

############################### Final Phylo-Data Object ###########################
### Check dataset: phylo + status trait
# includedile and comMatrix should match (all species lists were taxized)
dim(t(com.Alpes2))  # 1413, the communtiy matrix from analysisSkyIsl.R: merged species lists, pre-edit to remove ambiguous spcies, etc. 

##### Use Cody's code to extract the names used to search GenBank (inculdefile = taxized names) with output tip labels
###### Extract names and IDs == match taxonomy of PHLAWD output with includefile ########### 
# ./extract_names_and_ids_from_phlawd_db.py -d /Users/hannahmarx/Dropbox/Work/TankLab/Programs/PHLAWD/pln102314.db -n SpeciesPoolEcrinInclude > alpesincludefile.out20Apr2015.txt

### The table to lin input name (includefile), accepted name (GenBank == Phylogeny Tips) and NCBI_ID
extractID <- read.csv("output/0_IncludeFile/alpesincludefile.out20Apr2015.txt")
dim(extractID) #1375    4
head(extractID)
extractID$input_name <- sub("-", "_", extractID$input_name,)
extractID$input_name <- sub(" ", "_", extractID$input_name,)
extractID$accepted_name <- sub("-", "_", extractID$accepted_name,)
extractID$accepted_name <- sub(" ", "_", extractID$accepted_name,)
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = " var.*", replacement = "") 
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = " subsp.*", replacement = "") 
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = "x ", replacement = "") 
head(extractID)

length(includefile) #1376, the list that was used as the 'includefile'; -1 because ComplŽment flore
# == comMatrix$taxized == extractID$input_name
length(extractID$input_name) #1375
(extractID$input_name %in% rownames(t(com.Alpes2))) ## ALL TRUE: GOOD

### Pre-taxized tree:
temp <- treedata(AlpineFinalTree, t(com.Alpes2)) ##### NEED TO CHANGE NAMES OF DATA TO PHLAWD TAXONOMY TO MATCH TREE
phy <- temp$phy
phy
data <- as.data.frame(temp$data)
head(data)
dim(data) #1107
1107/1375 #0.8050909  % OF THE TIPS BEFORE TIP TAXONOMY CORRECTED TO MATCH COMMUNITY MATRIX

######### Re-label Tree Tips to match Taxonomy Harmony of other datasets (and includefile)
## The matching species between the phylogney tips and the input names 
length(AlpineFinalTree$tip.label) #1154
length(AlpineFinalTree$tip.label[AlpineFinalTree$tip.label %in% extractID$input_name]) #1107  
##### NEED TO CHANGE THE NAMES OF TIP LABELS TO MATCH THE TAXONOMY OF THE COMMUNITY MATRIX/TRAITS

### Relabel tips
phy.tmp <- AlpineFinalTree #Rename so don't overwrite

for (i in 1:length(phy.tmp$tip.label)){
  if (!phy.tmp$tip.label[i] %in% extractID$input_name){
    if (phy.tmp$tip.label[i] %in% extractID$accepted_name){
      x <- extractID[which(extractID$accepted_name==phy.tmp$tip.label[i]), "input_name"]
      y <- extractID[which(extractID$accepted_name==phy.tmp$tip.label[i]), "accepted_name"]
      print(paste(i, "input =", x, ", accepted =", y))
      phy.tmp$tip.label[i] <- x
      
    }  
}
}
# Output saved as: TipLabelsChanged.rtf N =32
#write.tree(phy.tmp, file="output/6_Visualize Tree/phy.Alpes.taxized.tre")
#write.tree(phy.tmp, file="data/AnalysesDatasets/phy.Alpes.taxized.tre")

### Prune tree to include just species in community data: don't actually need to do this (pez)
treedata.fin <- treedata(phy.tmp, t(com.Alpes2)) ##### CHANGED NAMES OF TREE TO MATCH DATASET TAXONOMY 
AlpineFinalTree.taxized <- treedata.fin$phy #1140
comMatrix.tredat <- as.data.frame(treedata.fin$data) 
1140/1375 #0.8290909 % Ecrin Species pool: of the species that were retrieved from GenBank using the PHLAWD approach

comMatrix.tredat["Summits"]
dim(filter(comMatrix.tredat, Summits > 0)) #231
231/322 # 0.7173913 % of summits retrived

##### Summaries: 
includefile <- extractID$input_name

### Check: The number of matching species between includefile and spcies in Community Matrix
length(includefile[includefile %in% colnames(CommunityMatrixFinal)]) #1375/1375 

### The matching species between include file and accepted names == GenBank & Phylogney tips
length(includefile[includefile %in% extractID$accepted_name]) #1112/1375

## The matching species in the phylogney tips and the accepted names 
length(AlpineFinalTree.taxized$tip.label[AlpineFinalTree.taxized$tip.label %in% extractID$input_name]) #1140/1375 #0.8283636 % of speices retrived using PHLAWD

### Those without [, 2:3] don't have info in genbank
dim(extractID[is.na(extractID$ncbi_id), ]) #222  
222/1375 # 0.1614545 % don't have sequences in GenBank

## check by match names of inputfile not in PHLAWD output
dim(extractID[!extractID$accepted_name %in% extractID$input_name, ]) #255 names do not mach, 222 of those are just not in GenBank
#1139 + 255 = ; the include file length

# List of names with confilcting taxonomy
conflict.names <- na.omit(extractID[!extractID$accepted_name %in% extractID$input_name, ])
dim(conflict.names) #33...should equal TipLabelsChange, N= 33 GOOD!!

#write.tree(AlpineFinalTree.taxized, file="output/6_Visualize Tree/phy.comTaxized.Alpes.tre")
#write.tree(AlpineFinalTree.taxized, file="data/Extra/phy.comTaxized.Alpes.tre")

########################### Plot Tree ###########################################
###### Color taxa Native / Invasive: Code modified from picante color.plot.phylo
##Plot phylo with calibrated nodes IDed and nodes lables with taxonomy (dated with PATHd8)

######## Nodes Object ################
#### Need to re-run Congruifier to get cal object (calibrated nodes)
genetree=pezAlpes$phy

tax=read.csv(file="output/5_Scaling/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) 
head(tax)
tips=sapply(genetree$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
head(SJ_tax)
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""
atol=read.tree("output/5_Scaling/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
ftax=tax[match(atol$tip.label, rownames(tax)),]
head(ftax)
ftax[,2]="Spermatophyta"
fatol=subset(atol, ftax, "family")
phy=genetree
swaptips=paste(1:length(tips),tips,sep="_")
phy$tip.label=swaptips
tax=SJ_tax
rownames(tax)=swaptips
res=congruify.phylo(fatol, phy, tax, tol=0, scale="PATHd8") # need to use PATHd8 to get res$phy
res
cal=res$calibrations
#### Vector of congruified nodes
vec=mrcaID(phy, cal) ## get vector for calibrated nodes
sum(vec)==nrow(cal) ## check on whether the function is working appropriately
# plot.phylo(phy, type="fan", show.tip=FALSE, edge.width=0.1)
## plot box at node only if calibrated
# nodelabels(text=NULL, cex=ifelse(vec==0, NA, 2), frame="n", bg="lightskyblue", col="lightgray", pch=22)


##### Vector of BS supports 
p2 <- character(length(pezAlpes$phy$node.label)) # create a vector of colors that correspond to different support ranges 
p2[] <- "#0000ff00" ## Transparent color: "#RRGGBBAA" and the AA portion is the opacity/trasparency.
p2[pezAlpes$phy$node.label >=  95] <- "black"
p2[pezAlpes$phy$node.label ==  100] <- "black"
p2[pezAlpes$phy$node.label < 95 & pezAlpes$phy$node.label >= 75] <- "slategray4"
#p2[treePLboots$node.label < 95 & treePLboots$node.label >= 75] <- "slategray"
#p2[treePLboots$node.label < 75] <- "red1"
p2[pezAlpes$phy$node.label == ""] <- "black" # node = 353 = 100
p2 

#### From SanJuans_TreeVisFinal.R 
colorsWes <- c("#F21A00", "#E1AF00")
colorsIsl <- c("cyan3", "chocolate1")

#pdf(file="output/6_Visualize Tree/VascularFloraEcrin.pdf") 
color.plot.phylo3NoTrait(phy.Alpes.drop, data.Alpes.drop.pa, taxa.names = "phyloTips", trait = "Summits", col.names = colorsIsl, label.offset=5.2, phy.cex=.05, cut.labs =c( "Ecrin NP", "Alpine Summits"), leg.cex=.8, main="Vascular Flora of the Le Parc national des Ecrins")
#dev.off()


pdf(file="output/4_Trees/phy.dots.pdf")
plot(pezAlpes$phy, show.tip.label = FALSE, type="fan", )
tiplabels(tip = which(pezAlpes$comm["Summits",] > 0), pch = 19, cex = .25, col = "darkorange2")
dev.off()

dat.plot[dat.plot > 0] <- 1
dat.plot <- as.data.frame(dat.plot)
head(dat.plot)

#pdf(file="output/4_Trees/phy.dots.v2.pdf")
trait.plot(pezAlpes$phy, dat.plot, cols = list(dat.plot = c("white", "darkorange2")), font = 0, cex.lab = .0001)
nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
#dev.off()

