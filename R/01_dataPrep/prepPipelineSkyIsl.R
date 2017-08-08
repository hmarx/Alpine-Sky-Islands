#####################################################################################################################
############# Prep datasets for phylogenetic diveristy of alpine summits in the Ecrins National Park, France ########
############# Hannah E. Marx, 25 April 2017 ###########################################################################
#####################################################################################################################

source("R/functions/treeFunctions.R")
source("R/functions/TaxonomyHarmony.R")
source("R/functions/annotate.fasta.R")

####################################### Read in Raw datasets, homogonize taxonomy ####################################### 
## Process occurence data from Ecrins NP, species on summits, and species that currently occur in areas that were exposed through LGM
## Uses iPlant taxonomy to make taxonomy consistent between datasets 

############# ALL RELEVES : from Ecrin National Park, 2009-2014 == Species Pool for French Alpes project (allows for null models of entire flora)
rel09 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2009.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel10 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2010.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel11 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2011.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel12 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2012.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel13 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2013.csv", as.is=T, header=T, stringsAsFactors=FALSE)
rel14 <- read.csv("data/RAW/Dentant_releves/fs_2014-11-14_2014.csv", as.is=T, header=T, stringsAsFactors=FALSE)
head(rel09)

## Combine all releves into one dataframe
releves <- rbind(rel09, rel10, rel11, rel12, rel13, rel14)
head(releves)
#write.csv(releves, file="data/RAW/CombinedReleveEcrin.csv")

## add a column that matches iPlant TNRS taxonomy
combinedname.releve.tnrs <- add.taxized.column(df = as.data.frame(releves), colnum = 4, spliton = " ", sepas = " ", source="iPlant_TNRS")
head(combinedname.releve.tnrs)
combinedname.releve.tnrs[16,] ## take a look to make sure it works with a known synonym 

## Count the number of times each  species is in a releve (~abundance)
releve.count <- combinedname.releve.tnrs %>% group_by(Id_station) %>% count(taxized)
dim(releve.count) #1442...1413 with taxonomy fixed ... 1390 after infraspecific taxa combined

## Species list for entire Ecrins releves: counts = number of times each species was found in a releve from 2009-2014
ecrins.releve.count <- as.data.frame(releve.count)
dim(ecrins.releve.count) #1390 species
#write.csv(combinedname.releve.tnrs, file="data/OutputDatasets/ecrins.releve.tnrs.csv")

############## Summit Matrix (from Julien, Jan 19 2015) 
#### Summit island "communities": GPS points from each releve assigned to summit
############# filtered above tree-line (from forestry national inventory layers) and elevation (>2500)
summitMatrix <- read.csv("data/RAW/final_matrix.csv", sep=";", row.names=1, as.is=T, header = T, stringsAsFactors=FALSE)
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
dim(summitMatrix.count) #322 species
#write.csv(summitMatrix.count, file="data/OutputDatasets/summit.count.tnrs.csv")

############## LGM count lists (from Julien) 
underIce.tmp <- read.csv('data/RAW/sp_pool_under_ice.csv')
dim(underIce.tmp) #165 species
persistent.tmp <- read.csv('data/RAW/sp_pool_persistent.csv')
dim(persistent.tmp) #261 species

LGM.tmp <- merge(persistent.tmp, underIce.tmp,  by=1, all=T) # number of releves with each species above or below LGM layer limit
LGM <- as.data.frame(LGM.tmp, strings.as.factors=F)
rownames(LGM) <- LGM.tmp[,1]
LGM[is.na(LGM)] <- 0
colnames(LGM) <- c("species", "persistent", "underIce")
LGM <- LGM[sort(rownames(LGM)),]
dim(LGM) #298
head(LGM)
length(which(LGM$persistent != 0)) #261
length(which(LGM$underIce != 0)) #165

LGM.taxize <- add.taxized.column(df = LGM, colnum = 1, spliton = "_", sepas = " ", source = "iPlant_TNRS")
head(LGM.taxize)
#write.csv(LGM.taxize, file="data/OutputDatasets/LGM.count.tnrs.csv")


################################################# Prepare Clean Summit Matrix  ################################################# 

#######  merge a column of all species from Ecrins releves counts to pres/abs summit matrix from Julien
comMerge <- as.data.frame(merge(as.data.frame(ecrins.releve.count), summitMatrix.count[c(1, 4:ncol(summitMatrix.count))], by.x=1, by.y=1, all=T))
comMerge[is.na(comMerge)] <- 0
head(comMerge)
dim(comMerge) #1413   21
comMerge$taxized
tail(comMerge) 
comMerge$n
comMerge[1] ### SPECIES LIST

# Combine with LGM species lists
comMerge2 <- as.data.frame(merge(LGM.taxize[c(1, 4:ncol(LGM.taxize))], comMerge, by=1, all=T))
comMerge2[is.na(comMerge2)] <- 0
head(comMerge2)
tail(comMerge2)
dim(comMerge2) #1438
comMerge2$taxized

## Colapse counts into each species
comMerge2.count <- comMerge2 %>% group_by(taxized) %>% summarise_each(funs(max))
head(comMerge2.count)
dim(comMerge2.count) #1414   23
which(duplicated(comMerge2.count[1]) == T) ## Should be no duplicated names (= 0)
#write.csv(comMerge2.count, file = "data/OutputDatasets/comMerge.csv")
#Community matrix: 1413 species,   22 “communities” (persistent, under ice, ecrins NP, 19 summits)
#file = comMerge.csv
#—Removed the genera not resolved to species from community matrix (N=37) for PHLAWD search
#—removed from includefile: ComplŽment flore

#PHLAWD species list:  species names from comMerge.csv
#EcrinTotalSpeciesList.csv —> SpeciesPoolEcrinInclude (N=1376)


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
length(comMatrix$taxized) #1377

#write.csv(comMatrix[1], file="output/0_IncludeFile/EcrinTotalSpeciesList.csv", fileEncoding = "UTF-8") ### SPECIES LIST
## saved in textwrangelr as "SpeciesPoolEcrinInclude"
#read.csv(file = "output/0_IncludeFile/EcrinTotalSpeciesList.csv")

includefile <- comMatrix$taxized 
length(includefile) #1376 the list that was used as the 'includefile'


#################################################  1_PHLAWD ################################################# 
####### Spermatophyta ####### 
## Restrict search of includefile and .keep files to just spermatophyta 
## for each alignment output for each gene region from PHLAWD, collapse intra-spacific taxa to species level, then keep longest sequence --> concatenate

source("R/ParsePHLAWD.R")
parsePHLAWD("output/1_PHLAWD/atpB/atpB.FINAL.aln.full") #119...118
parsePHLAWD("output/1_PHLAWD/ITS/ITS.FINAL.aln.full") #1202...910
parsePHLAWD("output/1_PHLAWD/matK/matK.FINAL.aln.full") #964...807
parsePHLAWD("output/1_PHLAWD/rbcL/rbcL.FINAL.aln.full") #962...821
parsePHLAWD("output/1_PHLAWD/trnTLF/trnTLF.FINAL.aln.full") #1045...875

#################################################  2_Remove 
## look at each gene tree, removed incorrect/outiler sequences, then re-clean and concatenate 

parseREMOVED("output/2_Remove/atpB.unique.GB.rem2.fasta") #114
parseREMOVED("output/2_Remove/ITS.unique.GB.fasta.rem") #908
parseREMOVED("output/2_Remove/trnTLF.unique.GB.fasta.rem") #870
parseREMOVED("output/2_Remove/matK.unique.GB.fasta.rem") #806
parseREMOVED("output/2_Remove/rbcL.unique.GB.fasta.rem") #820

## save GenBank Identifiers for uniques sequences kept and used for each species  
source("R/parseGBinfo.R")
ITS <- parseGBinfo(gene = "ITS", unique.GB.fasta = "output/2_Remove/ITS.unique.GB.fasta.rem") #908
tnrtlfGB <- parseGBinfo(gene = "trnTLF", unique.GB.fasta = "output/2_Remove/trnTLF.unique.GB.fasta.rem") #870
rbclGB <- parseGBinfo(gene = "rbcL", unique.GB.fasta = "output/2_Remove/rbcL.unique.GB.fasta.rem") #820
matkGB <- parseGBinfo(gene = "matK", unique.GB.fasta = "output/2_Remove/matK.unique.GB.fasta.rem") #806
atpbGB <- parseGBinfo(gene = "atpB", unique.GB.fasta = "output/2_Remove/atpB.unique.GB.rem2.fasta") #114

GB_accessions <- mergeGBinfo(ITS, tnrtlfGB, rbclGB, matkGB, atpbGB)
dim(GB_accessions) #1259 taxa, 5 gene regions
str(GB_accessions)
#write.csv(GB_accessions, "output/2_Remove/GB_accessions.csv")

#################################################  5_Trees
### Tried both vasular plants (Tracheophyta) and seed plants (Spermatophyta)
### Used Spermatophyta for analyses 

######## Use Congruifier (geiger) to assess quality of phylogeny estimation 
### Modified from code by Jon Eastman

## Concatenated Trees:
genetree=read.nexus("output/5_Trees/Concatenated/RAxML_bipartitions.align.concat.EcrinSpPool.121115.phy.1000.nex") # saved as rooted in FigTree
tax=read.csv(file="data/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

tips=sapply(genetree$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips, rownames(tax))
SJ_tax=tax[ll,]
rownames(SJ_tax)=names(tips)
SJ_tax=as.matrix(SJ_tax)
SJ_tax[is.na(SJ_tax)]=""
atol=read.tree("data/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
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

#pdf("output/5_Trees/Concatenated/Figures/EcrinSpPool.121115.1000.nex.congruifyTaxonomy.fan.pdf", width=10, height=10) 
#tree <- ladderize(out, right=F)
#plot.phylo(tree, type="fan", cex=0.05, label.offset = .05) 
#nodelabels(out$node.label, frame="n", col="red", cex=0.2)
#dev.off()

#pdf("output/5_Trees/Concatenated/Figures/EcrinSpPool.121115.1000.nex.congruifyTaxonomy.pdf") 
#tree <- ladderize(out, right=F)
#plot.phylo(edge.width = 0.25, tree, cex=0.05, label.offset = .05) 
#nodelabels(out$node.label, frame="n", col="red", cex=0.2)
#dev.off()

###help you evaluate where the tree is inconsistent with the taxonomy 
#see how the best node for lineage differs from the clade definition:
#out$FUN("Poaceae")
#missing from the clade in your tree OR unexpected within clade (but is found here in the tree)
out$FUN("Orobanchaceae")[[1]][[1]] #unexpected
# [1] "739_Phelipanche"

out$FUN("Orobanchaceae")
# * = Congruified Nodes
# Noded labels with ""  implies some inconsistency between the tips expected as defined in fleshed_genera.csv and the subtended tips in the tree at the best matching node


############ Label bootstrap support on unscaled ML tree
treeBS <- read.tree(file="output/5_Trees/Concatenated/RAxML_bipartitions.align.concat.EcrinSpPool.121115.phy.1000") 
treeBS #1084
## Find the node to root on
getMRCA(treeBS, tip=c("Taxus_baccata", "Pinus_nigra")) #2062

treeBSroot <- root(treeBS, node=2062, resolve.root=T)  
is.rooted(treeBSroot)
treeBSroot
#write.tree(treeBSroot, file="output/5_Trees/Concatenated/EcrinSpPool.121115.1000.unscaled.root.tre")

p1 <- character(length(treeBSroot$node.label)) # create a vector of colors that correspond to different support ranges
p1[treeBSroot$node.label >= 95] <- "turquoise3"
p1[treeBSroot$node.label < 95 & treeBSroot$node.label >= 75] <- "slategray"
p1[treeBSroot$node.label < 75] <- "red1"
p1[treeBSroot$node.label == 100] <- "green"
p1 
co = c("green", "turquoise3", "slategray", "red1")

pdf(file="output/5_Trees/Concatenated/Figures/EcrinSpPool.unscaled.BS.pdf") 
#EcrinSpPool.unscaled.BS.pdf
par(mar = rep(2, 4))
plot.phylo(ladderize(treeBSroot), type="fan", cex=0.2, label.offset=.01)
nodelabels(pch = 21, cex = 0.3, col=p1, bg = p1) #p
#title(main="Unmasked")
# Printing a legend:
legend("bottomleft", legend = expression(BS == 100, BS >= 95, 75 <= BS * " < 95", BS < 75), pch = 21, pt.bg = co, pt.cex = 1, cex=.5, bty = "n")
dev.off()


########################################################## 6_Scaling  ############################################################

##### WRITE treePL for TreeScaling: Re-run congruifier, changing scale = "NA", and writing out .treePL files

## Saved rooted in FigTree: rooted on all ferns & alies
genetree 

res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 

## WRITE OUT TREEPL FILE -- you'll need to do more than just run the exported file, but this gives you a start
nsites=5144 #SHOULD be number of nucleotides in alignment 
write.treePL(res1$target, res1$calibrations, nsites, base="output/6_Scaling/RAxML_bipartitions.align.concat.EcrinSpPool.121115.1000.nex.treePL", opts=list(prime=TRUE))
## modify .infile to adjust location of .intree and outfile (local)
## Run on fortytwo: treePL EcrinSpPool.120515.1000.treePL.infile

#### After prime, changed end of block to estimated parameters
#PLACE THE LINES BELOW IN THE CONFIGURATION FILE
#opt = 5
#optad = 5
#optcvad = 3
#moredetailcvad

### Read in treePL output to change tip labels back
treePL <- read.tree("output/6_Scaling/RAxML_bipartitions.align.concat.EcrinSpPool.121115.1000.nex.treePL.dated.tre")
treePL$tip.label=genetree$tip.label[match(swaptips, res1$target$tip.label)] ### CHANGE tip labels back

#Write .tre file
#write.tree(treePL, file="output/6_Scaling/EcrinSpPool.121115.1000.treePL.dated.rename.tre")


########### Match bootstrap node labels from ML tree to treePL scaled tree: Save as final Tree ###############
#Modified from: http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
#library(ape)

treeBSroot #<- read.tree(file="output/5_Trees/Concatenated/EcrinSpPool.121115.1000.unscaled.root.tre") 
treePL #<- read.tree("output/6_Scaling/RAxML_bipartitions.align.concat.EcrinSpPool.121115.1000.nex.treePL.dated.tre")
treePL #1084

ltt(treePL)
max(treePL$edge.length) #227.2309
     
##### Match bootstrap node labels from ML tree to treePL scaled tree ####
## Takes a while for large trees...
plotTreePLBoot(treePL=treePL, bootTree=treeBSroot, file="output/6_Scaling/EcrinSpPool.12112015.dated.bootstrap.tre") 
## Deleted "Root" label from nexus in TextWrangler 

########################################################### 6_VisualizeTree ######################################################### 

AlpineFinalTree <- read.tree(file="output/6_Scaling/EcrinSpPool.12112015.dated.bootstrap.tre")
AlpineFinalTree$tip.label #note: no hyphens
#1084 tips

############################### Final Phylo-Data Object ###########################
# Run previous analysisSkyIsl.R which used the TNRS taxonomy 
source("R/01_dataPrep/analysisSkyIsl_TNRS.R")

### Check dataset: phylo + status trait
# includedile and comMatrix should match (all species lists were taxized)
dim(t(alps.TNRS.sites))  # 1414, the communtiy matrix from analysisSkyIsl.R: merged species lists, pre-edit to remove ambiguous spcies, etc. 

##### Use Cody's code to extract the names used to search GenBank (inculdefile = taxized names) with output tip labels
###### Extract names and IDs == match taxonomy of PHLAWD output with includefile ########### 
# ./extract_names_and_ids_from_phlawd_db.py -d /Users/hannahmarx/Dropbox/Work/Programs/PHLAWD/pln170915.db -n SpeciesPoolEcrinInclude > alpesincludefile.out14Nov2015.txt

### The table to link input name (includefile), accepted name (GenBank == Phylogeny Tips) and NCBI_ID
extractID <- read.csv("output/0_IncludeFile/alpesincludefile.out14Nov2015.txt")
dim(extractID) #1376    4
head(extractID)
extractID$input_name <- sub("-", "_", extractID$input_name,)
extractID$input_name <- sub(" ", "_", extractID$input_name,)
extractID$accepted_name <- sub("-", "_", extractID$accepted_name,)
extractID$accepted_name <- sub(" ", "_", extractID$accepted_name,)
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = " var.*", replacement = "") 
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = " subsp.*", replacement = "") 
extractID$accepted_name <- gsub(extractID$accepted_name, pattern = "x ", replacement = "") 
head(extractID)

length(includefile) #1377, the list that was used as the 'includefile'; -1 because ComplŽment flore
# == comMatrix$taxized == extractID$input_name
length(extractID$input_name) #1376
(extractID$input_name %in% rownames(t(alps.sites))) ## ALL TRUE: GOOD
extractID[!(extractID$input_name %in% rownames(t(alps.sites))),] ## ALL TRUE: GOOD

tmp <- as.data.frame(t(alps.sites), stringsAsFactors = F)
head(tmp)
tax=read.csv(file="data/Congruify/fleshed_genera.csv", as.is=TRUE, row=1, stringsAsFactors = F) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

tips.ecrins.sperma=sapply(extractID$input_name, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips.ecrins.sperma, rownames(tax))
ecrins_tax.sperma=tax[ll,]
rownames(ecrins_tax.sperma)=names(tips.ecrins.sperma)
ecrins_tax.sperma=as.matrix(ecrins_tax.sperma)
ecrins_tax.sperma[is.na(ecrins_tax.sperma)]=""
head(ecrins_tax.sperma)
length(which(ecrins_tax.sperma[,"Spermatophyta"] == "Spermatophyta")) # 1288 species are in Spermatophyta
ecrins.include.sperma <- names(which(ecrins_tax.sperma[,"Spermatophyta"] == "Spermatophyta"))
ecrins.include.sperma.extract <- (extractID[extractID$input_name %in% ecrins.include.sperma,])

tips.sperma=sapply(rownames(tmp[tmp["Summits"] > 0,]), function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips.sperma, rownames(tax))
sperma_tax=tax[ll,]
rownames(sperma_tax)=names(tips.sperma)
sperma_tax=as.matrix(sperma_tax)
sperma_tax[is.na(sperma_tax)]=""
(sperma_tax)
length(which(sperma_tax[,"Spermatophyta"] == "Spermatophyta")) ## 298 summits species are in Spermatophyta

tips.persist=sapply(rownames(tmp[tmp["Persistent"] > 0,]), function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips.persist, rownames(tax))
persis_tax=tax[ll,]
rownames(persis_tax)=names(tips.persist)
persis_tax=as.matrix(persis_tax)
persis_tax[is.na(persis_tax)]=""
(persis_tax)
length(which(persis_tax[,"Spermatophyta"] == "Spermatophyta")) ## 222 persistent species are in Spermatophyta


### Pre-taxized tree:
temp <- treedata(AlpineFinalTree, t(alps.TNRS.sites)) ##### NEED TO CHANGE NAMES OF DATA TO PHLAWD TAXONOMY TO MATCH TREE
phy <- temp$phy
phy
data <- as.data.frame(temp$data)
head(data)
dim(data) #1049
1049/1288 #0.814441  % OF THE TIPS BEFORE TIP TAXONOMY CORRECTED TO MATCH COMMUNITY MATRIX

######### Re-label Tree Tips to match Taxonomy Harmony of other datasets (and includefile)
## The matching species between the phylogney tips and the input names 
length(AlpineFinalTree$tip.label) #1084
length(AlpineFinalTree$tip.label[AlpineFinalTree$tip.label %in% extractID$input_name]) #1049  
##### NEED TO CHANGE THE NAMES OF TIP LABELS TO MATCH THE TAXONOMY OF THE COMMUNITY MATRIX/TRAITS

### Relabel tips
phy.tmp <- AlpineFinalTree #Rename so don't overwrite

for (i in 1:length(phy.tmp$tip.label)){
  if (!phy.tmp$tip.label[i] %in% extractID$input_name){ #if the tip label has been corrected to accepted in NCBI (PHLAWD)
    if (phy.tmp$tip.label[i] %in% extractID$accepted_name){
      x <- extractID[which(extractID$accepted_name==phy.tmp$tip.label[i]), "input_name"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
      y <- extractID[which(extractID$accepted_name==phy.tmp$tip.label[i]), "accepted_name"] #name corected by NCBI
      print(paste(i, "input =", x, ", accepted =", y))
      phy.tmp$tip.label[i] <- x #replace tip with searched name
      
    }  
}
}

#"877 input = Deschampsia_flexuosa , accepted = Avenella_flexuosa"
phy.tmp #1084 
"Avenella_flexuosa" %in% phy.tmp.taxized$tip.label
"Capsella_bursa_pastoris" %in% phy.tmp.taxized$tip.label

#phy.tmp.taxized <- add.taxized.tips(phy = phy.tmp, spliton = "_", sepas = " ", taxonomy.source = "iPlant_TNRS")
#actually don't do this -- just use code above b/c taxonomy was changed to iPlant in community matrix first

# Output saved as: output/6_Visualize Tree/TipLabelsChanged.rtf N =32
#write.tree(phy.tmp, file="output/6_VisualizeTree/phy.Alpes.taxized.tre")
#write.tree(phy.tmp, file="data/AnalysesDatasets/phy.Alpes.taxized.tre")

### Prune tree to include just species in community data: don't actually need to do this (pez)
treedata.fin <- treedata(phy.tmp, t(alps.TNRS.sites)) ##### CHANGED NAMES OF TREE TO MATCH DATASET TAXONOMY 
#The following tips were not found in 'data' and were dropped from 'phy':
#Malus_x
#Phelipanche_purpurea
#Picris_sp

AlpineFinalTree.taxized <- treedata.fin$phy #1081
comMatrix.tredat <- as.data.frame(treedata.fin$data) 
dim(comMatrix.tredat) # 1081 species 19 communtiies
1081/1288 #0.8392857 % Ecrin Species pool: of the species that were retrieved from GenBank using the PHLAWD approach

dim(filter(comMatrix.tredat, comMatrix.tredat["Summits"] > 0)) #215 species on summits
215/298 # 0.7214765 % of summits retrieved


##### Summaries: 
length(ecrins.include.sperma) # species in the includefile for the PHLAWD serach that are seedplants 
# 1288
ecrins.include.sperma.extract

### The matching species between include file and accepted names == GenBank & Phylogney tips
length(ecrins.include.sperma[ecrins.include.sperma %in% extractID$accepted_name]) #1074/1288

## The matching species in the phylogney tips and the accepted names 
length(AlpineFinalTree.taxized$tip.label[AlpineFinalTree.taxized$tip.label %in% extractID$input_name]) #1081/1288 #0.8392857 % of speices retrived using PHLAWD

### Those without [, 2:3] don't have info in genbank
dim(ecrins.include.sperma.extract[is.na(ecrins.include.sperma.extract$ncbi_id), ]) #174  
174/1288 # 0.1350932 % don't have sequences in GenBank

## check by match names of inputfile not in PHLAWD output
dim(ecrins.include.sperma.extract[!ecrins.include.sperma.extract$accepted_name %in% ecrins.include.sperma.extract$input_name, ]) #206 names do not mach, 174 of those are just not in GenBank
#1074 + 206 = 1280; should be the include file length
# 206 - 32 (taxonomy updates on phylo) = 174 == true missing species 

# List of names with confilcting taxonomy, which were changed on phy.comTaxized.Alpes.tre
conflict.names <- na.omit(ecrins.include.sperma.extract[!ecrins.include.sperma.extract$accepted_name %in% ecrins.include.sperma.extract$input_name, ])
dim(conflict.names) #32...should equal TipLabelsChange, N= 32 GOOD!!



###############################################################################################################
###############################################################################################################
###################### Edit taxonomy of species in tree, community matrix, gene trees, and alignments to match TaxRef v.9


################################ Match Final Community Matrix Taxonomy with TaxRef v.9 ######################################################### 

# Last taxonomy table commuity matrix that was used for all analyses
v1tax <- com.endemic.TNRS.tmp
head(v1tax)
dim(v1tax) #1414

# Updated taxonomy to TaxRef V.9
taxref <- read.csv("data/TaxRef/SkyIslandsSPList_8May2017.csv", row.names = 1, header=T) #SkyIslandsSPList_24Apr2017.csv
head(taxref)
dim(taxref)

# Merge and remove NA's from both datasets (mostly genera without species ep.)
v2tax <- merge(taxref[4], v1tax, by.x=0, by.y=1)
dim(v2tax) #1376
#write.csv(v2tax, "data/TaxRef/v2tax.csv")

head(v2tax)

# Merge phylogeny to see 
length(pezalps.TNRS$phy$tip.label) #1081
v3tax <- merge(cbind(as.data.frame(pezalps.TNRS$phy$tip.label), rep(1, length(pezalps.TNRS$phy$tip.label))), v2tax, by.x=1, by.y=1, all.x=T, all.y=T)
dim(v3tax) 
tail(v3tax)
duplicated(v3tax[1])
#write.csv(v3tax, "data/TaxRef/v3tax.csv")

### Do some manual editing to remove blank cells
### Removed: Chara_Linnaeus; Compl_ment_flore; Zz_Lichens
### Added TaxRef from SkyIslandsSPList_CDcorrectedvs2.dbf, and copied origin of introduced species
# Galium_commune??
# Hieracium_inuloides
# Hieracium_pilosum
# Leucanthemum_atratum

v4tax <- read.csv("data/TaxRef/v3tax_edit.csv", row.names = 1, header=T) #v3tax_edit.csv
dim(v4tax) #1373
dim(v4tax[v4tax[2] == 1,]) #1081 == the tips in TNRS phylo

# Exptract dupliacted taxonomies
taxref.dup <- v4tax[duplicated(v4tax$Valid_LIB.C.150),]
taxref.dup <- taxref.dup[order(taxref.dup$Valid_LIB.C.150),]
dim(taxref.dup) #26
#write.csv(taxref.dup, "data/TaxRef/taxref.dup.csv")

# Sum duplicated speices rows -> remove duplicated
str(v4tax)
v4tax$rep.1..length.pezalps.TNRS.phy.tip.label..[is.na(v4tax$rep.1..length.pezalps.TNRS.phy.tip.label..)] <- 0
v5tax <- v4tax[-c(1)] %>% group_by(Valid_LIB.C.150, Alien_origin) %>% summarise_each(funs(max))
dim(v5tax) #1345
#write.csv(v5tax, "data/TaxRef/v5tax.csv") ## This will be the new community matrix

#1347	308	235
#1072	217	168
#0.795842613	0.704545455	0.714893617

## remove duplicates 
uniquTax <- v4tax[!duplicated(v4tax$Valid_LIB.C.150),]
uniquTax[duplicated(uniquTax$tip.label)]
v6tax <- merge(uniquTax[c(1,3)], v5tax, by.x=2, by.y=1)
dim(v6tax) #1345
duplicated(v6tax$Valid_LIB.C.150)
v6tax[duplicated(v6tax$Valid_LIB.C.150)]
#write.csv(v6tax, "data/TaxRef/v6tax.csv")



######## Relabel Tips to match TaxRef v9

### Relabel tips on phylogeny to match TaxRef v.9; drop duplicate tips
phy.tmp.taxref <- phy.tmp #Rename so don't overwrite

for (i in 1:length(phy.tmp.taxref$tip.label)){
  #if (!phy.tmp.taxref$tip.label[i] %in% v6tax$Valid_LIB.C.150){ #if the tip labe has been corrected TaxRef v9
    if (phy.tmp.taxref$tip.label[i] %in% v6tax$pezalps.TNRS.phy.tip.label){
      x <- v6tax[which(v6tax$pezalps.TNRS.phy.tip.label==phy.tmp.taxref$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
      y <- v6tax[which(v6tax$pezalps.TNRS.phy.tip.label==phy.tmp.taxref$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
      #print(paste(i, "input =", x, ", accepted =", y))
      phy.tmp.taxref$tip.label[i] <- as.character(y)#replace tip with taxref name
    } else {
    print(phy.tmp.taxref$tip.label[i])
  }
}

[1] "Picris_sp"
[1] "Erigeron_acris"
[1] "Scorzoneroides_autumnalis"
[1] "Campanula_cenisia"
[1] "Clinopodium_alpinum"
[1] "Phelipanche_purpurea"
[1] "Gentiana_verna"
[1] "Rosa_dumalis"
[1] "Malus_x"
[1] "Ulmus_minor"
[1] "Trifolium_badium"
[1] "Trifolium_campestre"
[1] "Poa_bulbosa"
[1] "Carex_sylvatica"
[1] "Carex_montana"
[1] "Carex_digitata"
[1] "Carex_echinata"
[1] "Luzula_sylvatica"

### edited v6tax with names changed from NCBI alpesincludefile.out14Nov2015.txt
v6tax_edit <- read.csv("data/TaxRef/v6tax_edit.csv")
phy.tmp.taxref <- phy.tmp #Rename so don't overwrite

for (i in 1:length(phy.tmp.taxref$tip.label)){
  #if (!phy.tmp.taxref$tip.label[i] %in% v6tax$Valid_LIB.C.150){ #if the tip labe has been corrected TaxRef v9
  if (phy.tmp.taxref$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.taxref$tip.label[i]), "pezAlpes.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.taxref$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.taxref$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}

phy.tmp.taxref$tip.label
phy.tmp.taxref$tip.label[duplicated(phy.tmp.taxref$tip.label)] #0
"Malus_x" %in% phy.tmp.taxref$tip.label

phy.tmp.taxref.dropdup <- drop.tip(phy = phy.tmp.taxref, tip = phy.tmp.taxref$tip.label[duplicated(phy.tmp.taxref$tip.label)])
#write.tree(phy.tmp.taxref.dropdup, file="data/TaxRef/phy.Alpes.taxized.taxref.tre")


################################################# 
extractID <- read.csv("output/0_IncludeFile/alpesincludefile.out5May2017.txt")
dim(extractID) #1376    4
extractID <- cbind(phylo.tip = extractID$accepted_name, extractID)
extractID$input_name <- sub("-", "_", extractID$input_name)
extractID$input_name <- sub(" ", "_", extractID$input_name)
extractID$accepted_name <- str_replace_all(extractID$accepted_name, fixed(" "), "_")
extractID$accepted_name <- sub("[.]", "_", extractID$accepted_name)

extractID$phylo.tip <- sub("-", "_", extractID$phylo.tip)
extractID$phylo.tip <- sub(" ", "_", extractID$phylo.tip)
extractID$phylo.tip <- gsub(extractID$phylo.tip, pattern = " var.*", replacement = "") 
extractID$phylo.tip <- gsub(extractID$phylo.tip, pattern = " subsp.*", replacement = "") 
extractID$phylo.tip <- gsub(extractID$phylo.tip, pattern = "x ", replacement = "") 
head(extractID)
### The table to link input name (includefile), accepted name (GenBank == Phylogeny Tips) and NCBI_ID

### Create a table to link GenBank accepted name (GenBank == Phylogeny Tips before trimming)
(GB_accessions[!rownames(GB_accessions) %in% extractID$accepted_name,])
(GB_accessions[!GB_accessions$species %in% extractID$accepted_name,])
GB_accessions <- cbind(GenBank_taxonomy = rownames(GB_accessions), GB_accessions)
head(GB_accessions)

taxon.dump <- merge(extractID, GB_accessions, by.x="accepted_name", by.y="species", all.y=T)
head(taxon.dump)
write.csv(taxon.dump, "data/TaxRef/taxon.dump.csv")
#tail(taxon.dump$phylo.tip)
#tail(v6tax_edit$pezAlpes.phy.tip.label)

head(v4tax) # the community matrix based on TaxRef v.9
taxon.dump.merge <- merge(taxon.dump, v4tax, by.x="input_name", by.y="pezalps.TNRS.phy.tip.label", all=T)
head(taxon.dump.merge)
dim(taxon.dump.merge) #1553 = the lenght of the include file + some duplicates
taxon.dump.merge.order <- taxon.dump.merge[order(taxon.dump.merge$input_name),]

#taxon.dump.merge.order[!taxon.dump.merge.order$phylo.tip %in% taxon.dump.merge.order$input_name, 3] # this column is redundant
taxon.dump.merge.order <- taxon.dump.merge.order[-c(which(colnames(taxon.dump.merge.order) %in% c("accepted_name", "name_class","Under.Ice")))] #remove "input_name",
head(taxon.dump.merge.order)
taxon.dump.merge.order2 <- taxon.dump.merge.order[order(taxon.dump.merge.order$Valid_LIB.C.150),]
head(taxon.dump.merge.order2)

taxon.dump.merge.order3 <- taxon.dump.merge.order2[c(11,12, 5:9, 13:ncol(taxon.dump.merge.order2))]

head(taxon.dump.merge.order3)
duplicated(taxon.dump.merge.order3$Valid_LIB.C.150)

str(taxon.dump.merge.order3)
taxon.dump.merge.order3[3:7] <- lapply(taxon.dump.merge.order3[3:7], function(x) as.numeric(as.character(x)))
taxon.dump.merge.order3[is.na(taxon.dump.merge.order3)] <- 0
taxon.dump.merge.order4 <- taxon.dump.merge.order3 %>% group_by(Valid_LIB.C.150, Alien_origin) %>% summarise_each(funs(max))
head(taxon.dump.merge.order4)
dim(taxon.dump.merge.order4) #1346
colnames(taxon.dump.merge.order4)
write.csv(taxon.dump.merge.order4, "output/Dryad1_EcrinsCommunityMatrix_27April2017April2017v3.csv")

colnames(taxon.dump.merge.order2) <- c("input_name" = "PHLAWD Input Name", "phylo.tip" = "TNRS iPlant Taxonomy", "ncbi_id" = "GenBank Taxonomy ID", 
                                       "GenBank_taxonomy" = "GenBank Taxonomy", "ITS", "trnTLF", "rbcL", "matK", "atpB", "rep.1..length.pezalps.TNRS.phy.tip.label.." = "Alps Phylogeny", 
                                       "Valid_LIB.C.150" = "TaxRef v.9", "Alien_origin" = "Alien Origin", 
                                       "Endemic" = "Endemic", "Ecrins.NP" = "Regional", "Summits" = "All Summits", "Persistent" = "LGM", 
                                       "Plat.de.la.Selle" ="Aiguille du plat de la Selle", "Barre.des.Ecrins" = "Barre des Ecrins", "La.Meije"= "la Meije",
                                       "Sirac" = "le Sirac", "Rouies" = "les Rouies", "Olan" = "Olan", "Pelvoux" = "Mont Pelvoux", 
                                       "Occidentale.Ailefroide", "Brevoort" = "Point Brevoort", "Burlan" = "Pointes de Burlan",
                                       "Rocher.de.la.Selle" = "Rocher de la Selle", "Choisy" = "Tour Choisy", "Muraillette" = "Tete de la Muraillette", 
                                       "Rateau" = "le Rateau", "Lauvitel" = "Signal de Lauvitel")
write.csv(taxon.dump.merge.order2, "output/Dryad1_EcrinsCommunityMatrix_27April2017April2017v2.csv")


#### Annotate taxonomy of each alignment:

atpB <- annotate.fasta(fastaFile= "output/3_Clean/atpB.unique.rem.name.fasta.clean", taxonomy = v6tax_edit[c(2,3)])
ITS <- annotate.fasta(fastaFile= "output/3_Clean/ITS.unique.rem.name.fasta.clean", taxonomy = v6tax_edit[c(2,3)])
matK <- annotate.fasta(fastaFile= "output/3_Clean/matK.unique.rem.name.fasta.clean", taxonomy = v6tax_edit[c(2,3)])
rbcL <- annotate.fasta(fastaFile= "output/3_Clean/rbcL.unique.rem.name.fasta.clean", taxonomy = v6tax_edit[c(2,3)])
trnTLF <- annotate.fasta(fastaFile= "output/3_Clean/trnTLF.unique.rem.name.fasta.clean", taxonomy = v6tax_edit[c(2,3)])
concat <- annotate.fasta(fastaFile= "output/4_Concatenate/align.concat.EcrinSpPool.121115.fasta", taxonomy = v6tax_edit[c(2,3)])


#### Annotate taxonomy of each gene tree:

phy.tmp.atpB <- read.tree("output/5_Trees/GeneTrees/atpB/RAxML_bipartitions.atpBrem.100.210915.sperma") 
for (i in 1:length(phy.tmp.atpB$tip.label)){
  if (phy.tmp.atpB$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.atpB$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.atpB$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.atpB$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}
"Arctostaphylos_uva-ursi" %in% phy.tmp.atpB$tip.label
write.tree(phy.tmp.atpB, file="output/5_Trees/GeneTrees/atpB/atpB.taxized.taxref.tre")

phy.tmp.ITS <- read.tree("output/5_Trees/GeneTrees/ITS/RAxML_bipartitions.ITSrem.100.210915.sperma") 
for (i in 1:length(phy.tmp.ITS$tip.label)){
  if (phy.tmp.ITS$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.ITS$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.ITS$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.ITS$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}
"Arctostaphylos_uva-ursi" %in% phy.tmp.ITS$tip.label
phy.tmp.ITS.dropdup <- ape::drop.tip(phy = phy.tmp.ITS, tip = c(phy.tmp.ITS$tip.label[duplicated(phy.tmp.ITS$tip.label)]))
write.tree(phy.tmp.ITS.dropdup, file="output/5_Trees/GeneTrees/ITS/ITS.taxized.taxref.tre")

phy.tmp.matK <- read.tree("output/5_Trees/GeneTrees/matK/RAxML_bipartitions.matK.100.180915.sperma") 
for (i in 1:length(phy.tmp.matK$tip.label)){
  if (phy.tmp.matK$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.matK$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.matK$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.matK$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}
"Arctostaphylos_uva-ursi" %in% phy.tmp.matK$tip.label
phy.tmp.matK.dropdup <- ape::drop.tip(phy = phy.tmp.matK, tip = c(phy.tmp.matK$tip.label[duplicated(phy.tmp.matK$tip.label)]))
write.tree(phy.tmp.matK.dropdup, file="output/5_Trees/GeneTrees/matK/matK.taxized.taxref.tre")

phy.tmp.rbcL <- read.tree("output/5_Trees/GeneTrees/rbcL/RAxML_bipartitions.rbcL.100.180915.sperma") 
for (i in 1:length(phy.tmp.rbcL$tip.label)){
  if (phy.tmp.rbcL$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.rbcL$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.rbcL$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.rbcL$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}
"Arctostaphylos_uva-ursi" %in% phy.tmp.rbcL$tip.label
phy.tmp.rbcL.dropdup <- ape::drop.tip(phy = phy.tmp.rbcL, tip = c(phy.tmp.rbcL$tip.label[duplicated(phy.tmp.rbcL$tip.label)]))
write.tree(phy.tmp.rbcL.dropdup, file="output/5_Trees/GeneTrees/rbcL/rbcL.taxized.taxref.tre")

phy.tmp.trnTLF <- read.tree("output/5_Trees/GeneTrees/trnTLF/RAxML_bipartitions.trnTLFrem.100.210915.sperma") 
for (i in 1:length(phy.tmp.trnTLF$tip.label)){
  if (phy.tmp.trnTLF$tip.label[i] %in% v6tax_edit$pezalps.TNRS.phy.tip.label){
    x <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.trnTLF$tip.label[i]), "pezalps.TNRS.phy.tip.label"] #the name that was used to serach (i.e. corrected to iPlant taxonomy in community matirx --> includefile)
    y <- v6tax_edit[which(v6tax_edit$pezalps.TNRS.phy.tip.label==phy.tmp.trnTLF$tip.label[i]), "Valid_LIB.C.150"] #name in TaxRef v9
    #print(paste(i, "input =", x, ", accepted =", y))
    phy.tmp.trnTLF$tip.label[i] <- as.character(y)#replace tip with taxref name
  } 
}
"Arctostaphylos_uva-ursi" %in% phy.tmp.trnTLF$tip.label
phy.tmp.trnTLF.dropdup <- ape::drop.tip(phy = phy.tmp.trnTLF, tip = c(phy.tmp.trnTLF$tip.label[duplicated(phy.tmp.trnTLF$tip.label)]))
write.tree(phy.tmp.trnTLF.dropdup, file="output/5_Trees/GeneTrees/trnTLF/trnTLF.taxized.taxref.tre")





