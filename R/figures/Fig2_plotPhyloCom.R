#####################################################################################################################
############# Functions to plot and visualize megaphylogenies ####################################################### 
############# Hannah E. Marx, 25 April 2017 #########################################################################
#####################################################################################################################

source("R/functions/treeFunctions.R")
source("R/functions/treeVisFunctions.R")

########################### Plot Commuity Phylogeny of Ecrins Alpine "Sky Islands" ###########################################
###### Color taxa on summit: Code modified from picante color.plot.phylo
###### Plot phylo with calibrated nodes IDed and nodes lables with taxonomy (dated with PATHd8)

######## Nodes Object ################
#### Need to re-run Congruifier to get cal object (calibrated nodes) to plot congruified nodes 

genetree=pezAlpes$phy
tax=read.csv(file="data/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) 
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
atol=read.tree("data/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011
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


##### Vector of bootstrap support values
p2 <- character(length(pezAlpes$phy$node.label)) # create a vector of colors that correspond to different support ranges 
p2[] <- "#0000ff00" ## Transparent color: "#RRGGBBAA" and the AA portion is the opacity/trasparency.
p2[pezAlpes$phy$node.label >=  95] <- "black"
p2[pezAlpes$phy$node.label ==  100] <- "black"
p2[pezAlpes$phy$node.label < 95 & pezAlpes$phy$node.label >= 75] <- "slategray4"
#p2[treePLboots$node.label < 95 & treePLboots$node.label >= 75] <- "slategray"
#p2[treePLboots$node.label < 75] <- "red1"
p2[pezAlpes$phy$node.label == ""] <- "black" # node = 353 = 100
p2 


######## Plot phylo with the presence/absence of each taxa on a summit color coded across tips
rownames(pezAlpes$comm)
## Add a row of species names
com.data.plot <- as.data.frame(cbind(species = names(as.data.frame(pezAlpes$comm)), 
                                     t(pezAlpes$comm[c(1:4, 6:nrow(pezAlpes$comm)),])))
# convert to numerics
com.data.plot[2:ncol(com.data.plot)] <- sapply(com.data.plot[2:ncol(com.data.plot)] , 
                                               function(x) as.numeric(as.character(x)))

# Add column to indicate endemic species 
colnames(com.data.plot)
com.data.plot.tmp <- merge(com.data.plot, endemics[3], by=0, all.x=T)
colnames(com.data.plot.tmp)[19] <- "Summits"
colnames(com.data.plot.tmp)[20] <- "Endemic"
com.data.plot.tmp$Endemic[!is.na(com.data.plot.tmp$Endemic)] <- 1
com.data.plot.tmp$Endemic[is.na(com.data.plot.tmp$Endemic)] <- 0
length(which(com.data.plot.tmp$Endemic == 1)) #58
head(com.data.plot.tmp)
com.data.plot <- com.data.plot.tmp[-1]
rownames(com.data.plot) <- com.data.plot$species

## Re-order column names (increasing species richness)
com.data.plot <- com.data.plot[,c("species", "Brevoort",
                                  "Pelvoux",
                                  "Choisy",
                                  "Occidentale Ailefroide",
                                  "Plat de la Selle",
                                  "Burlan",
                                  "La Meije",
                                  "Barre des Ecrins",
                                  "Rouies",
                                  "Muraillette",
                                  "Olan", 
                                  "Sirac",
                                  "Lauvitel",    
                                  "Rateau", 
                                  "Rocher de la Selle",
                                  "Persistent",   
                                  "Summits",
                                  "Endemic")]

# plot
pdf("figs/Figure2_Phylo.pdf", 10,10)
trait.plot.colorCom(tree = pezAlpes$phy,
                    dat = com.data.plot,
                    shape = "Endemic",
                    pch = 8,
                    cols = list("Brevoort" = c("white", "red4"),
                                "Pelvoux" = c("white", "red2"),
                                "Choisy" = c("white", "tomato2"),
                                "Occidentale Ailefroide" = c("white", "darkorange1"),
                                "Plat de la Selle" = c("white", "goldenrod"),
                                "Burlan"= c("white", "yellow"),
                                "La Meije"= c("white", "chartreuse"),
                                "Barre des Ecrins"= c("white", "limegreen"),
                                "Rouies"= c("white", "darkgreen"),
                                "Muraillette"= c("white", "darkslategrey"),
                                "Olan"= c("white", "darkturquoise"), 
                                "Sirac"= c("white", "dodgerblue"),
                                "Lauvitel"= c("white", "blue"),    
                                "Rateau"= c("white", "purple2"), 
                                "Rocher de la Selle"= c("white", "hotpink"),
                                "Persistent" = c("white", "deeppink"),   
                                "Summits" = c("white", "deeppink4"),
                                "Endemic" = c("white", "black")),
                    trait = "Summits",
                    datTr = com.data.plot,
                    taxa.names = "species",
                    col.names = c("grey", "black"),
                    cex.lab = 0.15,
                    font.lab = 0.5,
                    cex.shape=.5)
dev.off()
