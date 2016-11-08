
source("R/treeFunctions.R")
source("R/trait.plot.colorTips.R")
source("R/trait.plot.colorCom.R")


########################### Plot Tree ###########################################
###### Color taxa on summit: Code modified from picante color.plot.phylo
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

dat.plot <- pezAlpes$comm["Summits",]
dat.plot[dat.plot > 0] <- 1
dat.plot <- as.data.frame(dat.plot)
head(dat.plot)

#pdf(file="output/6_Visualize Tree/Spermatophyta/FloraEcrin.pdf") 
color.plot.phylo3NoTrait(pezAlpes$phy, as.data.frame(cbind("tiplabels" = rownames(dat.plot), "summitsPA" = dat.plot[,1])), 
                         taxa.names = "tiplabels", trait = "summitsPA", col.names = colorsIsl, label.offset=5.2, phy.cex=.05, 
                         cut.labs =c("Ecrin NP", "Alpine Summits"), leg.cex=.8, main="Flora of the Le Parc national des Ecrins")
#dev.off()


#pdf(file="output/6_Visualize Tree/Spermatophyta/FloraEcrin.dots.pdf")
plot(pezAlpes$phy, show.tip.label = FALSE, type="fan")
tiplabels(tip = which(pezAlpes$comm["Summits",] > 0), pch = 8, cex = .5, col = "darkorange2")
#dev.off()


#pdf(file="output/6_Visualize Tree/Spermatophyta/FloraEcrin.dots.v2.pdf")
trait.plot(pezAlpes$phy, dat.plot, cols = list(dat.plot = c("white", "darkorange2")), font = 0, cex.lab = .2) #.0001
nodelabels(text=NULL, cex=ifelse(vec==0, NA, 1), frame="n", col="black", pch=21, bg="white") # Plot congruified nodes
nodelabels(pch = 21, cex = 0.3, col=p2, bg = p2) # Plot bootstrap support
#dev.off()


######## Plot phylo with p/a on summits plotted across tips
rownames(pezAlpes$comm)
## Add a row of species names
com.data.plot <- as.data.frame(cbind(species = names(as.data.frame(pezAlpes$comm)), 
                                     t(pezAlpes$comm[c(1:4, 6:nrow(pezAlpes$comm)),])))
# convert to numerics
com.data.plot[2:ncol(com.data.plot)] <- sapply(com.data.plot[2:ncol(com.data.plot)] , 
                                               function(x) as.numeric(as.character(x)))
# change coounts to 1/0
com.data.plot[, 2:ncol(com.data.plot)][com.data.plot[, 2:ncol(com.data.plot)] > 0] <- 1
str(com.data.plot)
head(com.data.plot)
names(com.data.plot)

# remove "Under Ice"
com.data.plot <- com.data.plot[,-c(which(colnames(com.data.plot)%in%"Under Ice"))]
colnames(com.data.plot)

# Add endemics
head(endemics)
endemics.df <- endemics[1]
endemics.df


ncol(com.data.plot.tmp)
com.data.plot.tmp <- merge(com.data.plot, endemics.df, by=0, all.x=T)
colnames(com.data.plot.tmp)[20] <- "Endemic"
com.data.plot.tmp$Endemic[!is.na(com.data.plot.tmp$Endemic)] <- 1
com.data.plot.tmp$Endemic[is.na(com.data.plot.tmp$Endemic)] <- 0
length(which(com.data.plot.tmp$Endemic == 1)) #59
head(com.data.plot.tmp)
com.data.plot <- com.data.plot.tmp[-1]
rownames(com.data.plot) <- com.data.plot$species

# choose color palate
n=19
pie(rep(1,n), col=FALSE)
pie(rep(1,n), col=rainbow(n))

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
pdf("output/6_Visualize Tree/Spermatophyta/EcrinsCommunitiesPhyloEndemic.pdf", 10,10)
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

# plot
pdf("figs/EcrinsCommunitiesPhylo.pdf", 10,10)
trait.plot.colorCom(tree = pezAlpes$phy,
                    dat = com.data.plot,
                    shape = "Endemic",
                    pch = "*",
                    cex.shape = 1, 
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
                    font.lab = 0.5)
dev.off()