
########################### Beta diveristy (turnover) of Ecrins Alpine Summits ########################### 
#### 07 March, 2016 H. Marx
##########################################################################################################

################ spacodiR: PIst for each node in summit community phylogeny ################ 
### check branching times 
phy.nodetimes(pezAlpes$phy, time.range = c(0, max(alps.phy$edge.length)), proportion = TRUE)

########### community diversity statistics of Hardy and Senterre (2007): tree-based
### Ecrins Communtiy Phylogeny
ecrin.beta <- spacodi.calc(sp.plot = t(pezAlpes$comm), phy = pezAlpes$phy, pairwise = T)
ecrin.beta$PIst #0.002860691
ecrin.beta$pairwise.PIst
head(ecrin.beta$sp.plot)
(ecrin.beta$sp.tree)

### Prune community and phylogeny to just summits (for SES randomization from summit species pool)
spac.com=as.spacodi(pezAlpes$comm[-c(5, 12, 18, 19),])
head(spac.com)
colSums(spac.com)
spac.com.summits <- spac.com[!rowSums(spac.com) == 0,]
spac.phy <- treedata(pezAlpes$phy, spac.com.summits)[[1]]
com.beta <- spacodi.calc(sp.plot = spac.com.summits, phy = spac.phy, pairwise = T)
com.beta$PIst #0.004090516
com.beta$pairwise.PIst 

###spatial clustering: species within plots are more phylogenetically related on average
#than species from distinct plots where Pst > Ist, Bst > 0, or PIst > 0. Species are
#functionally more similar locally than those from distinct plots where Tst > Ist, Ust > 0,
#or TAUst > 0
###spatial overdispersion: species within plots are less phylogenetically related on average
#than species from distinct plots where Pst < Ist, Bst < 0, or PIst < 0. Species
#are functionally less similar locally than are species from distinct plots where Tst < Ist,
#Ust < 0, or TAUst < 0


############ Plotting observed and expected community structure across branching times of a phylogeny
# generate a plot of observed and expected PIst on summit community phylogeny 
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
sp.permut=spacodi.by.nodes(sp.plot=com.beta$sp.plot, phy=com.beta$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut$observed.PIst)
head(sp.permut$randomization.test)

pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/PIst.branching.time.pdf")
spacodi.permutplot(sp.permut, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut$observed.PIst, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/observed.PIst.csv")
write.csv(sp.permut$expected.PIst, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/expected.PIst.csv")
write.csv(sp.permut$randomization.test, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/randomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()


# PIst: 10 permutations to perform on the dataset (default)
PI=spacodi.by.nodes(sp.plot=com.beta$sp.plot, sp.parm="PIst", phy=com.beta$sp.tree, return.all=TRUE, method="1s")
pdf(file = "output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/10reps/PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(PI, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()
pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/10reps/PIst.branching.time.pdf")
spacodi.permutplot(PI, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(PI$observed.PIst, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/10reps/observed.PIst.csv")
write.csv(PI$expected.PIst, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/10reps/expected.PIst.csv")
write.csv(PI$randomization.test, file="output/9_PhyoDiversity/Spermatophyta/beta/PIst_summitPool/10reps/randomization.test.csv")




####################### Decompose beta diveristy (PhyloSor) to get 'true' turnover (independent od Species Richness) ####################### 

#source("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/ecrins_beta_diversity.R") ## run on terminal 

betadiv.ecrins <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/decompo_beta_betadiv3.csv")
head(betadiv.ecrins)

### Convert all beta metrics into matrix 
betadiv.ecrins.tmp <- cbind(com1 = sapply(strsplit(as.character(betadiv.ecrins$X), split = "-", fixed = T), "[", 1L), 
      com2 = sapply(strsplit(as.character(betadiv.ecrins$X), split = "-", fixed = T), "[", 2L), 
      betadiv.ecrins)
head(betadiv.ecrins.tmp)

### Plot in heatmap
betadiv.ecrins.tmp.melt <- melt(betadiv.ecrins.tmp)
head(betadiv.ecrins.tmp.melt)
pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/18Jan2016/all.beta.metics.pdf")
ggplot(betadiv.ecrins.tmp.melt, aes(variable, X)) +
  geom_tile(aes(fill=value), colour ="white") + #, alpha = .75
  scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white") +
  theme_grey(base_size = 6) + labs(x = "",  y = "") + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
dev.off()

### NOTE: values for SES are different betweeen decompo_beta_betadiv3.csv and the distances from output.dist = T
## The rest of the anlyses are done using the output distance matrices

betadiv.ecrins.dist.PhyloSor.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_PhyloSor.csv")
betadiv.ecrins.dist.UniFrac.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_UniFrac.csv")

betadiv.ecrins.dist.PhyloSorTurn.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_PhyloSor_turn.csv")
betadiv.ecrins.dist.UniFracTurn.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_UniFrac_turn.csv")

betadiv.ecrins.dist.PhyloSorPD.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_PhyloSor_PD.csv")
betadiv.ecrins.dist.UniFracPD.tmp <- read.csv("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/distance_SES/SES_UniFrac_PD.csv")


### Heatmap of summits PhyloSor 
col.order <- paste(alps.env.sprich$Row.names)
betadiv.ecrins.dist.PhyloSor <- betadiv.ecrins.dist.PhyloSor.tmp[,-1]
colnames(betadiv.ecrins.dist.PhyloSor) <- betadiv.ecrins.dist.PhyloSor.tmp$X
rownames(betadiv.ecrins.dist.PhyloSor) <- betadiv.ecrins.dist.PhyloSor.tmp$X
betadiv.ecrins.dist.PhyloSor <- betadiv.ecrins.dist.PhyloSor[match(col.order, rownames(betadiv.ecrins.dist.PhyloSor)),]
betadiv.ecrins.dist.PhyloSor <- betadiv.ecrins.dist.PhyloSor[,match(col.order, colnames(betadiv.ecrins.dist.PhyloSor))]
betadiv.ecrins.dist.PhyloSor
betadiv.ecrins.dist.PhyloSor.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSor)
betadiv.ecrins.dist.PhyloSor.matrix.summits <- betadiv.ecrins.dist.PhyloSor.matrix[!rownames(betadiv.ecrins.dist.PhyloSor.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                   !colnames(betadiv.ecrins.dist.PhyloSor.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.PhyloSor.matrix.summits,  na.rm = T) 
  
### Heatmap of summits UniFrac 
betadiv.ecrins.dist.UniFrac <- betadiv.ecrins.dist.UniFrac.tmp[,-1]
colnames(betadiv.ecrins.dist.UniFrac) <- betadiv.ecrins.dist.UniFrac.tmp$X
rownames(betadiv.ecrins.dist.UniFrac) <- betadiv.ecrins.dist.UniFrac.tmp$X
betadiv.ecrins.dist.UniFrac <- betadiv.ecrins.dist.UniFrac[match(col.order, rownames(betadiv.ecrins.dist.UniFrac)),]
betadiv.ecrins.dist.UniFrac <- betadiv.ecrins.dist.UniFrac[,match(col.order, colnames(betadiv.ecrins.dist.UniFrac))]
betadiv.ecrins.dist.UniFrac
betadiv.ecrins.dist.UniFrac.matrix <- data.matrix(betadiv.ecrins.dist.UniFrac)
betadiv.ecrins.dist.UniFrac.matrix.summits <- betadiv.ecrins.dist.UniFrac.matrix[!rownames(betadiv.ecrins.dist.UniFrac.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                  !colnames(betadiv.ecrins.dist.UniFrac.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.UniFrac.matrix.summits,  na.rm = T) 

### Heatmap of summits PhyloSor Turn
betadiv.ecrins.dist.PhyloSorTurn <- betadiv.ecrins.dist.PhyloSorTurn.tmp[,-1]
colnames(betadiv.ecrins.dist.PhyloSorTurn) <- betadiv.ecrins.dist.PhyloSorTurn.tmp$X
rownames(betadiv.ecrins.dist.PhyloSorTurn) <- betadiv.ecrins.dist.PhyloSorTurn.tmp$X
betadiv.ecrins.dist.PhyloSorTurn <- betadiv.ecrins.dist.PhyloSorTurn[match(col.order, rownames(betadiv.ecrins.dist.PhyloSorTurn)),]
betadiv.ecrins.dist.PhyloSorTurn <- betadiv.ecrins.dist.PhyloSorTurn[,match(col.order, colnames(betadiv.ecrins.dist.PhyloSorTurn))]
betadiv.ecrins.dist.PhyloSorTurn
betadiv.ecrins.dist.PhyloSorTurn.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSorTurn)
betadiv.ecrins.dist.PhyloSorTurn.matrix.summits <- betadiv.ecrins.dist.PhyloSorTurn.matrix[!rownames(betadiv.ecrins.dist.PhyloSorTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                      !colnames(betadiv.ecrins.dist.PhyloSorTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits,  na.rm = T) 

### Heatmap of summits UniFrac Turn
betadiv.ecrins.dist.UniFracTurn <- betadiv.ecrins.dist.UniFracTurn.tmp[,-1]
colnames(betadiv.ecrins.dist.UniFracTurn) <- betadiv.ecrins.dist.UniFracTurn.tmp$X
rownames(betadiv.ecrins.dist.UniFracTurn) <- betadiv.ecrins.dist.UniFracTurn.tmp$X
betadiv.ecrins.dist.UniFracTurn <- betadiv.ecrins.dist.UniFracTurn[match(col.order, rownames(betadiv.ecrins.dist.UniFracTurn)),]
betadiv.ecrins.dist.UniFracTurn <- betadiv.ecrins.dist.UniFracTurn[,match(col.order, colnames(betadiv.ecrins.dist.UniFracTurn))]
betadiv.ecrins.dist.UniFracTurn
betadiv.ecrins.dist.UniFracTurn.matrix <- data.matrix(betadiv.ecrins.dist.UniFracTurn)
betadiv.ecrins.dist.UniFracTurn.matrix.summits <- betadiv.ecrins.dist.UniFracTurn.matrix[!rownames(betadiv.ecrins.dist.UniFracTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                          !colnames(betadiv.ecrins.dist.UniFracTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.UniFracTurn.matrix.summits,  na.rm = T) 

### Heatmap of summits PhyloSor PD
betadiv.ecrins.dist.PhyloSorPD <- betadiv.ecrins.dist.PhyloSorPD.tmp[,-1]
colnames(betadiv.ecrins.dist.PhyloSorPD) <- betadiv.ecrins.dist.PhyloSorPD.tmp$X
rownames(betadiv.ecrins.dist.PhyloSorPD) <- betadiv.ecrins.dist.PhyloSorPD.tmp$X
betadiv.ecrins.dist.PhyloSorPD <- betadiv.ecrins.dist.PhyloSorPD[match(col.order, rownames(betadiv.ecrins.dist.PhyloSorPD)),]
betadiv.ecrins.dist.PhyloSorPD <- betadiv.ecrins.dist.PhyloSorPD[,match(col.order, colnames(betadiv.ecrins.dist.PhyloSorPD))]
betadiv.ecrins.dist.PhyloSorPD
betadiv.ecrins.dist.PhyloSorPD.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSorPD)
betadiv.ecrins.dist.PhyloSorPD.matrix.summits <- betadiv.ecrins.dist.PhyloSorPD.matrix[!rownames(betadiv.ecrins.dist.PhyloSorPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                         !colnames(betadiv.ecrins.dist.PhyloSorPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.PhyloSorPD.matrix.summits,  na.rm = T) 

### Heatmap of summits UniFrac PD
betadiv.ecrins.dist.UniFracPD <- betadiv.ecrins.dist.UniFracPD.tmp[,-1]
colnames(betadiv.ecrins.dist.UniFracPD) <- betadiv.ecrins.dist.UniFracPD.tmp$X
rownames(betadiv.ecrins.dist.UniFracPD) <- betadiv.ecrins.dist.UniFracPD.tmp$X
betadiv.ecrins.dist.UniFracPD <- betadiv.ecrins.dist.UniFracPD[match(col.order, rownames(betadiv.ecrins.dist.UniFracPD)),]
betadiv.ecrins.dist.UniFracPD <- betadiv.ecrins.dist.UniFracPD[,match(col.order, colnames(betadiv.ecrins.dist.UniFracPD))]
betadiv.ecrins.dist.UniFracPD
betadiv.ecrins.dist.UniFracPD.matrix <- data.matrix(betadiv.ecrins.dist.UniFracPD)
betadiv.ecrins.dist.UniFracPD.matrix.summits <- betadiv.ecrins.dist.UniFracPD.matrix[!rownames(betadiv.ecrins.dist.UniFracPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                        !colnames(betadiv.ecrins.dist.UniFracPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
heatmap.2(betadiv.ecrins.dist.UniFracPD.matrix.summits,  na.rm = T) 

## aestitic values
#scale = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(98)
fontsize_row = 10 - nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits) / 15
#bk2 = unique(seq(-3,3, length=100))

#create the breaks
bk2 = unique(c(seq(-6, -.0001, length=100), 0, seq(.0001, 6, length=100)))

#set different color vectors for each interval
col1 = colorRampPalette(rev(brewer.pal(n = 9, name ="Blues")))(99) #set the order of greys
col2 <- rep("white", 1)
col3 = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(99)
colors2 <- c(col1, col2, col3)

######## Undecomposed
pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_PhyloSor.pdf", onefile=FALSE)
### Matrix with significant SES values
mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSor.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSor.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.PhyloSor.matrix.summits, col=colors2, main="SES_PhyloSor", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, breaks=bk2, scale="none")
dev.off()

### PhyloSor Significant
length(which(betadiv.ecrins.dist.PhyloSor.matrix.summits > 1.96))/2 #0
length(which(betadiv.ecrins.dist.PhyloSor.matrix.summits < -1.96))/2 #24
### TOTAL pairwise comparisons
(nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits)*(nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits) - 1))/ 2
## 24 / 105

pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_UniFrac.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.UniFrac.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFrac.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.UniFrac.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.UniFrac.matrix.summits, col=colors2, main="SES_UniFrac", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.UniFrac.matrix.summits > 1.96))/2 #0
length(which(betadiv.ecrins.dist.UniFrac.matrix.summits < -1.96))/2 #24

########## TURN
pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_PhyloSor_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits, main="SES_PhyloSor_turn", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96))/2 #8
length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96))/2 #7


pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_UniFrac_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.UniFracTurn.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFracTurn.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.UniFracTurn.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.UniFracTurn.matrix.summits, main="SES_UniFrac_turn", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.UniFracTurn.matrix.summits > 1.96))/2 #7
length(which(betadiv.ecrins.dist.UniFracTurn.matrix.summits < -1.96))/2 #8

######## PD
pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_PhyloSor_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.PhyloSorPD.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.PhyloSorPD.matrix.summits, main="SES_PhyloSor_PD", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96))/2 #8
length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96))/2 #16


pdf(file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.UniFracPD.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.UniFracPD.matrix.summits, main="SES_UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96))/2 #10
length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96))/2 #14




#### Alternate methods for heatmap visulaization
heatmap(betadiv.ecrins.dist.UniFrac.matrix.summits, Colv = F, scale='none', col=scale)
heatmap(betadiv.ecrins.dist.PhyloSor.matrix.summits, Colv = F, scale='none', col=scale)
heatmap(betadiv.ecrins.dist.PhyloSor.matrix.summits, Colv = F, scale='none', col=scale)


annotation_cols = list(Summit = c("Brevoort"="red2", "Pelvoux" = "tomato", "Choisy" = "orangered1", 
                                  "Occidentale Ailefroide" = "tan1", "Plat de la Selle" = "orange", "Burlan"="yellow", "La Meije"= "limegreen", 
                                  "Barre des Ecrins"= "green3", "Rouies"= "green4", "Muraillette"= "darkslategrey", "Olan"= "darkturquoise",
                                  "Sirac"= "aquamarine2", "Lauvitel"= "cyan", "Rateau"= "dodgerblue", "Rocher de la Selle"= "blue"))


length(annotation_cols$Summit)
annotation_row <- data.frame(Summit = factor(names(annotation_cols$Summit)))
rownames(annotation_row) <- (rownames(betadiv.ecrins.dist.PhyloSor.matrix.summits))
annotation_col <- data.frame(Summit = factor(names(annotation_cols$Summit)))
rownames(annotation_col) <- (rownames(betadiv.ecrins.dist.PhyloSor.matrix.summits))

pheatmap(betadiv.ecrins.dist.PhyloSor.matrix.summits, col=scale, main="SES_PhyloSor_turn", cluster_cols=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, annotation_col = annotation_col,
         annotation_row = annotation_row, annotation_colors = annotation_cols[1])





