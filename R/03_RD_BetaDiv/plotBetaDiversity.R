#####################################################################################################################
############# Plotting phylogenetic turnover between alpine summits (beta) ##########################################
############# Decomposed beta diveristy (PhyloSor) to get 'true' turnover (independent of Species Richness) #########
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################

source("R/functions/chooseClade.R")
source("R/functions/plotBetaPair.R")

########### Spermatophyta ########### 

betadiv.ecrins <- read.csv(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/Dryad_decompo_beta_betadiv.csv")
head(betadiv.ecrins)

### Convert all beta metrics into matrix 
betadiv.ecrins.tmp <- cbind(com1 = sapply(strsplit(as.character(betadiv.ecrins$X), split = "-", fixed = T), "[", 1L), 
      com2 = sapply(strsplit(as.character(betadiv.ecrins$X), split = "-", fixed = T), "[", 2L), 
      betadiv.ecrins)
head(betadiv.ecrins.tmp)

### Plot all beta diversity output from both metrics in heatmap
betadiv.ecrins.tmp.melt <- melt(betadiv.ecrins.tmp)
head(betadiv.ecrins.tmp.melt)
pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/all.beta.metics.pdf")
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

betadiv.ecrins.dist.PhyloSor.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_PhyloSor.csv")
betadiv.ecrins.dist.UniFrac.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_UniFrac.csv")

betadiv.ecrins.dist.PhyloSorTurn.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_PhyloSor_turn.csv")
betadiv.ecrins.dist.UniFracTurn.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_UniFrac_turn.csv")

betadiv.ecrins.dist.PhyloSorPD.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_PhyloSor_PD.csv")
betadiv.ecrins.dist.UniFracPD.tmp <- read.csv("output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/distance_SES/Spermatophyta_SES_UniFrac_PD.csv")


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
pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_PhyloSor.pdf", onefile=FALSE)
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

pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_UniFrac.pdf", onefile=FALSE)
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
pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_PhyloSor_turn.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits, main="SES_PhyloSor_turn", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96))/2 #8
length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96))/2 #7


pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_UniFrac_turn.pdf", onefile=FALSE)
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
pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_PhyloSor_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.PhyloSorPD.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.PhyloSorPD.matrix.summits, main="SES_PhyloSor_PD", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96))/2 #8
length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96))/2 #16


pdf(file="output/8_PhyoDiversity/beta/decomposedBeta/seedPlants/figs/SES_UniFrac_PD.pdf", onefile=FALSE)
mat <- matrix(ifelse(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96, "*", ""), 
              nrow(betadiv.ecrins.dist.UniFracPD.matrix.summits))
mat[is.na(mat)] <-  ""
pheatmap(betadiv.ecrins.dist.UniFracPD.matrix.summits, main="SES_UniFrac_PD", cluster_rows=T, 
         fontsize_row=fontsize_row, border_color=NA,  
         display_numbers = mat, col=colors2, breaks=bk2)
dev.off()
length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96))/2 #10
length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96))/2 #14


########### Each clade ########### 

beta.Sperma <- plotBetaPair(clade = "Spermatophyta", 
                            outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/18Jan2016/")

beta.Aster <- plotBetaPair(clade = "Asterales", 
                           outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Asterales/")

beta.Caryo <-plotBetaPair(clade = "Caryophyllales", 
                          outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Caryophyllales/")

beta.Lam <-plotBetaPair(clade = "Lamiales", 
                        outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Lamiales/")

beta.Poa <-plotBetaPair(clade = "Poales", 
                        outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Poales/")

beta.Ros <-plotBetaPair(clade = "Rosales", 
                        outputFile = "output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/Rosales/")

## Summarize the number of signigicant comparisons
write.csv(rbind(beta.Sperma, beta.Aster, beta.Caryo, beta.Lam, beta.Poa, beta.Ros), file="output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades/decomposedBetaClades.csv")






