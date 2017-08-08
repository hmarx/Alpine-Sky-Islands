#####################################################################################################################
############# Function to read output from decomposed beta-diversity and plot pairwise heatmaps #####################
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################


plotBetaPair <- function(clade, outputFile){
  betadiv.ecrins.dist.PhyloSor.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_PhyloSor.csv", sep = "_"), sep = ""))
  betadiv.ecrins.dist.UniFrac.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_UniFrac.csv", sep = "_"), sep = ""))
  
  betadiv.ecrins.dist.PhyloSorTurn.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_PhyloSor_turn.csv", sep = "_"), sep = ""))
  betadiv.ecrins.dist.UniFracTurn.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_UniFrac_turn.csv", sep = "_"), sep = ""))
  
  betadiv.ecrins.dist.PhyloSorPD.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_PhyloSor_PD.csv", sep = "_"), sep = ""))
  betadiv.ecrins.dist.UniFracPD.tmp <- read.csv(paste(paste("output/9_PhyoDiversity/Spermatophyta/beta/decomposedBeta/clades", clade, "distance_SES/", sep = "/"), paste(clade, "SES_UniFrac_PD.csv", sep = "_"), sep = ""))
  
  
  ### Heatmap of summits PhyloSor 
  betadiv.ecrins.dist.PhyloSor <- betadiv.ecrins.dist.PhyloSor.tmp[,-1]
  colnames(betadiv.ecrins.dist.PhyloSor) <- betadiv.ecrins.dist.PhyloSor.tmp$X
  rownames(betadiv.ecrins.dist.PhyloSor) <- betadiv.ecrins.dist.PhyloSor.tmp$X
  
  betadiv.ecrins.dist.PhyloSor.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSor)
  betadiv.ecrins.dist.PhyloSor.matrix.summits <- betadiv.ecrins.dist.PhyloSor.matrix[!rownames(betadiv.ecrins.dist.PhyloSor.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                     !colnames(betadiv.ecrins.dist.PhyloSor.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  
  ### Heatmap of summits UniFrac 
  betadiv.ecrins.dist.UniFrac <- betadiv.ecrins.dist.UniFrac.tmp[,-1]
  colnames(betadiv.ecrins.dist.UniFrac) <- betadiv.ecrins.dist.UniFrac.tmp$X
  rownames(betadiv.ecrins.dist.UniFrac) <- betadiv.ecrins.dist.UniFrac.tmp$X
  betadiv.ecrins.dist.UniFrac.matrix <- data.matrix(betadiv.ecrins.dist.UniFrac)
  betadiv.ecrins.dist.UniFrac.matrix.summits <- betadiv.ecrins.dist.UniFrac.matrix[!rownames(betadiv.ecrins.dist.UniFrac.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                   !colnames(betadiv.ecrins.dist.UniFrac.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  ### Heatmap of summits PhyloSor Turn
  betadiv.ecrins.dist.PhyloSorTurn <- betadiv.ecrins.dist.PhyloSorTurn.tmp[,-1]
  colnames(betadiv.ecrins.dist.PhyloSorTurn) <- betadiv.ecrins.dist.PhyloSorTurn.tmp$X
  rownames(betadiv.ecrins.dist.PhyloSorTurn) <- betadiv.ecrins.dist.PhyloSorTurn.tmp$X
  betadiv.ecrins.dist.PhyloSorTurn.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSorTurn)
  betadiv.ecrins.dist.PhyloSorTurn.matrix.summits <- betadiv.ecrins.dist.PhyloSorTurn.matrix[!rownames(betadiv.ecrins.dist.PhyloSorTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                             !colnames(betadiv.ecrins.dist.PhyloSorTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  
  ### Heatmap of summits UniFrac Turn
  betadiv.ecrins.dist.UniFracTurn <- betadiv.ecrins.dist.UniFracTurn.tmp[,-1]
  colnames(betadiv.ecrins.dist.UniFracTurn) <- betadiv.ecrins.dist.UniFracTurn.tmp$X
  rownames(betadiv.ecrins.dist.UniFracTurn) <- betadiv.ecrins.dist.UniFracTurn.tmp$X
  betadiv.ecrins.dist.UniFracTurn.matrix <- data.matrix(betadiv.ecrins.dist.UniFracTurn)
  betadiv.ecrins.dist.UniFracTurn.matrix.summits <- betadiv.ecrins.dist.UniFracTurn.matrix[!rownames(betadiv.ecrins.dist.UniFracTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                           !colnames(betadiv.ecrins.dist.UniFracTurn.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  
  ### Heatmap of summits PhyloSor PD
  betadiv.ecrins.dist.PhyloSorPD <- betadiv.ecrins.dist.PhyloSorPD.tmp[,-1]
  colnames(betadiv.ecrins.dist.PhyloSorPD) <- betadiv.ecrins.dist.PhyloSorPD.tmp$X
  rownames(betadiv.ecrins.dist.PhyloSorPD) <- betadiv.ecrins.dist.PhyloSorPD.tmp$X
  betadiv.ecrins.dist.PhyloSorPD.matrix <- data.matrix(betadiv.ecrins.dist.PhyloSorPD)
  betadiv.ecrins.dist.PhyloSorPD.matrix.summits <- betadiv.ecrins.dist.PhyloSorPD.matrix[!rownames(betadiv.ecrins.dist.PhyloSorPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                         !colnames(betadiv.ecrins.dist.PhyloSorPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  
  ### Heatmap of summits UniFrac PD
  betadiv.ecrins.dist.UniFracPD <- betadiv.ecrins.dist.UniFracPD.tmp[,-1]
  colnames(betadiv.ecrins.dist.UniFracPD) <- betadiv.ecrins.dist.UniFracPD.tmp$X
  rownames(betadiv.ecrins.dist.UniFracPD) <- betadiv.ecrins.dist.UniFracPD.tmp$X
  betadiv.ecrins.dist.UniFracPD.matrix <- data.matrix(betadiv.ecrins.dist.UniFracPD)
  betadiv.ecrins.dist.UniFracPD.matrix.summits <- betadiv.ecrins.dist.UniFracPD.matrix[!rownames(betadiv.ecrins.dist.UniFracPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice"), 
                                                                                       !colnames(betadiv.ecrins.dist.UniFracPD.matrix) %in% c("Summits", "Ecrins NP", "Persistent", "Under Ice")]
  
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
  pdf(file=paste(outputFile, "SES_PhyloSor.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.PhyloSor.matrix.summits <- betadiv.ecrins.dist.PhyloSor.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.PhyloSor.matrix.summits))!=nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits)-1, 
                                                                                                     colSums(is.na(betadiv.ecrins.dist.PhyloSor.matrix.summits))!=ncol(betadiv.ecrins.dist.PhyloSor.matrix.summits)-1]
  ### Matrix with significant SES values
  mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSor.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSor.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.PhyloSor.matrix.summits, col=colors2, main=paste(clade, "SES_PhyloSor", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, breaks=bk2, scale="none")
  dev.off()
  
  ### PhyloSor Significant
  PhyloSor.sig.high <- length(which(betadiv.ecrins.dist.PhyloSor.matrix.summits > 1.96))/2 #0
  PhyloSor.sig.low <- length(which(betadiv.ecrins.dist.PhyloSor.matrix.summits < -1.96))/2 #24
  ### TOTAL pairwise comparisons
  total.pairwise <- (nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits)*(nrow(betadiv.ecrins.dist.PhyloSor.matrix.summits) - 1))/ 2
  ## 24 / 105
  
  pdf(file=paste(outputFile, "SES_UniFrac.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.UniFrac.matrix.summits <- betadiv.ecrins.dist.UniFrac.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.UniFrac.matrix.summits))!=nrow(betadiv.ecrins.dist.UniFrac.matrix.summits)-1, 
                                                                                             colSums(is.na(betadiv.ecrins.dist.UniFrac.matrix.summits))!=ncol(betadiv.ecrins.dist.UniFrac.matrix.summits)-1]
  mat <- matrix(ifelse(betadiv.ecrins.dist.UniFrac.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFrac.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.UniFrac.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.UniFrac.matrix.summits, col=colors2, main=paste(clade, "SES_UniFrac", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, breaks=bk2)
  dev.off()
  UniFrac.sig.high <- length(which(betadiv.ecrins.dist.UniFrac.matrix.summits > 1.96))/2 #0
  UniFrac.sig.low <- length(which(betadiv.ecrins.dist.UniFrac.matrix.summits < -1.96))/2 #24
  
  ########## TURN
  pdf(file=paste(outputFile, "SES_PhyloSorTurn.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.PhyloSorTurn.matrix.summits <- betadiv.ecrins.dist.PhyloSorTurn.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits))!=nrow(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits)-1, 
                                                                                                     colSums(is.na(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits))!=ncol(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits)-1]
  mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits, main=paste(clade, "SES_PhyloSor_turn", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, col=colors2, breaks=bk2)
  dev.off()
  PhyloSorTurn.sig.high <- length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits > 1.96))/2 #8
  PhyloSorTurn.sig.low <- length(which(betadiv.ecrins.dist.PhyloSorTurn.matrix.summits < -1.96))/2 #7
  
  
  pdf(file=paste(outputFile, "SES_UniFracTurn.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.UniFracTurn.matrix.summits <- betadiv.ecrins.dist.UniFracTurn.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.UniFracTurn.matrix.summits))!=nrow(betadiv.ecrins.dist.UniFracTurn.matrix.summits)-1, 
                                                                                                     colSums(is.na(betadiv.ecrins.dist.UniFracTurn.matrix.summits))!=ncol(betadiv.ecrins.dist.UniFracTurn.matrix.summits)-1]
  mat <- matrix(ifelse(betadiv.ecrins.dist.UniFracTurn.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFracTurn.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.UniFracTurn.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.UniFracTurn.matrix.summits, main=paste(clade, "SES_UniFrac_turn", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, col=colors2, breaks=bk2)
  dev.off()
  UniFracTurn.sig.high <- length(which(betadiv.ecrins.dist.UniFracTurn.matrix.summits > 1.96))/2 #7
  UniFracTurn.sig.low <- length(which(betadiv.ecrins.dist.UniFracTurn.matrix.summits < -1.96))/2 #8
  
  ######## PD
  pdf(file=paste(outputFile, "SES_PhyloSorPD.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.PhyloSorPD.matrix.summits <- betadiv.ecrins.dist.PhyloSorPD.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.PhyloSorPD.matrix.summits))!=nrow(betadiv.ecrins.dist.PhyloSorPD.matrix.summits)-1, 
                                                                                                   colSums(is.na(betadiv.ecrins.dist.PhyloSorPD.matrix.summits))!=ncol(betadiv.ecrins.dist.PhyloSorPD.matrix.summits)-1]
  mat <- matrix(ifelse(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96 | betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.PhyloSorPD.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.PhyloSorPD.matrix.summits, main=paste(clade, "SES_PhyloSor_PD", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, col=colors2, breaks=bk2)
  dev.off()
  PhyloSorPD.sig.high <- length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits > 1.96))/2 #8
  PhyloSorPD.sig.low <- length(which(betadiv.ecrins.dist.PhyloSorPD.matrix.summits < -1.96))/2 #16
  
  
  pdf(file=paste(outputFile, "SES_UniFracPD.pdf", sep = ""), onefile=FALSE)
  betadiv.ecrins.dist.UniFracPD.matrix.summits <- betadiv.ecrins.dist.UniFracPD.matrix.summits[rowSums(is.na(betadiv.ecrins.dist.UniFracPD.matrix.summits))!=nrow(betadiv.ecrins.dist.UniFracPD.matrix.summits)-1, 
                                                                                                 colSums(is.na(betadiv.ecrins.dist.UniFracPD.matrix.summits))!=ncol(betadiv.ecrins.dist.UniFracPD.matrix.summits)-1]
  mat <- matrix(ifelse(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96 | betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96, "*", ""), 
                nrow(betadiv.ecrins.dist.UniFracPD.matrix.summits))
  mat[is.na(mat)] <-  ""
  pheatmap(betadiv.ecrins.dist.UniFracPD.matrix.summits, main=paste(clade, "SES_UniFrac_PD", sep="_"), cluster_rows=T, 
           fontsize_row=fontsize_row, border_color=NA,  
           display_numbers = mat, col=colors2, breaks=bk2)
  dev.off()
  UniFracPD.sig.high <- length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits > 1.96))/2 #10
  UniFracPD.sig.low <- length(which(betadiv.ecrins.dist.UniFracPD.matrix.summits < -1.96))/2 #14
  
  
  beta.pair <- rbind(cbind(clade, total.pairwise, metric = "PhyloSor", sig.high = PhyloSor.sig.high, sig.low = PhyloSor.sig.low, 
                           PD.sig.high = PhyloSorPD.sig.high, PD.sig.low = PhyloSorPD.sig.low, Turn.sig.high = PhyloSorTurn.sig.high, Turn.sig.low = PhyloSorTurn.sig.low), 
                     cbind(clade, total.pairwise, metric = "UniFrac", sig.high = UniFrac.sig.high, sig.low = UniFrac.sig.low, 
                           PD.sig.high = UniFracPD.sig.high, PD.sig.low = UniFracPD.sig.low, Turn.sig.high = UniFracTurn.sig.high, Turn.sig.low = UniFracTurn.sig.low))
  print(beta.pair)
  
}


