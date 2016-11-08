

betaDivSummary <- function(betaOutput, geog.dist, nperm, clade, dependent, independent){
  phylosor.alps <- read.csv(betaOutput, row.names=1)
  
  #convert to distance matric
  phylosor.alps.summits <- as.dist(phylosor.alps[-c(which(rownames(phylosor.alps)%in%c("Ecrins NP","Summits", "Under Ice", "Persistent"))), 
                                                 -c(which(colnames(phylosor.alps)%in%c("Ecrins.NP","Summits", "Under.Ice", "Persistent")))])
  ordering <- labels(phylosor.alps.summits) #make sure distance matrices are in the same order for mantel, etc.
  spatial.dist.mat <- as.matrix(geog.dist)[ordering, ordering]
  geog.dist.dis <- as.dist(spatial.dist.mat)
  
  MRM_geog_dist <- MRM(phylosor.alps.summits ~ geog.dist.dis, nperm = nperm)
  MRM_geog_dist
  
  MRM_out <- cbind(dependent, independent, slope = signif(unlist(MRM_geog_dist)[2], digits = 3), 
                   intercept = signif(unlist(MRM_geog_dist)[1], digits = 3), r.squared.R2 = signif(unlist(MRM_geog_dist)[5], digits = 3),
                   r.squared.pval = signif(unlist(MRM_geog_dist)[6], digits = 3), F.test.F = signif(unlist(MRM_geog_dist)[7], digits = 3), 
                   F.test.F.pval = signif(unlist(MRM_geog_dist)[8], digits = 3))
  rownames(MRM_out) <- clade
  print(MRM_out)
}

  
  
  
