##############################
## reduce communtiy phylogeny to clades of interest

alps.phy <- read.tree(file ="data/AnalysesDatasets/phy.Alpes.taxized.tre")

alps.sites <- read.csv(file="data/AnalysesDatasets/alps.sites.csv", row.names=1, header = T)
alps.sites <- data.matrix(alps.sites)

pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites) #traits = alps.traits

## Change community matrix to presence / absence 
tmp <- pezAlpes$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa <- data.matrix(cbind("taxa"=rownames(t(tmp)), t(tmp[4:nrow(tmp),])))
alps.phy <- pezAlpes$phy

####
tax=read.csv(file="data/AnalysesDatasets/fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

pruneCladeTaxonomyLookup <- function(tip.labels, tax, level, taxonomy){
  tips.ecrins=sapply(tip.labels, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips.ecrins, rownames(tax))
  ecrins_tax=tax[ll,]
  rownames(ecrins_tax)=names(tips.ecrins)
  ecrins_tax=as.matrix(ecrins_tax)
  ecrins_tax[is.na(ecrins_tax)]=""
  head(ecrins_tax)
  #length(which(ecrins_tax[,"Angiospermae"] == "Angiospermae")) # 1064 species are in Spermatophyta
  ecrins.clade <- names(which(ecrins_tax[,level] == taxonomy))
  #return(as.data.frame(ecrins_tax))
  return(ecrins.clade)
}

lookupTabPool <- function(pool, tip.labels, tax, level, taxonomy){
  tips.ecrins=sapply(tip.labels, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips.ecrins, rownames(tax))
  ecrins_tax=tax[ll,]
  rownames(ecrins_tax)=names(tips.ecrins)
  ecrins_tax=as.matrix(ecrins_tax)
  ecrins_tax[is.na(ecrins_tax)]=""
  ecrins_tax <- cbind(Pool = rep(pool, times = nrow(ecrins_tax)), Species = rownames(ecrins_tax), ecrins_tax)
  rownames(ecrins_tax) <- NULL
  return(ecrins_tax)
  
}


