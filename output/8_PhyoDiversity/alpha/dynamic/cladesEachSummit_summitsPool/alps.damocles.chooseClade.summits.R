##############################
## reduce communtiy phylogeny to clades of interest
require(dplyr)

alps.phy <- read.tree(file ="phy.Alpes.taxized.tre")

alps.sites <- read.csv(file="alps.sites.csv", row.names=1, header = T)
alps.sites <- data.matrix(alps.sites)

pezAlpes <- comparative.comm(phy = alps.phy, comm = alps.sites) #traits = alps.traits

## Change community matrix to presence / absence 
tmp <- pezAlpes$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa <- data.matrix(cbind("taxa"=rownames(t(tmp)), t(tmp[4:nrow(tmp),])))
alps.phy <- pezAlpes$phy

##### SUMMITS
alps.sites.df <- as.data.frame(alps.sites)
## Contemporary species pool = summits 
summits.sites.tmp <- as.data.frame(cbind("taxa" = names(alps.sites.df), as.data.frame(t(alps.sites.df))))
summits.sites <- filter(summits.sites.tmp, Summits > 0)
head(summits.sites)
rownames(summits.sites) <- summits.sites$taxa
dim(summits.sites)
summits.sites <- t(summits.sites[-1])
summits.sites <- data.matrix(summits.sites)
tmp2 <- summits.sites
tmp2[which(tmp2 != 0)] <- 1
alps.summits.pa <- data.matrix(cbind("taxa"=rownames(t(tmp2)), t(tmp2[5:nrow(tmp2),])))


####
tax=read.csv(file="fleshed_genera.csv", as.is=TRUE, row=1) ##"linkage table", From Jon (NESCent working group on plant rates and traits)

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

