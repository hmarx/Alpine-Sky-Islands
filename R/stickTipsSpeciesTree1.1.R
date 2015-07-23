###################### Incorporate missing taxa into community phylo ######################

library(geiger)
library(picante)
library(spacodiR)
library(phytools)
library(phylobase)
library(phangorn)
library(apTreeshape)
library(distory)
library(taxize)

## For each missing species (n);
#### pick a branch (x) at random from a subtree that includes the genus + sister species
################# * might costrain within sup clade here...eg. endemics**
#### randomnly pick a time (t) from a uniform distibution of the branch length (x)
#### split branch x at time t 
#### add the unsampled taxon with branch length t 
#### recurisevly call the tree  with newly added branch until N new taxa have been incorporated 

###### Repeat the function X times to get a distibution of trees
## -> calculate PD metrics across the distribution of trees 

#### FUNCTION #####
# table: a dataframe of a taxonomic correspondance table (can include taxa that are not in tree); 
#   colnames = taxonomic levels (i.e. genus, family, order): 3 levels max
#   rownames = taxonomy names (i.e. genus1, family3...)
# genephy:  object of class phylo, ultrametric
# splist: a character string of species to be added to the tree, in format "Genus_species"

stickTipsSpeciesTree <- function(table, genephy, splist) {
  if(!is.ultrametric(genephy)){
    stop("this code requries an ultrametric tree")
  }
  # convert to lookuptable to characters
  table[,1] <- as.character(table[,1])
  table[,2] <- as.character(table[,2])
  table[,3] <- as.character(table[,3])
  
  ## Check to make sure table is in taxonomy order
  nochanges <- which(table[,1]==table[,2])
  if( sum(duplicated(c(unique(table[-nochanges,1]), unique(table[-nochanges,2])))) != 0 ){
    stop("two tips have the same name!\n")}   # this will avoid the creation of a weird object
  
  ## Unlist the tips of the phylogeny; get lookup table for the tips in phylogeny
  tree_tips_genus <- sapply(genephy$tip.label, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  tree_tips_genus
  ll_tree=match(tree_tips_genus, (table[,1]))
  looktree=table[ll_tree,]
  head(looktree)
  rownames(looktree)=names(tree_tips_genus)
  looktree <- as.data.frame(looktree) # taxonomy table for tree tips
  looktree[looktree == ""] <- NA
    
  # species names in to add list that are not in the tree tips
  outers.sp <- splist[!(splist %in% genephy$tip.label)]
  outers.sp 
  
  # Match outlier species to lookup table == species to stick
  tips_out=sapply(outers.sp, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll_out=match(tips_out, (table[,1]))
  lookout=table[ll_out,]
  rownames(lookout)=names(tips_out)
  lookout <- as.data.frame(lookout) ## table for species to stick into tree
  lookout[lookout == ""] <- NA
  
  #select tips names that are not in the table
  outers.tree <- unique(tree_tips_genus[!(tree_tips_genus %in% table[,1])])
  #warning if the tree contains more names than the table
  if (length(outers.tree)!=0) { warning(paste(outers.tree,"is/are not in your tab !\n"),call.=FALSE)    }
  
  genephy.tmp <- (genephy) 
  #plot(genephy.tmp, label.offset = 10)
  #nodelabels()
  #tiplabels()
  
  # sticking missing tips
  for(i in 1:nrow(lookout)){ #For each missing species (n)
    print(i)
    n <- lookout[i,] #species to add with taxnomy info

    # Find highest taxonomy it shares with spcies in tree
    if (!is.na(n[,1]) && n[,1] %in% looktree[,1]){
      print(paste("Species", rownames(n), "matches a", colnames(looktree[1]), sep=" "))
      cladeinclude <- looktree[looktree[,1] == as.character(n[,1]),] #the clade it should be included in
    } else {
      if (!is.na(n[,2]) && n[,2] %in% looktree[,2]) {
        print(paste("Species", rownames(n), "matches a", colnames(looktree[2]), sep=" "))
        cladeinclude <- looktree[looktree[,2] == as.character(n[,2]),]
      } else {
        if (!is.na(n[,3]) && n[,3] %in% looktree[,3]) {
          print(paste("Species", rownames(n), "matches an", colnames(looktree[3]), sep=" "))
          cladeinclude <- looktree[looktree[,3] == as.character(n[,3]),]   
        } else {
          print(paste("Species", rownames(n), "is not in taxonomy table", sep=" "))
          next
        }   
      }
    }
    
    ### prune subtree that includes the genus + sister species (so taxa can be sister to clade)
    # Node that subtends clade to place missing taxa 
    mrca.cladeinclude <- getMRCA(genephy.tmp, tip=c(rownames(cladeinclude))) #17
    
    if (is.null(mrca.cladeinclude)){ #if the subtending node is the root
      # List of all decendant nodes from root and node of clade to place taxa
      decendent.tip.nodes <- genephy.tmp$edge[which(genephy.tmp$edge[,2] == which(genephy.tmp$tip.label %in% rownames(cladeinclude))),]   
      nodeadd <- sample(decendent.tip.nodes, 1) 

            
    } else {
      mrca.sister <- getSisters(genephy.tmp, mrca.cladeinclude, mode="number")
      sister.decendents <- getDescendants(genephy.tmp, mrca.sister)
      sister.tips <- tips(genephy.tmp, mrca.sister)
      mrca.sisterandtarget <- getMRCA(genephy.tmp, c(c(rownames(cladeinclude)), sister.tips))
      
      # List of all decendants from node subtending clade to include and sister clade
      branchlist <- getDescendants(genephy.tmp, mrca.sisterandtarget)
      ## # pick a node (x) at random from node subtending target clade
      nodeadd <- sample(c(branchlist[which(!branchlist %in% c(mrca.sister, sister.decendents))]), 1) #mrca.sisterandtarget, 
      
    }
      
      # randomnly pick a time (t) from a uniform distibution of the branch length (x)
      randomLength <- runif(1, min=0.0001, max=0.9999) #min > 0, max < 1
      
      # add the unsampled taxon with branch length t
      position.length <- randomLength * genephy.tmp$edge.length[which(genephy.tmp$edge[,2]==nodeadd)] #branchadd
      
      # The missing taxa is added below the node (nodeadd), at lenght (position.length)
      tree.added <- bind.tip(tree = genephy.tmp, tip.label = rownames(n), where = nodeadd, position = position.length)
      #plot(tree.added)
      genephy.tmp <- tree.added # #recurisevly call the tree  with newly added branch until N new taxa have been incorporated 

  }
  
  return(genephy.tmp) # return the tree with the added tips

}
    

