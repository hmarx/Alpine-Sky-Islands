#####################################################################################################################
############# Functions to manipulate megaphylogenies ####################################### 
############# Hannah E. Marx, 16 Nov 2015 ####################################### 
#####################################################################################################################

# ID subtrees for plotting support values on nodes of phylogeny 
# Code borrowed from Rich Glor (http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html)
## "The getAllSubTrees function below is a necessary subfunction that atomizes 
## a tree into each individual subclade and was provided compliments of Luke Harmon."
getAllSubtrees<-function(phy, minSize=2) {
  res<-list()
  count=1
  ntip<-length(phy$tip.label)
  for(i in 1:phy$Nnode) {
    l<-tips(phy, ntip+i)
    bt<-match(phy$tip.label, l)
    if(sum(is.na(bt))==0) {
      st<-phy} 
    else st<-drop.tip(phy, phy$tip.label[is.na(bt)])
    if(length(st$tip.label)>=minSize) {
      res[[count]]<-st
      count<-count+1
    }
  }
  res
}

## Modified version of plotBayesBoot from Rich Glor (http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html)
## Instead of Bayesian posterior probabilies, plots bootstrap support from ML analyses onto 
# time-scaled phylogeny from treePL output
### "The plotBayesBoot function below plots both posterior probability and bootstrap values on each node 
### of the consensus tree obtained from your Bayesian analysis. Bootstrap values will appear in bold text 
### immediately below and to the left of the node they support, whereas Bayesian posterior probabilies 
### will appear in regular face above and to the left of the node."
plotTreePLBoot <- function(treePL, bootTree, file) {
  getAllSubtrees(treePL)->treePLSub
  getAllSubtrees(bootTree)->bootSub
  #bootList<-matrix("<50",Nnode(treePL),1)
  bootList<-matrix("", Nnode(treePL),1)
  #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
  for(i in 1:Nnode(treePL)) {
    for(j in 1:Nnode(bootTree)) {
      match(treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)], bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)])->shared
      match(bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)], treePLSub[[i]]$tip.label[order(treePLSub[[i]]$tip.label)])->shared2
      if(sum(is.na(c(shared,shared2)))==0) {
        bootTree$node.label[j]->bootList[i]
      }}}
  treePLBS <- treePL
  treePLBS$node.label <- bootList
  plot(ladderize(treePLBS, right=F), cex=.10, lwd=0.1) #Plots your Bayesian consensus tree
  nodelabels(treePLBS$node.label, adj=c(1.2, -0.3), frame="n", cex=.2, font=2) #Adds bootstrap values.
  tree <- as.phylo(treePLBS)
  write.tree(tree, file=file) #SAVE as tree
}


#### Rerurn a vector of nodes in phylogeny that were matched to external time calibraitons with the Congruification approach
mrcaID=function(phy, cal){
  cal=as.matrix(cal)
  res=sapply(1:nrow(cal), function(idx){ ## loop over rows
    tips=cal[idx, c("taxonA", "taxonB")] ## fetch spanning taxa
    return(geiger:::.mrca(tips, phy)) ## MRCA of spanning taxa (node ID)
  })
  N=Ntip(phy)
  n=Nnode(phy)
  nn=integer(N+n) ## create empty vector of same length as branches in tree
  nn[res]=1 ## identify nodes that appear within calibrations
  nn=nn[-c(1:N)] ## exclude tip branches
  return(nn) ## return vector (ordered from first to last internal node in the tree)
}
