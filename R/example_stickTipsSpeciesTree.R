###################### Incorporate missing taxa into community phylo ######################
## For each missing species (n);
#### pick a branch (x) at random from a subtree that includes the genus + sister species
################# * might costrain within sup clade here...eg. endemics**
#### randomnly pick a time (t) from a uniform distibution of the branch length (x)
#### split branch x at time t 
#### add the unsampled taxon with branch length t 
#### recurisevly call the tree  with newly added branch until N new taxa have been incorporated 

###### Repeat the function X times to get a distibution of trees
## -> calculate PD metrics across the distribution of trees 

source("R/stickTipsSpeciesTree1.1.R")

## EXAMPLE :
sal=data(caudata) 
saltree <- caudata$phy #197
plot(saltree)
length(sal$dat)

res=congruify.phylo(sal$fam, sal$phy, sal$tax, tol=0, scale="PATHd8") 
print(res$calibrations)
plot(ladderize(res$phy,right=FALSE), cex=0.35, type="fan", label.offset=2.5) 
is.ultrametric(res$phy)

treePrune <- drop.random(phy = res$phy, n = 170) #27
plot(treePrune)
nodelabels()
tiplabels()
treedata(phy = treeDrop, data = caudata$dat)

treeDrop <- drop.random(phy = treePrune, n = 5) #22

#speciesAdd <- names(caudata$dat[which(!names(caudata$dat) %in% treeDrop$tip.label)])
speciesAdd <- treePrune$tip.label[which(!treePrune$tip.label %in% treeDrop$tip.label)]

tmp1 <- stickTipsSpeciesTree(table = caudata$tax, genephy = treeDrop, splist = speciesAdd)
plot(tmp1)

caudata$tax[rownames(caudata$tax) %in% c(treeDrop$tip.label, speciesAdd), ]

par(mfrow=c(1,2))
plot(treeDrop, cex=.5) ; title("Original")
plot(tmp1, cex=.5) ; title("Add Species")


## Simulate tree
# Ultrametric tree, tips.labes 
sim_tre <- read.tree(text="((((genus1_sp1:0.3333333333,(genus1_sp2:0.2222222222,(genus1_sp3:0.1111111111,genus1_sp4:0.1111111111):0.1111111111):0.1111111111):0.2222222222,(genus2_sp1:0.1111111111,genus2_sp2:0.1111111111):0.4444444444):0.1111111111,genus3:0.6666666667):0.3333333333,((genus4_sp1:0.1111111111,genus4_sp2:0.1111111111):0.1111111111,genus4_sp3:0.2222222222):0.7777777778);")
is.ultrametric(sim_tre)
plot.phylo(sim_tre); title("genus tree")
summary(sim_tre)
sum(sim_tre$edge.length)
sim_tre$edge
sim_tre$edge.length
nodelabels()
tiplabels()
my.subtrees <- subtrees(sim_tre)
plot(my.subtrees[[1]])
nodelabels()

## make correspondance table: can include taxa that are not in tree
## E.g. congrify taxonomy lookup table: genus, family, order, rownames = genus_species (or the tip label)
#
sim_tmp <- data.frame(genus = sort(paste("genus", rep(1:10, 5), sep="")), species = (paste("sp", sample(5:60, 50, replace=F),sep="")), 
                      family = sort(paste("fam",rep(1:5, 10),sep="")), order = sort(paste("order",rep(1:2, 25),sep="")))
sim_ctab <- unique(sim_tmp[,-2])
#sim_ctab <- sim_ctabtmp[,-1]
#rownames(sim_ctab) <- sim_ctabtmp[,1]
sim_ctab

## Species List: species to be added to the tree
sim_splist <- c(sample(apply(sim_tmp[ , 1:2 ] , 1 , paste , collapse = "_" ), 10), sim_tre$tip.label)
sim_splist


### Function: modified from stickTips_1.1.r
table <- sim_ctab
genephy <- sim_tre
splist <- sim_splist
outgroupList <- c("genus4_sp3", "genus4_sp2", "genus4_sp1")

treeRAND <-add.random(genephy,n=10)
plot(treeRAND)
summary(treeRAND)

test.stick <- stickTipsSpeciesTree(table = sim_ctab, genephy = sim_tre, splist = sim_splist)
plot(test.stick)
summary(test.stick)

trees <- replicate(20, stickTipsSpeciesTree(table = sim_ctab, genephy = sim_tre, splist = sim_splist), simplify=FALSE)
trees 
class(trees)<-"multiPhylo"

plot.multiPhylo(trees, 1)
plot(trees[[3]])

summary(sim_tre)
summary((trees[[9]]))

# Sum of branch lenghts will be higher for addes tips
sum(sim_tre$edge.length)
max((trees[[9]]$edge.length))

#Check to see id the depth of the tips to the root is the same:
distRoot(sim_tre, method = "patristic")
distRoot(trees[[3]], method = "patristic")

## Compare plots
plot(trees[[1]])
plot(trees[[9]])
plot.phylo(sim_tre); title("genus tree")


tree.summary <- lapply(1:length(trees), function(x) summary(trees[[x]])) 

is.ultrametric(sim_tre)
is.rooted(sim_tre)
is.ultrametric(trees[[1]])
is.rooted(trees[[1]])

untangle1<-function(x) reorder(reorder(x,"pruningwise"))
plot(untangle1(trees[[4]]))


par(mfrow=c(1,2))
plot(sim_tre) ; title("Original")
plot(untangle1(trees[[1]])) ; title("Add Species")


############### TEST OF TREE TOPOLOGY #####################

#[R-sig-phylo] Comparing the topology of two trees: Summary of responses
# treedist() and RF.dist() [phangorn] - Steel&Penny (1993)
# SH.test() [phangorn] - Shimodaira, H.&Hasegawa, M. (1999)
# dist.topo() [ape] - Penny&Hendy (1985), Kuhner&Felsenstein (1994)
# Geodesic distance between phylogenetic trees and associated functions
#[distory]
# Icong index (website:
#http://www.ese.u-psud.fr/utilisateurs/devienne/index.html) - de Vienne et al. (2007)


### Kishino-Hasegawa (KH) Test
# Null hypothesis: both topologies are equally well supported by the data (i.e. the difference in log- likelihoods, ??, is expected to be 0.0)



###### Robinson-Foulds Distance #######
RF.dist(trees[[1]], trees[[4]])
#The following R code read a newick file with multiple trees, and calculate the topological distance between all pairs
# The topological distance is defined as twice the number of internal branches defining different bipartitions of the tips (Penny and Hendy 1985).
dist.topo( trees[[1]], trees[[10]])
d.tree.PH85 = matrix( nrow=length(trees), ncol=length(trees) )
for ( i in 1:length(trees) ) {
  for (j in i:length(trees)) {
    print ( c(i, j))
    d.tree.PH85[i,j] = dist.topo( trees[[i]], trees[[j]], method="PH85")
  }
}
d.tree.PH85




#The branch length score may be seen as similar to the previous distance but taking branch lengths into account. 
#Kuhner and Felsenstein (1994) proposed to calculate the square root of the sum of the squared differences 
#of the (internal) branch lengths defining similar bipartitions (or splits) in both trees.
d.tree.score = matrix( nrow=length(trees), ncol=length(trees) )
for ( i in 1:length(trees) ) {
  for (j in i:length(trees)) {
    print ( c(i, j))
    d.tree.score[i,j] = dist.topo( trees[[i]], trees[[j]], method="score")
  }
}
d.tree.score

### Geodesic Distance Between Phylogenetic Trees
tree.dists <- dist.multiPhylo(trees)





#####  Sakin's Index
#convert multiphylo objec to treeshape
trees.tmp <- lapply(1:length(trees), function(x) as.treeshape(trees[[x]])) 

#The Sackin???s index is computed as the sum of the number of ancestors for each tips of the tree. 
#The less balanced a tree is and the larger its Sackin???s index. It can be normalized in order to obtain a statistic 
#that does not depend on the tree size, and so compare trees with different sizes. The normalization depends on the reference model (Yule or PDA). 
hist(sapply(trees.tmp,FUN=sackin), freq=FALSE,main="Histogram of Sakin's Indices",xlab="Sakin's Index")


## colless function: unnormalized index (default)
#The Colless??? index Ic computes the sum of absolute values |L ??? R| at each node of the tree where L (resp. R) 
#is the size of the left (resp. right) daughter clade at the node.
colless(as.treeshape(trees[[1]])) 
hist(sapply(trees.tmp,FUN=colless), freq=FALSE,main="Histogram of Colless' Indices",xlab="Colless' Index")

# Histogram of random clade sizes
hist(sapply(trees.tmp,FUN=cladesize),
     freq=FALSE,main="Random clade sizes for random generated trees",xlab="clade size")



