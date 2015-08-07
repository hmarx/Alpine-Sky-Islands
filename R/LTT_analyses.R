
library(phytools)
library(paleotree)

## plotting the numbers of lineages through time
max(branching.times(alps.phy))
pdf(file="output/9_PhyoDiversity/ltt.pdf")
ltt.plot(phy = alps.phy, log="y")
dev.off()
#ltt.coplot(phy = treePL, log="y",tiplabels=F)



### Trait Disparity tthrough Time
pic.tree <- comparative.comm(phy = treePL, comm = alps.sites, 
                 env = alps.env, traits = (alps.data.complete)) #599

traits <- pic.tree$data
head(traits)
pic.shuffle<- function(x){
  x$tip.label <- sample(x$tip.label, length(x$tip.label), replace =F)
  pic(traits[x$tip.label, 1], x)
}

node.null <- replicate(999, pic.shuffle(x = pic.tree$phy))

combined.oss.null <- cbind(pic(traits[pic.tree$phy$tip.label, 1], pic.tree$phy), node.null)



my.ranks <- apply(combined.oss.null, MARGIN =1, rank)[1,]

times <- branching.times(pic.tree$phy)
plot(-1*times, my.ranks, pch =16)

my.dtt <- dtt(pic.tree$phy, traits[1:2], index="avg.sq", nsim=99)







