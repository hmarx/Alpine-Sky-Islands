
library(phytools)
library(paleotree)

## plotting the numbers of lineages through time
max(branching.times(pezAlpes$phy)) #352.2347
#pdf(file="output/9_PhyoDiversity/Spermatophyta/ltt.blank.pdf")
ltt.plot(phy = pezAlpes$phy, log="y", xaxt="n", yaxt = "n", ann=F)
#dev.off()
#ltt.coplot(phy = treePL, log="y",tiplabels=F)

## LTT for 5 clades
source("alps.damocles.chooseClade.R")
ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.rosales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Rosales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")

pruned.tree.alps.aster <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.asterales])
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.sites.pa)

pruned.tree.alps.poales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.poales])
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.sites.pa)

pruned.tree.alps.rosales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.rosales])
alps.damocles.rosales <- treedata(pruned.tree.alps.rosales, alps.sites.pa)

pruned.tree.alps.lamiales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.lamiales])
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.sites.pa)

pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.sites.pa)

pdf(file="output/9_PhyoDiversity/Spermatophyta/ltt.clades.pdf")
mltt.plot(pezAlpes$phy, alps.damocles.aster$phy, alps.damocles.poales$phy, alps.damocles.rosales$phy, 
          alps.damocles.lamiales$phy, alps.damocles.caryophyllales$phy, log="y", legend = F,  ylab = "ln(species number)", xlab = "Time (millions of years)")
dev.off()

ltt.plot(phy = alps.damocles.aster$phy, log="y")
ltt.plot(phy = pezAlpes$phy, log="y", xaxt="n", yaxt = "n", ann=F)
ltt.plot(phy = pezAlpes$phy, log="y", xaxt="n", yaxt = "n", ann=F)
ltt.plot(phy = pezAlpes$phy, log="y", xaxt="n", yaxt = "n", ann=F)
ltt.plot(phy = pezAlpes$phy, log="y", xaxt="n", yaxt = "n", ann=F)



### Trait Disparity through Time
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

plot(my.dtt)





