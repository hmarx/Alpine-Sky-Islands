#####################################################################################################################
############# Historical drivers of phylogenetic beta diversity patterns ############################################ 
############# beta-diversity at node times (spacodi) ################################################################
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################

source("R/spacodiWrap.R")
source("R/plotBetaPair.R")
source("R/chooseClade.R")

################ spacodiR: PIst for each node in summit community phylogeny ################ 
### check branching times 
phy.nodetimes(pezAlpes$phy, time.range = c(0, max(alps.phy$edge.length)), proportion = TRUE)

########### community diversity statistics of Hardy and Senterre (2007): tree-based
### Ecrins Communtiy Phylogeny
ecrin.beta <- spacodi.calc(sp.plot = t(pezAlpes$comm), phy = pezAlpes$phy, pairwise = T)
ecrin.beta$PIst #0.002860691
ecrin.beta$pairwise.PIst
head(ecrin.beta$sp.plot)
(ecrin.beta$sp.tree)

### Prune community and phylogeny to just summits (for SES randomization from summit species pool)
spac.com=as.spacodi(pezAlpes$comm[-c(5, 12, 18, 19),])
head(spac.com)
colSums(spac.com)
spac.com.summits <- spac.com[!rowSums(spac.com) == 0,]
spac.phy <- treedata(pezAlpes$phy, spac.com.summits)[[1]]
com.beta <- spacodi.calc(sp.plot = spac.com.summits, phy = spac.phy, pairwise = T)
com.beta$PIst #0.004090516
com.beta$pairwise.PIst 

###spatial clustering: species within plots are more phylogenetically related on average
#than species from distinct plots where Pst > Ist, Bst > 0, or PIst > 0. Species are
#functionally more similar locally than those from distinct plots where Tst > Ist, Ust > 0,
#or TAUst > 0
###spatial overdispersion: species within plots are less phylogenetically related on average
#than species from distinct plots where Pst < Ist, Bst < 0, or PIst < 0. Species
#are functionally less similar locally than are species from distinct plots where Tst < Ist,
#Ust < 0, or TAUst < 0


############ Plotting observed and expected community structure across branching times of a phylogeny
# generate a plot of observed and expected PIst on summit community phylogeny 
# PIst: Bst analogue for presence/absence data, expressing phylogenetic turnover (independently of species turnover/richness).
sp.permut=spacodi.by.nodes(sp.plot=com.beta$sp.plot, phy=com.beta$sp.tree, sp.parm="PIst", n.rep=999, method = "1s") #shuffling of abundances across entire dataset
head(sp.permut$observed.PIst)
head(sp.permut$randomization.test)

pdf(file="output/8_PhyoDiversity/beta/Time/seedPlants/PIst.branching.time.pdf")
spacodi.permutplot(sp.permut, bty="n", add.id = F, sig.plot=TRUE)
dev.off()
write.csv(sp.permut$observed.PIst, file="output/8_PhyoDiversity/beta/Time/seedPlants/observed.PIst.csv")
write.csv(sp.permut$expected.PIst, file="output/8_PhyoDiversity/beta/Time/seedPlants/expected.PIst.csv")
write.csv(sp.permut$randomization.test, file="output/8_PhyoDiversity/beta/Time/seedPlants/randomization.test.csv")

#### plotting diversity turnover on trees
pdf(file = "output/8_PhyoDiversity/beta/Time/seedPlants/PIst.phylo.nodes.pdf", 10, 10)
spacodi.treeplot(sp.permut, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
dev.off()


################# Each Clade
spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Asterales", 
            output = "output/8_PhyoDiversity/beta/Time/clades/")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Caryophyllales", 
            output = "output/8_PhyoDiversity/beta/Time/clades/")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Lamiales", 
            output = "output/8_PhyoDiversity/beta/Time/clades/")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Poales", 
            output = "output/8_PhyoDiversity/beta/Time/clades/")

spacodiWrap(phy= alps.phy, sites = alps.sites.pa, clade = "Rosales", 
            output = "output/8_PhyoDiversity/beta/Time/clades/")
