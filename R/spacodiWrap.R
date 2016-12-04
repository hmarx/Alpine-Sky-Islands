#####################################################################################################################
############# Calculate phylogenetic diversity at each node in phylogeny, and plot ##################################
############# beta-diversity at node times (spacodi) ################################################################
############# Hannah E. Marx, 7 Mar 2016 ############################################################################
#####################################################################################################################

spacodiWrap <- function(phy, sites, clade, output){
  ecrins.clade <- pruneCladeTaxonomyLookup(tip.labels = phy$tip.label, tax, level = "order", taxonomy = clade)
  pruned.tree.alps.clade <-drop.tip(phy, phy$tip.label[!phy$tip.label %in% ecrins.clade])
  alps.damocles.clade <- treedata(pruned.tree.alps.clade, sites)
  pruned.sites.alps.clade <- as.data.frame(alps.damocles.clade$data[,-1], strings.as.factors=F)
  
  ########### community diversity statistics of Hardy and Senterre (2007): tree-based  
  ### Prune community and phylogeny to just summits (for SES randomization from summit species pool)
  spac.com=as.data.frame(pruned.sites.alps.clade[-c(which(colnames(pruned.sites.alps.clade)%in%c("Ecrins NP", "Summits", "Under Ice", "Persistent")))])
  head(spac.com)
  spac.com[] <- apply(spac.com, 2, function(x) as.numeric(as.character(x)))
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
  
  pdf(file=paste(output, paste(clade, "PIst.branching.time.pdf", sep = "."), sep = "/"))
  spacodi.permutplot(sp.permut, bty="n", add.id = F, sig.plot=TRUE)
  dev.off()
  write.csv(sp.permut$observed.PIst, file=paste(output, paste(clade, "observed.PIst.csv", sep = "."), sep = "/"))
  write.csv(sp.permut$expected.PIst, file=paste(output, paste(clade, "expected.PIst.csv", sep = "."), sep = "/"))
  write.csv(sp.permut$randomization.test, file=paste(output, paste(clade, "randomization.test.csv", sep = "."), sep = "/"))
  
  #### plotting diversity turnover on trees
  pdf(file=paste(output, paste(clade, "PIst.phylo.nodes.pdf", sep = "."), sep = "/"), 10, 10)
  spacodi.treeplot(sp.permut, com.beta$sp.tree, sig.plot=TRUE, add.id=FALSE, type="fan", tip=1, pch=.5)
  dev.off()
  
}

