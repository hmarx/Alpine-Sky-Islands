 
##############################
## Prune Ecrins phylogeny and community to clades of interest 

#source("alps.damocles.chooseClade.R")

ecrins.asterales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Asterales")
ecrins.poales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Poales")
ecrins.lamiales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")

############ just on asterales
pruned.tree.alps.aster <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.asterales]) #144 species in Ecrins
alps.damocles.aster <- treedata(pruned.tree.alps.aster, alps.summits.pa) # 39 species on summits

############ just on poales
pruned.tree.alps.poales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.poales]) # 153 species in Ecrins
alps.damocles.poales <- treedata(pruned.tree.alps.poales, alps.summits.pa) # 37 tips on summits

############ just on lamiales
pruned.tree.alps.lamiales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.lamiales]) #109 species in Ecrins
alps.damocles.lamiales <- treedata(pruned.tree.alps.lamiales, alps.summits.pa) # 21 tips on summits

############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(alps.phy, alps.phy$tip.label[!alps.phy$tip.label %in% ecrins.caryophyllales]) # 74 species in Ecrins
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.summits.pa) # 19 tips on summits

