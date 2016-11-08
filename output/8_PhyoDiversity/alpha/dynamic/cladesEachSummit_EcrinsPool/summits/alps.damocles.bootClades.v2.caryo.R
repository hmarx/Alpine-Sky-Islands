##############################
## Dynamic null model of communtiy assembly in the Ecrin NP, France

require(DAMOCLES)
require(ape)
require(pez)
require(geiger)
require(picante)

source("alps.damocles.chooseClade.R")

ecrins.caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = alps.phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")

#set parameter values
idparsopt = 1:2 #optimize extinction and colonization 

# inital mutation = max edge length * 1/10 (10 colonization or extinction per million years )
# == 227.2309 * 1/10 = 22 / mil yrs
initparsopt = c(0.1,0.1)

############ just on caryophyllales
pruned.tree.alps.caryophyllales <-drop.tip(pezAlpes$phy, pezAlpes$phy$tip.label[!pezAlpes$phy$tip.label %in% ecrins.caryophyllales])
alps.damocles.caryophyllales <- treedata(pruned.tree.alps.caryophyllales, alps.sites.pa)

obs.mntd.caryophyllales <- ses.mntd(samp =t(alps.damocles.caryophyllales$data)[c(3,16),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
obs.mpd.caryophyllales <- ses.mpd(samp =t(alps.damocles.caryophyllales$data)[c(3,16),], dis = cophenetic(alps.damocles.caryophyllales$phy), null.model = "phylogeny.pool", runs = 1000)
write.csv(obs.mntd.caryophyllales, file="obs.mntd.caryophyllales.csv")
write.csv(obs.mpd.caryophyllales, file="obs.mpd.caryophyllales.csv")

damo.bs.caryophyllales <- DAMOCLES_bootstrap(
  phy = alps.damocles.caryophyllales$phy,
  pa = alps.damocles.caryophyllales$data[,c(1, 16)],  #summits P/A,
  idparsopt = 1:length(initparsopt),
  parsfix = 0,
  idparsfix = (1:3)[-idparsopt],
  pars2 = c(1E-3,1E-4,1E-5,1000),
  pchoice = 0,
  runs = 1000,
  estimate_pars = FALSE,
  conf.int = 0.95
)

write.csv(damo.bs.caryophyllales$null_community_data, file="caryophyllales.null_community_data.csv")
write.csv(damo.bs.caryophyllales$summary_table, file="caryophyllales.summary_table.csv")



