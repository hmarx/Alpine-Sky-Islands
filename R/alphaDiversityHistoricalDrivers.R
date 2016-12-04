#####################################################################################################################
############# Historical drivers of phylogenetic alpha diversity patterns ###########################################
############# Both Static & Dynamic Null Models #####################################################################
############# Hannah E. Marx, 6 June 2016 ###########################################################################
#####################################################################################################################

source("R/chooseClade.R")

############################################ Priority Effects among clades ##########################################
############################################ Alpha diveristy at time intervals ######################################

max(branching.times(pezAlpes$phy)) ##352.2347
350 / 50 # 7    

trees.slice1 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 200)
plot(trees.slice1)
max(branching.times(trees.slice1))

trees.slice2 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 150)
plot(trees.slice2)
max(branching.times(trees.slice2))

trees.slice3 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 100)
plot(trees.slice3)
max(branching.times(trees.slice3))

trees.slice4 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 90)
plot(trees.slice4)
max(branching.times(trees.slice4))

trees.slice5 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 80)
plot(trees.slice5)
max(branching.times(trees.slice5))

trees.slice6 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 70)
plot(trees.slice6)
max(branching.times(trees.slice6))

trees.slice7 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 60)
plot(trees.slice7)
max(branching.times(trees.slice7))

trees.slice8 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 50)
plot(trees.slice8)
max(branching.times(trees.slice8))

trees.slice9 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 40)
plot(trees.slice9)
max(branching.times(trees.slice9))

trees.slice10 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 30)
plot(trees.slice10)
max(branching.times(trees.slice10))

trees.slice11 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 20)
plot(trees.slice11)
max(branching.times(trees.slice11))

trees.slice12 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 10)
plot(trees.slice12)
max(branching.times(trees.slice12))

trees.slice13 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 5)
plot(trees.slice13)
max(branching.times(trees.slice13))

trees.slice14 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 4)
plot(trees.slice14)
max(branching.times(trees.slice14))

trees.slice15 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 3)
plot(trees.slice15)
max(branching.times(trees.slice15))

trees.slice16 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 2)
plot(trees.slice16)
max(branching.times(trees.slice16))

trees.slice17 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 1)
plot(trees.slice17) #1012 / 1081 = 94%
max(branching.times(trees.slice17))

trees.slice18 <- timeSliceTree(ttree = pezAlpes$phy, sliceTime = 0)

multi.tree.slices <- c(trees.slice1, trees.slice2, trees.slice3, trees.slice4, trees.slice5, trees.slice6, trees.slice7, trees.slice8,
                       trees.slice9, trees.slice10, trees.slice11, trees.slice12, trees.slice13, trees.slice14, trees.slice15, trees.slice16,
                       trees.slice17, trees.slice18)

multi.tree.slices.mpd <- lapply(1:18, function(x) ses.mpd(samp = comparative.comm(phy = multi.tree.slices[[x]], 
                                                                                  com = pezAlpes$comm)$comm, dis = cophenetic(multi.tree.slices[[x]]), 
                                                          null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999))

multi.tree.slices.mntd <- lapply(1:18, function(x) ses.mntd(samp = comparative.comm(phy = multi.tree.slices[[x]], 
                                                                                    com = pezAlpes$comm)$comm, dis = cophenetic(multi.tree.slices[[x]]), 
                                                            null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999))

time.list <- (c(200, 150, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 4, 3, 2, 1, 0))

for (i in 1:length(multi.tree.slices.mpd)){
  multi.tree.slices.mpd[[i]] <- cbind(summit = rownames(multi.tree.slices.mpd[[i]]), multi.tree.slices.mpd[[i]], slice = rep(time.list[[i]], times = nrow(multi.tree.slices.mpd[[i]])))
}

for (i in 1:length(multi.tree.slices.mntd)){
  multi.tree.slices.mntd[[i]] <- cbind(summit = rownames(multi.tree.slices.mntd[[i]]), multi.tree.slices.mntd[[i]], slice = rep(time.list[[i]], times = nrow(multi.tree.slices.mntd[[i]])))
}


head(multi.tree.slices.mpd[[1]])

multi.tree.slices.mpd.all <- cbind(do.call(rbind, multi.tree.slices.mpd))
rownames(multi.tree.slices.mpd.all) <- NULL
multi.tree.slices.mpd.all <- cbind(multi.tree.slices.mpd.all, metric = rep("mpd", times= nrow(multi.tree.slices.mpd.all)))
colnames(multi.tree.slices.mpd.all) <- c("summit", "ntaxa",  "obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p", "runs", "slice", "metric")

multi.tree.slices.mntd.all <- do.call(rbind, multi.tree.slices.mntd)
rownames(multi.tree.slices.mntd.all) <- NULL
multi.tree.slices.mntd.all <- cbind(multi.tree.slices.mntd.all, metric = rep("mntd", times= nrow(multi.tree.slices.mntd.all)))
colnames(multi.tree.slices.mntd.all) <- c("summit", "ntaxa",  "obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p", "runs", "slice", "metric")

multi.tree.slices.alpha <- rbind(multi.tree.slices.mntd.all, multi.tree.slices.mpd.all)
head(multi.tree.slices.alpha)
#write.csv(multi.tree.slices.alpha, "output/8_PhyoDiversity/alpha/static/Time/multi.tree.slices.alpha.csv")



############################################ Recent and rapid diversification #######################################
############################################ Alpha diveristy without endemics #######################################

############################### Contemporary Source Pool : Regional Ecrins NP #################
###################################################################################### 
# Convert to pres/abs
tmp <- pezAlpes.NOendmics$comm
tmp[which(tmp != 0)] <- 1
alps.sites.pa.NOendemics <- t(tmp)

# Pruend phylogeny and community matrix for each clade 
ecrins.NOendemics.Spermatophyta <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "Spermatophyta", taxonomy = "Spermatophyta")
pruned.tree.ecrins.NOendemics.Spermatophyta <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Spermatophyta])
Spermatophyta.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Spermatophyta, alps.sites.pa.NOendemics)
Spermatophyta.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Spermatophyta.NOendemics.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Spermatophyta.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Spermatophyta.NOendemics.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Asterales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.ecrins.NOendemics.Asterales <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Asterales])
Asterales.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Asterales, alps.sites.pa.NOendemics)
Asterales.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Asterales.NOendemics.treedata$data), cophenetic.phylo(Asterales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Asterales.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Asterales.NOendemics.treedata$data), cophenetic.phylo(Asterales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Poales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.ecrins.NOendemics.Poales <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Poales])
Poales.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Poales, alps.sites.pa.NOendemics)
Poales.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Poales.NOendemics.treedata$data), cophenetic.phylo(Poales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Poales.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Poales.NOendemics.treedata$data), cophenetic.phylo(Poales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Rosales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "order", taxonomy = "Rosales")
pruned.tree.ecrins.NOendemics.Rosales <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Rosales])
Rosales.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Rosales, alps.sites.pa.NOendemics)
Rosales.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Rosales.NOendemics.treedata$data), cophenetic.phylo(Rosales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Rosales.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Rosales.NOendemics.treedata$data), cophenetic.phylo(Rosales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Lamiales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.ecrins.NOendemics.Lamiales <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Lamiales])
Lamiales.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Lamiales, alps.sites.pa.NOendemics)
Lamiales.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Lamiales.NOendemics.treedata$data), cophenetic.phylo(Lamiales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Lamiales.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Lamiales.NOendemics.treedata$data), cophenetic.phylo(Lamiales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.ecrins.NOendemics.Caryophyllales <-drop.tip(pezAlpes.NOendmics$phy, pezAlpes.NOendmics$phy$tip.label[!pezAlpes.NOendmics$phy$tip.label %in% ecrins.NOendemics.Caryophyllales])
Caryophyllales.NOendemics.treedata <- treedata(pruned.tree.ecrins.NOendemics.Caryophyllales, alps.sites.pa.NOendemics)
Caryophyllales.NOendemics.treedata.sesmpd.phylonull <- ses.mpd(t(Caryophyllales.NOendemics.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Caryophyllales.NOendemics.treedata.sesmntd.phylonull <- ses.mntd(t(Caryophyllales.NOendemics.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)


##### Combine all output of randomization from phylogeny pool
phylogeny.pool.NOendemic.SESmpd <- rbind((Spermatophyta.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                         (Asterales.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Asterales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                                         (Poales.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(Poales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Poales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                         (Rosales.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Rosales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                         (Lamiales.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Lamiales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                         (Caryophyllales.NOendemics.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.treedata.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.NOendemic.SESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.pool.NOendemic.SESmntd <- rbind((Spermatophyta.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                          (Asterales.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Asterales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                                          (Poales.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(Poales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Poales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                          (Rosales.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Rosales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                          (Lamiales.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Lamiales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                          (Caryophyllales.NOendemics.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.treedata.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(phylogeny.pool.NOendemic.SESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

phylogeny.pool.NOendemic.SES <- rbind(phylogeny.pool.NOendemic.SESmntd, phylogeny.pool.NOendemic.SESmpd)                               
head(phylogeny.pool.NOendemic.SES)     
phylogeny.pool.NOendemic.SES <- phylogeny.pool.NOendemic.SES[!phylogeny.pool.NOendemic.SES$summits == "Ecrins NP",]
phylogeny.pool.NOendemic.SES$clade <- factor(phylogeny.pool.NOendemic.SES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
phylogeny.pool.NOendemic.SES <- cbind(phylogeny.pool.NOendemic.SES, pool = rep(x = "Ecrins NP", times = nrow(phylogeny.pool.NOendemic.SES)))
#write.csv(phylogeny.pool.NOendemic.SES, file="output/8_PhyoDiversity/alpha/static/Endemics/phylogeny.pool.NO.endemic.SES.csv")



############################### Contemporary Source Pool : All Summits #################
######################################################################################
# Pres/abs for summits pool
tmp2 <- pezAlpes.NOendmics.summits$comm
tmp2[which(tmp2 != 0)] <- 1
alps.summits.pa.NOendemics <- t(tmp2)

# Pruend phylogeny and community matrix for each clade 
ecrins.NOendemics.summits.Spermatophyta <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "Spermatophyta", taxonomy = "Spermatophyta")
pruned.tree.ecrins.NOendemics.summits.Spermatophyta <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Spermatophyta])
Spermatophyta.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Spermatophyta, alps.summits.pa.NOendemics)
Spermatophyta.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Spermatophyta.NOendemics.summits.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Spermatophyta.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Spermatophyta.NOendemics.summits.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.summits.Asterales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.ecrins.NOendemics.summits.Asterales <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Asterales])
Asterales.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Asterales, alps.summits.pa.NOendemics)
Asterales.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Asterales.NOendemics.summits.treedata$data), cophenetic.phylo(Asterales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Asterales.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Asterales.NOendemics.summits.treedata$data), cophenetic.phylo(Asterales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.summits.Poales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.ecrins.NOendemics.summits.Poales <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Poales])
Poales.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Poales, alps.summits.pa.NOendemics)
Poales.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Poales.NOendemics.summits.treedata$data), cophenetic.phylo(Poales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Poales.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Poales.NOendemics.summits.treedata$data), cophenetic.phylo(Poales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.summits.Rosales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "order", taxonomy = "Rosales")
pruned.tree.ecrins.NOendemics.summits.Rosales <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Rosales])
Rosales.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Rosales, alps.summits.pa.NOendemics)
Rosales.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Rosales.NOendemics.summits.treedata$data), cophenetic.phylo(Rosales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Rosales.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Rosales.NOendemics.summits.treedata$data), cophenetic.phylo(Rosales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.summits.Lamiales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.ecrins.NOendemics.summits.Lamiales <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Lamiales])
Lamiales.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Lamiales, alps.summits.pa.NOendemics)
Lamiales.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Lamiales.NOendemics.summits.treedata$data), cophenetic.phylo(Lamiales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Lamiales.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Lamiales.NOendemics.summits.treedata$data), cophenetic.phylo(Lamiales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.summits.Caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.NOendmics.summits$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.ecrins.NOendemics.summits.Caryophyllales <-drop.tip(pezAlpes.NOendmics.summits$phy, pezAlpes.NOendmics.summits$phy$tip.label[!pezAlpes.NOendmics.summits$phy$tip.label %in% ecrins.NOendemics.summits.Caryophyllales])
Caryophyllales.NOendemics.summits.treedata <- treedata(pruned.tree.ecrins.NOendemics.summits.Caryophyllales, alps.summits.pa.NOendemics)
Caryophyllales.NOendemics.summits.treedata.sesmpd.phylonull <- ses.mpd(t(Caryophyllales.NOendemics.summits.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Caryophyllales.NOendemics.summits.treedata.sesmntd.phylonull <- ses.mntd(t(Caryophyllales.NOendemics.summits.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.summits.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)


##### Combine all output of randomization from phylogeny pool
summit.pool.NOendemic.SESmpd <- rbind((Spermatophyta.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                      (Asterales.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Asterales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                                      (Poales.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(Poales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Poales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                      (Rosales.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Rosales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                      (Lamiales.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Lamiales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                      (Caryophyllales.NOendemics.summits.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.summits.treedata.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.summits.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))))
names(summit.pool.NOendemic.SESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

summit.pool.NOendemic.SESmntd <- rbind((Spermatophyta.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                       (Asterales.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Asterales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                                       (Poales.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(Poales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Poales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                       (Rosales.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Rosales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                       (Lamiales.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Lamiales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                       (Caryophyllales.NOendemics.summits.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.summits.treedata.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.summits.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(summit.pool.NOendemic.SESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

summit.pool.NOendemic.SES <- rbind(summit.pool.NOendemic.SESmntd, summit.pool.NOendemic.SESmpd)                               
head(summit.pool.NOendemic.SES)     
summit.pool.NOendemic.SES <- summit.pool.NOendemic.SES[!summit.pool.NOendemic.SES$summits == "Ecrins NP",]
summit.pool.NOendemic.SES$clade <- factor(summit.pool.NOendemic.SES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
summit.pool.NOendemic.SES <- cbind(summit.pool.NOendemic.SES, pool = rep(x = "Summits", times = nrow(summit.pool.NOendemic.SES)))
write.csv(summit.pool.NOendemic.SES, file="output/9_PhyoDiversity/Endemics/summit.pool.NO.endemic.SES.csv")



############################### Contemporary Source Pool : LGM #################
######################################################################################
# Pres/abs for Persistent pool
tmp2 <- pezAlpes.persistent.NOendemics$comm
tmp2[which(tmp2 != 0)] <- 1
alps.Persistent.pa.NOendemics <- t(tmp2)

# Pruend phylogeny and community matrix for each clade 
ecrins.NOendemics.Persistent.Spermatophyta <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "Spermatophyta", taxonomy = "Spermatophyta")
pruned.tree.ecrins.NOendemics.Persistent.Spermatophyta <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Spermatophyta])
Spermatophyta.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Spermatophyta, alps.Persistent.pa.NOendemics)
Spermatophyta.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Spermatophyta.NOendemics.Persistent.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Spermatophyta.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Spermatophyta.NOendemics.Persistent.treedata$data), cophenetic.phylo(Spermatophyta.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Persistent.Asterales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "order", taxonomy = "Asterales")
pruned.tree.ecrins.NOendemics.Persistent.Asterales <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Asterales])
Asterales.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Asterales, alps.Persistent.pa.NOendemics)
Asterales.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Asterales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Asterales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Asterales.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Asterales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Asterales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Persistent.Poales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "order", taxonomy = "Poales")
pruned.tree.ecrins.NOendemics.Persistent.Poales <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Poales])
Poales.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Poales, alps.Persistent.pa.NOendemics)
Poales.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Poales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Poales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Poales.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Poales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Poales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Persistent.Rosales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "order", taxonomy = "Rosales")
pruned.tree.ecrins.NOendemics.Persistent.Rosales <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Rosales])
Rosales.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Rosales, alps.Persistent.pa.NOendemics)
Rosales.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Rosales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Rosales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Rosales.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Rosales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Rosales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Persistent.Lamiales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "order", taxonomy = "Lamiales")
pruned.tree.ecrins.NOendemics.Persistent.Lamiales <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Lamiales])
Lamiales.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Lamiales, alps.Persistent.pa.NOendemics)
Lamiales.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Lamiales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Lamiales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Lamiales.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Lamiales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Lamiales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)

ecrins.NOendemics.Persistent.Caryophyllales <- pruneCladeTaxonomyLookup(tip.labels = pezAlpes.persistent.NOendemics$phy$tip.label, tax, level = "order", taxonomy = "Caryophyllales")
pruned.tree.ecrins.NOendemics.Persistent.Caryophyllales <-drop.tip(pezAlpes.persistent.NOendemics$phy, pezAlpes.persistent.NOendemics$phy$tip.label[!pezAlpes.persistent.NOendemics$phy$tip.label %in% ecrins.NOendemics.Persistent.Caryophyllales])
Caryophyllales.NOendemics.Persistent.treedata <- treedata(pruned.tree.ecrins.NOendemics.Persistent.Caryophyllales, alps.Persistent.pa.NOendemics)
Caryophyllales.NOendemics.Persistent.treedata.sesmpd.phylonull <- ses.mpd(t(Caryophyllales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)
Caryophyllales.NOendemics.Persistent.treedata.sesmntd.phylonull <- ses.mntd(t(Caryophyllales.NOendemics.Persistent.treedata$data), cophenetic.phylo(Caryophyllales.NOendemics.Persistent.treedata$phy), null.model = "phylogeny.pool", abundance.weighted = FALSE, runs = 999)


##### Combine all output of randomization from phylogeny pool
persistent.pool.NOendemic.SESmpd <- rbind((Spermatophyta.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                          (Asterales.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Asterales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))), 
                                          (Poales.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Poales", sig=ifelse(Poales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Poales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                          (Rosales.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Rosales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                          (Lamiales.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Lamiales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))),
                                          (Caryophyllales.NOendemics.Persistent.treedata.sesmpd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.Persistent.treedata.sesmpd.phylonull), metric = "mpd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.Persistent.treedata.sesmpd.phylonull[,"mpd.obs.p"] > 0.95, 1,0))))
names(persistent.pool.NOendemic.SESmpd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

persistent.pool.NOendemic.SESmntd <- rbind((Spermatophyta.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Spermatophyta.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Spermatophyta", sig=ifelse(Spermatophyta.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Spermatophyta.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                           (Asterales.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Asterales.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Asterales", sig=ifelse(Asterales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Asterales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))), 
                                           (Poales.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Poales.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Poales", sig=ifelse(Poales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Poales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                           (Rosales.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Rosales.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Rosales", sig=ifelse(Rosales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Rosales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                           (Lamiales.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Lamiales.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Lamiales", sig=ifelse(Lamiales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Lamiales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))),
                                           (Caryophyllales.NOendemics.Persistent.treedata.sesmntd.phylonull %>%  mutate(summits = rownames(Caryophyllales.NOendemics.Persistent.treedata.sesmntd.phylonull), metric = "mntd", clade = "Caryophyllales", sig=ifelse(Caryophyllales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] <= 0.05 | Caryophyllales.NOendemics.Persistent.treedata.sesmntd.phylonull[,"mntd.obs.p"] > 0.95, 1,0))))
names(persistent.pool.NOendemic.SESmntd)[2:7] <- c("obs", "rand.mean", "rand.sd", "obs.rank", "obs.z", "obs.p")

persistent.pool.NOendemic.SES <- rbind(persistent.pool.NOendemic.SESmntd, persistent.pool.NOendemic.SESmpd)                               
head(persistent.pool.NOendemic.SES)     
persistent.pool.NOendemic.SES$clade <- factor(persistent.pool.NOendemic.SES$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
persistent.pool.NOendemic.SES <- cbind(persistent.pool.NOendemic.SES, pool = rep(x = "Persistent LGM", times = nrow(persistent.pool.NOendemic.SES)))
write.csv(persistent.pool.NOendemic.SES, file="output/9_PhyoDiversity/Endemics/persistent.pool.NO.endemic.SES.csv")

master.SES.NOendemic <- rbind(phylogeny.pool.NOendemic.SES, summit.pool.NOendemic.SES, persistent.pool.NOendemic.SES)
#write.csv(master.SES.NOendemic, file="output/8_PhyoDiversity/alpha/static/Endemics/master.SES.NOendemic.csv")



##################################  Statistical comparison of SES alpha diveristy with and without endemics ##################################  
#Paired t-test (use shapiro to test normal distribution of residuals)

master.ses.static.NOendemics <- read.csv(file="output/8_PhyoDiversity/alpha/static/Endemics/master.SES.NOendemic.csv", row.names=1)
head(master.ses.static.NOendemics)

master.ses.static.NOendemics <- master.ses.static.NOendemics[master.ses.static.NOendemics$summits %in% rownames(alps.env.sprich.summits),]
unique(master.ses.static.NOendemics$summits)
unique(master.ses.static.summtis$summits)

## All spermatophyta

t.test.sperm.endemic.no.ecrins.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                               master.ses.static.summtis$pool == "Ecrins NP" & 
                                                               master.ses.static.summtis$metric == "mntd", "obs.z"],
                                   master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                  master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                  master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.ecrins.mntd.out <- cbind(test = "Endemics_No", pool = "Ecrins", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.ecrins.mntd$estimate, t.value = t.test.sperm.endemic.no.ecrins.mntd$statistic, df = t.test.sperm.endemic.no.ecrins.mntd$parameter,
                                      p.value= t.test.sperm.endemic.no.ecrins.mntd$p.value, conf.low = t.test.sperm.endemic.no.ecrins.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.ecrins.mntd$conf.int[2])


t.test.sperm.endemic.no.ecrins.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                          master.ses.static.summtis$metric == "mpd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                             master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.ecrins.mpd.out <- cbind(test = "Endemics_No", pool = "Ecrins", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.ecrins.mpd$estimate, t.value = t.test.sperm.endemic.no.ecrins.mpd$statistic, df = t.test.sperm.endemic.no.ecrins.mpd$parameter,
                                                 p.value= t.test.sperm.endemic.no.ecrins.mpd$p.value, conf.low = t.test.sperm.endemic.no.ecrins.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.ecrins.mpd$conf.int[2])


t.test.sperm.endemic.no.summit.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Summits" & 
                                                                          master.ses.static.summtis$metric == "mntd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Summits" &
                                                                             master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.summit.mntd.out <- cbind(test = "Endemics_No", pool = "summit", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.summit.mntd$estimate, t.value = t.test.sperm.endemic.no.summit.mntd$statistic, df = t.test.sperm.endemic.no.summit.mntd$parameter,
                                                 p.value= t.test.sperm.endemic.no.summit.mntd$p.value, conf.low = t.test.sperm.endemic.no.summit.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.summit.mntd$conf.int[2])


t.test.sperm.endemic.no.summit.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                         master.ses.static.summtis$pool == "Summits" & 
                                                                         master.ses.static.summtis$metric == "mpd", "obs.z"],
                                             master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                            master.ses.static.NOendemics$pool == "Summits" &
                                                                            master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.summit.mpd.out <- cbind(test = "Endemics_No", pool = "summit", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.summit.mpd$estimate, t.value = t.test.sperm.endemic.no.summit.mpd$statistic, df = t.test.sperm.endemic.no.summit.mpd$parameter,
                                                p.value= t.test.sperm.endemic.no.summit.mpd$p.value, conf.low = t.test.sperm.endemic.no.summit.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.summit.mpd$conf.int[2])


t.test.sperm.endemic.no.persistent.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                          master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                          master.ses.static.summtis$metric == "mntd", "obs.z"],
                                              master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                             master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                             master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.sperm.endemic.no.persistent.mntd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Spermatophyta", metric="mntd", mean.of.diff = t.test.sperm.endemic.no.persistent.mntd$estimate, t.value = t.test.sperm.endemic.no.persistent.mntd$statistic, df = t.test.sperm.endemic.no.persistent.mntd$parameter,
                                                 p.value= t.test.sperm.endemic.no.persistent.mntd$p.value, conf.low = t.test.sperm.endemic.no.persistent.mntd$conf.int[1], conf.high =t.test.sperm.endemic.no.persistent.mntd$conf.int[2])


t.test.sperm.endemic.no.persistent.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade == "Spermatophyta" & 
                                                                         master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                         master.ses.static.summtis$metric == "mpd", "obs.z"],
                                             master.ses.static.NOendemics[master.ses.static.NOendemics$clade == "Spermatophyta" & 
                                                                            master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                            master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.sperm.endemic.no.persistent.mpd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Spermatophyta", metric="mpd", mean.of.diff = t.test.sperm.endemic.no.persistent.mpd$estimate, t.value = t.test.sperm.endemic.no.persistent.mpd$statistic, df = t.test.sperm.endemic.no.persistent.mpd$parameter,
                                                p.value= t.test.sperm.endemic.no.persistent.mpd$p.value, conf.low = t.test.sperm.endemic.no.persistent.mpd$conf.int[1], conf.high =t.test.sperm.endemic.no.persistent.mpd$conf.int[2])


### Each major clade separately

t.test.clades.endemic.no.ecrins.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                            master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                            master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                               master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                               master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.ecrins.mntd.out <- cbind(test = "Endemics_No", pool = "Ecrins NP", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.ecrins.mntd$estimate, t.value = t.test.clades.endemic.no.ecrins.mntd$statistic, df = t.test.clades.endemic.no.ecrins.mntd$parameter,
                                                   p.value= t.test.clades.endemic.no.ecrins.mntd$p.value, conf.low = t.test.clades.endemic.no.ecrins.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.ecrins.mntd$conf.int[2])


t.test.clades.endemic.no.ecrins.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                           master.ses.static.summtis$pool == "Ecrins NP" & 
                                                                           master.ses.static.summtis$metric == "mpd", "obs.z"],
                                               master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                              master.ses.static.NOendemics$pool == "Ecrins NP" &
                                                                              master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.ecrins.mpd.out <- cbind(test = "Endemics_No", pool = "Ecrins NP", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.ecrins.mpd$estimate, t.value = t.test.clades.endemic.no.ecrins.mpd$statistic, df = t.test.clades.endemic.no.ecrins.mpd$parameter,
                                                  p.value= t.test.clades.endemic.no.ecrins.mpd$p.value, conf.low = t.test.clades.endemic.no.ecrins.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.ecrins.mpd$conf.int[2])


t.test.clades.endemic.no.summits.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                               master.ses.static.summtis$pool == "Summits" & 
                                                                               master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                   master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                  master.ses.static.NOendemics$pool == "Summits" &
                                                                                  master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.summits.mntd.out <- cbind(test = "Endemics_No", pool = "Summits", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.summits.mntd$estimate, t.value = t.test.clades.endemic.no.summits.mntd$statistic, df = t.test.clades.endemic.no.summits.mntd$parameter,
                                                      p.value= t.test.clades.endemic.no.summits.mntd$p.value, conf.low = t.test.clades.endemic.no.summits.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.summits.mntd$conf.int[2])


t.test.clades.endemic.no.summits.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                              master.ses.static.summtis$pool == "Summits" & 
                                                                              master.ses.static.summtis$metric == "mpd", "obs.z"],
                                                  master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                 master.ses.static.NOendemics$pool == "Summits" &
                                                                                 master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.summits.mpd.out <- cbind(test = "Endemics_No", pool = "Summits", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.summits.mpd$estimate, t.value = t.test.clades.endemic.no.summits.mpd$statistic, df = t.test.clades.endemic.no.summits.mpd$parameter,
                                                     p.value= t.test.clades.endemic.no.summits.mpd$p.value, conf.low = t.test.clades.endemic.no.summits.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.summits.mpd$conf.int[2])


t.test.clades.endemic.no.persistent.mntd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                              master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                              master.ses.static.summtis$metric == "mntd", "obs.z"],
                                                  master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                 master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                                 master.ses.static.NOendemics$metric == "mntd", "obs.z"], paired = T)

t.test.clades.endemic.no.persistent.mntd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Clades", metric="mntd", mean.of.diff = t.test.clades.endemic.no.persistent.mntd$estimate, t.value = t.test.clades.endemic.no.persistent.mntd$statistic, df = t.test.clades.endemic.no.persistent.mntd$parameter,
                                                     p.value= t.test.clades.endemic.no.persistent.mntd$p.value, conf.low = t.test.clades.endemic.no.persistent.mntd$conf.int[1], conf.high =t.test.clades.endemic.no.persistent.mntd$conf.int[2])


t.test.clades.endemic.no.persistent.mpd <- t.test(master.ses.static.summtis[master.ses.static.summtis$clade != "Spermatophyta" & 
                                                                             master.ses.static.summtis$pool == "Persistent LGM" & 
                                                                             master.ses.static.summtis$metric == "mpd", "obs.z"],
                                                 master.ses.static.NOendemics[master.ses.static.NOendemics$clade != "Spermatophyta" & 
                                                                                master.ses.static.NOendemics$pool == "Persistent LGM" &
                                                                                master.ses.static.NOendemics$metric == "mpd", "obs.z"], paired = T)

t.test.clades.endemic.no.persistent.mpd.out <- cbind(test = "Endemics_No", pool = "Persistent LGM", clade = "Clades", metric="mpd", mean.of.diff = t.test.clades.endemic.no.persistent.mpd$estimate, t.value = t.test.clades.endemic.no.persistent.mpd$statistic, df = t.test.clades.endemic.no.persistent.mpd$parameter,
                                                    p.value= t.test.clades.endemic.no.persistent.mpd$p.value, conf.low = t.test.clades.endemic.no.persistent.mpd$conf.int[1], conf.high =t.test.clades.endemic.no.persistent.mpd$conf.int[2])


################
#The sign of a t-value tells us the direction of the difference in sample means

endemic.stats <- as.data.frame(rbind(t.test.sperm.endemic.no.ecrins.mntd.out, t.test.sperm.endemic.no.ecrins.mpd.out, t.test.sperm.endemic.no.summit.mntd.out,
      t.test.sperm.endemic.no.summit.mpd.out, t.test.sperm.endemic.no.persistent.mntd.out, t.test.sperm.endemic.no.persistent.mpd.out,
      t.test.clades.endemic.no.ecrins.mntd.out, t.test.clades.endemic.no.ecrins.mpd.out, t.test.clades.endemic.no.summits.mntd.out,
      t.test.clades.endemic.no.summits.mpd.out, t.test.clades.endemic.no.persistent.mntd.out, t.test.clades.endemic.no.persistent.mpd.out))
rownames(endemic.stats) <- NULL
head(endemic.stats)
str(biotic.corr.stats)

endemic.stats[5:9] <- apply(endemic.stats[5:9], 2, as.character)
endemic.stats[5:9] <- apply(endemic.stats[5:9], 2, as.numeric)

endemic.stats$p.value <- round(endemic.stats$p.value, 9)

#write.csv(endemic.stats, file="output/8_PhyoDiversity/alpha/static/Endemics/diff.SES.endemics.csv")


