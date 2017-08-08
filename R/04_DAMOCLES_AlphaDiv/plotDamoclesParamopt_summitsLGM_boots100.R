#####################################################################################################################
############# Plot optimization of ML parameter estimates for rates of ##############################################
############# colonization (gamma) and local extinction (mu) for dynamic DAMOCLES null model ########################
############# Regional species pool to LGM community  ###############################################################
############# Hannah E. Marx, 12 Feb 2017 ###########################################################################
#####################################################################################################################

## Observed:
### read in ML estimates of parameters 
sperma.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.sperma.summary.csv", stringsAsFactors=F)
aster.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.asterales.summary.csv", stringsAsFactors=F)
caryo.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.caryoph.summary.csv", stringsAsFactors=F)
lamial.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.lamiales.summary.csv", stringsAsFactors=F)
rosales.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.rosales.summary.csv", stringsAsFactors=F)
poales.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/best/output/opt.poales.summary.csv", stringsAsFactors=F)

sperma.opt[] <- lapply(sperma.opt, function(x){as.numeric(x)})
aster.opt[] <- lapply(aster.opt, function(x){as.numeric(x)})
caryo.opt[] <- lapply(caryo.opt, function(x){as.numeric(x)})
lamial.opt[] <- lapply(lamial.opt, function(x){as.numeric(x)})
poales.opt[] <- lapply(poales.opt, function(x){as.numeric(x)})
rosales.opt[] <- lapply(rosales.opt, function(x){as.numeric(x)})

## summary statitisics of parameter optimization
obs.param.summary <- rbind(cbind(summarize_each(na.omit(sperma.opt), funs(mean)), clade = "Spermatophyta"), 
             cbind(summarize_each(na.omit(aster.opt), funs(mean)), clade = "Asterales"), 
             cbind(summarize_each(na.omit(poales.opt), funs(mean)), clade = "Poales"), 
             cbind(summarize_each(na.omit(rosales.opt), funs(mean)), clade = "Rosales"),
             cbind(summarize_each(na.omit(lamial.opt), funs(mean)), clade = "Lamiales"),
             cbind(summarize_each(na.omit(caryo.opt), funs(mean)), clade = "Caryophyllales"))

## Organize output
# for i in {1..100}; do cat scale.$i.opt.params.asterales.* > scale.$i.aster.summary.csv; done
# for i in {1..100}; do cat scale.$i.opt.params.spermatophyta.* > scale.$i.sperma.summary.csv; done
# for i in {1..100}; do cat scale.$i.opt.params.lamiales.* > scale.$i.lam.summary.csv; done
# for i in {1..100}; do cat scale.$i.opt.params.caryophyllales.* > scale.$i.cary.summary.csv; done
# for i in {1..100}; do cat scale.$i.opt.params.rosales.* > scale.$i.ros.summary.csv; done
# for i in {1..100}; do cat scale.$i.opt.params.poales.* > scale.$i.poa.summary.csv; done

#scp hmarx@128.196.193.176:/home/hmarx/Documents/Ecrins/alpha_dynamic/summits_LGM/boots100/output/*summary.csv

null.param.df <- data.frame()
## read in summary for each tree
for (i in 1:100){
  ### read in ML estimates of parameters 
  sperma.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "sperma.summary.csv", sep="."), stringsAsFactors=F)
  aster.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "aster.summary.csv", sep="."), stringsAsFactors=F)
  caryo.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "cary.summary.csv", sep="."), stringsAsFactors=F)
  lamial.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "lam.summary.csv", sep="."), stringsAsFactors=F)
  rosales.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "ros.summary.csv", sep="."), stringsAsFactors=F)
  poales.LGM.opt <- read.csv(paste("output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/output/scale", i, "poa.summary.csv", sep="."), stringsAsFactors=F)
  
  sperma.LGM.opt[] <- lapply(sperma.LGM.opt, function(x){as.numeric(x)})
  aster.LGM.opt[] <- lapply(aster.LGM.opt, function(x){as.numeric(x)})
  caryo.LGM.opt[] <- lapply(caryo.LGM.opt, function(x){as.numeric(x)})
  lamial.LGM.opt[] <- lapply(lamial.LGM.opt, function(x){as.numeric(x)})
  poales.LGM.opt[] <- lapply(poales.LGM.opt, function(x){as.numeric(x)})
  rosales.LGM.opt[] <- lapply(rosales.LGM.opt, function(x){as.numeric(x)})
  
  ## summary statitisics of parameter optimization
  tmp <- rbind(cbind(summarize_each(na.omit(sperma.LGM.opt), funs(mean)), clade = "Spermatophyta", tree = i), 
        cbind(summarize_each(na.omit(aster.LGM.opt), funs(mean)), clade = "Asterales", tree = i), 
        cbind(summarize_each(na.omit(poales.LGM.opt), funs(mean)), clade = "Poales", tree = i), 
        cbind(summarize_each(na.omit(rosales.LGM.opt), funs(mean)), clade = "Rosales", tree = i),
        cbind(summarize_each(na.omit(lamial.LGM.opt), funs(mean)), clade = "Lamiales", tree = i),
        cbind(summarize_each(na.omit(caryo.LGM.opt), funs(mean)), clade = "Caryophyllales", tree = i))
  null.param.df <- rbind(null.param.df, tmp)
  
}

head(null.param.df)

# plot null mean mu agains observed from ML tree (dot)
mean_mu <- ggplot(data = null.param.df) +
  geom_boxplot(aes(x=clade, y=mu)) +
  geom_point(data = obs.param.summary, aes(x = clade, y = mu))
ggsave(mean_mu, file="output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/figs/mean_mu.pdf")

# plot null mean gamma_0 agains observed from ML tree (dot)
mean_gamma <- ggplot(data = null.param.df) +
  geom_boxplot(aes(x=clade, y=gamma_0)) +
  geom_point(data = obs.param.summary, aes(x = clade, y = gamma_0))
ggsave(mean_gamma, file="output/8_PhyoDiversity/alpha/dynamic/paramOpt/2_AllSummits/boots100/figs/mean_gamma.pdf")



######### Extras: plot mu and gamma for each rep
opt.summary.comb <- rbind(cbind(sperma.opt, clade = c(rep("Spermatophyta", n = 1:nrow(sperma.opt)))), 
                          cbind(aster.opt, clade = c(rep("Asterales", n = 1:nrow(aster.opt)))),
                          cbind(poales.opt, clade = c(rep("Poales", n = 1:nrow(poales.opt)))),
                          cbind(rosales.opt, clade = c(rep("Rosales", n = 1:nrow(rosales.opt)))),
                          cbind(lamial.opt, clade = c(rep("Lamiales", n = 1:nrow(lamial.opt)))),
                          cbind(caryo.opt, clade = c(rep("Caryophyllalaes", n = 1:nrow(caryo.opt)))))

opt.plot.comb <- ggplot(opt.summary.comb, aes(x=as.numeric(as.character(gamma_0)), y=as.numeric(as.character(mu))))
opt.plot.comb <- opt.plot.comb + geom_point()
opt.plot.comb <- opt.plot.comb + facet_grid(clade ~ .)
opt.plot.comb <- opt.plot.comb + ylab(expression(mu))
opt.plot.comb <- opt.plot.comb + xlab(expression(gamma))
opt.plot.comb <- opt.plot.comb + theme_bw()
opt.plot.comb
pdf(file="output/8_PhyoDiversity/alpha/dynamic/figs/parameter_allsummits_LGM.opt.pdf")
opt.plot.comb
dev.off()
