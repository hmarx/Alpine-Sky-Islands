#####################################################################################################################
############# Plot optimization of ML parameter estimates for rates of ##############################################
############# colonization (gamma) and local extinction (mu) for dynamic DAMOCLES null model ########################
############# Regional species pool to LGM community  ###############################################################
############# Hannah E. Marx, 12 Feb 2017 ###########################################################################
#####################################################################################################################

### read in ML estimates of parameters 
sperma.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.sperma.summary.csv", stringsAsFactors=F)
aster.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.asterales.summary.csv", stringsAsFactors=F)
caryo.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.caryoph.summary.csv", stringsAsFactors=F)
lamial.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.lamiales.summary.csv", stringsAsFactors=F)
rosales.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.rosales.summary.csv", stringsAsFactors=F)
poales.LGM.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/output/opt.poales.summary.csv", stringsAsFactors=F)

sperma.LGM.opt[] <- lapply(sperma.LGM.opt, function(x){as.numeric(x)})
aster.LGM.opt[] <- lapply(aster.LGM.opt, function(x){as.numeric(x)})
caryo.LGM.opt[] <- lapply(caryo.LGM.opt, function(x){as.numeric(x)})
lamial.LGM.opt[] <- lapply(lamial.LGM.opt, function(x){as.numeric(x)})
poales.LGM.opt[] <- lapply(poales.LGM.opt, function(x){as.numeric(x)})
rosales.LGM.opt[] <- lapply(rosales.LGM.opt, function(x){as.numeric(x)})

## summary statitisics of parameter optimization
write.csv(rbind(Spermatophyta = summary(sperma.LGM.opt)[4,], 
                Asterales = summary(aster.LGM.opt)[4,],
                Poales = summary(poales.LGM.opt)[4,],
                Rosales = summary(rosales.LGM.opt)[4,],
                Lamiales = summary(lamial.LGM.opt)[4,],
                Caryophyllales = summary(caryo.LGM.opt)[4,]), file="output/8_PhyoDiversity/alpha/dynamic/paramOpt/1_regionalEcrins/LGM/clade.LGM.opt.summary.mean.csv")


plot(as.numeric(as.character(sperma.LGM.opt$gamma_0)), as.numeric(as.character((sperma.LGM.opt$mu))))
plot(as.numeric(as.character(aster.LGM.opt$gamma_0)), as.numeric(as.character(aster.LGM.opt$mu)))
plot(as.numeric(as.character(caryo.LGM.opt$gamma_0)), as.numeric(as.character(caryo.LGM.opt$mu)))
plot(as.numeric(as.character(lamial.LGM.opt$gamma_0)), as.numeric(as.character(lamial.LGM.opt$mu)))
plot(as.numeric(as.character(poales.LGM.opt$gamma_0)), as.numeric(as.character(poales.LGM.opt$mu)))
plot(as.numeric(as.character(rosales.LGM.opt$gamma_0)), as.numeric(as.character(rosales.LGM.opt$mu)))


opt.summary.comb <- rbind(cbind(sperma.LGM.opt, clade = c(rep("Spermatophyta", n = 1:nrow(sperma.LGM.opt)))), 
                          cbind(aster.LGM.opt, clade = c(rep("Asterales", n = 1:nrow(aster.LGM.opt)))),
                          cbind(poales.LGM.opt, clade = c(rep("Poales", n = 1:nrow(poales.LGM.opt)))),
                          cbind(rosales.LGM.opt, clade = c(rep("Rosales", n = 1:nrow(rosales.LGM.opt)))),
                          cbind(lamial.LGM.opt, clade = c(rep("Lamiales", n = 1:nrow(lamial.LGM.opt)))),
                          cbind(caryo.LGM.opt, clade = c(rep("Caryophyllalaes", n = 1:nrow(caryo.LGM.opt)))))

opt.plot.comb <- ggplot(opt.summary.comb, aes(x=as.numeric(as.character(gamma_0)), y=as.numeric(as.character(mu))))
opt.plot.comb <- opt.plot.comb + geom_point()
opt.plot.comb <- opt.plot.comb + facet_grid(clade ~ .)
opt.plot.comb <- opt.plot.comb + ylab(expression(mu))
opt.plot.comb <- opt.plot.comb + xlab(expression(gamma))
opt.plot.comb <- opt.plot.comb + theme_bw()
opt.plot.comb
pdf(file="output/8_PhyoDiversity/alpha/dynamic/figs/parameter_regional_LGM.opt.pdf")
opt.plot.comb
dev.off()
