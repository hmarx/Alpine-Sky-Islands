#####################################################################################################################
############# Plot optimization of ML parameter estimates for rates of ##############################################
############# colonization (gamma) and local extinction (mu) for dynamic DAMOCLES null model ########################
############# Hannah E. Marx, 6 June 2016 ###########################################################################
#####################################################################################################################

### read in ML estimates of parameters 
sperma.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.sperma.summary.csv", stringsAsFactors=F)
aster.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.asterales.summary.csv", stringsAsFactors=F)
caryo.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.caryoph.summary.csv", stringsAsFactors=F)
lamial.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.lamiales.summary.csv", stringsAsFactors=F)
poales.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.poales.summary.csv", stringsAsFactors=F)
rosales.opt <- read.csv("output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/opt.rosales.summary.csv", stringsAsFactors=F)

sperma.opt[] <- lapply(sperma.opt, function(x){as.numeric(x)})
aster.opt[] <- lapply(aster.opt, function(x){as.numeric(x)})
caryo.opt[] <- lapply(caryo.opt, function(x){as.numeric(x)})
lamial.opt[] <- lapply(lamial.opt, function(x){as.numeric(x)})
poales.opt[] <- lapply(poales.opt, function(x){as.numeric(x)})
rosales.opt[] <- lapply(rosales.opt, function(x){as.numeric(x)})

## summary statitisics of parameter optimization
write.csv(rbind(Spermatophyta = summary(sperma.opt)[4,], 
                Asterales = summary(aster.opt)[4,],
                Poales = summary(poales.opt)[4,],
                Rosales = summary(rosales.opt)[4,],
                Lamiales = summary(lamial.opt)[4,],
                Caryophyllales = summary(caryo.opt)[4,]), file="output/8_PhyoDiversity/alpha/dynamic/paramOpt/clades/clade.opt.summary.mean.csv")


lot(as.numeric(as.character(sperma.opt$gamma_0)), as.numeric(as.character((sperma.opt$mu))))
plot(as.numeric(as.character(aster.opt$gamma_0)), as.numeric(as.character(aster.opt$mu)))
plot(as.numeric(as.character(caryo.opt$gamma_0)), as.numeric(as.character(caryo.opt$mu)))
plot(as.numeric(as.character(lamial.opt$gamma_0)), as.numeric(as.character(lamial.opt$mu)))
plot(as.numeric(as.character(poales.opt$gamma_0)), as.numeric(as.character(poales.opt$mu)))
plot(as.numeric(as.character(rosales.opt$gamma_0)), as.numeric(as.character(rosales.opt$mu)))


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
pdf(file="output/8_PhyoDiversity/alpha/dynamic/figs/parameter.opt.pdf")
opt.plot.comb
dev.off()
