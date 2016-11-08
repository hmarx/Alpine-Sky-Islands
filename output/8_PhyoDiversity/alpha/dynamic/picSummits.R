
library(ape)
library(caper)
library(geiger)
### Histogram of PIC of presence / absence to see clades with signal for many absences in short terminal branches

contrast <- sort(pic(x = as.numeric(alps.damocles$data[,16]), phy = alps.damocles$phy))

plot.phylo(alps.damocles$phy, show.node.label = T, show.tip.label = F)
c(1983, 1923, 1630, 1528, 1290, 1652)
tips(alps.damocles$phy, 1652)

hist(contrast, breaks = 10)
plot(contrast)

contrast.cap <- comparative.data(alps.damocles$phy, as.data.frame(alps.damocles$data[,c(1,16)]), taxa)

head(contrast.cap$data)
crunchMod <- crunch(formula = Summits ~ Summits, data=contrast.cap, factor.action = "allow")
crunchModRobust <- caic.robust(crunchMod)
caic.diagnostics(crunchModRobust, outlier=2)

### plots of parameter optimization tests:
sperma.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.sperma.summary.csv", stringsAsFactors=F)
aster.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.asterales.summary.csv", stringsAsFactors=F)
caryo.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.caryoph.summary.csv", stringsAsFactors=F)
lamial.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.lamiales.summary.csv", stringsAsFactors=F)
poales.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.poales.summary.csv", stringsAsFactors=F)
rosales.opt <- read.csv("output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_2/opt.rosales.summary.csv", stringsAsFactors=F)

sperma.opt[] <- lapply(sperma.opt, function(x){as.numeric(x)})
aster.opt[] <- lapply(aster.opt, function(x){as.numeric(x)})
caryo.opt[] <- lapply(caryo.opt, function(x){as.numeric(x)})
lamial.opt[] <- lapply(lamial.opt, function(x){as.numeric(x)})
poales.opt[] <- lapply(poales.opt, function(x){as.numeric(x)})
rosales.opt[] <- lapply(rosales.opt, function(x){as.numeric(x)})

write.csv(rbind(Spermatophyta = summary(sperma.opt)[4,], 
      Asterales = summary(aster.opt)[4,],
      Poales = summary(poales.opt)[4,],
      Rosales = summary(rosales.opt)[4,],
      Lamiales = summary(lamial.opt)[4,],
      Caryophyllales = summary(caryo.opt)[4,]), file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/clade.opt.summary.mean.csv")


plot(as.numeric(as.character(sperma.opt$gamma_0)), as.numeric(as.character((sperma.opt$mu))), )
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
                          
str(opt.summary.comb)
opt.plot.comb <- ggplot(opt.summary.comb, aes(x=as.numeric(as.character(gamma_0)), y=as.numeric(as.character(mu))))
opt.plot.comb <- opt.plot.comb + geom_point()
opt.plot.comb <- opt.plot.comb + facet_grid(clade ~ .)
opt.plot.comb <- opt.plot.comb + ylab(expression(mu))
opt.plot.comb <- opt.plot.comb + xlab(expression(gamma))
opt.plot.comb <- opt.plot.comb + theme_bw()
opt.plot.comb
pdf(file="output/9_PhyoDiversity/Spermatophyta/dynamic/parameter.opt.pdf")
opt.plot.comb
dev.off()
