
## Plot parameter optimization for DAMOCLES
param.opt.aster <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.asterales.summary.csv")
param.opt.caryo <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.caryoph.summary.csv")
param.opt.poa <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.poales.summary.csv")
param.opt.ros <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.rosales.summary.csv")
param.opt.sperm <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.sperma.summary.csv")
param.opt.lam <- read.csv(file="output/9_PhyoDiversity/Spermatophyta/dynamic/paramOptClades/out_1/summary/opt.lamiales.summary.csv")


plot(as.numeric(as.character(param.opt.sperm$gamma_0)), as.numeric(as.character((param.opt.sperm$mu))))
plot(as.numeric(as.character(param.opt.aster$gamma_0)), as.numeric(as.character(param.opt.aster$mu)))
plot(as.numeric(as.character(param.opt.caryo$gamma_0)), as.numeric(as.character(param.opt.caryo$mu)))
plot(as.numeric(as.character(param.opt.lam$gamma_0)), as.numeric(as.character(param.opt.lam$mu)))
plot(as.numeric(as.character(param.opt.poa$gamma_0)), as.numeric(as.character(param.opt.poa$mu)))
plot(as.numeric(as.character(param.opt.ros$gamma_0)), as.numeric(as.character(param.opt.ros$mu)))
