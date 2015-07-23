
#### from phylomorphoPlotsFINAL.R

require(mvtnorm)
require(snowfall)
require(nlme)
require(qpcR)
require(reshape)
require(calibrate)
library(ggbiplot)
require(adephylo)
library(devtools)
library(ggbiplot)
require(phytools)
require(phylobase)
library(picante)
library(cati)


############################ Prepare Data ########################################################  

head(pezAlpes$data) #1100
summary(pezAlpes$data)

alps.data.complete <- (pezAlpes$data[complete.cases(pezAlpes$data[c(3,5,7,9,10)]),c(3,5,7,9,10)])
head(alps.data.complete)
str(alps.data.complete)
levels(alps.data.complete$WOODY)

dim(alps.data.complete)
comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = (alps.data.complete)) #599
comparative.comm(phy = alps.phy, comm = summits.sites[-c(1:3),], env = alps.env, traits = (alps.data.complete)) #103
comparative.comm(phy = alps.phy, comm = persistent.sites[-3,], env = alps.env, traits = (alps.data.complete)) #82
comparative.comm(phy = alps.phy, comm = UnderIce.sites[-3,], env = alps.env, traits = (alps.data.complete)) #50

rownames(pezAlpes$comm)
head(as.data.frame(t(pezAlpes$comm))["Summits"])

tmp <- merge(alps.data.complete, as.data.frame(t(pezAlpes$comm))["Summits"],  by=0)
tmp[tmp["Summits"] > 0, "Summits"] <- 1
alps.data.complete <- tmp
rownames(alps.data.complete) <- tmp$Row.names

pezAlpes.Traits <- comparative.comm(phy = alps.phy, comm = alps.sites, env = alps.env, traits = (alps.data.complete)) #599
pezAlpes.Traits$data$PL_H_flora_MAX <- log(pezAlpes.Traits$data$PL_H_flora_MAX)
pezAlpes.Traits$data$SEEDM <- log(pezAlpes.Traits$data$SEEDM)
pezAlpes.Traits$data$PLOIDY <- as.numeric(pezAlpes.Traits$data$PLOIDY)
pezAlpes.Traits$data$VEG_DISP <- as.numeric(pezAlpes.Traits$data$VEG_DISP)
pezAlpes.Traits$data$WOODY <- as.numeric(pezAlpes.Traits$data$WOODY)


############################ Hill Smith coordinate of variation ############################ 

head(alps.data.complete)
alps.hs <- dudi.hillsmith(alps.data.complete[2:7], scannf = FALSE, nf =2)
scatter(alps.hs, cex.label=0.05)

##https://web.stanford.edu/class/bios221/labs/multivariate/lab_5_multivariate.html
ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")
# make the scree plot in a viewport
myscree <- function(eigs, x=0.8, y=0.1, just=c("right","bottom")){
  vp <- viewport(x=x, y=y, width=0.2, height=0.2, just=just)
  sp <- qplot(factor(1:length(eigs)), eigs, 
              geom="bar", stat="identity") +  
    labs(x = NULL, y = NULL)
  print(sp, vp=vp)
}

alpine <- factor(alps.data.complete$Summits, levels=c("Subalpine"= 0 , "Alpine" = 1))
pca1.dfs <- data.frame(alps.hs$li, alpine)

# multiply the loadings by 5 so they are more spread out
subject <- names(alps.data.complete[2:7])
pca1.dfl <- data.frame(6*alps.hs$co[,1:2])

#pdf(file="output/8_TraitDiveristy/pca.hillsmith.pdf") 
ppp + geom_point(data=pca1.dfs, aes(x=Axis1, y=Axis2, col=alpine)) 
myscree(alps.hs$eig / sum(alps.hs$eig))
#dev.off()

#pdf(file="output/8_TraitDiveristy/pca.hillsmith.labels.pdf") 
ppp + geom_point(data=pca1.dfs, aes(x=Axis1, y=Axis2, col=alpine)) +
 geom_text(data=pca1.dfl, aes(x=Comp1, y=Comp2, label=rownames(pca1.dfl)), cex=2)
myscree(alps.hs$eig / sum(alps.hs$eig))
#dev.off()

############################ Decomposition of traits variances using nested factors ############################ 
data(finch.ind)

dim(traits.finch)
#the trait matrix contains 2513 individuals values for 4 traits
table(sp.finch)
#the species names vector contains 2513 individuals belonging to 12 species
table(ind.plot.finch)
#the sites names vector contains 2513 individuals belonging to 6 sites (Islands)

old.par<-par(no.readonly = TRUE)
vec<- seq(1,length(sp.finch)*2, by = 2)
genus<-as.vector(unlist(strsplit(as.vector(sp.finch),"_"))[vec])
fact<-cbind(genus = as.factor(genus),
            species = as.factor(as.vector(sp.finch)),
            sites = as.factor(as.vector(ind.plot.finch)))
ptm <- proc.time()
res.partvar.finch<-partvar(traits = traits.finch, factors = fact)
## The partvar function decompose the variance accross nested scales. Thus choose
the order of the factors very carefully!
  proc.time_partvar <- proc.time() - ptm
res.partvar.finch
par(mfrow = c(2,2), mai = c(0.2,0.2,0.2,0.2)) #save graphical parameters
colors<-c(rgb(102,167,0, maxColorValue = 255),
          rgb(185,210,0, maxColorValue = 255),
          rgb(98,174,255, maxColorValue = 255),
          rgb(158,30,240, maxColorValue = 255))
piePartvar(res.partvar.finch, col = colors)

str(fact)


melt.alps <- AbToInd(pezAlpes.Traits$data, t(pezAlpes.Traits$comm), type.sp.val = "count")
melt.alps$traits
alps.traits.melt <- as.data.frame(data.matrix(melt.alps$traits))
str(alps.traits.melt)
colnames(alps.traits.melt) <- c("LogMaxHeight", "LogSeedMass", "Ploidy", "VegDisp", "Woody")

vec<- seq(1,length(melt.alps$sp)*2, by = 2)
genus<-as.vector(unlist(strsplit(as.vector(melt.alps$sp),"_"))[vec])
fact<-cbind(genus = as.factor(genus),
            sites = as.factor(as.vector((melt.alps$ind.plot))))
ptm <- proc.time()
res.partvar.finch <- partvar(traits = alps.traits.melt, factors = fact)
## The partvar function decompose the variance accross nested scales. Thus choose
# the order of the factors very carefully!
proc.time_partvar <- proc.time() - ptm
res.partvar.finch
par(mfrow = c(4,4), mai = c(0.2,0.2,0.2,0.2)) #save graphical parameters
colors<-c(rgb(102,167,0, maxColorValue = 255),
          rgb(185,210,0, maxColorValue = 255),
          rgb(98,174,255, maxColorValue = 255),
          rgb(158,30,240, maxColorValue = 255))
piePartvar(res.partvar.finch, col = colors)

par(new=TRUE)
par(mfrow = c(5,1), cex = 0.5)

pdf(file="output/8_TraitDiveristy/valueplots.pdf")
plotDistri(alps.traits.melt, rep("all_sp", times = dim(alps.traits.melt)[1]),
           melt.alps$ind.plot, ylim.cex = 3, plot.ask = F, cex.leg = 0.5)
dev.off()

#plotDistri(alps.traits.melt, rep("region", times = dim(alps.traits.melt)[1]),
 #          melt.alps$sp, ylim.cex = 6, plot.ask = F, leg = F)


res.alps <- Tstats(alps.traits.melt, ind.plot = melt.alps$ind.plot, sp = melt.alps$sp,
                    nperm = 9, print = FALSE)

pdf(file="output/8_TraitDiveristy/res.pdf")
plot(res.alps)
dev.off()




