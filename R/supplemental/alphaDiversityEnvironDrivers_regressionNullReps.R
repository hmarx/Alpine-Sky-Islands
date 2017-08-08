#####################################################################################################################
############# Environmental predictors of phylogenetic diversity  ###################################################
############# Hannah E. Marx, 11 Feb 2017 ###########################################################################
#####################################################################################################################

############# Regression of SES on environmental variables ##########################################################
############# Observed SES compared to distibution from 999 replicates of SES on each summit ~ environmental variables

library(sjPlot)
library(sjmisc)
library(lme4)
#awk '
#    FNR==1 && NR!=1 { while (/^<header>/) getline; }
#    1 {print}
#' ses.static.nullCom* >ses.null.test.csv
#' 

#### observed 
head(alpha.div.master.env.pca.pca)

## Repeat SES n times, with 999 iterations each 
nullComOut <- read.csv("output/8_PhyoDiversity/alpha/static/nullCom/ses.null.test.csv", stringsAsFactors = F)
head(nullComOut)
str(nullComOut)

### Get environmental dataset cleaned up
geog.var <- merge(spatial.dist, alps.env.sprich.summits[, c("elevation", "S.area", "area", "max_slope", "simpson.d.geol")], by=0)
## PCA categorical variables
geog.env <- merge(geog.var, env.pca, by.x=1,  by.y=0)
  
nullComOut.env.pca.geog <- merge(nullComOut, geog.env, by.x=10, by.y=1)
#head(nullComOut.env.pca.geog)
#dim(nullComOut.env.pca.geog)
#str(nullComOut.env.pca.geog)
nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
#unique(nullComOut.env.pca.geog.filter$summits)
head(nullComOut.env.pca.geog.filter)

#################################### Regression ###################################################################

################################################  Geographic: Maximum elevation (elevation)  ###############################################################  

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.elevation <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$elevation <- log(nullComOut.env.pca.geog$elevation)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]

  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ elevation, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("elevation", nrow(j_coefs.tmp)))
  out.file.elevation <- rbind(out.file.elevation, j_coefs.tmp)
}
head(out.file.elevation)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ elevation, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
  
}
j_coefs.obs.elevation <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.elevation <- cbind(j_coefs.obs.elevation, variable = rep("elevation", nrow(j_coefs.tmp)))

plot(obs.z ~ elevation, data =  alpha.div.master.env.pca %>% filter(metric == "mntd" & clade == "Spermatophyta" & pool == "Ecrins NP"))
plot(obs.z ~ elevation, data =  nullComOut.env.pca.geog %>% filter(metric == "mntd" & clade == "Spermatophyta" & pool == "Ecrins NP"))


out.file.elevation$pool <- factor(out.file.elevation$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.elevation$clade <- factor(out.file.elevation$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.elevation$pool <- factor(j_coefs.obs.elevation$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.elevation$clade <- factor(j_coefs.obs.elevation$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

elevation.slope <- ggplot(filter(out.file.elevation, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.elevation, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ elevation: slope")
ggsave(elevation.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/elevation.mntd.slope.pdf")

elevation.intercept <- ggplot(filter(out.file.elevation, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.elevation, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ elevation: intercept")
ggsave(elevation.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/elevation.mntd.intercept.pdf")

elevation.slope <- ggplot(filter(out.file.elevation, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.elevation, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ elevation: slope")
ggsave(elevation.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/elevation.mpd.slope.pdf")

elevation.intercept <- ggplot(filter(out.file.elevation, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.elevation, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ elevation: elevation.intercept")
ggsave(elevation.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/elevation.mpd.intercept.pdf")




############################################### Geographic : Total summit area (S.area) ############################################### 

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.sarea <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$S.area <- log(nullComOut.env.pca.geog$S.area)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ S.area, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("S.area", nrow(j_coefs.tmp)))
  out.file.sarea <- rbind(out.file.sarea, j_coefs.tmp)
}
head(out.file.sarea)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ S.area, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.sarea <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.sarea <- cbind(j_coefs.obs.sarea, variable = rep("S.area", nrow(j_coefs.tmp)))

out.file.sarea$pool <- factor(out.file.sarea$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.sarea$clade <- factor(out.file.sarea$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.sarea$pool <- factor(j_coefs.obs.sarea$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.sarea$clade <- factor(j_coefs.obs.sarea$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

S.area.slope <- ggplot(filter(out.file.sarea, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ S.area: slope")
ggsave(S.area.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/S.area.mntd.slope.pdf")

S.area.intercept <- ggplot(filter(out.file.sarea, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ S.area: intercept")
ggsave(S.area.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/S.area.mntd.intercept.pdf")

S.area.slope <- ggplot(filter(out.file.sarea, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ S.area: slope")
ggsave(S.area.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/S.area.mpd.slope.pdf")

S.area.intercept <- ggplot(filter(out.file.sarea, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ S.area: S.area.intercept")
ggsave(S.area.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/S.area.mpd.intercept.pdf")


################################################  Geographic: Refugial area (area) ###############################################################  

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.refarea <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$area <- log(nullComOut.env.pca.geog$area)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ area, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("area", nrow(j_coefs.tmp)))
  out.file.refarea <- rbind(out.file.refarea, j_coefs.tmp)
}
head(out.file.refarea)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ area, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.refarea <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.refarea <- cbind(j_coefs.obs.refarea, variable = rep("area", nrow(j_coefs.tmp)))

out.file.refarea$pool <- factor(out.file.refarea$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.refarea$clade <- factor(out.file.refarea$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.refarea$pool <- factor(j_coefs.obs.refarea$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.refarea$clade <- factor(j_coefs.obs.refarea$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

area.slope <- ggplot(filter(out.file.refarea, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.refarea, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ area: slope")
ggsave(area.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/area.mntd.slope.pdf")

area.intercept <- ggplot(filter(out.file.refarea, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.refarea, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ area: intercept")
ggsave(area.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/area.mntd.intercept.pdf")

area.slope <- ggplot(filter(out.file.refarea, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.refarea, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ area: slope")
ggsave(area.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/area.mpd.slope.pdf")

area.intercept <- ggplot(filter(out.file.refarea, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.refarea, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ area: area.intercept")
ggsave(area.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/area.mpd.intercept.pdf")


################################################  Geographic: Maximum slope  (max_slope) ###############################################################  

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.max_slope <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$max_slope <- log(nullComOut.env.pca.geog$max_slope)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ max_slope, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("max_slope", nrow(j_coefs.tmp)))
  out.file.max_slope <- rbind(out.file.max_slope, j_coefs.tmp)
}
head(out.file.max_slope)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ max_slope, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.max_slope <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.max_slope <- cbind(j_coefs.obs.max_slope, variable = rep("max_slope", nrow(j_coefs.tmp)))

out.file.max_slope$pool <- factor(out.file.max_slope$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.max_slope$clade <- factor(out.file.max_slope$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.max_slope$pool <- factor(j_coefs.obs.max_slope$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.max_slope$clade <- factor(j_coefs.obs.max_slope$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

max_slope.slope <- ggplot(filter(out.file.max_slope, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.max_slope, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ max_slope: slope")
ggsave(max_slope.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/max_slope.mntd.slope.pdf")

max_slope.intercept <- ggplot(filter(out.file.max_slope, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.max_slope, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ max_slope: intercept")
ggsave(max_slope.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/max_slope.mntd.intercept.pdf")

max_slope.slope <- ggplot(filter(out.file.max_slope, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.max_slope, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ max_slope: slope")
ggsave(max_slope.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/max_slope.mpd.slope.pdf")

max_slope.intercept <- ggplot(filter(out.file.max_slope, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.max_slope, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ max_slope: max_slope.intercept")
ggsave(max_slope.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/max_slope.mpd.intercept.pdf")

################################################  Geographic: Heterogeneity Topography (pc.hetero.score) ###############################################################  

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.pc.hetero.score <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$pc.hetero.score <- log(nullComOut.env.pca.geog$pc.hetero.score)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ pc.hetero.score, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("pc.hetero.score", nrow(j_coefs.tmp)))
  out.file.pc.hetero.score <- rbind(out.file.pc.hetero.score, j_coefs.tmp)
}
head(out.file.pc.hetero.score)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ pc.hetero.score, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.pc.hetero.score <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.pc.hetero.score <- cbind(j_coefs.obs.pc.hetero.score, variable = rep("pc.hetero.score", nrow(j_coefs.tmp)))

out.file.pc.hetero.score$pool <- factor(out.file.pc.hetero.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.pc.hetero.score$clade <- factor(out.file.pc.hetero.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.pc.hetero.score$pool <- factor(j_coefs.obs.pc.hetero.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.pc.hetero.score$clade <- factor(j_coefs.obs.pc.hetero.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

pc.hetero.score.slope <- ggplot(filter(out.file.pc.hetero.score, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.hetero.score, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.hetero.score: slope")
ggsave(pc.hetero.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.hetero.score.mntd.slope.pdf")

pc.hetero.score.intercept <- ggplot(filter(out.file.pc.hetero.score, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.hetero.score, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.hetero.score: intercept")
ggsave(pc.hetero.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.hetero.score.mntd.intercept.pdf")

pc.hetero.score.slope <- ggplot(filter(out.file.pc.hetero.score, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.hetero.score, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.hetero.score: slope")
ggsave(pc.hetero.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.hetero.score.mpd.slope.pdf")

pc.hetero.score.intercept <- ggplot(filter(out.file.pc.hetero.score, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.hetero.score, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.hetero.score: pc.hetero.score.intercept")
ggsave(pc.hetero.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.hetero.score.mpd.intercept.pdf")

################################################  Geographic: Heterogeneity of lithology (simpson.d.geol)  ###############################################################  

path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.simpson.d.geol <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$simpson.d.geol <- log(nullComOut.env.pca.geog$simpson.d.geol)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ simpson.d.geol, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("simpson.d.geol", nrow(j_coefs.tmp)))
  out.file.simpson.d.geol <- rbind(out.file.simpson.d.geol, j_coefs.tmp)
}
head(out.file.simpson.d.geol)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ simpson.d.geol, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.simpson.d.geol <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.simpson.d.geol <- cbind(j_coefs.obs.simpson.d.geol, variable = rep("simpson.d.geol", nrow(j_coefs.tmp)))

out.file.simpson.d.geol$pool <- factor(out.file.simpson.d.geol$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.simpson.d.geol$clade <- factor(out.file.simpson.d.geol$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.simpson.d.geol$pool <- factor(j_coefs.obs.simpson.d.geol$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.simpson.d.geol$clade <- factor(j_coefs.obs.simpson.d.geol$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

simpson.d.geol.slope <- ggplot(filter(out.file.simpson.d.geol, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.simpson.d.geol, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ simpson.d.geol: slope")
ggsave(simpson.d.geol.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/simpson.d.geol.mntd.slope.pdf")

simpson.d.geol.intercept <- ggplot(filter(out.file.simpson.d.geol, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.simpson.d.geol, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ simpson.d.geol: intercept")
ggsave(simpson.d.geol.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/simpson.d.geol.mntd.intercept.pdf")

simpson.d.geol.slope <- ggplot(filter(out.file.simpson.d.geol, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.simpson.d.geol, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ simpson.d.geol: slope")
ggsave(simpson.d.geol.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/simpson.d.geol.mpd.slope.pdf")

simpson.d.geol.intercept <- ggplot(filter(out.file.simpson.d.geol, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.simpson.d.geol, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ simpson.d.geol: simpson.d.geol.intercept")
ggsave(simpson.d.geol.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/simpson.d.geol.mpd.intercept.pdf")

############################################# Environment: Available Energy (pc.availener.score) ############################################# 
path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.pc.availener.score <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$pc.availener.score <- log(nullComOut.env.pca.geog$pc.availener.score)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ pc.availener.score, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("pc.availener.score", nrow(j_coefs.tmp)))
  out.file.pc.availener.score <- rbind(out.file.pc.availener.score, j_coefs.tmp)
}
head(out.file.pc.availener.score)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ pc.availener.score, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.pc.availener.score <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.pc.availener.score <- cbind(j_coefs.obs.pc.availener.score, variable = rep("pc.availener.score", nrow(j_coefs.tmp)))


out.file.pc.availener.score$pool <- factor(out.file.pc.availener.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.pc.availener.score$clade <- factor(out.file.pc.availener.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.pc.availener.score$pool <- factor(j_coefs.obs.pc.availener.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.pc.availener.score$clade <- factor(j_coefs.obs.pc.availener.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

pc.availener.score.slope <- ggplot(filter(out.file.pc.availener.score, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.availener.score, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.availener.score: slope")
ggsave(pc.availener.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.availener.score.mntd.slope.pdf")

pc.availener.score.intercept <- ggplot(filter(out.file.pc.availener.score, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.availener.score, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.availener.score: intercept")
ggsave(pc.availener.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.availener.score.mntd.intercept.pdf")

pc.availener.score.slope <- ggplot(filter(out.file.pc.availener.score, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.availener.score, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.availener.score: slope")
ggsave(pc.availener.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.availener.score.mpd.slope.pdf")

pc.availener.score.intercept <- ggplot(filter(out.file.pc.availener.score, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.availener.score, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.availener.score: pc.availener.score.intercept")
ggsave(pc.availener.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.availener.score.mpd.intercept.pdf")

############################################# Environment: Stress (pc.stress.score) #############################################
path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.pc.stress.score <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$pc.stress.score <- log(nullComOut.env.pca.geog$pc.stress.score)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ pc.stress.score, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("pc.stress.score", nrow(j_coefs.tmp)))
  out.file.pc.stress.score <- rbind(out.file.pc.stress.score, j_coefs.tmp)
}
head(out.file.pc.stress.score)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ pc.stress.score, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.pc.stress.score <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.pc.stress.score <- cbind(j_coefs.obs.pc.stress.score, variable = rep("pc.stress.score", nrow(j_coefs.tmp)))

out.file.pc.stress.score$pool <- factor(out.file.pc.stress.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.pc.stress.score$clade <- factor(out.file.pc.stress.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.pc.stress.score$pool <- factor(j_coefs.obs.pc.stress.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.pc.stress.score$clade <- factor(j_coefs.obs.pc.stress.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

pc.stress.score.slope <- ggplot(filter(out.file.pc.stress.score, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stress.score, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.stress.score: slope")
ggsave(pc.stress.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stress.score.mntd.slope.pdf")

pc.stress.score.intercept <- ggplot(filter(out.file.pc.stress.score, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stress.score, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.stress.score: intercept")
ggsave(pc.stress.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stress.score.mntd.intercept.pdf")

pc.stress.score.slope <- ggplot(filter(out.file.pc.stress.score, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stress.score, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.stress.score: slope")
ggsave(pc.stress.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stress.score.mpd.slope.pdf")

pc.stress.score.intercept <- ggplot(filter(out.file.pc.stress.score, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stress.score, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.stress.score: pc.stress.score.intercept")
ggsave(pc.stress.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stress.score.mpd.intercept.pdf")

############################################# Environment: Present Stability (pc.stabil.score) #############################################
path = "output/8_PhyoDiversity/alpha/static/nullCom/output/"
out.file.pc.stabil.score <- data_frame()
file.names <- dir(path, pattern =".csv")
# loop over each RD replicate 
for (i in 1:length(file.names)){
  file <- read.csv(paste(path, file.names[i], sep=""), stringsAsFactors = T)
  
  ### Merge with environmental dataset 
  nullComOut.env.pca.geog <- merge(file, geog.env, by.x=10, by.y=1)
  nullComOut.env.pca.geog$obs <- as.numeric(nullComOut.env.pca.geog$obs)
  nullComOut.env.pca.geog$obs.z <- as.numeric(nullComOut.env.pca.geog$obs.z)
  nullComOut.env.pca.geog$rand.mean <- as.numeric(nullComOut.env.pca.geog$rand.mean)
  nullComOut.env.pca.geog$pc.stabil.score <- log(nullComOut.env.pca.geog$pc.stabil.score)
  nullComOut.env.pca.geog <- nullComOut.env.pca.geog[-15]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog[!nullComOut.env.pca.geog$summits == "Persistent",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Summits",]
  nullComOut.env.pca.geog.filter <- nullComOut.env.pca.geog.filter[!nullComOut.env.pca.geog.filter$summits == "Under Ice",]
  
  ###### fit a linear model to each predictor variable | metric, clade, and species pool
  le_lin_fit <- function(dat) {
    the_fit <- lm(obs.z ~ pc.stabil.score, dat)
    setNames(coef(the_fit), c("intercept", "slope"))
  }
  j_coefs.tmp <- ddply(nullComOut.env.pca.geog.filter, .(metric, clade, pool), le_lin_fit)
  j_coefs.tmp <- cbind(rep = rep(i, nrow(j_coefs.tmp)), j_coefs.tmp, variable = rep("pc.stabil.score", nrow(j_coefs.tmp)))
  out.file.pc.stabil.score <- rbind(out.file.pc.stabil.score, j_coefs.tmp)
}
head(out.file.pc.stabil.score)

###### fit a linear model to each predictor variable | metric, clade, and species pool
le_lin_fit <- function(dat) {
  the_fit <- lm(obs.z ~ pc.stabil.score, dat)
  setNames(coef(the_fit), c("intercept", "slope"))
}
j_coefs.obs.pc.stabil.score <- ddply(alpha.div.master.env.pca, .(metric, clade, pool), le_lin_fit)
j_coefs.obs.pc.stabil.score <- cbind(j_coefs.obs.pc.stabil.score, variable = rep("pc.stabil.score", nrow(j_coefs.tmp)))

plot(obs.z ~ pc.stabil.score, data =  alpha.div.master.env.pca %>% filter(metric == "mntd" & clade == "Spermatophyta" & pool == "Ecrins NP"))
plot(obs.z ~ pc.stabil.score, data =  nullComOut.env.pca.geog %>% filter(metric == "mntd" & clade == "Spermatophyta" & pool == "Ecrins NP"))


out.file.pc.stabil.score$pool <- factor(out.file.pc.stabil.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
out.file.pc.stabil.score$clade <- factor(out.file.pc.stabil.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))
j_coefs.obs.pc.stabil.score$pool <- factor(j_coefs.obs.pc.stabil.score$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
j_coefs.obs.pc.stabil.score$clade <- factor(j_coefs.obs.pc.stabil.score$clade, levels = c("Spermatophyta", "Asterales", "Poales","Rosales","Lamiales",  "Caryophyllales"))

pc.stabil.score.slope <- ggplot(filter(out.file.pc.stabil.score, metric == "mntd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stabil.score, metric == "mntd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.stabil.score: slope")
ggsave(pc.stabil.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stabil.score.mntd.slope.pdf")

pc.stabil.score.intercept <- ggplot(filter(out.file.pc.stabil.score, metric == "mntd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stabil.score, metric == "mntd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MNTD ~ pc.stabil.score: intercept")
ggsave(pc.stabil.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stabil.score.mntd.intercept.pdf")

pc.stabil.score.slope <- ggplot(filter(out.file.pc.stabil.score, metric == "mpd"), aes(x = slope)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stabil.score, metric == "mpd"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.stabil.score: slope")
ggsave(pc.stabil.score.slope, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stabil.score.mpd.slope.pdf")

pc.stabil.score.intercept <- ggplot(filter(out.file.pc.stabil.score, metric == "mpd"), aes(x = intercept)) + 
  geom_histogram(bins = 30) +
  geom_vline(data = filter(j_coefs.obs.pc.stabil.score, metric == "mpd"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(clade ~ pool) +
  labs(title = "MPD ~ pc.stabil.score: pc.stabil.score.intercept")
ggsave(pc.stabil.score.intercept, file = "output/8_PhyoDiversity/alpha/static/nullCom/figs/pc.stabil.score.mpd.intercept.pdf")





#########################################################################

enviro.regress.null <- rbind(out.file.elevation, out.file.sarea, out.file.refarea, out.file.max_slope,
                             out.file.max_slope, out.file.pc.hetero.score, out.file.simpson.d.geol, out.file.pc.availener.score,
                             out.file.pc.stress.score, out.file.pc.stabil.score)
  
enviro.regress.obs <- rbind(j_coefs.obs.elevation, j_coefs.obs.sarea, j_coefs.obs.refarea, j_coefs.obs.max_slope,
                             j_coefs.obs.pc.hetero.score, j_coefs.obs.simpson.d.geol, j_coefs.obs.pc.availener.score, 
                             j_coefs.obs.pc.stress.score, j_coefs.obs.pc.stabil.score)


enviro.regress.null.mpd.slope <- ggplot(filter(enviro.regress.null, metric == "mpd" & pool == "Regional" & clade == "Spermatophyta"), aes(x = slope)) + 
  geom_histogram(binwidth= .05) +
  geom_vline(data = filter(enviro.regress.obs, metric == "mpd" & pool == "Regional" & clade == "Spermatophyta"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(variable ~ .) +
  labs(title = "MPD ~ variable (slope)")
ggsave(enviro.regress.null.mpd.slope, file="output/9_Environment/alpha/regression/null/enviro.regress.null.mpd.slope.pdf", width = 5, height = 10)

enviro.regress.null.mpd.intercept <- ggplot(filter(enviro.regress.null, metric == "mpd" & pool == "Regional" & clade == "Spermatophyta"), aes(x = intercept)) + 
  geom_histogram(binwidth= .05) +
  geom_vline(data = filter(enviro.regress.obs, metric == "mpd" & pool == "Regional" & clade == "Spermatophyta"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(variable ~ .) +
  labs(title = "MPD ~ variable (intercept)")
ggsave(enviro.regress.null.mpd.intercept, file="output/9_Environment/alpha/regression/null/enviro.regress.null.mpd.intercept.pdf", width = 5, height = 10)

enviro.regress.null.mntd.slope <-ggplot(filter(enviro.regress.null, metric == "mntd" & pool == "Regional" & clade == "Spermatophyta"), aes(x = slope)) + 
  geom_histogram(binwidth= .05) +
  geom_vline(data = filter(enviro.regress.obs, metric == "mntd" & pool == "Regional" & clade == "Spermatophyta"), aes(xintercept=slope), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(variable ~ .) +
  labs(title = "MNTD ~ variable (slope)")
ggsave(enviro.regress.null.mntd.slope, file="output/9_Environment/alpha/regression/null/enviro.regress.null.mntd.slope.pdf", width = 5, height = 10)

enviro.regress.null.mntd.intercept <- ggplot(filter(enviro.regress.null, metric == "mntd" & pool == "Regional" & clade == "Spermatophyta"), aes(x = intercept)) + 
  geom_histogram(binwidth= .05) +
  geom_vline(data = filter(enviro.regress.obs, metric == "mntd" & pool == "Regional" & clade == "Spermatophyta"), aes(xintercept=intercept), color = "red", show.legend = F, lty=2) + #observed mean SESmntd from random draw null, DAMOCLES
  facet_grid(variable ~ .) +
  labs(title = "MNTD ~ variable (intercept)")
ggsave(enviro.regress.null.mntd.intercept, file="output/9_Environment/alpha/regression/null/enviro.regress.null.mntd.intercept.pdf", width = 5, height = 10)

fit1 <- lmer(obs.z ~ elevation + S.area + area + max_slope + pc.hetero.score + simpson.d.geol + (1 | clade), data = filter(nullComOut.env.pca.geog.filter, pool == "Ecrins NP"))
sjt.lmer(fit1)

fit1.mntd <- lm(obs.z ~ elevation + S.area + area + max_slope + pc.hetero.score + simpson.d.geol, data = filter(nullComOut.env.pca.geog.filter, metric =="mntd" & pool == "Ecrins NP" & clade=="Spermatophyta"))
fit1.mpd <- lm(obs.z ~ elevation + S.area + area + max_slope + pc.hetero.score + simpson.d.geol, data = filter(nullComOut.env.pca.geog.filter, metric =="mpd" & pool == "Ecrins NP" & clade=="Spermatophyta"))
sjt.lmer(fit1.mntd, fit1.mpd)



