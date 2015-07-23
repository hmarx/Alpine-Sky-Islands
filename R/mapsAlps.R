
library(maps)
library(mapdata)
library(sp)
library(rgdal)
library(rgeos)
library(ggmap)
library(mapplots)
#install.packages("ggmap")
#install.packages("mapplots")

map("worldHires","France", col="gray90", fill=TRUE)

france <- map_data("france")
head(france)
unique(france$region)
states <- subset(france, region %in% c("Hautes-Alpes", "Isere"))

names(envAlps)

#pdf(file="output/7_Maps/FranceSummits.pdf") 
p <- ggplot()
p <- p + geom_polygon(data=france, aes(x=long, y=lat, group=group), colour="grey10", fill="grey10")
p <- p + geom_polygon(data=states, aes(x=long, y=lat, group=group), colour="white", fill="grey10")
p <- p + geom_point(data = envAlps, aes(x=X_WGS84, y=Y_WGS84, color='corall'))
p <- p + geom_point(aes(x=2.3508, y=48.8567), color="red")
p <- p + geom_point(aes(x=5.7222, y=45.2002), color="red")
p <- p + annotate("text", x=5.7222, y=45.5, color="red", label="Grenoble")
p <- p + annotate("text", x=2.3508, y=49.5, color="red", label="Paris")
p <- p + theme(legend.position = "none", axis.title = element_blank())
p
#dev.off()

#pdf(file="output/7_Maps/FranceStatesSummits.pdf") 
p <- ggplot()
p <- p + geom_polygon(data=states, aes(x=long, y=lat, group=group), colour="white", fill="grey10")
p <- p + geom_point(data = envAlps, aes(x=X_WGS84, y=Y_WGS84, color='corall'))
p <- p + coord_equal() 
p
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummits.pdf") 
p <- ggplot()
p <- p + geom_polygon(data=states, aes(x=long, y=lat, group=group), colour="white", fill="grey10")
p <- p + geom_point(data = envAlps, aes(x=X_WGS84, y=Y_WGS84, color='corall')) ##annual temperature
p <- p+ coord_cartesian(xlim=c(5.8, 6.5), ylim=c(44.5, 45.5)) 
p
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummitsPrecip.pdf") 
p <- ggplot()
p <- p + geom_polygon(data=states, aes(x=long, y=lat, group=group), colour="white", fill="grey10")
p <- p + geom_point(data = envAlps, aes(x=X_WGS84, y=Y_WGS84, color='corall', size=present_bio_1)) ##annual temperature
p <- p + coord_cartesian(xlim=c(5.8, 6.5), ylim=c(44.5, 45.5)) 
p
#dev.off()


# determine the bounding box of the spatial object
min(states$long)
min(states$lat)
max(states$long)
max(states$lat)

# get and plot a map
ecrins <- get_map(location = c(5.9, 44.7, 6.7, 45.1), source="google", maptype = "satellite", color="bw", zoom=10)

dispersion.phylonull <- merge(comm.sesmntd, comm.sesmpd.phylonull, by=0)
rownames(dispersion.phylonull) <- dispersion.phylonull$Row.names
env.dispersion <- merge(envAlps, dispersion.phylonull[2:ncol(dispersion.phylonull)], by=0, all.x=F)

#pdf(file="output/7_Maps/EcrinsSummitsTopoMNTD.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color=mntd.obs.z, size=ntaxa.x)) +  
  scale_colour_gradientn(colours=rainbow(5)) +
  scale_size_continuous(range = c(3, 10)) +
  ggtitle("SES mntd")
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummitsTopoMPD.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color=mpd.obs.z, size=ntaxa.x)) +  
  scale_colour_gradientn(colours=rainbow(5)) +
  scale_size_continuous(range = c(3, 10)) +
  ggtitle("SES mpd")
#dev.off()


#pdf(file="output/7_Maps/EcrinsSummitsTopoTemp.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color="coral", size=present_bio_1)) +  ##annual mean temperature 
  ggtitle("Annual Mean Temperature")
#dev.off()

#pdf(file="output/7_Maps/EcrinsSummitsTopoPrecip.pdf")
ggmap(ecrins) + 
  geom_point(data = env.dispersion, aes(x=X_WGS84, y=Y_WGS84, color="coral", size=present_bio_12)) +  ##annual mean temperature 
  ggtitle("Annual Mean Precipitation")
#dev.off()
