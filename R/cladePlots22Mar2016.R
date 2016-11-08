
source("R/chooseClade.ecrinsPool.R")


## ECRINS pool: 
ecrins_tax <- lookupTabPool(pool = "Ecrins", tip.labels = alps.phy$tip.label, tax = tax, level = "Angiospermae", taxonomy = "Angiospermae")
ecrins_tax.df <- as.data.frame(ecrins_tax)
head(ecrins_tax.df)

## SUMMIT pool: 
summits_tax <- as.data.frame(lookupTabPool(pool = "All Summits", tip.labels = pezAlpes.summits$phy$tip.label, tax = tax, level = "Angiospermae", taxonomy = "Angiospermae"))
dim(summits_tax)

## PERSISTENT pool: 
lgm_tax <- as.data.frame(lookupTabPool(pool = "LGM", tip.labels = pezAlpes.persistent$phy$tip.label, tax = tax, level = "Angiospermae", taxonomy = "Angiospermae"))
dim(lgm_tax)

lookup_sum <- rbind(ecrins_tax.df, summits_tax, lgm_tax)
str(lookup_sum)

newOrder <- rev(names(sort(summary(lookup_sum$order))))

lookup_sum$order <- factor(lookup_sum$order, levels = newOrder)

pdf(file="output/9_PhyoDiversity/Spermatophyta/CladesPlots/clade_counts_bar.pdf")
ggplot(lookup_sum[(lookup_sum$order) != "", ], aes(x=order, fill = Pool)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ylab("number of species")
dev.off()


pdf(file="figs/supplemental/clade_counts_bar.pdf")
ggplot(lookup_sum[(lookup_sum$order) != "", ], aes(x=order, fill = Pool)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ylab("number of species")
dev.off()

Data <- as.data.frame((lookup_sum %>% group_by(Pool, order) %>% summarise(n = n_distinct(Species))))
head(Data)

Data.tab <- dcast(Data, Pool ~ order, value.var = "n")
write.csv(Data.tab, file="output/9_PhyoDiversity/Spermatophyta/CladesPlots/clade.counts.csv")




# http://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
# ggpie: draws a pie chart.
# give it:
# * `dat`: your dataframe  ggplot(dat, aes_st
# * `by` {character}: the name of the fill column (factor)
# * `totals` {character}: the name of the column that tracks
#    the time spent per level of `by` (percentages work too).
# returns: a plot object.

ggpie <- function (dat, by, totals, mycolours) {
  ggplot(dat, aes_string(x = factor(1), y=totals, fill=by)) +
    geom_bar(stat='identity', color='black') +
    guides(fill=guide_legend(override.aes=list(colour=NA))) + # removes black borders from legend
    scale_fill_manual(values = mycolours)  +
    coord_polar(theta='y') +
    theme(axis.ticks=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_text(colour='black'),
          axis.title=element_blank(),
          legend.position="none",
          panel.grid  = element_blank()) +
    scale_y_continuous(breaks=cumsum(dat[[totals]]) - dat[[totals]]/2, labels=dat[[by]]) 
}

summary(ecrins_tax.df)
dim(ecrins_tax.df)
df <- as.data.frame(summary(ecrins_tax.df$order))
df.tmp <- cbind(rownames(df), df[,1])
df.tmp <- as.data.frame(df.tmp)
colnames(df.tmp) <- c("order", "count")
df.tmp$order <- as.character(df.tmp$order)
df.tmp$count <- as.numeric(as.character(df.tmp$count))
df.tmp <- rbind(c("other", sum(as.data.frame(rbind(df.tmp[df.tmp$count < 24,], df.tmp[df.tmp$order == "",]))[2])), df.tmp[!df.tmp$count < 24,])
df.tmp <- df.tmp[!df.tmp$order == "",]
df.tmp$order <- as.character(df.tmp$order)
df.tmp$count <- as.numeric(as.character(df.tmp$count))
df.tmp <- df.tmp[order(df.tmp$count),]
sum(df.tmp$count)


mycolours1 <- c("Asterales" =  "#CC79A7", 
               "Poales" = "#D55E00", 
               "Rosales" = "#0072B2", 
               "Lamiales" = "#F0E442", 
               "Caryophyllales" =  "#009E73", 
               "other" = "#56B4E9", 
               "Brassicales" = "#E69F00",   
               "Malpighiales" = "#999999", 
               "Fabales" = "#CC79A7", 
               "Ericales" = "#D55E00", 
               "Saxifragales" = "#0072B2", 
               "Ranunculales" = "#F0E442", 
               "Gentianales" = "#009E73",  
               "Apiales" = "#56B4E9",
               "Asparagales" = "#E69F00", 
               "Dipsacales"= "#999999")

pdf(file="output/9_PhyoDiversity/Spermatophyta/CladesPlots/EcrinsCladeDiv.pdf", 11, 11)
ggpie(dat=df.tmp, by="order", totals="count",mycolours = mycolours1) +
  ggtitle("Ecrins NP") +
  theme(axis.ticks.margin=unit(0,"lines"),
        plot.margin=rep(unit(0, "lines"),4)) +
  theme(panel.background = element_rect(fill = "white")) 
dev.off()


p <- ggplot(theTable.fam, aes(family, fill = family))
p <- p+  geom_bar()
p <- p+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggplot(theTable.ord,
       aes(x = factor(""), fill = order) ) +
  geom_bar() +
  coord_polar(theta = "y") +
  scale_x_discrete("")




tips.ecrins.summits=sapply(pezAlpes.summits$phy$tip.label, function(x){
  unlist(strsplit(x,"_",fixed=TRUE))[1]
})
ll=match(tips.ecrins.summits, rownames(tax))
ecrins_tax.summits=tax[ll,]
rownames(ecrins_tax.summits)=names(tips.ecrins.summits)
ecrins_tax.summits=as.matrix(ecrins_tax.summits)
ecrins_tax.summits[is.na(ecrins_tax.summits)]=""
head(ecrins_tax.summits)
length(which(ecrins_tax.summits[,"Angiospermae"] == "Angiospermae")) # 211 species are in Spermatophyta
summits.angiosperm <- names(which(ecrins_tax.summits[,"Angiospermae"] == "Angiospermae"))
summits.aster <- names(which(ecrins_tax.summits[,"order"] == "Asterales"))

ecrins_tax.summits <- as.data.frame(ecrins_tax.summits)
summary(ecrins_tax.summits)
dim(ecrins_tax.summits)
df.sum <- as.data.frame(summary(ecrins_tax.summits$order))
df.tmp.summit <- cbind(rownames(df.sum), df.sum[,1])
df.tmp.summit <- as.data.frame(df.tmp.summit)
colnames(df.tmp.summit) <- c("order", "count")
df.tmp.summit$order <- as.character(df.tmp.summit$order)
df.tmp.summit$count <- as.numeric(as.character(df.tmp.summit$count))
df.tmp.summit <- rbind(c("other", sum(as.data.frame(rbind(df.tmp.summit[df.tmp.summit$count < 10,], df.tmp.summit[df.tmp.summit$order == "",]))[2])), df.tmp.summit[!df.tmp.summit$count < 10,])
df.tmp.summit <- df.tmp.summit[!df.tmp.summit$order == "",]
df.tmp.summit$order <- as.character(df.tmp.summit$order)
df.tmp.summit$count <- as.numeric(as.character(df.tmp.summit$count))
df.tmp.summit <- df.tmp.summit[order(df.tmp.summit$count),]
sum(df.tmp.summit$count)

pal <- wes_palette(14, name = "Zissou", type = "continuous")
image(volcano, col = pal)

wes_palettes$Darjeeling2

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


mycolours <- c("Asterales" =  "#CC79A7", 
               "Poales" = "#D55E00", 
               "Rosales" = "#0072B2", 
               "Lamiales" = "#F0E442", 
               "Caryophyllales" =  "#009E73", 
               "other" = "#56B4E9", 
               "Brassicales" = "#E69F00",   
               "Malpighiales" = "#999999", 
               "Fabales" = "#CC79A7")

"Ericales" = "#D55E00", 
"Saxifragales" = "#0072B2", 
"Ranunculales" = "#F0E442", 
"Gentianales" = "#009E73",  
"Apiales" = "#56B4E9"

pdf(file="output/9_PhyoDiversity/Spermatophyta/CladesPlots/EcrinsSumitsCladeDiv.pdf", 13, 13)
#par(mar=c(5,3,5,5)) #bottom, left, top and right, 
ggpie(dat=df.tmp.summit, by="order", totals="count", mycolours) +
  ggtitle("Ecrins Sumits") +
  theme(axis.ticks.margin=unit(0,"lines"),
        plot.margin=rep(unit(0, "lines"),4)) +
  theme(axis.text.x = element_text(size=15)) +
  theme(axis.title.x = element_text(hjust=5)) +
  theme(title = element_text(size=20)) +
  theme(panel.background = element_rect(fill = "white")) 
dev.off()

### dataframes of major lineage represnted on summits
df.tmp.summit
ecrins.angiosperm <- names(which(ecrins_tax[,"Angiospermae"] == "Angiospermae"))
ecrins.asterales <- names(which(ecrins_tax[,"order"] == "Asterales"))
ecrins.poales <- names(which(ecrins_tax[,"order"] == "Poales"))
ecrins.rosales <- names(which(ecrins_tax[,"order"] == "Rosales"))
ecrins.lamiales <- names(which(ecrins_tax[,"order"] == "Lamiales"))
ecrins.caryophyllales <- names(which(ecrins_tax[,"order"] == "Caryophyllales"))

39 + 37 + 22 +21 + 19 
138/215 #0.6418605


