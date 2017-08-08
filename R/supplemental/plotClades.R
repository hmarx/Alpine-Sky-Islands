#####################################################################################################################
############# Plot most species rich orders of plants within each species pool ######################################
############# Both Static & Dynamic Null Models #####################################################################
############# Hannah E. Marx, 22 March 2016 #########################################################################
#####################################################################################################################

source("R/chooseClade.R")

## ECRINS pool: 
ecrins_tax <- lookupTabPool(pool = "Regional", tip.labels = alps.phy$tip.label, tax = tax, level = "Angiospermae", taxonomy = "Angiospermae")
ecrins_tax.df <- as.data.frame(ecrins_tax)
head(ecrins_tax.df)
dim(ecrins_tax.df)

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

pdf(file="figs/supplemental/clade_counts_bar.pdf")
ggplot(lookup_sum[(lookup_sum$order) != "", ], aes(x=order, fill = Pool)) + 
  geom_bar() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  ylab("number of species")
dev.off()

Data <- as.data.frame((lookup_sum %>% group_by(Pool, order) %>% summarise(n = n_distinct(Species))))
head(Data)

Data.tab <- dcast(Data, Pool ~ order, value.var = "n")
write.csv(Data.tab, file="output/8_PhyoDiversity/cladePlots/clade.counts.csv")


