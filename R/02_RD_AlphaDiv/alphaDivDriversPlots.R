#####################################################################################################################
############# Plotting historical drivers of phylogenetic diversity patterns ########################################
############# Both Static & Dynamic Null Models #####################################################################
############# Hannah E. Marx, 6 June 2016 ###########################################################################
#####################################################################################################################

############################################ Priority Effects among clades ############################################ 
multi.tree.slices.alpha.df <- read.csv("output/8_PhyoDiversity/alpha/static/Time/multi.tree.slices.alpha.csv", as.is=T)
head(multi.tree.slices.alpha.df)

multi.tree.slices.alpha.df <- multi.tree.slices.alpha.df[!(multi.tree.slices.alpha.df$summit == "Under Ice"),]
multi.tree.slices.alpha.df <- multi.tree.slices.alpha.df[!(multi.tree.slices.alpha.df$summit == "Ecrins NP"),]
multi.tree.slices.alpha.df <- na.omit(multi.tree.slices.alpha.df)

multi.tree.slices.alpha.df$metric <- factor(multi.tree.slices.alpha.df$metric, labels = c( "mntd" = "MNTD", "mpd"="MPD"))


plot_alpha_time <- ggplot(multi.tree.slices.alpha.df, aes(x=slice, y=obs.z, colour=metric)) + 
  scale_color_manual(breaks = c("MNTD", "MPD"), values=c("MNTD" = "#AA7744", "MPD" = "#44AAAA"), name="Metric") +
  geom_point(aes(shape=ifelse(obs.p <= 0.05 | obs.p > 0.95, "dot", "no_dot"))) +
  scale_shape_manual(values=c(dot=19, no_dot=1), guide="none") +
  scale_x_continuous(limits = c(0, 100)) + 
  geom_smooth(method="lm") + 
  #geom_text(x = 66, y = 4.45, label = as.character(as.expression(eq1)), parse = TRUE, color = "black", family="Helvetica") +
  #geom_text(x = 67, y = 4.1, label = as.character(as.expression(eq2)), parse = TRUE, color = "black",  family="Helvetica") +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x = "time (mya)", y= "SES") 
#plot_alpha_time
pdf(file="output/8_PhyoDiversity/alpha/static/Time/plot_alpha_time.pdf")
plot_alpha_time
dev.off()
pdf(file="figs/Figure5a_plot_alpha_time.pdf")
plot_alpha_time
dev.off()

biotic.timeslice = multi.tree.slices.alpha.df %>% 
  group_by(metric) %>%
  do({model = lm(obs.z~slice, data=.)    # create your model
      data.frame(tidy(model),              # get coefficient info
                 glance(model))})

biotic.timeslice.df <- as.data.frame(biotic.timeslice %>% dplyr::select(metric, term, estimate, adj.r.squared, p.value.1)  %>% spread(term, estimate))
# column slice is slope (coefficient of variable area)

## P-value from t-test that compares slope to SE

#metric adj.r.squared    p.value.1 (Intercept)       slice
#1   mntd    0.21276500 1.229247e-14  -1.8525500 0.020583605
#2    mpd    0.01422435 3.391911e-02  -0.5985312 0.005138383

eq1 <- substitute(MNTD: italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p, 
                 list(a = format(biotic.timeslice.df[1,5], digits = 2), 
                      b = format(biotic.timeslice.df[1,4], digits = 2), 
                      r2 = format(biotic.timeslice.df[1,2], digits = 3),
                 p = signif(biotic.timeslice.df[1,3], digits = 3)))
as.character(as.expression(eq))

eq2 <- substitute(MPD: italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p, 
                  list(a = format(biotic.timeslice.df[2,5], digits = 2), 
                       b = format(biotic.timeslice.df[2,4], digits = 2), 
                       r2 = format(biotic.timeslice.df[2,2], digits = 3),
                       p = signif(biotic.timeslice.df[2,3], digits = 3)))
as.character(as.expression(eq))


############################################ Recent and rapid radiations ############################################ 
master.ses.static.NOendemics <- read.csv(file="output/8_PhyoDiversity/alpha/static/Endemics/master.SES.NOendemic.csv", row.names=1)
head(master.ses.static.NOendemics)
max(na.omit(master.ses.static.NOendemics$obs.z))


master.ses.static.NOendemics <- master.ses.static.NOendemics[!(master.ses.static.NOendemics$summits == "Under Ice"),]
unique(master.ses.static.NOendemics$summits)
master.ses.static.NOendemics$clade <- factor(master.ses.static.NOendemics$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
master.ses.static.NOendemics$pool <- factor(master.ses.static.NOendemics$pool, levels = c( "Ecrins NP", "Summits", "Persistent LGM"))

plot.mpd.NOendemics <- ggplot(master.ses.static.NOendemics, aes(y=clade, x=reorder(factor(summits), as.numeric(as.character(ntaxa))), fill=obs.z))
plot.mpd.NOendemics <- plot.mpd.NOendemics + geom_tile(colour="white", alpha=.75)
plot.mpd.NOendemics <- plot.mpd.NOendemics + scale_fill_gradient2(low="cyan2", high="red", mid="beige", na.value="white", limits=c(-5, 5))
plot.mpd.NOendemics <- plot.mpd.NOendemics + geom_point(aes(size=ifelse(sig, "dot", "no_dot")))
plot.mpd.NOendemics <- plot.mpd.NOendemics + scale_size_manual(values=c(dot=1, no_dot=NA), guide="none")
plot.mpd.NOendemics <- plot.mpd.NOendemics + theme_grey(base_size = 6) + labs(x = "",  y = "") 
plot.mpd.NOendemics <- plot.mpd.NOendemics + facet_grid(pool ~ metric)#,space="free",scales="free", as.table = F)   
plot.mpd.NOendemics <- plot.mpd.NOendemics + scale_x_discrete(expand = c(0, 0)) 
plot.mpd.NOendemics <- plot.mpd.NOendemics + scale_y_discrete(expand = c(0, 0)) 
plot.mpd.NOendemics <- plot.mpd.NOendemics + theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(size = base_size * 0.8, angle = -45, hjust = 0, colour = "black"))+
  guides(fill = guide_colorbar(barwidth = 1, barheight = 7,title.position = "top", title.hjust = 0.5, title = "SES"))
plot.mpd.NOendemics <- plot.mpd.NOendemics + coord_fixed(ratio=1)
plot.mpd.NOendemics <- plot.mpd.NOendemics + geom_point(aes(shape=ifelse(is.na(obs.z), "is_NA", "not_NA")))
plot.mpd.NOendemics <- plot.mpd.NOendemics + scale_shape_manual(values=c(is_NA=4, not_NA=NA), guide="none")
plot.mpd.NOendemics

#pdf(file="output/8_PhyoDiversity/alpha/static/Endemics/rd_SES_tile_NOendemics.pdf")
plot.mpd.NOendemics
#dev.off()

master.ses.static.NOendemics <- read.csv(file="output/8_PhyoDiversity/alpha/static/Endemics/master.SES.NOendemic.csv", row.names=1)
head(master.ses.static.NOendemics)
head(master.ses.alpha)
master.ses.static.NOendemics$pool <- factor(master.ses.static.NOendemics$pool, levels = c( "Ecrins NP",  "Summits", "Persistent LGM"))
master.ses.static.NOendemics$pool <- factor(master.ses.static.NOendemics$pool, labels = c( "Ecrins NP" = "Regional",  "Summits" = "All Summits", "Persistent LGM" = " LGM"))
master.ses.static.NOendemics$clade <- factor(master.ses.static.NOendemics$clade, levels = c( "Caryophyllales", "Lamiales", "Rosales", "Poales", "Asterales", "Spermatophyta"))
master.ses.static.NOendemics$metric <- factor(master.ses.static.NOendemics$metric, labels = c("mntd"="MNTD", "mpd"="MPD"))

alpha.endemics.df <- rbind(cbind(master.ses.alpha[,-1], subset = rep("All", times = nrow(master.ses.alpha))), 
                           cbind(master.ses.static.NOendemics, subset = rep("No Endemics", times = nrow(master.ses.static.NOendemics))))
head(alpha.endemics.df)
alpha.endemics.df <- alpha.endemics.df[!alpha.endemics.df$summits %in% c("Under Ice", "Summits", "Persistent"),] #reove this...not interesting
alpha.endemics.df$clade <- factor(alpha.endemics.df$clade, levels = c("Spermatophyta", "Asterales", "Poales", "Rosales",  "Lamiales", "Caryophyllales"))
head(alpha.endemics.df)

alpha.endemics.box <- ggplot(alpha.endemics.df, aes(x=clade, y=obs.z, fill=subset)) + 
  geom_boxplot() +
  facet_grid(metric ~ .) +
  labs(title = "Alpha Diveristy (RD) - Endemics") +
  labs(y = "SES") 
alpha.endemics.box
ggsave("figs/supplemental/alpha.endemics.box.pdf", alpha.endemics.box, width=8.5, height=8.5)


alpha.endemics.box.ecrins <- ggplot(alpha.endemics.df[alpha.endemics.df$pool == "Regional",], aes(x=clade, y=obs.z, fill=subset)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("All" = "#4477AA", "No Endemics"="#AA4455"), name="Species Subset", labels=c("All"="All species", "No Endemics" = "Without endemic species")) +
  facet_grid(metric ~ .) +
  labs(title = "Alpha Diveristy (RD) Regional Source Pool - Endemics") +
  labs(y = "SES")+
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    strip.text.y = element_text(size = 12, face="bold"),
    legend.position="bottom",
    axis.text.x = element_text(angle = -45, hjust = 0, colour = "black"))
alpha.endemics.box.ecrins
ggsave("figs/supplemental/alpha.endemics.box.ecrins.pdf", alpha.endemics.box.ecrins, width=8.5, height=8.5)
ggsave("figs/Figure5b_alpha.endemics.box.ecrins.pdf", alpha.endemics.box.ecrins, width=8.5, height=8.5)

alpha.endemics.box.summits <- ggplot(alpha.endemics.df[alpha.endemics.df$pool == "All Summits",], aes(x=clade, y=obs.z, fill=subset)) + 
  geom_boxplot() +
  facet_grid(metric ~ .) +
  labs(title = "Alpha Diveristy with and without Endemics \n(RD null - All Summits Species Pool)") +
  labs(y = "SES") 
alpha.endemics.box.summits
ggsave("figs/supplemental/alpha.endemics.box.summits.pdf", alpha.endemics.box.summits, width=8.5, height=8.5)

alpha.endemics.box.persistent <- ggplot(alpha.endemics.df[alpha.endemics.df$pool == " LGM",], aes(x=clade, y=obs.z, fill=subset)) + 
  geom_boxplot() +
  facet_grid(metric ~ .) +
  labs(title = "Alpha Diveristy without Endemics \n(RD null - LGM Species Pool)") +
  labs(y = "SES") 
alpha.endemics.box.persistent
ggsave("figs/supplemental/alpha.endemics.box.persistent.pdf", alpha.endemics.box.persistent, width=8.5, height=8.5)





