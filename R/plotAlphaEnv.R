#####################################################################################################################
############# Plot resression of SES alpha diveristy ~ PC of BioClim variables ######################################
############# Hannah E. Marx, 14 Mar 2016 ###########################################################################
#####################################################################################################################


plotAlphaEnv <- function(df, indepVar, filename, xx, dot=TRUE, eqn=TRUE){
  ## linear model (lm)
  lm_eqn = function(df){
    x <- df[,"obs.z"]
    y <- df[,indepVar]
    m = lm(x ~ y, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2*","~~italic(P)~"="~p, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)[8], digits = 3), #r.squared
                          p = signif(summary(m)[[4]][2,4], digits = 3))) #coef
    as.character(as.expression(eq));                 
  }
  
  eq.RD <- ddply(df %>% filter(model == "RD"), .(metric, pool, model, clade), lm_eqn)
  eq.RD.coord <- cbind(eq.RD, yy = rev(seq(length.out = length(unique(eq.RD$clade)), from = 3, by=.25))) #xx = rep(0, length(unique(eq.RD$clade))), 
  data.select = df %>% filter(model == "RD")
  
  p.RD <- ggplot(data = data.select, aes_string(x = indepVar, y = "obs.z", colour="clade")) +
    geom_smooth(method = "lm") + facet_grid(metric ~ pool + model)  
  if (dot & !eqn){
    p1.RD <- p.RD + geom_point(aes_string(colour = "clade")) 
    print(p1.RD)
    ggsave(filename = as.character(filename), p1.RD, width = 11, height = 8.5) 
  }  
  if (eqn & !dot){
    p2.RD <- p.RD + geom_text(data=eq.RD.coord, aes(x = xx, y = yy,label=V1), parse = TRUE, inherit.aes=T, size=3) 
    print(p2.RD)
    ggsave(filename = as.character(filename), p2.RD, width = 11, height = 8.5) 
  } else if (eqn & dot){
    p3.RD <- p.RD + geom_text(data=eq.RD.coord, aes(x = xx, y = yy, label=V1), parse = TRUE, inherit.aes=T, size=3) +
      geom_point(aes_string(colour = "clade")) 
    print(p3.RD)
    ggsave(filename = as.character(filename), p3.RD, width = 11, height = 8.5) 
  } else {
    print(p.RD)
    ggsave(filename = as.character(filename), p.RD, width = 11, height = 8.5) 
  }
}
