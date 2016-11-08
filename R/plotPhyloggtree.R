#http://guangchuangyu.github.io/2015/09/subsetting-data-in-ggtree/
#https://groups.google.com/forum/#!topic/bioc-ggtree/KUbqv3gFYW0
#http://guangchuangyu.github.io/2014/12/viewing-and-annotating-phylogenetic-tree-with-ggtree/

t(pezAlpes$comm)

p <- ggtree(pezAlpes$phy)
p <- p %<% endemics + geom_text(aes(shape=8, size=3, subset = (pezAlpes$phy %in% endemics))


plot.phylo(pezAlpes$phy, )


