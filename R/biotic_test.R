

max(branching.times(phy = pezAlpes$phy)) #352.2347

350 / 50

300-50

is.rooted(pezAlpes$phy)
phy.deresolve(phy = pezAlpes$phy, time.range = 2)

function (phy, time.range = c(0, 0), relative = TRUE) 
{
  if (is.null(phy$edge.length)) 
    stop("The tree does not appear to have branch lengths")
  if (class(phy) != "phylo") 
    stop("The tree does not appear to be a valid phylo object")
  if (length(time.range) > 2) 
    stop("Cannot interpret the range of time with more than two elements")
  if (length(time.range) == 1) 
    time.range = c(0, time.range)
  time.range = time.range[order(time.range)]
  bb <- branching.times(phy)
  if (relative) 
    bb = bb/max(bb)
  inr = sapply(bb, function(x) withinrange(x, time.range[1], 
                                           time.range[2]))
  ind <- as.numeric(names(bb[inr]))
  if (any(ind == Ntip(phy) + 1)) 
    ind = ind[-which(ind == Ntip(phy) + 1)]
  n <- length(ind)
  if (!n) {
    return(phy)
  }
  else {
    ind.tmp = match(ind, phy$edge[, 2])
    ind = ind.tmp[!is.na(ind.tmp)]
  }
  orig.edge = phy$edge
  orig.phy = phy
  ntips = Ntip(phy)
  reedge <- function(ancestor, des.to.drop) {
    wh <- which(phy$edge[, 1] == des.to.drop)
    dd <- which(orig.edge[, 2] == des.to.drop)
    dropped.branch <- phy$edge.length[dd]
    d.d <- c(get.desc.of.node(des.to.drop, orig.phy))
    if (length(d.d)) 
      phy$edge.length[match(d.d, orig.edge[, 2])] <<- phy$edge.length[match(d.d, 
                                                                            orig.edge[, 2])] + dropped.branch
    for (k in wh) {
      if (phy$edge[k, 2] %in% node.to.drop) {
        reedge(ancestor, phy$edge[k, 2])
      }
      else {
        phy$edge[k, 1] <<- ancestor
      }
    }
  }
  node.to.drop <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node.to.drop) 
      next
    reedge(anc[i], node.to.drop[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node.to.drop)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node.to.drop < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node.to.drop - length(phy$tip.label))]
  phy
}