source("beta_pd_decompo.r")

tree=read.tree("pone.0042760.s007.nwk")
plot(tree)

coms=read.csv("pone.0042760.s006.csv",h=T,row.names=1)
coms

samp = coms
decompo_beta <- beta.pd.decompo(com = coms,  tree=tree,type="both",output.dist=F, random=10)
decompo_beta
