#####################################################################################################################
############# Prep bootstrap distribution of trees ##############################################################
############# Hannah E. Marx, 16 Feb 2017 ###########################################################################
#####################################################################################################################

# load 1000 bootstrap trees
bstreesall <- read.tree(file = "output/5_Trees/Concatenated/RAxML_bootstrap.align.concat.EcrinSpPool.121115.phy.1000")
#write.tree(bstreesall[[1]], file="output/5_Trees/Concatenated/bootstrap100/bstreesall.1.tre")

### Actually don't want to root
# roottaxa <- c("Pinus_sylvestris", "Pinus_nigra", "Pinus_mugo", "Picea_abies", "Larix_decidua", "Abies_alba", "Cedrus_atlantica", "Juniperus_communis", "Juniperus_sabina", "Taxus_baccata")
# bstreesall.root.1 <- root(bstreesall[[1]], roottaxa, resolve.root = TRUE)
#write.tree(bstreesall.root.1, file="output/5_Trees/Concatenated/bootstrap100/bstreesall.root.1.tre")

# randomly draw 100
bstrees100 <- sample(bstreesall, size=100)
bstrees100[[1]]

tax.list=read.csv(file="data/Congruify/fleshed_genera.csv", as.is=TRUE, row=1) 
atol=read.tree("data/Congruify/out_dates.tre") #dataed reference tree, Soltis et al. 2011

for (i in 1:length(bstrees100)){
  
  ####### Congruify  ##########
  genetree=bstreesall[[i]]
  tips=sapply(genetree$tip.label, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips, rownames(tax.list))
  SJ_tax=tax.list[ll,]
  rownames(SJ_tax)=names(tips)
  SJ_tax=as.matrix(SJ_tax)
  SJ_tax[is.na(SJ_tax)]=""
  ftax=tax.list[match(atol$tip.label, rownames(tax.list)),]
  ftax[,2]="Spermatophyta"
  fatol=subset(atol, ftax, "family")
  phy=genetree
  swaptips=paste(1:length(tips),tips,sep="_")
  phy$tip.label=swaptips
  tax=SJ_tax
  rownames(tax)=swaptips
  
  ####### Get scaling set up  ##########
  res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 
  nsites=5144 #SHOULD be number of nucleotides in alignment 
  write.treePL(res1$target, res1$calibrations, nsites, base=paste("/Users/hannahmarx/Dropbox/Work/FranceLab/FranceProjects/AlpinePD_France/Alpine-Sky-Islands/output/5_Trees/Concatenated/bootstrap100/scaling/scale", i, sep="."), 
               opts=list(opt = 5, optad = 5, optcvad = 1))
  
}

####### Run scaling in shell: /Users/hannahmarx/Dropbox/Work/FranceLab/FranceProjects/AlpinePD_France/Alpine-Sky-Islands/output/5_Trees/Concatenated/bootstrap100/scaling/
  ## for i in *.infile; do treePL $i; done

####### Read in treePL output to change tip labels back
path = "/Users/hannahmarx/Dropbox/Work/FranceLab/FranceProjects/AlpinePD_France/Alpine-Sky-Islands/output/5_Trees/Concatenated/bootstrap100/scaling/"
file.names <- dir(path, pattern =".dated.tre")
for (i in 1:100){
  ####### Congruify  ##########
  genetree=bstreesall[[i]]
  tips=sapply(genetree$tip.label, function(x){
    unlist(strsplit(x,"_",fixed=TRUE))[1]
  })
  ll=match(tips, rownames(tax.list))
  SJ_tax=tax.list[ll,]
  rownames(SJ_tax)=names(tips)
  SJ_tax=as.matrix(SJ_tax)
  SJ_tax[is.na(SJ_tax)]=""
  ftax=tax.list[match(atol$tip.label, rownames(tax.list)),]
  ftax[,2]="Spermatophyta"
  fatol=subset(atol, ftax, "family")
  phy=genetree
  swaptips=paste(1:length(tips),tips,sep="_")
  phy$tip.label=swaptips
  tax=SJ_tax
  rownames(tax)=swaptips
  
  ####### Get scaling object  ##########
  res1=congruify.phylo(fatol, phy, tax, tol=0, scale="NA") 
  
  file <- read.tree(paste(path, file.names[i], sep="")) ### Read in treePL output
  file$tip.label=genetree$tip.label[match(swaptips, res1$target$tip.label)] ### CHANGE tip labels back
  write.tree(treePL, file=paste("output/5_Trees/Concatenated/bootstrap100/scale", i, "dated.rename.tre", sep=".")) #Write .tre file
  
}

