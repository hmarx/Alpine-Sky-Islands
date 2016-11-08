
######################################################################################################################
#### Hannah E. Marx
#### 21 July 2015
######################################################################################################################
########## Randomize occurence matrix of communites sampling from differnt source pools with defined probability
########## Calculate ses.mpd, ses.mntd

## Random: presence of a species in community | prob(drawn from source pool) 
##       : community matrix (species richness and occupancy rate (i.e. local community richness))
## Constant: distance matrix (phylogeny, trait, env, etc.) 

### resample only speceis in a defined pool, reshuffling = probaility of current occurence 
###### 1) Subset species in pool
###### 2) calulate resampling prob for species
###### 3) shuffle tips 
######################################################################################################################
############### Generate a list of random communities  
# drawProbs = list of resampling weights for each species
# dist = community phylogeny 
# com = community matrix, rownames = communtities (observed), colnames = species
# N = the number of null communities to generate

###### simulated phylogeny
#sim_tre <- read.tree(text="((((A:0.3333333333,(B:0.2222222222,(C:0.1111111111,D:0.1111111111):0.1111111111):0.1111111111):0.2222222222,(E:0.1111111111,F:0.1111111111):0.4444444444):0.1111111111,G:0.6666666667):0.3333333333,((H:0.1111111111,I:0.1111111111):0.1111111111,J:0.2222222222):0.7777777778);")
                          
###### simulated community matrix with abundances
#com.data.a <- matrix(c(3,0,1,0,5,2,1,2,0,0,4,3,1, 2,0,0,0,5,4,4,1,0,0,0,1,0,1), nrow=3)
#row.names(com.data.a) <- c("com1", "com2", "com3")
#colnames(com.data.a) <- c("A","B","C","D","E","F","G","H","I")

######  object to call in function 
#dist = cophenetic(sim_tre) # phylo object
#com = com.data.a # community data.frame (rownames = community)
#sourcePool = "com1" #the pool that you want to draw from to generate the null distribution; one row of the communtiy dataframe

#######  Test the community matrix randomization  
#randomizeSourcePool(dist=sim_tre, com = com.data.a, sourcePool = "com1", N =5)
#randomizeSourcePool(dist=sim_tre, com = com.data.a, sourcePool = "com2", N =5)

#######  Compare to picante output
#ses.mpd.sourcePool(dist=sim_tre, com = com.data.a, sourcePool = "com1", N =100)
#ses.mntd.sourcePool(dist=sim_tre, com = com.data.a, sourcePool = "com1", N =100)
#ses.mpd(samp = com.data.a, dis=cophenetic(sim_tre), null.model = "sample.pool", runs =100)

randomizeSourcePool <- function(dist, com, sourcePool, N){
  # Get the observed total occurence counts acrross all summits for each species
  orig.com <- t(com)
  orig.com.counts <- rowSums(orig.com)
  
  # Occurences in defined species pool
  sourcePoolCounts <- com[sourcePool,]
  
  # Proability to reshuffle each species = occurence in source pool / total occurence acrross communities
  samp.prob <- (sourcePoolCounts/orig.com.counts) #sum to one for each species across all communities
  samp.prob[is.na(samp.prob)] <- 0 #just in case there is missing data
  #print(samp.prob)
  
  # remove source pool from community to randomize
  focal.com <- as.data.frame(orig.com)[, names(as.data.frame(orig.com)) != sourcePool]
  
  community.matrix.list <- list() # create a list of randomized communities | sampling probability
  ### Randomization  
    rep <- 0
    while (rep < N){
      rep <- rep + 1    
      sim.comm.df <- data.frame()
      for (i in 1:nrow(focal.com)){ # for each species in community with source pool removed 
        sp.prob <- samp.prob[rownames(focal.com)[i]] # the sample probability for one species
        # draw occurence in community | probability of sampling in source pool
        sim.sp.occ <- (sample(c(0,1), size=length(focal.com[i,]), replace=T, 
                              prob = c(1-sp.prob, sp.prob)))
        #print(sim.sp.occ)
        sim.comm.df <- rbind(sim.comm.df, sim.sp.occ)
      }
      sim.comm.df
      colnames(sim.comm.df) = names(focal.com) 
      rownames(sim.comm.df) = rownames(focal.com) 
      
      community.matrix.list[[rep]] <- t(sim.comm.df)
  }
  return(community.matrix.list)
}

############### Calulate ses.mpd for source pool randomization

ses.mpd.sourcePool <- function(dist, com, sourcePool, N){

  mpd.obs.all <-  mpd(samp = com, dis=dist)
  names(mpd.obs.all) <- rownames(com)
  #rowSums(as.data.frame(!com==0))
  
  rand.com.tmp <- randomizeSourcePool(dist =dist, com = com, sourcePool = sourcePool, N =N)
  #colSums(as.data.frame((rand.com.tmp)[1]))
  #rowSums(as.data.frame((rand.com.tmp)[1]))
  null.mpd <- lapply(1:length(rand.com.tmp), function(x) mpd(rand.com.tmp[[x]], dist))
  null.mpd.c <- as.data.frame(do.call(rbind, null.mpd))
  colnames(null.mpd.c) <- rownames(as.data.frame((rand.com.tmp)[1]))
  
  output.df <- data.frame()
  ## For each community:
  for (i in 1:ncol(null.mpd.c)){
    ntaxa <- colSums(!t(com)==0)[i] #number of species
    mpd.obs <- mpd.obs.all[names(null.mpd.c[i])] #observed mean metric
    mpd.rand.mean <- mean(null.mpd.c[,i], na.rm = T) # mean of metrics from randomization
    mpd.rand.sd <- sd(null.mpd.c[,i], na.rm = T) # standard deviation of metrics from randomization
    mpd.obs.rank <- rank(c(mpd.obs, null.mpd.c[,i]))[1] # rank of observed vs. null communities
    mpd.obs.z <- (mpd.obs - mpd.rand.mean) / mpd.rand.sd #Standardized effect size of mpd vs. null communities (equivalent to -NRI)
    runs <- N #Number of randomizations
    mpd.obs.p <- mpd.obs.rank / (runs+1) #P-value (quantile) of observed mpd vs. null communities (= mpd.obs.rank / runs + 1)
    output.com <- cbind(ntaxa, mpd.obs, mpd.rand.mean, mpd.rand.sd, mpd.obs.rank, mpd.obs.z, mpd.obs.p, runs) 
    output.df <- rbind(output.df, output.com)
    
  }
  return(output.df)
  
}


############### Calulate ses.mntd for source pool randomization

ses.mntd.sourcePool <- function(dist, com, sourcePool, N){
  mntd.obs.all <-  mntd(samp = com, dis=dist)
  names(mntd.obs.all) <- rownames(com)
  
  rand.com.tmp <- randomizeSourcePool(dist =dist, com = com, sourcePool = sourcePool, N =N)
  null.mntd <- lapply(1:length(rand.com.tmp), function(x) mntd(rand.com.tmp[[x]], dist))
  null.mntd.c <- as.data.frame(do.call(rbind, null.mntd))
  colnames(null.mntd.c) <- rownames(as.data.frame((rand.com.tmp)[1]))
  
  output.df <- data.frame()
  ## For each community:
  for (i in 1:ncol(null.mntd.c)){
    
    ntaxa <- colSums(!t(com)==0)[i] #number of species
    mntd.obs <- mntd.obs.all[names(null.mntd.c[i])] #observed mean metric
    mntd.rand.mean <- mean(null.mntd.c[,i], na.rm = T) # mean of metrics from randomization
    mntd.rand.sd <- sd(null.mntd.c[,i], na.rm = T) # standard deviation of metrics from randomization
    mntd.obs.rank <- rank(c(mntd.obs, null.mntd.c[,i]))[1] # rank of observed vs. null communities
    mntd.obs.z <- (mntd.obs - mntd.rand.mean) / mntd.rand.sd #Standardized effect size of mntd vs. null communities (equivalent to -NRI)
    runs <- N #Number of randomizations
    mntd.obs.p <- mntd.obs.rank / (runs+1) #P-value (quantile) of observed mntd vs. null communities (= mntd.obs.rank / runs + 1)
    output.com <- cbind(ntaxa, mntd.obs, mntd.rand.mean, mntd.rand.sd, mntd.obs.rank, mntd.obs.z, mntd.obs.p, runs) 
    output.df <- rbind(output.df, output.com)
    
  }
  return(output.df)
  
}



