
#load('R/ContingencyTests/trends')
#load('R/ContingencyTests/groups')

underIce.tmp <- read.csv('data/sp_pool_under_ice.csv')
dim(underIce.tmp) #165
persistant.tmp <- read.csv('data/sp_pool_persistant.csv')
dim(persistant.tmp) #261

LGM.tmp <- merge(persistant.tmp, underIce.tmp,  by=1, all=T) # number of releves with each species above or below LGM layer limit
LGM <- as.data.frame(LGM.tmp, , strings.as.factors=F)
rownames(LGM) <- LGM.tmp[,1]
LGM[is.na(LGM)] <- 0
colnames(LGM) <- c("species", "persistant", "underIce")
LGM <- LGM[sort(rownames(LGM)),]
dim(LGM) #298
head(LGM)

# make the contingency table
# you wanna equate O/1 here to below/above the LGM glaciar for instance
#x <- table(groups, trends)
x <- LGM[,c(2:3)]

# compute expected frequencies under randomness assumption
exp <- round(rowSums(x)*(sum(x[,2])/sum(x)), digits=2) 

# compute the difference between the expected and the observed frequencies
# or the other way around, this is the statistics use by a classic chi-square test
diff.exp <- x[,2]-exp

# create a vector of p.values (empty so far)
p.values <- numeric(dim(x)[1])

# compute a fisher exact test for each group (each line against the rest of the table)
# this yields the probability of observing the odd-ratio pattern under the 
# assumption of random dispersion in the data #A / B freq above and below glacier 
for (i in 1:length(levels(groups))){
  newgroups <- groups
  levels(newgroups)[!(levels(newgroups)==levels(newgroups)[i])] <- "Others"
  x2 <- table(newgroups, trends)
  p.values[i] <- round(fisher.test(x2, simulate.p.value=FALSE)$p.value, digits=4)
}
results.tmp <- as.data.frame(cbind(x, exp, diff.exp, p.values))
dim(results.tmp) #298
results <- results.tmp[results.tmp$p.values<0.1, ] #print out results only interesting results
dim(results) #242

head(results)

# print out results to make you pretty result table
#write.csv(results, file='results.conting.RESCUE.csv')

resultsPersi <- (results[results$persistant > results$underIce,])
write.csv(resultsPersi, file='results.conting.persis.csv')


