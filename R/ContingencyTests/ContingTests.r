
load('trends')
load('groups')

# make the contingency table
# you wanna equate O/1 here to below/above the LGM glaciar for instance
x <- table(groups, trends)

# compute expected frequencies under randomness assumption
exp <- round(rowSums(x)*(sum(x[,2])/sum(x)), digits=2) 

# compute the difference between the expected and the observed frequencies
# or the other way around, this is the statistics use by a classic chi-square test
diff.exp <- x[,2]-exp

# create a vector of p.values (empty so far)
p.values <- numeric(dim(x)[1])

# compute a fisher exact test for each group (each line against the rest of the table)
# this yields the probability of observing the odd-ratio pattern under the 
# assumption of random dispersion in the data
for (i in 1:length(levels(groups))){
  newgroups <- groups
  levels(newgroups)[!(levels(newgroups)==levels(newgroups)[i])] <- "Others"
  x2 <- table(newgroups, trends)
  p.values[i] <- round(fisher.test(x2, simulate.p.value=FALSE)$p.value, digits=4)
}
results <- as.data.frame(cbind(x, exp, diff.exp, p.values))

results <- results[results$p.values<0.1, ] #print out results only interesting results
results

# print out results to make you pretty result table
write.csv(results, file='results.conting.RESCUE.csv')

