########## TAXONOMY ORGANIZATION ###########

# Maintain consistent taxon names across species lists and various datasets
# Rownames = species names 


library(taxize)
# Deal with French properly
options(encoding="latin1")

### Load dataset
#head(releves) 
#dim(releves)
#colnames(releves)


get.genus.species <- function(df, colnum, spliton, sepas){
  split <- strsplit(as.character(df[,colnum]), split=spliton, fixed=TRUE)
  genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
  species.name <- sapply(split, "[", 2L) #get just the species_var...
  combinedname <- paste(genus.name, species.name, sep=sepas) #get just genus_species
  return(combinedname) #30839
}



#df = the dataframe you would like to fix, species in rows, but names to change must also be in a column that is not 0
#colnum = the number of the column in the dataframe with the taxonomy you want to fix
#spliton = the punctuation separating names in taxonomy, e.g. "_" 
#sepas = the separation used for the format in the taxonomic lookup, e.g. " "
#source = source to match names against: iPlant_TNRS, NCBI, MSW3

### Get Just Genus_species...(remove infraspecific identifiers) 
add.taxized.column <- function(df, colnum, spliton, sepas, source){
  split <- strsplit(as.character(df[,colnum]), split=spliton, fixed=TRUE)
  genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
  species.name <- sapply(split, "[", 2L) #get just the species_var...
  ### Remove punctuation 
  combinedname <- paste(genus.name, species.name, sep=sepas) #get just genus_species
  combinedname <- gsub(combinedname, pattern = "\\.", replacement = sepas) 
  combinedname <- gsub(combinedname, pattern = "-", replacement = sepas)
  combinedname <- gsub(combinedname, pattern = " L ,$", replacement = "") 
  combinedname <- gsub(combinedname, pattern = " x$", replacement = "") 
  combinedname <- gsub(combinedname, pattern = " F.H.$", replacement = "") 
  combinedname <- gsub(combinedname, pattern = ",$", replacement = "") 
  df.nex <- as.data.frame(cbind(combinedname, df))
  df.nex[] <- lapply(df.nex, as.character)
  ### Make sure taxonomic names are spelled correctly, and are up to date
  ## Make a taxonomy lookup list using iPlant TNRS database
  tmp <-  tnrs(query = df.nex[!duplicated(df.nex[1]),1], source = source)[ , -c(5:7)] #accepted name, blank = no opinion
  #head(tmp)
  tmp.rm <- na.omit(tmp[tmp[2]!= "",]) #just those with accepted name not blank
  #write.csv(tmp, file="Sp.List.tmp")
  for (i in 1:nrow(df.nex)){
    #i=16# Euphorbia dulcis
    if (as.character(df.nex[i, 1]) %in% tmp.rm[[1]]){ #if the sp is in the taxonomy lookup list
      if (tmp.rm[tmp.rm[1] == as.character(df.nex[i,1]),1] == tmp.rm[tmp.rm[1] == as.character(df.nex[i,1]),2]){
        next
      } else {
        df.nex[i,1] <- tmp.rm[tmp.rm[1] == df.nex[i,1],2]
      }
    }
  }
  
  taxized <- get.genus.species(df = df.nex, colnum = 1, spliton = " ", sepas = " ")
  
  df.nex.tmp2 <- cbind(taxized, df.nex)
  #head(df.nex)
  #df.nex[16,]
  #dim(df.nex)
  #write.csv(df.nex, file="Sp.List.tmp2.csv")
  return(df.nex.tmp2) #30839
}




