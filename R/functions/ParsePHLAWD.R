## Last updated August 28, 2013
##With help from Matt Settles / Matt Pennell; Version of NEWfunction.R
##Modified May15,2013 by HEM to correct NCBI names, and add "_"

# Will have to install this package 
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")

library(Biostrings) # load package

#setwd("~/Dropbox/Hannah-Dave/SanJuans/HannahFINAL/2_SpeciesList/") # Navigate to the directory with PHLAWD output to be parsed

# This function will take the full alignment from the PHLAWD output and remove the NCBI ID, 
# and keep only the longest unique sequences if there are multiple hits for a single species
parsePHLAWD <- function(fasta.file){
  GBseqs <- readDNAStringSet(fasta.file) #read .aln.full
  namesGB <- names(GBseqs) #get the full NCBI names
  print(length(namesGB))
  split <- strsplit(namesGB, split="|", fixed=TRUE) #split names
  species.name <- sapply(split, "[", 2L) #get just the genus_species_var...
  genus.name <- sapply(strsplit(species.name, split="_"), "[", 1L) 
  species.name2 <- sapply(strsplit(species.name, split="_"), "[", 2L) 
  combinedname <- paste(genus.name, species.name2, sep="_") #get just genus_species
  sizes <- rowSums(alphabetFrequency(GBseqs)[,c("A","C","T","G")]) #get the nucleotide lenght of each sequence
  ord <- order(combinedname, -sizes)
  seqs <- GBseqs[ord] #order by lenght of sequence, longest first
  namesGBord <- names(seqs) #get the full NCBI names in correct order
  combinedname <- combinedname[ord]
  ID <- duplicated(combinedname) # identify duplicated combined names
  uniques <- seqs[!ID] #get only the unique sequences, choosing the longest since it is the first in the list
  uniquesnames <- combinedname[!ID]
  print(length(uniques))
  file.name <- strsplit(fasta.file, split=".", fixed=TRUE)[[1]][[1]]
  species_uniques <- uniques 
  names(species_uniques) <- uniquesnames
  writeXStringSet(species_uniques, file=paste(file.name, "unique", sep=".", format="fasta"))
  names(uniques) <- namesGBord[!ID] #full NCBI names
  writeXStringSet(uniques, file=paste(file.name, "unique.GB", sep=".", format="fasta"))
  #return(combinedname)
  #names(uniques) <- paste(genus.name[!ID], species.name2[!ID], sep="") # to get without space, eg for the Lamiales project b/c Nancy's seqs were like this, match for Mafft
 
}  

## To execute, run the above funtion, then call the file that you would like to parse. See the example atpB.FINAL.aln.full below:
#parsePHLAWD("atpB.FINAL.aln.full") 
## Output: *unique.fasta == the alignment trimed to just the longest sequences, i.e. the unique seqs
##        *unique.GB.fasta == same as the above, but with the ncbi info. and the species names

######## Fix names of files that were removed
#setwd("~/Documents/Idaho/Tank/Projects/SanJuans/FINAL/2b_Remove/")

parseREMOVED <- function(fasta.file){
  rem <- readDNAStringSet(fasta.file) #read .aln.full
  namesRem <- names(rem) #get the full NCBI names
  print(length(namesRem)) #62
  
  split <- strsplit(namesRem, split="|", fixed=TRUE) #split names
  
  species.name <- sapply(split, "[", 2L) #get just the genus_species_var...
  genus.name <- sapply(strsplit(species.name, split="_"), "[", 1L) 
  species.name2 <- sapply(strsplit(species.name, split="_"), "[", 2L) 
  combinedname <- paste(genus.name, species.name2, sep="_") #get just genus_species
  
  names(rem) <- combinedname
  file.name <- strsplit(fasta.file, split=".", fixed=TRUE)[[1]][[1]]
  writeXStringSet(rem, file=paste(file.name, "unique.rem.name", sep=".", format="fasta"))
}

#parseREMOVED("atpB.unique.GB.fasta.rem") #62

parseALIGNMENT <- function(fasta.file){
  GBseqs <- readDNAStringSet(fasta.file) #read .aln.full
  namesGB <- names(GBseqs) #get the full NCBI names
  print(length(namesGB))
  combinedname <- namesGB #get just genus_species
  file.name <- strsplit(fasta.file, split=".", fixed=TRUE)[[1]][[1]]
  write.csv(combinedname, file=paste(file.name, "species", ".csv", sep="."))
  return(combinedname)
  #names(uniques) <- paste(genus.name[!ID], species.name2[!ID], sep="") # to get without space, eg for the Lamiales project b/c Nancy's seqs were like this, match for Mafft
  
}

#fasta.file <- "~/Dropbox/Work/FranceLab/FranceProjects/IslandComparaive/WeigletWG/8Archepleagos/231014/4_Concatenate/align.concat.8arch.241014.fst"
#extractID <- read.csv("~/Dropbox/Work/FranceLab/FranceProjects/Islandcomparaive/WeigletWG/8Archepleagos/ExtractID/output241014.txt")

parseALIGNMENT.Input.to.acceptedName <- function(fasta.file, extractID, file.name){
  dim(extractID) #6477 = number of species used in include file == speices in island dataset

  ## Get just genus species for input names
  split <- strsplit(as.character(extractID$input_name), split=" ", fixed=TRUE) #split names
  species.name <- sapply(split, "[", 2L) #get just the genus_species_var...
  species.name2 <- strsplit(species.name, split="-", fixed=TRUE)
  species.name3 <- sapply(species.name2, "[", 1L) #get just the genus_species_var...
  genus.name <- sapply(split, "[", 1L) 
  combinedname.input <- paste(genus.name, species.name3, sep=" ") #get just genus_species
  
  
  extractID$input_name <- combinedname.input
  #write.csv(extractID.uniques$input_name, file="TEST.csv")
  
  extractID.uniques <- subset(extractID,!duplicated(extractID$input_name)) #remove duplicated input names
  dim(extractID.uniques) #4637    4
  
  GBseqs <- readDNAStringSet(fasta.file) #read concatenated alignmenzt
  print(length(GBseqs))  # 4383 = number of species in alignment 
  
  
  #print(dim(combinedname[(which(combinedname %in% extractID$input_name))])) # check to make sure all the names in alighment map to GenBank ID 
  #i = 4316
  matched.acceptedID.align <- DNAStringSet()
  matched.accepted.align <- DNAStringSet()
  matched.input.align <- DNAStringSet()
  for (i in 1:length(GBseqs)){ 
    
    namesGB <- names(GBseqs[i])
    
    splitnamesGB <- strsplit(as.character(namesGB), split="_", fixed=TRUE) #split names
    splitnamesGBspecies.name <- sapply(splitnamesGB, "[", 2L) #get just the genus_species_var...
    splitnamesGBgenus.name <- sapply(splitnamesGB, "[", 1L) 
    combinedname <- paste(splitnamesGBgenus.name, splitnamesGBspecies.name, sep=" ") #get just genus_species
    
    tmp <- extractID.uniques[extractID.uniques$input_name == combinedname,]
    if (nrow(tmp) == 0){
      print(paste(namesGB, "does not match PHLAWD includefile"))
      matched.acceptedID.align <- c(matched.acceptedID.align)
      matched.accepted.align <- c(matched.accepted.align)
    } else {
      
        acceptedID <- paste(as.character(tmp$accepted_name), tmp$ncbi_id, sep ="|" )
        seq.acceptedID <- GBseqs[i]
        names(seq.acceptedID) <- acceptedID
        matched.acceptedID.align <- c(matched.acceptedID.align, seq.acceptedID)
        
   
        accepted <- as.character(tmp$accepted_name)
        seq.accepted <- GBseqs[i]    
        names(seq.accepted) <- accepted
        matched.accepted.align <- c(matched.accepted.align, seq.accepted)    
        
     
        input <- as.character(tmp$input_name)
        seq.input <- GBseqs[i]
        names(seq.input) <- input
        matched.input.align <- c(matched.input.align, seq.input)    
    

    }

  }
  matched.acceptedID.align #4377
  matched.accepted.align #4377
  matched.input.align #4377
  
  writeXStringSet(matched.acceptedID.align, file=paste(file.name, "acceptedID", "fst", sep="."))
  writeXStringSet(matched.accepted.align, file=paste(file.name, "accepted", "fst", sep="."))
  writeXStringSet(matched.input.align, file=paste(file.name, "input", "fst", sep="."))
  
  
  #names(uniques) <- paste(genus.name[!ID], species.name2[!ID], sep="") # to get without space, eg for the Lamiales project b/c Nancy's seqs were like this, match for Mafft
  
}

#parseALIGNMENT.Input.to.acceptedName(fasta.file, extractID, file.name="align.concat.8arch.241014")

#[1] "Picris_sp does not match PHLAWD includefile"
#[1] "Chenopodium_glaucum does not match PHLAWD includefile"
#[1] "Polypodium_polypodioides does not match PHLAWD includefile"
#[1] "Hedyotis_corymbosa does not match PHLAWD includefile"
#[1] "Polygonum_chinense does not match PHLAWD includefile"
#[1] "Phelipanche_purpurea does not match PHLAWD includefile"












