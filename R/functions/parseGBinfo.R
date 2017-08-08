


#### This function will take the unique sequences from parsePHLAWD output and retrieve the species name and NCBI ID
#unique.GB.fasta <- ("output/2_Remove/atpB.unique.GB.rem2.fasta") #119...118
#gene <- "atpB"

parseGBinfo <- function(gene, unique.GB.fasta){
  GBseqs <- readDNAStringSet(unique.GB.fasta) #read .aln.full
  namesGB <- names(GBseqs) #get the full NCBI names
  print(length(namesGB))
  split <- strsplit(namesGB, split="|", fixed=TRUE) #split names
  species.name <- sapply(split, "[", 2L) #get just the genus_species_var...
  ncbi <- sapply(split, "[", 4L) #get NCBI identifier
  combinedname_ncbi.tmp <- as.list(ncbi)
  names(combinedname_ncbi.tmp) <- species.name
  combinedname_ncbi.df <- cbind(combinedname_ncbi.tmp)
  colnames(combinedname_ncbi.df) <- gene
  return(combinedname_ncbi.df)
} 


#### This function will merge the species names to record which NCBI ID was used for each gene region
#atpbGB <- parseGBinfo(gene = "atpB", unique.GB.fasta = "output/2_Remove/atpB.unique.GB.rem2.fasta")
#rbclGB <- parseGBinfo(gene = "rbcL", unique.GB.fasta = "output/2_Remove/rbcL.unique.GB.fasta.rem")
#tnrtlfGB <- parseGBinfo(gene = "trnTLF", unique.GB.fasta = "output/2_Remove/trnTLF.unique.GB.fasta.rem")

mergeGBinfo <- function(...){
  input_list <- list(...)
  joined <- join_all(lapply(input_list, function(x) data.frame(x, rn = row.names(x))), by = 'rn', type = 'full')
  rownames(joined) <- joined$rn
  joined <- joined[ , !(colnames(joined) == "rn")]
  joined_sort <- joined[order(row.names(joined)), ]
  joined_sort[joined_sort == "NULL"] <- "NA"
  dd  <-  as.data.frame(matrix(unlist(joined_sort), nrow=length(unlist(joined_sort[1]))), row.names = rownames(joined_sort))
  split <- strsplit(rownames(dd), split="_", fixed=TRUE) #split names
  species <- sapply(split, "[", 2L) #get just the genus_species_var...
  genus <- sapply(split, "[", 1L) #get just the genus_species_var...
  combinedname <- paste(genus, species, sep="_") #get just genus_species
  dd.df <- cbind(combinedname, dd)
  colnames(dd.df) <- c("species", colnames(joined_sort))
  head(dd.df)
  return(dd.df)
}

#mergeGBinfo(atpbGB, rbclGB, tnrtlfGB)


