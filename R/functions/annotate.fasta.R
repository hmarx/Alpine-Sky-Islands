# Annotates .fasta alignment to updated nomenclature

#fastaFile= "output/5_Trees/Concatenated/align.concat.EcrinSpPool.121115.fasta"
#taxonomy = v6tax_edit[c(2,3)]


annotate.fasta <- function(fastaFile, taxonomy){
  GBseqs <- readDNAStringSet(fastaFile) #read .aln.rn
  ## Make a tmp file so don't overwrite 
  tmp.fasta <- DNAStringSet(GBseqs)
  fasta.names <- names(tmp.fasta)
  
  for (i in 1:length(tmp.fasta)){
    #print(i)
    if (names(tmp.fasta[i]) %in% taxonomy[,2]){
    fasta.names[i] = as.character(taxonomy[which(names(tmp.fasta[i])==taxonomy[,2]), 1])
        
    }
    else {
      fasta.names[i] = fasta.names[i]
       
    }
      
  }
  names(tmp.fasta) <- fasta.names
  file.name <- strsplit(fastaFile, split=".", fixed=TRUE)[[1]][[1]]
  writeXStringSet(tmp.fasta, file=paste(file.name, "unique.rem.fasta.clean.annot", sep=".", format="fasta"))
  return(tmp.fasta)
}

