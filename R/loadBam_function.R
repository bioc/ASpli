loadBAM <- function( targets, cores = 1 ) {
  
  #datac <- mclapply( as.character(targets$bam) , mc.cores = cores, readGAlignments )
  datac <- lapply( as.character(targets$bam), function(x){
      r <- readGAlignments(x)
      #Normalize seqnames. If . present in name, changes it to _ and warns the user
      if(length(grep("[.]", seqlevels(r)) > 0)){
        seqlevels(r) <- gsub("[.]", "_", seqlevels(r))
        warning("Some seqnames had a '.' present in their names. ASpli had to normalize them using '_'.")
      }
      gc()
      return(r)
    })
  names( datac ) <- rownames(targets)
  
  return(datac)
  
}
