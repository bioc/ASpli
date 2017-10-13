loadBAM <- function( targets, cores = 1 ) {
  
  datac <- mclapply( as.character(targets$bam) , mc.cores = cores, 
      readGAlignments )
  names( datac ) = rownames(targets)
  
  return(datac)
  
}
