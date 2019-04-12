loadBAM <- function( targets, cores = 1 ) {
  
  #datac <- mclapply( as.character(targets$bam) , mc.cores = cores, readGAlignments )
  datac <- lapply( as.character(targets$bam), function(x){
      r <- readGAlignments(x)
      gc()
      return(r)
    })
  names( datac ) <- rownames(targets)
  
  return(datac)
  
}
