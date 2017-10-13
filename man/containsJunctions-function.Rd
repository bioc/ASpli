\name{ Examine ASpliDU objects }

\alias{containsJunctions}
\alias{containsGenesAndBins}

\title{ Examine ASpliDU objects }

\description{
  AspliDU object may contain results of differential expression of genes, 
  differential usage of bins and junctions, however not everything is 
  calculated at the same or even present. Calculations for genes and bins can 
  be done independently from junctions. Functions \code{containsJunctions} and
  \code{containsGenesAndBins} allow to interrogate an ASpliDU object about the 
  kind of results it contain. 
}

\usage{
  containsJunctions( du )
  containsGenesAndBins( du )
}

\arguments{
  \item{ du }{ 
    An ASpliDU object.
  }
}
\value{ 
  A logical value that indicates that results for genes and bins, or results for
  junctions are available in the object. }

\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Create an ASpliDU object.
  # 
  du       <- DUreport( counts, targets )
  
  # Verify if du contains results for genes and bins.
  containsGenesAndBins( du )
  containsJunctions( du )

  # Append junction results  
  du       <- junctionDUreport( counts, targets, appendTo = du )

  # Verify if du contains results for genes and bins.
  containsJunctions( du )
  
  
}

