\name{readCounts}

\alias{readCounts}

\title{
  Summarize read overlaps
}

\description{ Summarize read overlaps against all feature levels }

\usage{
  readCounts( features, bam, targets, cores = 1, readLength,
                 maxISize, minAnchor = NULL)
}

\arguments{
  
  \item{features}{ An object of class ASpliFeatures. It is a list of GRanges 
    at gene, bin and junction level
  }
  
  \item{bam}{ List of GAlignments objects corresponding to bam files of 
    samples. }
  
  \item{targets}{ A dataframe containing sample, bam and experimental factors 
    columns }
  
  \item{readLength}{ Read length of sequenced library. It is used for compute 
    E1I and IE2 read summarization }
    
  \item{maxISize}{maximum intron expected size. Junctions longer than this size
    will be dicarded }
  
  \item{cores}{ Number of cores to use. Default is 1 }
  
  \item{minAnchor}{ Minimum percentage of read that should be aligned to 
    an exon-intron boundary  }
}

\value{
  An object of class ASpliCounts. Each slot is a dataframe containing features 
  metadata  and read counts. Summarization is reported at gene, bin, junction 
  and intron flanking regions (E1I, IE2).
}
\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz 
}

\seealso{
  Accesors: \code{\link{countsg}}, 
            \code{\link{countsb}},
            \code{\link{countsj}},
            \code{\link{countse1i}},
            \code{\link{countsie2}},
            \code{\link{rdsg}},
            \code{\link{rdsb}},
  Export:   \code{\link{writeCounts}}
}

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
  # Export data
  writeCounts( counts, output.dir = "only_counts" )
  
}

