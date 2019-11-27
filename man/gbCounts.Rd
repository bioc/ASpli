\name{gbCounts}

\alias{gbCounts}
\alias{readCounts}

\title{
  Summarize read overlaps
}

\description{ Summarize read overlaps against all feature levels }

\usage{
  gbCounts( features, bam = NULL, targets, cores = 1, minReadLength,
                 maxISize, minAnchor = 10)
}

\arguments{
  
  \item{features}{ An object of class ASpliFeatures. It is a list of GRanges 
    at gene, bin and junction level
  }
  
  \item{bam}{ List of GAlignments objects corresponding to bam files of 
    samples. If NULL is passed (default), bam files are read and processed sequentially. This is usually preferred if machine has low resources. }
  
  \item{targets}{ A dataframe containing sample, bam and experimental factors 
    columns }
  
  \item{minReadLength}{ Minimum read length of sequenced library. It is used for computing 
    E1I and IE2 read summarization. Make sure this number is smaller than the maximum read length in every bam file, otherwise no E1I or IE2 will be found.}
    
  \item{maxISize}{Maximum intron expected size. Junctions longer than this size
    will be dicarded }
  
  \item{cores}{ Number of cores to use. Default is 1. If bam is NULL, only 1 core is used independently of this parameter. }
  
  \item{minAnchor}{ Minimum percentage of read that should be aligned to 
    an exon-intron boundary.}
}

\value{
  An object of class ASpliCounts. Each slot is a dataframe containing features 
  metadata and read counts. Summarization is reported at gene, bin, junction 
  and intron flanking regions (E1I, IE2).

}
\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz 
}

\seealso{
  Accesors: \code{\link{countsg}}, 
            \code{\link{countsb}},
            \code{\link{countsj}},
            \code{\link{countse1i}},
            \code{\link{countsie2}},
            \code{\link{rdsg}},
            \code{\link{rdsb}},
            \code{\link{condition.order}},
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
  
  # Read counts from bam files
  counts   <- readCounts( features = features, targets = targets, minReadLength = 100, maxISize = 50000 )
  
  # Access summary and gene and bin counts and display them
  counts
  countsg(counts)
  countsb(counts)

  # Export data
  writeCounts( counts, output.dir = "only_counts" )
  
}

