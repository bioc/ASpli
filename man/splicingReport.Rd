\name{splicingReport}

\alias{splicingReport}

\title{
  Splicing report
}
\description{
  This function integrates differential usage information from different sources.
  Receives bin usage and junction usage and merges them in three steps.
  First merges all bins with annotated junction information (jir, jes, jalt).
  Then finds all bins that overlap locale J3 junctions (the junctions that cover the entire locale region) and merges them.
  Finally it does the same with anchor.
}

\usage{
  splicingReport(bdu, jdu, counts)
}

\arguments{
  \item{bdu}{ 
    An object of class \code{ASpliDU}
  }
  \item{jdu}{ 
    An object of class \code{ASpliJDU}
  }
  \item{counts}{ 
    An object of class \code{ASpliCounts}
  }
}

\value{ 
  An \code{ASpliSplicingReport} object with junction differential usage report.
  \code{binbased(ASpliSplicingReport)} is a data.frame with merged information from bins and annotated junctions. Junction event statistics are located in columns prefixed with cluster. If intron retention event, cluster.dPIR has information, else, cluster.dPSI. It also includes counts for all the junctions for debugging purposes.
  \code{localebased(ASpliSplicingReport)} is a data.frame with merged information from bins and locale junctions. Statistics for each junction are located in columns prefixed with "junction" and statistics for the entire cluster are prefixed with "cluster". Inclusion statistics are called dparticipation for contrasts. If locale matched a bin, then bin statistics are also included.
  \code{anchorbased(ASpliSplicingReport)} is a data.frame similar to  \code{localebased(ASpliSplicingReport)}.
   All data.frames also include junction counts for debugging purposes.
}

\author{
  Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
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
               factor1 = c( 'C','C','C','D','D','D'),
               subject = c(0, 1, 2, 0, 1, 2))
  
  # Read counts from bam files
  gbcounts  <- readCounts( features = features, targets = targets, minReadLength = 100, maxISize = 50000 )
  jcounts   <- jCounts(counts = gbcounts, features = features, minReadLength = 100)
  
  # Test for factor1 controlling for paired subject
  gbPaired <- gbDUreport(gbcounts, formula = formula(~subject+factor1))
  jPaired  <- jDUreport(jcounts, formula = formula(~subject+factor1))
  
  # Generate a splicing report merging bins and junctions DU
  sr              <- splicingReport(gbPaired, jPaired, gbcounts)

  # Access splicing report elements
  sr  
  localebased(sr)
  anchorbased(sr)
  binbased(sr)
}

\seealso{
    \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}
}


