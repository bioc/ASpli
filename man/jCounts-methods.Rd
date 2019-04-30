\name{jCounts}
\alias{jCounts}
\alias{AsDiscover}

\title{ Report PIN, PJU and PIR using experimental junctions }

\description{ Given a bin, it is possible to calculate PIN/PIR metric using 
junctions to estimate changes in the use of it along different conditions. 
PSI or PIR metric is calculated for each bin and experimental condition. The 
selection of which metric is used is based on the kind of splicing event 
associated with each bin. 
}
 
\usage{ 
 jCounts( counts, features, bam = NULL, readLength, threshold , cores )
}
\arguments{
\item{counts}{
An object of class ASpliCounts.
}

\item{features}{ An object of class ASpliFeatures.}

\item{bam}{ A list with BAM files contents.}

\item{readLength}{ Read length of sequenced read. Default 100L }

\item{threshold}{ Minimum  number of reads supporting junctions. Default=5 }

\item{cores}{ Number of processing cores to use }
}

\value{ 
An object of class ASpliAS 
\item{irPIR}{
  reports: event, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by 
  condition. J1, J2, J3 sum of junctions (J1, J2, J3) by condition. 
}

\item{altPIN}{ 
  reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of 
  junctions (J1, J2, J3) by condition. 
}

\item{esPIN}{ reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2,
  J3 sum of junctions (J1, J2, J3) by condition. 
}

\item{junctionsPIR}{
  PIR metric for each experimental junction using e1i and ie2 counts. Exclusion 
  junction is the junction itself. This output helps to discover new introns as 
  well as new retention events.
}

\item{junctionsPJU}{
  Given a junction, it is possible to analyze if it shares start, end or both 
  with another junction. If so, is because there is more than one way for/of 
  splicing. Ratio between them along samples is reported.
}

}
\author{
  Estefania Mancini, Marcelo Yanovsky and Ariel Chernomoretz 
}

\examples{

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  #library(GenomicFeatures)
  #genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
  #                               package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  #features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  #bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
  #                   "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  #targets <- data.frame( 
  #             row.names = paste0('Sample_',c(1:6)),
  #             bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
  #             factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  #bams <- loadBAM( targets )
  
  # Read counts from bam files
  #counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
  #                        maxISize = 50000 )
  
  # Calculate differential usage of genes, bins and junctions 
  #du       <- DUreport.norm( counts, targets )

  # Calculate PSI / PIR for bins and junction.
  #as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
  #                        threshold = 5, cores = 1 )
  
  #writeAS( as = as, output.dir = "only_as" )

}

\seealso{

Accesors: \code{\link{irPIR}}, \code{\link{altPIN}}, \code{\link{esPIN}}, 
\code{\link{junctionsPIR}}, \code{\link{junctionsPJU} }

Export: \code{ \link{writeAS} }

}
