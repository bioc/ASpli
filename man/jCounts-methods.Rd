\name{jCounts}
\alias{jCounts}
\alias{AsDiscover}

\title{ Report PSI, PJU and PIR using experimental junctions }

\description{ Given a bin, it is possible to calculate PSI/PIR metrics using 
              junctions to estimate changes in the bin usage along different conditions. 
              PSI or PIR metrics are calculated for each bin and experimental condition. The 
              selection of which metric is used is based on the kind of splicing event 
              associated with each bin. 
}
 
\usage{ 
 jCounts( counts, features, bam = NULL, minReadLength, threshold
                 = 5, cores = 1, minAnchor = 10 )
}
\arguments{
\item{counts}{
An object of class ASpliCounts.
}

\item{features}{ An object of class ASpliFeatures.}

\item{bam}{ A list with BAM files contents or NULL to let the function read bam files. If NULL is passed (default), bam files are read and processed sequentially. This is usually preferred if machine has low resources.}

\item{minReadLength}{ Minimum read length of sequenced library. It is used for computing 
    E1I and IE2 read summarization. Make sure this number is smaller than the maximum read length in every bam file, otherwise no junctions will be found.}

\item{threshold}{ Minimum  number of reads supporting junctions. Default=5 }

\item{cores}{ Number of processing cores to use. If bam is NULL, only 1 core is used independently of this parameter.  }

\item{minAnchor}{ An intronic junction must overlap completely and at least an minAnchor\% into the exon region and the intron region.
                  The regions can be exon1-intron or intron-exon2.}
}

\value{ 
An object of class ASpliAS. 
\item{irPIR}{
  Reports event, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by 
  condition. J1, J2, J3 sum of junctions (J1, J2, J3) by condition. 
}

\item{altPSI}{ 
  Reports event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum of 
  junctions (J1, J2, J3) by condition. 
}

\item{esPSI}{ 
  Reports event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2,
  J3 sum of junctions (J1, J2, J3) by condition. 
}

\item{junctionsPIR}{
  PIR metric for each experimental junction using e1i and ie2 counts. Exclusion 
  junction is the junction itself. This output helps to discover new introns as 
  well as new retention events.
  Reports hitIntron, hitIntronEvent, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by 
  condition.
}

\item{junctionsPJU}{
  Given a junction, it is possible to analyze if it shares start, end or both 
  with another junction. If so, it is because there is alternative splicing.
  Reports for every junction, junction information (junction, gene, strand, multipleHit, symbol, gene_coordinates, bin_spanned, j_within_bin), StartHit with all the junctions sharing the start with this junction, $PJU_{J1}=J3/(J1+ J3)$ for each condition, EndHit with all the junctions sharing the end with this junction and $PJU_{J2}=J3/(J2+ J3)$ for each condition. 
}

}
\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky and Ariel Chernomoretz 
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

  # Calculate PSI / PIR for bins and junction.
  junctionCounts <- jCounts(counts = counts, features = features, minReadLength = 100)

  # Access summary and bin and junction counts and metrics and display them
  junctionCounts
  irPIR(junctionCounts)
  junctionsPJU(junctionCounts)
  
  # Export data
  writeAS( as = junctionCounts, output.dir = "only_as" )

}

\seealso{

Accesors: \code{\link{irPIR}}, \code{\link{altPSI}}, \code{\link{esPSI}}, 
\code{\link{junctionsPIR}}, \code{\link{junctionsPJU} }

Export: \code{ \link{writeAS} }

}
