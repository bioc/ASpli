\name{plotTopTags}
\alias{plotTopTags}

\title{
Coverage plots 
}
\description{
Using genomic coordinates and BAM files this function is useful for make coverage plots
}
\usage{
plotTopTags(auxdf, genome, targetsPlot, output.dir)
}

\arguments{
\item{auxdf}{ A data frame: row.naMes=bin names, gene coordinates, bin  coordinates and  event name columns
}
\item{genome}{
TxDb genome
}
\item{targetsPlot}{
A dataframe containting: bam files name, condition (y axis tag), color for each condition}
\item{output.dir}{ Name of directory where plots are supossed to be exported }
}

\value{
  Coverage plots in Portable Network Graphics (PNG) format of selected events.
}
\author{
Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
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
  
  # Calculate differential usage of genes, bins and junctions 
  du       <- DUreport( counts, targets )
  du       <- junctionDUreport( counts, targets, appendTo = du)
  
  # Calculate PSI / PIR for bins and junction.
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )
  
  # Select top tags
  binsdu <- binsDU(du)
  topTagsBins <- which( binsdu$bin.fdr <= 0.1 & 
                        abs(binsdu$logFC) >= log( 1.5, 2 )

  targetsPlot <- data.frame( bam=targets$bam, 
                             sample=targets$factor1, 
                             color=c(rep("blue", 3),rep("red", 3)), 
                             stringsAsFactors=FALSE )
                          
  auxdf<-binsdu[topTagsBins,]
  
  plotTopTags( auxdf["GENE02:E002",], genomeTxDb, targetsPlot, 
               output.dir = "testPlots" )

}


