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
  \item{binbased}{
  bin: name of bin
  feature: type of feature
  bin.event: type of event
  locus: bin locus
  locus_overlap: does this bin overlap other bins
  gene_coordinates: bin gene coordinates
  start: bin start
  end: bin end
  length: bin length
  bin.logFC: bin log fold change for differential usage
  bin.pvalue: bin pvalue for differential usage
  bin.fdr: bin fdr for differential usage
  J3: J3 junction/s across the region
  J3.multiplicity: do multiple junctions cross the region
  J3.logFC: J3 log fold change for differential usage
  J3.log.mean: J3 log of mean counts accross all conditions for this junction
  cluster.LR: likelihood ratio of cluster for differential usage
  cluster.pvalue: pvalue of cluster for differential usage
  cluster.fdr: fdr of cluster for differential usage
  dPIR: cluster delta PIR
  dPSI: cluster delta PSI
  cluster.NonUniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  From cluster.NonUniformity to the end, J1, J2 and J3 counts for all samples
  }     
  \item{localebased}{
  junction: junction name
  junction.annotated: annotated or new junction
  junction.log.mean: log of mean counts of junction
  junction.logFC: log fold change of junction for differential usage
  junction.pvalue: pvalue of junction for differential usage
  junction.fdr: fdr of junction for differential usage
  From junction.fdr to junction.cluster, junction counts for all samples
  junction.cluster: name of the cluster the junction belongs to
  junction.participation: the maximal participation value observed across contrasted condictions for the junction
  junction.dParticipation: delta participation of the maximal participation value observed across contrasted condictions for the junction
  cluster.locus: locus the cluster belongs to
  cluster.size: number of junctions belonging to the cluster.
  cluster.LR: likelihood ratio of cluster differential usage.
  cluster.pvalue: pvalue of cluster differential usage.
  cluster.FDR: fdr of cluster differential usage.
  cluster.range: cluster location.
  cluster.participation: participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster
  cluster.dParticipation: delta participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster 
  bin: if the junction matched with a bin, name of the bin
  bin.pvalue: pvalue of the bin the junction matched with
  bin.fdr: fdr of the bin the junction matched with
  }  
  \item{anchorbased}{
  junction: junction name
  junction.annotated: annotated or new junction
  junction.log.mean: log of mean counts of junction
  junction.logFC: log fold change of junction for differential usage
  junction.LR: likelihood ratio of junction for differential usage
  junction.pvalue: pvalue of junction for differential usage
  junction.fdr: fdr of junction for differential usage
  J1.pvalue: pvalue of J1 for differential usage
  J2.pvalue: pvalue of J2 for differential usage
  NonUniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  dPIR: junction delta PIR
  From junction.fdr to cluster.locus, junction counts for all samples
  cluster.locus: locus the cluster belongs to
  cluster.pvalue: pvalue of cluster differential usage.
  cluster.FDR: fdr of cluster differential usage.
  bin: if the junction matched with a bin, name of the bin
  bin.pvalue: pvalue of the bin the junction matched with
  bin.fdr: fdr of the bin the junction matched with
  }     
  \item{contrast}{
  Conditions contrasted by ASpli
  }  
}

\author{
  Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\seealso{
  Accesors: \code{\link{binbased}}, 
            \code{\link{localebased}},
            \code{\link{anchorbased}},
  Export:   \code{\link{writeSplicingReport}}
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}
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



