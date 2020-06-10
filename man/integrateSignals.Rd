\name{integrateSignals}

\alias{integrateSignals}

\title{
  Integrate signals
}
\description{
  Integrates differential usage signals from different sources using overlaping regions.
  Returns all overlaping signals present in the region filtered by different parameters.
  To construct this table, the region is marked with a 1 in the corresponding signal if the signal in that region was originated in that particular event.
  If other signals match the region, the other signals are also marked as 1. If other signals overlap but don't match, then those signals are marked with *.

  For example:
                      region b bjs ja jl
   1:       Chr1:45560-45645 1   1  0  0
   2:     Chr1:387673-388267 1   1  1  1
   3:     Chr1:406793-406902 1   1  0  *
   ...
  
   b: bin coverage signal / bjs: bin junction support signal / ja: junction anchor / jl: junction locale
   
   In region 3, the region was originaly a bin event, but a bjs event also matched the region so both were merge. Finally, a locale overlaped the region without matching, so an * was used for that signal.
Signals are marked as follows:   
For each region with bin coverage signal, tries to find a matching region in binbased(sr). If any match, asigns that event to the bin event. If no match, the   signal is * and no statistical information is available.  
For each region with bin junction support (bjs) signal, finds a matching J3 junction in binbased(sr). If any match, the junction statistics are retrieved. If   no match, the signal is * and no statistical information is available.  
For each region with anchor signal, finds a matching J3 junction in anchorbased(sr). If any match, the junction statistics are retrieved. If no match, finds a matching Io bin in binbased(sr) with the same J3 junction. If any match, the event is marked as IoR. If no match, finds a bin signal with a matching region. If no match, then the event is marked as Novel Alternative Splicing Pattern (NASP). If any match but the bin is marked as *, then the event is marked as Complex Splicing Pattern (CSP). Finally, if no anchorbased(sr) with a matching J3 region, then the signal is marked with * and no statistical information is available.
For each region with locale signal, finds a matching cluster with the same range in localebased(sr). If any match, that means that the signal was originaly from that cluster and all statistical information from that cluster is retrieved.  If no match, tries to find a cluster with a matching J3. If no match, then it is marked with *.
Finaly, it tries to label each locale event according to the following heuristic:
If all junctions share the same end or start, finds all bins in between. If they are all the same type of bins (ie exons or introns), then this is an Alt 5'/3' event. If not, the event is marked as CSP.
If the cluster is composed of three junctions and not all junctions share the same end or start. If there is a cluster with the same J3, then it's an ES event. If not, it's a CSP. Finally, if some of the junctions in the cluster are novel, then the event is classified as a Novel one.
For each of the other unclassified events, it checs if its an exonic feature. If it is, then the event is marked as ASCE. If the feature is intronic and no signal is *, then it's an IR. Finally, if it's an intronic feature an some signals are *, then the event is marked as CSP.

}

\usage{
  integrateSignals(sr = NULL, asd = NULL, bin.FC = 3, bin.fdr = 0.05,
                 nonunif = 1, usenonunif = FALSE, bin.inclussion = 0.2,
                 bjs.inclussion = 10.3, bjs.fdr = 0.01, a.inclussion =
                 0.3, a.fdr = 0.01, l.inclussion = 0.3, l.fdr = 0.01,
                 otherSources = NULL, overlapType = "any")
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{asd}{ 
    An object of class \code{ASpliDU}
  }  
  \item{bin.FC}{
    Filter bin signals by fold change. Actually, log2 fold change is return, so default would return only bin signlas with bin.fc > log2(3).
  }
  \item{bin.fdr}{
    Filter bin signals by fdr.
  }
  \item{nonunif}{
    Filter intronic bins with non uniform support (nonunif << 1 is uniform)
  }
  \item{usenonunif}{
    Use non uniformity as filter.
  }
  \item{bin.inclussion}{
    Filter bin signals by junction support with dPIR or dPSI accordingly.
  }
  \item{bjs.inclussion}{
    Filter annotated junction signals by junction inclussion with dPIR or dPSI accordingly.
  }
  \item{bjs.fdr}{
    Filter annotated junction signals by fdr.
  }
  \item{a.inclussion}{
    Filter anchor junction signals by junction inclussion with dPIR.
  }
  \item{a.fdr}{
    Filter anchor junction signals by fdr.
  }
  \item{l.inclussion}{
    Filter locale junction signals by junction inclussion with dPSI.
  }
  \item{l.fdr}{
    Filter locale junction signals by fdr.
  }
  \item{otherSources}{
    If user wants to compare ASpli results with results from other methods, otherSources must be a GenomicRange object with all the regions found with the other methods. It will be integrated with a new column next to signals information.
  }
  \item{overlapType}{
    Type of regions overlap matching between the different signals. Defaults to "any" and can be any of the following: "any", "start", "end", "within", "equal".
  }    
}

\value{ 
  A \code{ASpliIntegratedSignals} with integrated signals
  \item{signals}{
  region: region of event
  locus: locus of event
  bin.event: type of event
  b: 1 if coverage evidence of differential usage, 0 if not
  bjs: 1 if bin junction evidence of differential usage, 0 if not
  ja: 1 if junction anchor evidence of differential usage, 0 if not
  jl: 1 if junction locale evidence of differential usage, 0 if not
  bin: bin in the region
  feature: feature type
  J3: J3 junctions accross the region
  binreg: bin region without chromosome
  locus_overlap: other genes overlaping this locus
  b.fdr: coverage fdr evidence of differential usage
  b.logfc: coverage log fold change evidence of differential usage
  bjs.lr: bin junction likelihood ratio evidence of differential usage
  bjs.fdr: bin junction fdr evidence of differential usage
  bjs.logfc: bin junction log fold change evidence of differential usage
  bjs.nonuniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  bjs.inclussion: junction PIR or PSI
  a.lr: anchorage likelihood ratio evidence of differential usage
  a.fdr: anchorage fdr evidence of differential usage
  a.logfc: anchorage log fold change evidence of differential usage
  a.nonuniformity: if non uniformity test was performed, numbers closer to zero mean uniformity and closer to one mean non uniformity
  a.dpir: junction delta PIR
  l.lr: anchorage likelihood ratio evidence of differential usage
  l.fdr: anchorage fdr evidence of differential usage
  l.participation: participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster
  l.dparticipation: delta participation of the significant junction (FDR < maxFDRForParticipation) presenting maximal participation value inside the cluster 
  meanLogCPM: mean log of CPM
  }   
  \item{filters}{
  Filters used when integrating signals
  }
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
  sr       <- splicingReport(gbPaired, jPaired, gbcounts)
  is       <- integrateSignals(sr, jcounts)
  
  # Show integrate signals results and filters used
  is
  is@signals 
  is@filters

}

\seealso{
Accesors: \link{signals}, 
          \link{filters},
Export:   \link{exportIntegratedSignals}
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}, \
  code{\link{ASpliIntegratedSignals}}

}


