\name{jDUreport}

\alias{jDUreport}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level. 
  It makes separate estimates for different types of junctions.
  Estimation is made at junction level using diffSpliceDGE function from edgeR package. 
  The statistics works as follows. All junctions bellonging to the same event are put togheter and constitute the event "gene". 
  Then, each junction is tested against this "gene" in a similar fashion that bins are tested against their gene in diffSpliceDGE.
  Localec are clusters made of junctions that share an end with at least another junction in the cluster.
}

\usage{
  jDUreport(asd, 
            minAvgCounts                       = 5, 
            contrast                           = NULL,
            filterWithContrasted               = TRUE,
            runUniformityTest                  = FALSE,
            mergedBams                         = NULL,
            maxPValForUniformityCheck          = 0.2,
            strongFilter                       = TRUE,
            maxConditionsForDispersionEstimate = 24,
            formula                            = NULL,
            coef                               = NULL,
            maxFDRForParticipation             = 0.05,
            useSubset                          = FALSE)
}

\arguments{
  \item{asd}{ 
    An object of class \code{ASpliAS} with results of PSI and PIR using experimental junctions
  }
  \item{minAvgCounts}{ 
    Minimum average counts for filtering.
  }  
  \item{contrast}{
    Define the comparison between conditions to be tested.
    \code{contrast} should be a vector with length equal to the number of 
    experimental conditions defined by \code{targets}. The values of this vector
    are the coefficients that will be used to weight each condition, the order 
    of the values corresponds to the order given by \code{getConditions} 
    function. When \code{contrast} is NULL, defaults to a vector containing -1,
    as the first value, 1 as the second an zero for all the remaining values, 
    this corresponds to a pair comparison where the first condition is assumed 
    to be a control and the second condition is the treatment, all other 
    conditions are ignored. 
    If NULL must provide a formula.  
  }
  \item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is TRUE, it is strongly
    recommended to do not change this value.
  }
  \item{runUniformityTest}{
    Run uniformity test on Intron Retention. Sometimes Mutually Exclusive Exons (MEX) events can be confused with Intron Retention events. This test compares the standard deviation of the inner intron region (11 bases from both ends) to the mean of both intron ends. Numbers closer to 0 mean the event is more probably an Intron Retention event than an MEX event. The test takes some time to run so it defaults to FALSE. 
  }
  \item{mergedBams}{ 
    Path to merged bams for each testing condition. If no merged bams exist (for example, paired samples without replicates), use the same bams as targets.
  }
  \item{maxPValForUniformityCheck}{
    To speed up uniformity test only check junctions with pval < maxPValForUniformityCheck
  } 
  \item{strongFilter}{
    If strongFilter is TRUE, then we remove all events with at least one junction that doesn't pass the filter.
  }  
  \item{maxConditionsForDispersionEstimate}{  
    In order to reduce resource usage, estimate dispersion for statistics tests with a reduced number of conditions.
  }
  \item{formula}{
    Either a formula or a contrast can be tested.
    If formula is used, complex tests can be run.
    \code{formula} should be a formula specifying which experimental conditions defined by \code{targets} to test.
    If coef is specified, then that coefficient will be tested. If not, it defaults to the last term in the formula.
    }
  \item{coef}{
    For formula only. The coefficient to be tested. If null the test defaults to the last term in the formula}   
  \item{maxFDRForParticipation}{
   In order to calculate junctionPSI participation, only use significant junctions (ie junctions with FDR < maxFDRForParticipation).}   
  \item{useSubset}{
    Experimental. It is strongly recomended to leave the default, FALSE.}
}

\value{ 
  An \code{ASpliJDU} object with results of differential usage at \code{junctions} level.
}

\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
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
  jPaired <- jDUreport(jcounts, formula = formula(~subject+factor1))
  
  # Show junctions information
  jPaired
  localej(jPaired)
  localec(jPaired)
  anchorj(jPaired)
  anchorc(jPaired)
  jir(jPaired)
  jes(jPaired)
  jalt(jPaired)
  
  # Export results  
  writeJDU( du = jPaired, output.dir = "jPaired" )
}

\seealso{
  \code{\link{edgeR}}, \code{\link{ASpliAS}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


