\name{jDUreport}

\alias{jDUreport}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level.
}

\usage{
  jDUreport(asd, minAvgCounts = 5, contrast = NULL,
                 filterWithContrasted = TRUE, runUniformityTest =
                 FALSE, mergedBams = NULL, maxPValForUniformityCheck =
                 0.2, strongFilter = TRUE,
                 maxConditionsForDispersionEstimate = 24)
}

\arguments{
  \item{asd}{ 
    An object of class \code{ASpliAS}
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
    Default = NULL 
  }
  \item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is TRUE, it is strongly
    recommended to do not change this value.
  }
  \item{runUniformityTest}{
    Run uniformity test on Intron Retention
  }
  \item{mergedBams}{ 
    Path to merged bams for each testing condition
  }
  \item{maxPValForUniformityCheck}{
    To speed up uniformity test only check junctions with pval < maxPValForUniformityCheck
  } 
  \item{strongFilter}{
    Use only junctions that passed filters for clustering
  }  
  \item{maxConditionsForDispersionEstimate}{  
    In order to reduce resource usage, estimate dispersion with a reduced number of conditions.
  }
}

\value{ 
  An \code{ASpliJDU} object with results of differential usage of junctions }

\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
  #TODO
  
}

\seealso{
  \code{\link{edgeR}}, \code{\link{ASpliAS}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


