\name{jDUreport}

\alias{jDUreport}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level. COMPLETAR
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
    Contrasts 
  }
  \item{filterWithContrasted}{
    Apply filter on all conditions or only on contrasts
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


