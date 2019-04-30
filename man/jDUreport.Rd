\name{jDUreport}

\alias{jDUreport}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  jDUreport(asd, 
            minAvgCounts              = 5, 
            contrast                  = NULL,
            filterWithContrasted      = FALSE,
            runUniformityTest         = FALSE,
            mergedBams                = NULL,
            maxPValForUniformityCheck = 0.2,
            strongFilter              = FALSE)
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
}

\value{ 
  An \code{ASpliJDU} object with results of differential usage of junctions }

\author{
  Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{

  
}

\seealso{
  \code{\link{edgeR}}, \code{\link{ASpliAS}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


