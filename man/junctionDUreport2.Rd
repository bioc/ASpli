\name{junctionDUreportExt}

\alias{junctionDUreportExt}

\title{
  Differential junction usage estimation
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  junctionDUreportExt(asd, 
                    targets, 
                    minAvgCounts              = 5, 
                    contrast                  = NULL,
                    filterWithContrasted      = FALSE,
                    runUniformityTest         = FALSE,
                    mergedBams                = NULL,
                    maxPValForUniformityCheck = 0.2)
}

\arguments{
  \item{asd}{ 
    An object of class \code{ASpliAS}
  }
  \item{targets}{
    A dataframe containing sample, bam and experimental factor columns.
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
}

\value{ 
  An \code{ASpliJDU} object with results of differential usage of junctions }

\author{
  Andres Rabinovich, Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  
  #as <- new("ASpliAS")
  #targets <- 1
  #a <- junctionDUreportExt(as, targets)
  
}

\seealso{
  \code{\link{edgeR}}, \code{\link{ASpliAS}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


