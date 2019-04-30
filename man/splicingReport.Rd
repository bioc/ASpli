\name{splicingReport}

\alias{splicingReport}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  splicingReport(bdu, 
               jdu, maxBinFDR = 0.2, maxJunctionFDR = 0.2)
}

\arguments{
  \item{bdu}{ 
    An object of class \code{ASpliDU}
  }
  \item{jdu}{ 
    An object of class \code{ASpliJDU}
  }
  \item{maxBinFDR}{Description TODO}
  \item{maxJunctionFDR}{Description TODO}
}

\value{ 
  An \code{ASpliSplicingReport} object with results of differential usage of junctions }

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
  \code{\link{jDUreport}}, \code{\link{gbDUreport}}, \code{\link{ASpliSplicingReport}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


