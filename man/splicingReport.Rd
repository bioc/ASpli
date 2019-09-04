\name{splicingReport}

\alias{splicingReport}

\title{
  Differential junction usage estimation.
}
\description{
  Estimate differential usage at junction level.
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
  An \code{ASpliSplicingReport} object with results of differential usage of junctions }

\author{
  Estefania Mancini, Andres Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
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
}


