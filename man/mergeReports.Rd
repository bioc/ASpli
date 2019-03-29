\name{mergeReports}

\alias{mergeReports}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  mergeReports(bdu, 
               jdu)
}

\arguments{
  \item{bdu}{ 
    An object of class \code{ASpliDU}
  }
  \item{jdu}{ 
    An object of class \code{ASpliJDU}
  }

}

\value{ 
  An \code{ASpliMergedReports} object with results of differential usage of junctions }

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
  \code{\link{junctionDUreport}}, \code{\link{binDUreport}}, \code{\link{ASpliMergedReports}}
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


