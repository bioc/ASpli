\name{exportSplicingReports}

\alias{exportSplicingReports}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
exportSplicingReports( sr, output.dir="sr" , openInBrowser = F, maxBinFDR = 0.2, maxJunctionFDR = 0.2 )
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{output.dir}{Description TODO}
  \item{openInBrowser}{Description TODO}
  \item{maxBinFDR}{Description TODO}
  \item{maxJunctionFDR}{Description TODO}
}

\value{ 
  Produces html reports
}
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
  Accessors: \code{\link{junctionsDU}}
  Export: \code{\link{writeDU}}
}


