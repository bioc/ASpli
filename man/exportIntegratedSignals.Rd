\name{exportIntegratedSignals}

\alias{exportIntegratedSignals}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
exportIntegratedSignals( is, output.dir="is", 
                         sr, counts, features, asd,
                         mergedBams, 
                         jCompletelyIncluded = FALSE, zoomRegion = 1.5, 
                         useLog = FALSE, tcex = 1, ntop = NULL, openInBrowser = F, 
                         makeGraphs = T, bforce=FALSE
                        )
}

\arguments{
  \item{is}{ 
    An object of class \code{ASpliIntegratedSignals}
  }
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{counts}{ 
    An object of class \code{ASpliCounts}
  }
  \item{features}{ 
    An object of class \code{ASpliFeatures}
  }
  \item{asd}{ 
    An object of class \code{ASpliAS}
  }  
  \item{output.dir}{Description TODO}
  \item{mergedBams}{Description TODO}
  \item{jCompletelyIncluded}{Description TODO}
  \item{zoomRegion}{Description TODO}
  \item{useLog}{Description TODO}
  \item{tcex}{Description TODO}
  \item{ntop}{Description TODO}
  \item{openInBrowser}{Description TODO}
  \item{makeGraphs}{Description TODO}
  \item{bforce}{Description TODO}
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


