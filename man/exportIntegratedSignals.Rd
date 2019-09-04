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
  \item{output.dir}{HTML reports output directory}
  \item{mergedBams}{Dataframe with path to merged bams for each condition }
  \item{jCompletelyIncluded}{If TRUE only plot junctions completely included in plot region. Else plot any overlapping junction in the region}
  \item{zoomRegion}{ Magnify plot region by this factor}
  \item{useLog}{Plot counts log}
  \item{tcex}{Text size}
  \item{ntop}{Only show n top signals}
  \item{openInBrowser}{Open reports in browser when done}
  \item{makeGraphs}{Generate graphs in reports}
  \item{bforce}{Force plot generation even if plot already exists}
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
}


