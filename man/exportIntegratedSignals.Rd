\name{exportIntegratedSignals}

\alias{exportIntegratedSignals}

\title{
  Export integrated signals.
}
\description{
  Export integrated signals in an easy to analyze HTML table.
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
  ##Calculates junction differential usage from junction counts 
  #j           <- jDUreport(asd, contrast, mergedBams = mergedBams)
  
  ##Calculates bin differential usage from bin counts 
  #gb          <- gbDUreport(counts, contrast)
  
  ##Generate splicig reports from different sources
  #sr          <- splicingReport(gb, j, counts)
  
  ##Integrate signals
  #is          <- integrateSignals(sr, asd)
  
  ##Export integrated signals
  #exportIntegratedSignals(is, output.dir = "is", sr, counts, features, mergedBams)
  
}

\seealso{
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}, \
  code{\link{ASpliIntegratedSignals}}
}


