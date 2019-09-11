\name{exportSplicingReports}

\alias{exportSplicingReports}

\title{
  Export splicing reports 
}
\description{
  Export splicing reports in easy to analyze HTML tables.
}

\usage{
exportSplicingReports( sr, output.dir="sr" , openInBrowser = F, maxBinFDR = 0.2, maxJunctionFDR = 0.2 )
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{output.dir}{HTML reports output directory}
  \item{openInBrowser}{Open reports in browser when done}
  \item{maxBinFDR}{Only show bins with FDR < maxBinFDR}
  \item{maxJunctionFDR}{Only show junctions with FDR < maxJunctionFDR}
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
  
  ##Export splicing report
  #exportSplicingReports(output.dir = "sr", sr)
  
}

\seealso{
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{splicingReport}}, \code{\link{ASpliSplicingReport}}
}


