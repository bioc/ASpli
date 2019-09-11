\name{splicingReport}

\alias{splicingReport}

\title{
  Splicing report
}
\description{
  This function integrates differential usage information from different sources.
  Receives bin usage and junction usage and merges them in three steps.
  First merges all bin based junction information (jir, jes, jalt).
  Then merges locale information with coverage information in the regions.
  Finally merges anchor information with coverage information in the regions.
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
  An \code{ASpliSplicingReport} object with junction differential usage report }

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
    
}

\seealso{
    \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}
}


