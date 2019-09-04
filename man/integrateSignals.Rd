\name{integrateSignals}

\alias{integrateSignals}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  integrateSignals(sr = NULL, asd = NULL, bin.FC = 3, bin.fdr = 0.05,
                 nonunif = 1, usenonunif = FALSE, bin.inclussion = 0.1,
                 bjs.inclussion = 0.2, bjs.fdr = 0.1, a.inclussion =
                 0.3, a.fdr = 0.05, l.inclussion = 0.3, l.fdr = 0.05,
                 usepvalBJS = FALSE, otherSources = NULL, overlapType =
                 "any")
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{asd}{ 
    An object of class \code{ASpliDU}
  }  
  \item{bin.FC}{Description TODO}
  \item{bin.fdr}{Description TODO}
  \item{nonunif}{Description TODO}
  \item{usenonunif}{Description TODO}
  \item{bin.inclussion}{Description TODO}
  \item{bjs.inclussion}{Description TODO}
  \item{bjs.fdr}{Description TODO}
  \item{a.inclussion}{Description TODO}
  \item{a.fdr}{Description TODO}
  \item{l.inclussion}{Description TODO}
  \item{l.fdr}{Description TODO}
  \item{usepvalBJS}{Description TODO}
  \item{otherSources}{Description TODO}
  \item{overlapType}{Description TODO}    
}

\value{ 
  A \code{ASpliIntegratedSignals} with integrated signals }

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
  \code{\link{splicingReport}}
}


