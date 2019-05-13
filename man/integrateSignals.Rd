\name{integrateSignals}

\alias{integrateSignals}

\title{
  Differential junction usage estimation. COMPLETAR
}
\description{
  Estimate differential usage at junction level. COMPLETAR
}

\usage{
  integrateSignals(sr=NULL, asd = NULL, bin.fdr=0.05,unif=0.1,dPIN=0.05,dPIR=0.05,j.fdr=0.05,j.particip=0.1,usepvalBJS=FALSE,bjs.fdr=0.1, otherSources = NULL)
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{asd}{ 
    An object of class \code{ASpliDU}
  }  
  \item{bin.fdr}{Description TODO}
  \item{maxBinFDR}{Description TODO}
  \item{unif}{Description TODO}
  \item{dPIN}{Description TODO}
  \item{dPIR}{Description TODO}
  \item{j.fdr}{Description TODO}
  \item{j.particip}{Description TODO}
  \item{usepvalBJS}{Description TODO}
  \item{bjs.fdr}{Description TODO}
  \item{otherSources}{Description TODO}
}

\value{ 
  A \code{data.frame} with integrated signals }

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
  \code{\link{splicingReport}}
}


