\name{integrateSignals}

\alias{integrateSignals}

\title{
  Integrate signals
}
\description{
  Integrates differential usage signals from different sources using overlaping regions.
  Returns all overlaping signals present in the region filtered by different parameters. 
  For example:
                      region b bjs ja jl
   1:       Chr1:45560-45645 1   1  0  0
   2:     Chr1:387673-388267 1   1  1  1
   3:     Chr1:406793-406902 1   1  0  0
   ...
  
   b: bin coverage signal / bjs: bin junction support signal / ja: junction anchor / jl: junction locale
}

\usage{
  integrateSignals(sr = NULL, asd = NULL, bin.FC = 3, bin.fdr = 0.05,
                 nonunif = 1, usenonunif = FALSE, bin.inclussion = 0.1,
                 bjs.inclussion = 0.2, bjs.fdr = 0.1, a.inclussion =
                 0.3, a.fdr = 0.05, l.inclussion = 0.3, l.fdr = 0.05,
                 otherSources = NULL, overlapType =
                 "any")
}

\arguments{
  \item{sr}{ 
    An object of class \code{ASpliSplicingReport}
  }
  \item{asd}{ 
    An object of class \code{ASpliDU}
  }  
  \item{bin.FC}{Filter bin signals by fold change}
  \item{bin.fdr}{Filter bin signals by fdr}
  \item{nonunif}{Filter intronic bins with non uniform support (nonunif < 1 is uniform)}
  \item{usenonunif}{Use non uniformity as filter}
  \item{bin.inclussion}{Filter bin signals by junction support with dPIR or dPSI accordingly}
  \item{bjs.inclussion}{Filter annotated junction signals by junction inclussion with dPIR or dPSI accordingly}
  \item{bjs.fdr}{Filter annotated junction signals by fdr}
  \item{a.inclussion}{Filter anchor junction signals by junction inclussion with dPIR}
  \item{a.fdr}{Filter anchor junction signals by fdr}
  \item{l.inclussion}{Filter locale junction signals by junction inclussion with dPSI}
  \item{l.fdr}{Filter locale junction signals by fdr}
  \item{otherSources}{Integrate ASpli signals with signals from other methods.}
  \item{overlapType}{Type of region overlaps. Defaults to "any" and can be any of the following: "any", "start", "end", "within", "equal".}    
}

\value{ 
  A \code{ASpliIntegratedSignals} with integrated signals }

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

}

\seealso{
  \code{\link{gbDUreport}}, \code{\link{jDUreport}}, \code{\link{ASpliSplicingReport}},  \code{\link{splicingReport}}, \
  code{\link{ASpliIntegratedSignals}}
}


