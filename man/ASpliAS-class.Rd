\name{ASpliAS-class}

\Rdversion{1.1}

\docType{class}

\alias{ASpliAS}
\alias{ASpliAS-class}
\alias{altPIN,ASpliAS-method}
\alias{esPIN,ASpliAS-method}
\alias{irPIR,ASpliAS-method}
\alias{joint,ASpliAS-method}
\alias{junctionsPIR,ASpliAS-method}
\alias{junctionsPJU,ASpliAS-method}
\alias{altPIN<-,ASpliAS,data.frame-method}
\alias{esPIN<-,ASpliAS,data.frame-method}
\alias{irPIR<-,ASpliAS,data.frame-method}
\alias{joint<-,ASpliAS,data.frame-method}
\alias{junctionsPIR<-,ASpliAS,data.frame-method}
\alias{junctionsPJU<-,ASpliAS,data.frame-method}
\alias{subset,ASpliAS-method}
\alias{writeAS,ASpliAS-method}
\alias{jDUreport,ASpliAS-method}

\title{Class \code{"ASpliAS"}}

\description{ Results of PIN and PIR using experimental junctions }

\section{ Slots }{
  \describe{
    \item{\code{irPIR}:}{
      Reports: event, e1i counts (J1), ie1 counts (J2), j_within (J3), PIR by 
      condition. J1, J2, J3 sum of junctions (J1, J2, J3) by condition.}
    \item{\code{altPIN}:}{
      Reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum 
      of junctions (J1, J2, J3) by condition.}
    \item{\code{esPIN}:}{
      Reports: event, J1 (start), J2 (end), J3 (exclusion), PSI. J1, J2, J3 sum 
      of junctions (J1, J2, J3) by condition.}
    \item{\code{join}:}{
      It is a combination of irPIR, altPIN and esPIN tables}
    \item{\code{junctionsPIR}:}{
      PIR metric for each experimental junction using e1i and ie2 counts. 
      Exclusion junction is the junction itself. This output helps to discover 
      new introns as well as new retention events}
    \item{\code{junctionsPJU}:}{
      Given a junction, it is possible to analyze if it shares start, end or 
      both with another junction. If so, is because there is more than one way 
      for/of splicing. Ratio between them along samples is reported.}
    \item{\code{targets}:}{
      TODO}      
  }
}

\author{ Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz}

\seealso{ 
  Methods: \code{\link{AsDiscover}}
  
  Accesors: \code{\link{altPIN}},
            \code{\link{irPIR}},
            \code{\link{esPIN}},
            \code{\link{joint}},
            \code{\link{junctionsPIR}},
            \code{\link{junctionsPJU}}
}