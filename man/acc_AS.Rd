\name{AS accessors}

\alias{altPIN}
\alias{esPIN}
\alias{irPIR}
\alias{joint}
\alias{junctionsPIR}
\alias{junctionsPJU}
\alias{altPIN<-}
\alias{esPIN<-}
\alias{irPIR<-}
\alias{joint<-}
\alias{junctionsPIR<-}
\alias{junctionsPJU<-}

\title{ Accessors for ASpliAS object }

\description{ Methods to retrieve and set data in ASpliAS object. Setting data 
into an ASpliAS object is not a typical task and must be done with care, because
it can affect the integrity of the object. }

\usage{
altPIN( x )
esPIN( x )
irPIR( x )
joint( x )
junctionsPIR( x )
junctionsPJU( x )
}

\arguments{
  \item{ x }{ An ASpliAS object }
}
\value{ Returns dataframes with genomic metadata and  PSI and PIR metrics }

\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }

\examples{

# Accessing data tables from an ASpliAS object

#as <- aspliASexample()

#ap <- altPIN(as)
#ep <- esPIN(as)
#ip <- irPIR(as)
#j  <- joint(as)
#jpi <- junctionsPIR(as)
#jps <- junctionsPJU(as)

# Setting data tables to an ASpliAS object

#as2 <- new( 'ASpliAS' )

#altPIN( as2 )  <- ap
#esPIN( as2 )  <- ep
#irPIR( as2 )  <- ip
#joint( as2 )  <- j
#junctionsPIR( as2 )  <- jpi
#junctionsPJU( as2 )  <- jps

}

