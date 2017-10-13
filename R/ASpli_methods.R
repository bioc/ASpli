# ---------------------------------------------------------------------------- #
# Class definitions
setClass( Class = "ASpliFeatures",
    representation = representation(
        genes = "GRangesList",
        bins = "GRanges",
        junctions = "GRanges"))

setClass( Class = "ASpliCounts",
    representation = representation(
        gene.counts = "data.frame", 
        exon.intron.counts = "data.frame",
        junction.counts = "data.frame",
        e1i.counts = "data.frame", 
        ie2.counts = "data.frame",
        gene.rd = "data.frame",
        bin.rd = "data.frame"))

setClass( Class="ASpliAS",
    representation = representation(
        irPIR = "data.frame",
        altPSI = "data.frame",
        esPSI = "data.frame",
        junctionsPIR = "data.frame",
        junctionsPSI = "data.frame",
        join = "data.frame") )

setClass( Class = "ASpliDU",
    representation = representation(
        genes = "data.frame",
        bins = "data.frame",
        junctions = "data.frame" ))
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Set methods
setGeneric ( name = "binGenome", 
    def = function( genome, geneSymbols = NULL, 
        logTo = "ASpli_binFeatures.log" ) standardGeneric( "binGenome" ) )

setMethod(
    f = "binGenome",
    signature = "TxDb",
    definition = function ( genome, geneSymbols = NULL, logTo = "ASpli_binFeatures.log") {
      
      features <- new( Class = "ASpliFeatures" )
      
      if ( is.null( geneSymbols ) ) {
        # Recupera los nombres de los genes
        geneSymbols <- data.frame( names( transcriptsBy(genome) ), stringsAsFactors = FALSE)
        row.names(geneSymbols) <- names( transcriptsBy(genome) )
        colnames(geneSymbols)  <- "symbol" 
      }
      if ( ! is.null ( logTo ) ) {
        sink( file = logTo )
      }
      genes.by.exons <- .createGRangesGenes( genome, geneSymbols ) 
      
      msg <- paste( "* Number of extracted Genes =" , length( genes.by.exons ) )
      message( msg )
      if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
      
      
      exon.bins <- .createGRangesExons(genome, geneSymbols)
      #add locus_overlap
      index <- match(exon.bins@elementMetadata$locus, names(genes.by.exons))
      locus_overlap <- rep("-", length(exon.bins))
      locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
      mcols(exon.bins) <- append(mcols(exon.bins), DataFrame(locus_overlap=locus_overlap))
      
      msg <- paste( "* Number of extracted Exon Bins =", length( exon.bins ) )
      message( msg )
      if ( sink.number() > 0 ) cat( msg, "\n" ) 
      
      
      intron.tot <- .createGRangesIntrons(genome, geneSymbols)
      #add locus_overlap
      index <- match(intron.tot@elementMetadata$locus, names(genes.by.exons))
      locus_overlap <- rep("-", length(intron.tot))
      locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[index]
      mcols(intron.tot) <- append( mcols(intron.tot), 
          DataFrame(locus_overlap=locus_overlap) )
      
      msg <- paste( "* Number of extracted intron bins =", length( intron.tot ) )
      message( msg )
      if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
      
      transcripts <- .createGRangesTranscripts(genome)
      
      msg <- paste( "* Number of extracted trascripts =" , length( unlist( transcripts ) ) )
      message( msg )
      if ( sink.number() > 0 ) cat( paste( msg,"\n") )
      
      junctions <- .createGRangesJunctions( genome ) 
      
      #add locus_overlap
      index <- match( junctions@elementMetadata$locus, names( genes.by.exons ) )
      locus_overlap <- rep( "-", length( junctions ) )
      locus_overlap <- genes.by.exons@elementMetadata$locus_overlap[ index ]
      mcols(junctions) <- append( mcols( junctions ), 
          DataFrame( locus_overlap = locus_overlap ) )
      
      msg <- paste("* Number of extracted junctions =", length(junctions) )
      message( msg)
      if ( sink.number() > 0 ) cat( paste( msg, "\n" ) )
      
      intron.bins <- intron.tot[ intron.tot@elementMetadata$feature == "I"  ]
      intron.orig <- intron.tot[ intron.tot@elementMetadata$feature == "Io" ]
      
      class  <- rep("fullI", length( intron.orig ))
      mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( class=class ))
      event  <- rep("-", length( intron.orig ))
      eventJ <- rep("-", length( intron.orig ))
      mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( event=event ))
      mcols( intron.orig ) <- append( mcols( intron.orig ), DataFrame( eventJ=eventJ ))
      
      exons <- exons( genome )
      
      exons.introns <- .findAsBin( exons, exon.bins, intron.bins, transcripts, 
          junctions, logTo )
      fullT <- c(exons.introns,intron.orig)
      
      features@genes <- genes.by.exons
      features@bins <- fullT
      features@junctions <- junctions
      
      return( features ) 
    })
    
    
setGeneric( name = "rds",
  def = function( counts, targets ) standardGeneric("rds") )

# TODO: Las densidades de reads de genes y bins se calculan dos veces. Una 
# vez aca y otra vez cuando se hace el filtrado para hacer DU.
setMethod(
  f= "rds",
  signature = "ASpliCounts",
  definition = function(counts, targets) {
    geneStart <- ncol(countsg(counts))-nrow(targets)+1
    gene.rd <- cbind( countsg(counts)[,1:geneStart-1], 
                      countsg(counts)[,geneStart:ncol(countsg(counts))] / 
                          countsg(counts)$effective_length
                   )
    binStart <- ncol(countsb(counts))-nrow(targets)+1
    bin.rd <- cbind(countsb(counts)[, 1:binStart-1], 
                  countsb(counts)[,binStart:ncol(countsb(counts))]
                  /countsb(counts)$length)
                  
    tb <- match(bin.rd$locus, rownames(gene.rd))
          rdfinalb=cbind(bin.rd, bin.rd[,binStart:ncol(bin.rd)]
                   /gene.rd[tb,geneStart:ncol(countsg(counts))])
                            
    counts@gene.rd <- gene.rd
    counts@bin.rd <- rdfinalb
    return(counts)
  }
)

# readCounts
setGeneric (
  name = "readCounts",
  def = function( features, bam,  targets, cores = 1, readLength, maxISize, 
      minAnchor = NULL)
  standardGeneric("readCounts") )

setMethod(
  f = "readCounts",
  signature = "ASpliFeatures",
  definition = function( features, bam,  targets, cores = 1, readLength,  
      maxISize, minAnchor = 10) {

    counts <- new(Class="ASpliCounts")
    minA <- round( minAnchor * readLength / 100 )
    
    # Count Genes
    gene.hits <- .counterGenes( bam, featuresg( features ), cores )
    message("Read summarization by gene completed")
    counts@gene.counts <- gene.hits

    # Count exons 
    bins <- featuresb( features )
    exons.hits <- .counterBin( bam, bins, gene.hits, cores )
    message( "Read summarization by bin completed" )
    
    # Count introns
    introns <- c( bins[ mcols(bins)$feature == "I" ], 
                  bins[ mcols(bins)$feature == "Io"],
                  bins[ mcols(bins)$eventJ  == "IR"])

    # Count exon1 - intron regions
    e1i <- introns
    start( e1i ) <- start( introns ) - ( readLength - minA )
    end( e1i )   <- start( introns ) + ( readLength - minA )
    e1i.hits     <- .counterJbin(bam, e1i, gene.hits, cores, readLength)
    message("Read summarization by ei1 region completed")
    
    # Count intron - exon2 regions
    ie2 <- introns
    start( ie2 ) <- end( introns ) - ( readLength - minA )
    end( ie2 )   <- end( introns ) + ( readLength - minA )
    ie2.hits     <- .counterJbin( bam, ie2, gene.hits, cores, readLength )
    message("Read summarization by ie2 region completed")
    
    # Count junctions
    junction.hits = .counterJunctions( features, bam, cores, maxISize )
    message("Junction summarization completed")

    # Create result object
    counts@gene.counts <- gene.hits
    counts@exon.intron.counts <- exons.hits
    counts@junction.counts <- junction.hits 
    counts@e1i.counts <- e1i.hits  
    counts@ie2.counts <- ie2.hits
    counts <- rds( counts, targets )
    
    return(counts)
    
  }
)

setGeneric (
  name= "AsDiscover",
  def = function( counts, 
      targets, 
      features, 
      bam, 
      readLength, 
      threshold = 5, 
      cores = 1 ) standardGeneric("AsDiscover") )

setMethod(
    f = "AsDiscover",
    signature = "ASpliCounts",
    definition = function( counts, 
        targets, 
        features, 
        bam, 
        readLength, 
        threshold = 5, 
        cores = 1 ) {
      as <- new(Class = "ASpliAS")
      
      df0 <- countsj(counts)[ countsj(counts)$multipleHit == "-", ]
      df0 <- df0[ df0$gene != "noHit" , ]
      
      targets <- .condenseTargetsConditions( targets )
      
      jcounts <- .filterJunctionBySample( df0=df0,
          targets=targets,
          threshold=threshold )
      
      # Junctions PSI:
      junctionsPSI    <- .junctionsPSI_SUM( df0, targets )
      as@junctionsPSI <- junctionsPSI
      message("Junctions PSI completed")
      
      # Junctions PIR:
      junctionsPIR <- .junctionsDiscover( df=jcounts, 
          bam, 
          cores = cores , 
          readLength, 
          targets, 
          features )
      
      as@junctionsPIR <- junctionsPIR
      message("Junctions PIR completed")
      
      jranges <- .createGRangesExpJunctions( rownames( jcounts ) )
      
      # TODO: refactor this code to other functions 
      # ---------------------------------------------------------------------- #
      # Get all bins that are intronic or are associated to a Intron retention 
      # event
      ic <- rbind( countsb(counts)[countsb(counts)$feature == "I",], 
          countsb(counts)[countsb(counts)$feature == "Io",], 
          countsb(counts)[countsb(counts)$event   == "IR*",],
          countsb(counts)[countsb(counts)$event   == "IR",])
      # Get A GRanges object for intron bins, ordered by ic
      intranges <- featuresb(features)[ rownames(ic) ]
      
      # get exclusion junction counts, and make and index to ordered by ic
      dfe1e2 <- .e1e2JPIR( intranges, jcounts, targets )
      indexOrder <- match( dfe1e2$jbin, rownames( ic ) )
      
      # Get counts of inclusion junctions
      e1i <- .extractCountColumns( countse1i( counts ), targets )[ rownames(ic) ,]
      ie2 <- .extractCountColumns( countsie2( counts ), targets )[ rownames(ic) ,]
      
      j3 <- data.frame( matrix( NA, 
              nrow =  nrow( e1i ), 
              ncol =  length( targets$condition ) ), 
          stringsAsFactors = FALSE )
      colnames( j3 ) <- colnames( e1i )  
      
      j3bin <- rep( NA , nrow( j3 ) )
      j3bin[ indexOrder ] <- rownames( dfe1e2 )
      j3[ indexOrder, ] <- .extractCountColumns( dfe1e2, targets )
      
      # Sum exclusion and inclusion counts by condition
      sumE1i <- .sumByCond( e1i, targets )
      sumIe2 <- .sumByCond( ie2, targets )
      sumJ3  <- .sumByCond( j3,  targets )
      
      # Calculates pir
      pirValues <- ( sumE1i + sumIe2 ) / ( sumE1i + sumIe2 + 2 * sumJ3 )  
      
      # Creates result object
      result <- cbind( 
          data.frame( event = ic$event ), 
          data.frame( J1 = paste( rownames( e1i ), "E1I", sep="_") ), 
          e1i, 
          data.frame( J2 = paste( rownames( ie2 ), "IE2", sep="_") ), 
          ie2,
          data.frame( J3 = j3bin ),
          j3, 
          pirValues ) 
      
      
      message("Junctions IR PIR completed")
      
      as@irPIR <- result
      # ---------------------------------------------------------------------- #
      
      # ---------------------------------------------------------------------- #
      # Get all exons, except those that are associated to a intron retention
      # event
      ec <- countsb(counts)[countsb(counts)$feature == "E",]
      ec <- ec[ec$event != "IR",]
      ec <- ec[ec$event != "IR*",]
      
      exranges <- featuresb( features )[ rownames( ec ) ]
      
      fillAndReorderBy <- function( df , orderNames ) {
        indexOrder <- match( rownames( df ) , orderNames )
        result <- data.frame( 
            matrix( 
                NA,
                nrow = length( orderNames ),
                ncol = ncol( df ) ) )
        result[ indexOrder, ] <- df
        colnames( result ) <- colnames( df )
        rownames( result ) <- orderNames
        return( result )
      }
      
      dfstart  <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'start' )
      dfstart  <- fillAndReorderBy( dfstart , rownames( ec ) )
      dfend    <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'end' )   
      dfend    <- fillAndReorderBy( dfend , rownames( ec ) )
      dfwithin <- .getJPSIByOverlap( jranges, exranges, jcounts, targets, 'within' )
      dfwithin <- fillAndReorderBy( dfwithin , rownames( ec ) )
      
      events   <- mcols( exranges ) $ event
      # ---------------------------------------------------------------------- #
      
      # ---------------------------------------------------------------------- #
      # Get the subset of previosly selected exons and gets only those associated
      # with an alternative splicing site usage event 
      getAlternativeSS <- function( df, events ) {
        rbind(
            df[ events == "Alt3ss", ],
            df[ events == "Alt5ss", ],
            df[ events == "Alt3ss*", ],
            df[ events == "Alt5ss*", ] )
      }
      
      altJ1 <- getAlternativeSS( dfstart , events )
      altJ2 <- getAlternativeSS( dfend , events )
      altJ3 <- getAlternativeSS( dfwithin , events )
      
      sumAltJ1 <- .sumByCond( .extractCountColumns( altJ1, targets ), targets )
      sumAltJ1[is.na(sumAltJ1)] <- 0 
      sumAltJ2 <- .sumByCond( .extractCountColumns( altJ2, targets ), targets )
      sumAltJ2[is.na(sumAltJ2)] <- 0
      sumAltJ3 <- .sumByCond( .extractCountColumns( altJ3, targets ), targets )
      sumAltJ3[is.na(sumAltJ3)] <- 0
      
      altPsiValues <- ( sumAltJ1 + sumAltJ2 ) / ( sumAltJ1 + sumAltJ2 + sumAltJ3 )
      
      result <- cbind( 
          data.frame( event = mcols( exranges[ rownames( altJ1) ] )$ event ), 
          data.frame( J1 = altJ1$overlappedSubjectNames ), 
          .extractCountColumns( altJ1, targets ),
          data.frame( J2 = altJ2$overlappedSubjectNames ), 
          .extractCountColumns( altJ2, targets ),
          data.frame( J3 = altJ3$overlappedSubjectNames ),
          .extractCountColumns( altJ3, targets ), 
          altPsiValues )
      
      message("Junctions AltSS PSI completed")
      altPSI( as ) <- result
      # ---------------------------------------------------------------------- #
      
      # ---------------------------------------------------------------------- #
      # Get the subset of previosly selected exons and gets only those associated
      # with an exon skipping event and those not assigned to any splice event.  
      getES <- function( df, events ) {
        rbind(
            df[ events == "ES", ],
            df[ events == "-", ],
            df[ events == "ES*", ] )
      }
      
      esJ1 <- getES( dfstart , events )
      esJ2 <- getES( dfend , events )
      esJ3 <- getES( dfwithin , events )

      sumEsJ1 <- .sumByCond( .extractCountColumns( esJ1, targets ), targets )
      sumEsJ1[is.na(sumEsJ1)] <- 0 
      sumEsJ2 <- .sumByCond( .extractCountColumns( esJ2, targets ), targets )
      sumEsJ2[is.na(sumEsJ2)] <- 0
      sumEsJ3 <- .sumByCond( .extractCountColumns( esJ3, targets ), targets )
      sumEsJ3[is.na(sumEsJ3)] <- 0
      
      esPsiValues <- ( sumEsJ1 + sumEsJ2 ) / ( sumEsJ1 + sumEsJ2 + 2 * sumEsJ3 )
      
      result <- cbind( 
          data.frame( event = mcols( exranges[ rownames( esJ1) ] )$ event ), 
          data.frame( J1 = esJ1$overlappedSubjectNames ), 
          .extractCountColumns( esJ1, targets ), 
          data.frame( J2 = esJ2$overlappedSubjectNames ), 
          .extractCountColumns( esJ2, targets ),
          data.frame( J3 = esJ3$overlappedSubjectNames ),
          .extractCountColumns( esJ3, targets ),
          esPsiValues )
      
      message("Junctions ES PSI completed")
      
      esPSI( as ) <- result
      # ---------------------------------------------------------------------- #
      
      # TODO: joint podria ser un getter, pero no es necesario mantener toda
      # esta data repetida
      joint( as ) <- rbind( altPSI( as ), esPSI( as ), irPIR( as ) )
      
      return( as )

    })

setMethod( 
    f = 'subset',
    signature = 'ASpliAS',
    def = function( x, targets, select) .subset.ASpliAS( x, targets, select ) )
    
# ---------------------------------------------------------------------------- #
# writeAS
setGeneric (
  name = "writeAS",
  def  = function(as, output.dir = "as" )
  standardGeneric( "writeAS" ) )

setMethod(
  f = "writeAS",
  signature = "ASpliAS",
  definition = function( as, output.dir = "as" ) {
    
    # Creates output folder structure
    exonsFilePSI     <- file.path( output.dir, "exons", "exon.altPSI.tab" )       
    exonsFileES      <- file.path( output.dir, "exons", "exon.altES.tab" )       
    intronsFile      <- file.path( output.dir, "introns", "intron.irPIR.tab" )
    junctionsFilePIR <- file.path( output.dir, "junctions", "junction.PIR.tab" )                     	     	     	     	    
    junctionsFilePSI <- file.path( output.dir, "junctions", "junction.PSI.tab" )
    asDiscoverFile   <- file.path( output.dir, "as_discovery.tab" )
    
    file.exists( output.dir ) || dir.create( output.dir )
    for ( folder in unique( lapply( c( exonsFilePSI, exonsFileES, intronsFile, 
                junctionsFilePIR, junctionsFilePSI ), dirname ) ) ) {
      dir.create( folder )
    }

    # Export exons
    write.table( altPSI(as), exonsFilePSI, sep="\t", quote=FALSE, col.names=NA)
    write.table( esPSI(as), exonsFileES, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Introns
    write.table( irPIR(as), intronsFile, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Junctions
    write.table(junctionsPIR(as), junctionsFilePIR, sep="\t", quote=FALSE, col.names=NA)
    write.table(junctionsPSI(as), junctionsFilePSI, sep="\t", quote=FALSE, col.names=NA)

    # Export AS discovery table
    write.table( joint(as), asDiscoverFile, sep="\t", quote=FALSE, col.names=NA)
  }
)

# TODO:  Es necesario agregar todos los parametros con valores por default en
# la firma del metodo ? 
setGeneric (
    name = "DUreport",
    def = function( counts, 
        targets, 
        minGenReads  = 10,
        minBinReads  = 5,
        minRds = 0.05,
        offset = FALSE,
        offsetAggregateMode = c( "geneMode", "binMode" )[1],
        offsetUseFitGeneX = TRUE,
        contrast = NULL,
        forceGLM = FALSE,
        ignoreExternal = TRUE,
        ignoreIo = TRUE, 
        ignoreI = FALSE,
        filterWithContrasted = FALSE,
        verbose = FALSE
      ) standardGeneric("DUreport") )

#setGeneric (
#  name = "DUreport_DEXSeq",
#  def = function ( counts, ... ) standardGeneric("DUreport_DEXSeq") )

setMethod(
  f = "DUreport",
  signature = "ASpliCounts",
  definition = function( counts, 
      targets, 
      minGenReads  = 10,
      minBinReads  = 5,
      minRds = 0.05,
      offset = FALSE,
      offsetAggregateMode = c( "geneMode", "binMode" )[1],
      offsetUseFitGeneX = TRUE,
      contrast = NULL,
      forceGLM = FALSE,
      ignoreExternal = TRUE,
      ignoreIo = TRUE, 
      ignoreI = FALSE,
      filterWithContrasted = FALSE,
      verbose = FALSE
    ) { 
      .DUreport( counts, targets, minGenReads, minBinReads, minRds, offset, 
          offsetAggregateMode, offsetUseFitGeneX, contrast, forceGLM,
        ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose  )
  }
)

setGeneric( name = 'DUreportBinSplice',
    def = function( counts, targets, minGenReads  = 10, minBinReads = 5, 
        minRds = 0.05, contrast = NULL, forceGLM = FALSE,  
        ignoreExternal = TRUE, ignoreIo = TRUE, ignoreI = FALSE, 
        filterWithContrasted = FALSE, verbose = TRUE ) 
      standardGeneric( 'DUreportBinSplice'))

setMethod( 
    f = 'DUreportBinSplice',
    signature = 'ASpliCounts',
    definition = function( counts, 
			targets, 
			minGenReads  = 10, 
			minBinReads = 5,
			minRds = 0.05, 
			contrast = NULL, 
			forceGLM = FALSE, 
			ignoreExternal = TRUE,
			ignoreIo = TRUE, 
			ignoreI = FALSE, 
			filterWithContrasted = FALSE,
			verbose = TRUE ) {
      .DUreportBinSplice( counts, targets, minGenReads, minBinReads, minRds, 
          contrast, forceGLM, ignoreExternal, ignoreIo, ignoreI, 
          filterWithContrasted, verbose = TRUE ) 
    })

setGeneric( name = "junctionDUreport",
    def = function (  counts, 
        targets, 
        appendTo = NULL, 
        minGenReads = 10,
        minRds = 0.05,
        threshold = 5,
        offset   = FALSE,
        offsetUseFitGeneX = TRUE,
        contrast = NULL,
        forceGLM = FALSE 
        ) standardGeneric("junctionDUreport") )

  
setMethod(
    f = "junctionDUreport",
    signature = "ASpliCounts",
    definition = function ( 
        counts, 
        targets, 
        appendTo = NULL,
        minGenReads = 10,
        minRds = 0.05,
        threshold = 5,
        offset = FALSE,
        offsetUseFitGeneX = TRUE,
        contrast = NULL,
        forceGLM = FALSE 
        # -------------------------------------------------------------------- #
        # Comment to disable priorcounts usage in bin normalization 
        # , priorCounts = 0 
        # -------------------------------------------------------------------- #
        ) {
      
      .junctionDUreport( counts, targets, appendTo,  minGenReads,  minRds, 
          threshold, offset, offsetUseFitGeneX, contrast, 
          forceGLM ) 
    }
)
    
setMethod( f = 'subset',
    signature = 'ASpliCounts',
    def = function( x, targets, select ) { .subset.ASpliCounts( x, targets, select ) }  )

setGeneric( 
    name = 'filterDU',
    def = function( du, what = c( 'genes','bins','junctions'), fdr = 1, 
        logFC = 0, absLogFC = TRUE, logFCgreater = TRUE ) standardGeneric('filterDU') )

setMethod(
    f = 'filterDU',
    signature = "ASpliDU",
    definition = function( du, what = c( 'genes','bins','junctions'), fdr = 1, 
        logFC = 0, absLogFC = TRUE, logFCgreater = TRUE ) {
      .filter.ASpliDU( du, what, fdr, logFC, absLogFC, logFCgreater ) } )


# ---------------------------------------------------------------------------- #
# writeDU

setGeneric( name = 'containsGenesAndBins', 
    def = function ( du ) standardGeneric("containsGenesAndBins") )

setMethod( f = 'containsGenesAndBins', 
    signature = "ASpliDU",
    definition = function ( du ) {
      nrow( genesDE( du ) ) > 0  & nrow( binsDU( du) ) > 0 
    } )

setGeneric( name = 'containsJunctions', 
    def = function ( du ) standardGeneric("containsJunctions") )

setMethod( f = 'containsJunctions', 
    signature = "ASpliDU",
    definition = function ( du ) {
      nrow( junctionsDU( du ) ) > 0
    } )

setGeneric( name = "writeDU", 
    def = function ( du, output.dir="du"  ) standardGeneric( "writeDU" ) )

setMethod(
    f = "writeDU",
    signature = "ASpliDU",
    definition = function( du, output.dir="du" ) {
      
      paths <- list()
      # Creates output folder structure
      if ( containsGenesAndBins( du ) ) {
        genesFile     <- file.path( output.dir, "genes","gene.de.tab" )       
        exonsFile     <- file.path( output.dir, "exons", "exon.du.tab")       
        intronsFile   <- file.path( output.dir, "introns", "intron.du.tab")
        paths <- append( paths, list( genesFile, exonsFile, intronsFile  ))
      }
      
      if ( containsJunctions( du ) ) {
        junctionsFile <- file.path( output.dir, "junctions", "junction.du.tab")
        paths <- append( paths , list( junctionsFile ) )
      }
      
      file.exists( output.dir ) || dir.create( output.dir )
      for (filename in paths ) {
        dir.create( dirname( filename ) )
      }
      
      if ( containsGenesAndBins( du ) ) {
        # Export Genes  
        write.table( genesDE( du ), genesFile, sep = "\t", quote = FALSE, 
            col.names = NA )
        
        # Export Exons
        exonBins <- binsDU(du)[binsDU(du)$feature == "E",]
        exonBins <- exonBins[exonBins$event !="IR",]
        write.table( exonBins, exonsFile, sep="\t", quote=FALSE, col.names=NA)
        
        
        # Export Introns 
        intronBins <- rbind( 
            binsDU(du)[binsDU(du)$feature == "I" ,], 
            binsDU(du)[binsDU(du)$feature == "Io",],
            binsDU(du)[binsDU(du)$event   == "IR",])
        write.table( intronBins, intronsFile, sep = "\t", quote = FALSE, 
            col.names = NA )
      }
      # Export Junctions
      if ( containsJunctions( du ) ) {
        write.table( junctionsDU( du ), junctionsFile, sep = "\t", quote = FALSE, 
            col.names=NA )
      }
    }
)

setGeneric( name = 'mergeBinDUAS',
    def = function( du, as, targets, contrast = NULL  ) 
      standardGeneric( 'mergeBinDUAS' ))

setMethod( f = 'mergeBinDUAS',
    signature = c( 'ASpliDU', 'ASpliAS' ),
    definition = function( du, as, targets, contrast = NULL  ) {
      .mergeBinDUAS( du, as, targets, contrast ) } )
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Filter ASpliCounts by reads counts and read densisty
#setGeneric( name = 'filterReadCounts',
#    def = function (counts) standardGeneric( 'filterReadCounts') )
#
#
#
#setMethod( f = 'filterReadCounts',
#    signature = 'ASpliCounts',
#    definition = function( counts, targets, minGenRead = 10, minRdReads= 0.05,
#        types = c( 'minByGeneSet', 'minByGeneCondition', 'avgByGeneCondition') ) )

# filtering parameters
# ------------------------------------------------ #
# quantifier | grouping | what   | whom  | filter  #
# ------------------------------------------------ #
# min        | set      | count  | gene  | all     #
# avg        | cond.    | rd     | bin   | any     #
#            |          |        | junc. |         #
# -------------------------------------------------#

# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# plotBins
setGeneric( name = "plotBins",
    def = function(
        counts, 
        as,
        bin, 
        factorsAndValues, 
        targets,
        main            = NULL,
        colors          = c( '#2F7955', '#79552F', '#465579', '#A04935', '#752020', 
            '#A07C35') ,
        panelTitleColors = '#000000',
        panelTitleCex   = 1,
        innerMargins    = c( 2.1, 3.1, 1.1, 1.1 ),
        outerMargins    = c( 0, 0, 2.4, 0 ), 
        useBarplots     = NULL,
        barWidth        = 0.9,
        barSpacer       = 0.4,
        las.x           = 2,
        useHCColors     = FALSE,
        legendAtSide    = TRUE,
        outfolder       = NULL,
        outfileType     = c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')[1],
        deviceOpt       = NULL ) standardGeneric( 'plotBins' ) )

setMethod( 
    f = "plotBins",
    signature = 'ASpliCounts',
    definition = function( 
        counts, 
        as,
        bin, 
        factorsAndValues, 
        targets,
        main            = NULL,
        colors          = c( '#2F7955', '#79552F', '#465579', '#A04935', 
            '#752020', '#A07C35') ,
        panelTitleColors = '#000000',
        panelTitleCex   = 1,
        innerMargins    = c( 2.1, 3.1, 1.1, 1.1 ),
        outerMargins    = c( 0, 0, 2.4, 0 ), 
        useBarplots     = NULL,
        barWidth        = 0.9,
        barSpacer       = 0.4,
        las.x           = 2,
        useHCColors     = FALSE,
        legendAtSide    = TRUE,
        outfolder       = NULL,
        outfileType     = c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')[1],
        deviceOpt       = NULL ) {
      
      .plotBins( counts, as, bin, factorsAndValues, targets, main, colors, 
          panelTitleColors, panelTitleCex, innerMargins, outerMargins, 
          useBarplots, barWidth, barSpacer, las.x, useHCColors, legendAtSide,
          outfolder, outfileType, deviceOpt )
    } 
)
# ---------------------------------------------------------------------------- # 




# ---------------------------------------------------------------------------- #
# plotGenomicRegions
setGeneric( 
    name = "plotGenomicRegions", 
    def = function( 
        #counts,
        features,
        x, 
        genomeTxDb, 
        targets, 
        xIsBin = TRUE, 
        layout = 'auto', 
        colors = 'auto', 
        plotTitles = 'auto', 
        sashimi = FALSE, 
        zoomOnBins= FALSE, 
        deviceOpt = NULL, 
        highLightBin = TRUE, 
        outfolder = NULL, 
        outfileType = 'png', 
        mainFontSize = 24, 
        annotationHeight = 0.2,
        annotationCol = 'black', 
        annotationFill = 'gray', 
        annotationColTitle = 'black',
        preMergedBAMs = NULL,
        useTransparency = FALSE,
        tempFolder = 'tmp',
        avoidReMergeBams = FALSE,
        verbose = TRUE ) standardGeneric( "plotGenomicRegions" ) )

setMethod(
    f = "plotGenomicRegions",
    signature = "ASpliFeatures",
    definition = function ( 
#        counts,
        features, 
        x, 
        genomeTxDb, 
        targets, 
        xIsBin = TRUE, 
        layout = 'auto', 
        colors = 'auto', 
        plotTitles = 'auto', 
        sashimi = FALSE, 
        zoomOnBins= FALSE, 
        deviceOpt = NULL, 
        highLightBin = TRUE, 
        outfolder = NULL, 
        outfileType = 'png', 
        mainFontSize = 24, 
        annotationHeight = 0.2, 
        annotationCol = 'black', 
        annotationFill = 'gray', 
        annotationColTitle = 'black',
        preMergedBAMs = NULL,
        useTransparency = FALSE,
        tempFolder = 'tmp',
        avoidReMergeBams = FALSE,
        verbose = TRUE ) {
      
          .plotGenomicRegions(
              x, 
              genomeTxDb, 
#              counts,
              features,
              targets, 
              xIsBin, 
              layout, 
              colors, 
              plotTitles, 
              sashimi, 
              zoomOnBins, 
              deviceOpt, 
              highLightBin, 
              outfolder, 
              outfileType,
              mainFontSize, 
              annotationHeight, 
              annotationCol, 
              annotationFill, 
              annotationColTitle,
              preMergedBAMs,
              useTransparency,
              tempFolder,
              avoidReMergeBams ,
              verbose )
    }
)
        
# ---------------------------------------------------------------------------- #

