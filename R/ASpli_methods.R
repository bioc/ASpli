# ---------------------------------------------------------------------------- #
# Class definitions
setClass( Class = "ASpliFeatures",
          representation = representation(
            genes = "GRangesList",
            bins = "GRanges",
            junctions = "GRanges",
            transcriptExons = "GRangesList"))

setClass( Class = "ASpliCounts",
          representation = representation(
            gene.counts = "data.frame", 
            exon.intron.counts = "data.frame",
            junction.counts = "data.frame",
            e1i.counts = "data.frame", 
            ie2.counts = "data.frame",
            gene.rd = "data.frame",
            bin.rd = "data.frame", 
            targets = "data.frame",
            condition.order = "character"))

setClass( Class="ASpliAS",
          representation = representation(
            irPIR = "data.frame",
            altPIN = "data.frame",
            esPIN = "data.frame",
            junctionsPIR = "data.frame",
            junctionsPJU = "data.frame",
            join = "data.frame", 
            targets = "data.frame") )

setClass( Class = "ASpliDU",
          representation = representation(
            genes = "data.frame",
            bins = "data.frame",
            junctions = "data.frame",
            contrast = "numeric" ))

setClass( Class = "ASpliJDU",
          representation = representation(
            localec = "data.frame",
            localej = "data.frame",
            anchorc = "data.frame",
            anchorj = "data.frame",
            jir     = "data.frame",
            jes     = "data.frame",
            jalt    = "data.frame",
            contrast= "numeric"))

setClass( Class = "ASpliSplicingReport",
          representation = representation(
            binbased    = "data.frame",
            localebased = "data.frame",
            anchorbased = "data.frame",
            contrast    = "numeric"))

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
    
    transcriptExons <- exonsBy(genome, use.names=TRUE)
    
    features@genes <- genes.by.exons
    features@bins <- fullT
    features@junctions <- junctions
    features@transcriptExons <- transcriptExons
    
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
    counts@targets <- targets
    group                  <- .condenseTargetsConditions(targets)$condition
    counts@condition.order <- levels(factor( group, unique( group ), ordered = TRUE ))
    return(counts)
  }
)

# gbCounts is a wrapper around readCounts for improved legibility
setGeneric (
  name = "gbCounts",
  def = function( features, bam=NULL, targets, cores = 1, minReadLength, maxISize, 
                  minAnchor = 10)
    standardGeneric("gbCounts") )

setMethod(
  f = "gbCounts",
  signature = "ASpliFeatures",
  definition = function( features, bam=NULL, targets, cores = 1, minReadLength,  
                         maxISize, minAnchor = 10) {
    return(readCounts( features, bam=bam, targets, cores = cores, minReadLength, maxISize, minAnchor = minAnchor))
  }
)

# jCounts is a wrapper around AsDiscover for improved legibility
setGeneric (
  name= "jCounts",
  def = function( counts, 
                  features, 
                  bam = NULL, 
                  minReadLength, 
                  threshold = 5, 
                  cores = 1,
                  minAnchor = 10) standardGeneric("jCounts") )

setMethod(
  f = "jCounts",
  signature = "ASpliCounts",
  definition = function( counts, 
                         features, 
                         bam = NULL, 
                         minReadLength, 
                         threshold = 5, 
                         cores = 1,
                         minAnchor = 10) {
    return(AsDiscover( counts, features, bam = bam, minReadLength, threshold, cores, minAnchor ))
  }
)

# readCounts
setGeneric (
  name = "readCounts",
  def = function( features, bam=NULL, targets, cores = 1, minReadLength, maxISize, 
                  minAnchor = 10)
    standardGeneric("readCounts") )

setMethod(
  f = "readCounts",
  signature = "ASpliFeatures",
  definition = function( features, bam=NULL, targets, cores = 1, minReadLength,  
                         maxISize, minAnchor = 10) {
    
    
    #Create result object
    counts <- new(Class="ASpliCounts")
    
    #Minimal anchors
    minAnchor <- if ( ! is.null(minAnchor) ) minAnchor else 10
    minA <- round( minAnchor * minReadLength / 100 )
    ptm <- proc.time()
    if(is.null(bam)) {
      ntargets <- nrow(targets)
    }else{
      ntargets <- 1
    }
    
    #Generates sample names in case there arent any
    targets <- .generateSamplesNames(targets)
    
    for(target in 1:ntargets){
      
      if(ntargets > 1){
        #Verbose
        message(paste("Summarizing", rownames(targets)[target]))
        
        #Load bam from current target
        bam <- loadBAM(targets[target, ])
      }      
      
      # Count Genes
      gene.hits <- .counterGenes( bam, featuresg( features ), cores )
      if(ncol(counts@gene.counts) == 0){
        counts@gene.counts <- gene.hits
      }else{
        counts@gene.counts <- cbind(counts@gene.counts, .extractCountColumns(gene.hits, targets[target, ]))
        colnames(counts@gene.counts)[ncol(counts@gene.counts)] <- rownames(targets)[target]
      }
      if(ntargets == 1) message("Read summarization by gene completed")
      
      # Count exons 
      bins <- featuresb( features )
      exons.hits <- .counterBin( bam, bins, gene.hits, cores )
      if(ncol(counts@exon.intron.counts) == 0){
        counts@exon.intron.counts <- exons.hits
      }else{
        counts@exon.intron.counts <- cbind(counts@exon.intron.counts, .extractCountColumns(exons.hits, targets[target, ]))
        colnames(counts@exon.intron.counts)[ncol(counts@exon.intron.counts)] <- rownames(targets)[target]      
      }
      if(ntargets == 1) message( "Read summarization by bin completed" )
      
      # Count introns
      introns <- c( bins[ mcols(bins)$feature == "I" ], 
                    bins[ mcols(bins)$feature == "Io"],
                    bins[ mcols(bins)$eventJ  == "IR"])
      
      # Count exon1 - intron regions
      e1i <- introns
      start( e1i ) <- start( introns ) - ( minReadLength - minA )
      end( e1i )   <- start( introns ) + ( minReadLength - minA )
      e1i.hits     <- .counterJbin(bam, e1i, gene.hits, cores, minReadLength)
      if(ncol(counts@e1i.counts) == 0){
        counts@e1i.counts <- e1i.hits
      }else{
        counts@e1i.counts <- cbind(counts@e1i.counts, .extractCountColumns(e1i.hits, targets[target, ]))
        colnames(counts@e1i.counts)[ncol(counts@e1i.counts)] <- rownames(targets)[target]         
      }
      if(ntargets == 1) message("Read summarization by ei1 region completed")
      
      # Count intron - exon2 regions
      ie2 <- introns
      start( ie2 ) <- end( introns ) - ( minReadLength - minA )
      end( ie2 )   <- end( introns ) + ( minReadLength - minA )
      ie2.hits     <- .counterJbin( bam, ie2, gene.hits, cores, minReadLength )
      if(ncol(counts@ie2.counts) == 0){
        counts@ie2.counts <- ie2.hits
      }else{
        counts@ie2.counts <- cbind(counts@ie2.counts, .extractCountColumns(ie2.hits, targets[target, ]))
        colnames(counts@ie2.counts)[ncol(counts@ie2.counts)] <- rownames(targets)[target]   
      }
      if(ntargets == 1) message("Read summarization by ie2 region completed")
      
      # Count junctions
      junction.hits    <- .counterJunctions( features, bam, cores, maxISize )
      if(ncol(counts@junction.counts) == 0){
        counts@junction.counts <- junction.hits
      }else{
        dt1                    <- data.table(counts@junction.counts, keep.rownames = T)
        dt2                    <- data.table(.extractCountColumns(junction.hits, targets[target, ]), keep.rownames = T)
        dt3                    <- data.frame(merge(dt1, dt2, by="rn", all.x=T, all.y=T))
        for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
          dt3[, s]           <- as.character(dt3[, s])
          junction.hits[, s] <- as.character(junction.hits[, s])
        }
        rownames(dt3)          <- dt3[, "rn"]
        dt3                    <- dt3[, -1]
        dt3[dt2$rn, 1:8]       <- .extractDataColumns(junction.hits, targets[target, ])
        counts@junction.counts <- dt3
        counts@junction.counts[is.na(counts@junction.counts)] <- 0
      }
      if(ntargets == 1) message("Junction summarization completed")
      
      gc()
      if(ntargets > 1){     
        sptm <- (proc.time() - ptm)[3]/target/60        
        message(paste("ETA:", round(sptm*(nrow(targets) - target)), "min"))      
      }
    }
    
    for(s in c("junction", "gene", "strand", "multipleHit", "symbol", "gene_coordinates", "bin_spanned", "j_within_bin")){
      counts@junction.counts[, s] <- as.factor(counts@junction.counts[, s])
    }
    colnames(counts@junction.counts)[9:ncol(counts@junction.counts)] <- rownames(targets)
    junctions.order <- sort(rownames(counts@junction.counts))
    junctions.order <- strsplit2(junctions.order, "[.]")
    junctions.order <- GRanges(seqnames=junctions.order[, 1], IRanges(start=as.numeric(junctions.order[, 2]), end=as.numeric(junctions.order[, 3])))
    junctions.order <- sort(junctions.order)
    junctions.order <- paste(junctions.order@seqnames, junctions.order@ranges@start, (junctions.order@ranges@start+junctions.order@ranges@width-1), sep=".")
    counts@junction.counts <- counts@junction.counts[junctions.order, ]
    
    # Create result object
    counts <- rds( counts, targets )
    gc()
    return(counts)
    
  }
)

setGeneric (
  name= "AsDiscover",
  def = function( counts, 
                  features, 
                  bam = NULL, 
                  minReadLength, 
                  threshold = 5, 
                  cores = 1,
                  minAnchor = 10) standardGeneric("AsDiscover") )

setMethod(
  f = "AsDiscover",
  signature = "ASpliCounts",
  definition = function( counts, 
                         features, 
                         bam = NULL, 
                         minReadLength, 
                         threshold = 5, 
                         cores = 1,
                         minAnchor = 10) {
    
    targets <- counts@targets
    
    as  <- new(Class = "ASpliAS")
    as@targets <- targets
    
    df0 <- countsj(counts)[ countsj(counts)$multipleHit == "-", ]
    df0 <- df0[ df0$gene != "noHit" , ]
    
    targets <- .condenseTargetsConditions( targets )
    
    jcounts <- .filterJunctionBySample( df0=df0,
                                        targets=targets,
                                        threshold=threshold )
    
    # Junctions PSI:
    junctionsPSI    <- .junctionsPSI_SUM( df0, targets )
    as@junctionsPJU <- junctionsPSI
    message("Junctions PJU completed")
    
    # Junctions PIR:
    if(is.null(bam)) {
      ntargets <- nrow(targets)
      for(target in 1:ntargets){
        
        if(ntargets > 1){
          #Load bam from current target
          bam <- loadBAM(targets[target, ])
          junctionsPIR <- .junctionsDiscover( df=jcounts, 
                                              bam, 
                                              cores = cores , 
                                              minReadLength, 
                                              targets[target, ], 
                                              features,
                                              minAnchor = minAnchor)        
        }
  
        if(ncol(as@junctionsPIR) == 0){
          as@junctionsPIR <- junctionsPIR
        }else{
          as@junctionsPIR <- cbind(as@junctionsPIR, junctionsPIR[, 3:6])
        }
      }
      junctions.order <- c(1, 2, 
                           seq(from=3, to=ncol(as@junctionsPIR), by=4),
                           seq(from=4, to=ncol(as@junctionsPIR), by=4),
                           seq(from=5, to=ncol(as@junctionsPIR), by=4))
      as@junctionsPIR <- as@junctionsPIR[, junctions.order]
      colnames(as@junctionsPIR)[c(-2, -1)] <- rep(rownames(targets), times=3)
      inicio_j1 <- 3
      inicio_j2 <- inicio_j1+nrow(targets)
      inicio_j3 <- inicio_j2+nrow(targets)
      j1 <- .sumByCond( as@junctionsPIR[, inicio_j1:(inicio_j1+nrow(targets)-1)],     targets )
      j2 <- .sumByCond( as@junctionsPIR[, inicio_j2:(inicio_j2+nrow(targets)-1)],     targets )
      j3 <- .sumByCond( as@junctionsPIR[, inicio_j3:(inicio_j3+nrow(targets)-1)],     targets )
      pirValues <- ( j1 + j2 ) / ( j1 + j2 + 2 * j3 )
      as@junctionsPIR <- cbind(as@junctionsPIR, pirValues)
    }else{
      
      junctionsPIR <- .junctionsDiscover( df=jcounts, 
                                          bam, 
                                          cores = cores , 
                                          minReadLength, 
                                          targets, 
                                          features,
                                          minAnchor = minAnchor ) 
      as@junctionsPIR <- junctionsPIR
    }
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
    colnames(dfe1e2)[c(-1,-2)] <- rownames(targets)
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
    colnames(altJ1)[-ncol(altJ1)] <- rownames(targets)
    colnames(altJ2)[-ncol(altJ2)] <- rownames(targets)
    colnames(altJ3)[-ncol(altJ3)] <- rownames(targets)
      
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
    
    message("Junctions AltSS PIN completed")
    altPIN( as ) <- result
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
    colnames(esJ1)[-ncol(esJ1)] <- rownames(targets)
    colnames(esJ2)[-ncol(esJ2)] <- rownames(targets)
    colnames(esJ3)[-ncol(esJ3)] <- rownames(targets)
    
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
    
    message("Junctions ES PIN completed")
    
    esPIN( as ) <- result
    # ---------------------------------------------------------------------- #
    
    # TODO: joint podria ser un getter, pero no es necesario mantener toda
    # esta data repetida
    joint( as ) <- rbind( altPIN( as ), esPIN( as ), irPIR( as ) )
    
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
    write.table( altPIN(as), exonsFilePSI, sep="\t", quote=FALSE, col.names=NA)
    write.table( esPIN(as), exonsFileES, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Introns
    write.table( irPIR(as), intronsFile, sep="\t", quote=FALSE, col.names=NA)
    
    # Export Junctions
    write.table(junctionsPIR(as), junctionsFilePIR, sep="\t", quote=FALSE, col.names=NA)
    write.table(junctionsPJU(as), junctionsFilePSI, sep="\t", quote=FALSE, col.names=NA)
    
    # Export AS discovery table
    write.table( joint(as), asDiscoverFile, sep="\t", quote=FALSE, col.names=NA)
  }
)

# TODO:  Es necesario agregar todos los parametros con valores por default en
# la firma del metodo ? 
setGeneric (
  name = "DUreport.norm",
  def = function( counts, 
                  targets, 
                  minGenReads  = 10,
                  minBinReads  = 5,
                  minRds = 0.05,
                  contrast = NULL,
                  forceGLM = FALSE,
                  ignoreExternal = TRUE,
                  ignoreIo = TRUE, 
                  ignoreI = FALSE,
                  filterWithContrasted = FALSE,
                  verbose = FALSE,
                  threshold = 5
  ) standardGeneric("DUreport.norm") )

#setGeneric (
#  name = "DUreport_DEXSeq",
#  def = function ( counts, ... ) standardGeneric("DUreport_DEXSeq") )

setMethod(
  f = "DUreport.norm",
  signature = "ASpliCounts",
  definition = function( counts, 
                         targets, 
                         minGenReads  = 10,
                         minBinReads  = 5,
                         minRds = 0.05,
                         contrast = NULL,
                         forceGLM = FALSE,
                         ignoreExternal = TRUE,
                         ignoreIo = TRUE, 
                         ignoreI = FALSE,
                         filterWithContrasted = FALSE,
                         verbose = FALSE,
                         threshold = 5
  ) { 
    offset = FALSE
    offsetAggregateMode = c( "geneMode", "binMode" )[1]
    offsetUseFitGeneX = TRUE    
    .DUreport( counts, targets, minGenReads, minBinReads, minRds, offset, 
               offsetAggregateMode, offsetUseFitGeneX, contrast, forceGLM,
               ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose, threshold  )
  }
)

setGeneric (
  name = "DUreport.offset",
  def = function( counts, 
                  targets, 
                  minGenReads  = 10,
                  minBinReads  = 5,
                  minRds = 0.05,
                  offsetAggregateMode = c( "geneMode", "binMode" )[1],
                  offsetUseFitGeneX = TRUE,
                  contrast = NULL,
                  forceGLM = FALSE,
                  ignoreExternal = TRUE,
                  ignoreIo = TRUE, 
                  ignoreI = FALSE,
                  filterWithContrasted = FALSE,
                  verbose = FALSE
  ) standardGeneric("DUreport.offset") )

setMethod(
  f = "DUreport.offset",
  signature = "ASpliCounts",
  definition = function( counts, 
                         targets, 
                         minGenReads  = 10,
                         minBinReads  = 5,
                         minRds = 0.05,
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
    offset = TRUE
    .DUreport( counts, targets, minGenReads, minBinReads, minRds, offset, 
               offsetAggregateMode, offsetUseFitGeneX, contrast, forceGLM,
               ignoreExternal, ignoreIo, ignoreI, filterWithContrasted, verbose  )
  }
)

setGeneric( name = 'gbDUreport',
            def = function( counts, minGenReads  = 10, minBinReads = 5, 
                            minRds = 0.05, contrast = NULL, forceGLM = FALSE,  
                            ignoreExternal = TRUE, ignoreIo = TRUE, ignoreI = FALSE, 
                            filterWithContrasted = FALSE, verbose = TRUE ) 
              standardGeneric( 'gbDUreport'))

setMethod( 
  f = 'gbDUreport',
  signature = 'ASpliCounts',
  definition = function( counts, 
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
    .DUreportBinSplice( counts, minGenReads, minBinReads, minRds, 
                        contrast, forceGLM, ignoreExternal, ignoreIo, ignoreI, 
                        filterWithContrasted, verbose = TRUE ) 
  })


setGeneric( name = "jDUreport",
            def = function (asd, 
                            minAvgCounts              = 5, 
                            contrast                  = NULL,
                            filterWithContrasted      = FALSE,
                            runUniformityTest         = FALSE,
                            mergedBams                = NULL,
                            maxPValForUniformityCheck = 0.2,
                            strongFilter              = TRUE,
                            maxConditionsForDispersionEstimate = 24
            ) standardGeneric("jDUreport") )


setMethod(
  f = "jDUreport",
  signature = "ASpliAS",
  definition = function (
    asd,
    minAvgCounts              = 5, 
    contrast                  = NULL,
    filterWithContrasted      = FALSE,
    runUniformityTest         = FALSE,
    mergedBams                = NULL,
    maxPValForUniformityCheck = 0.2,
    strongFilter              = TRUE,
    maxConditionsForDispersionEstimate = 24
  ) {
    
    .junctionDUreportExt( asd, minAvgCounts, contrast, 
                          filterWithContrasted, runUniformityTest, mergedBams, maxPValForUniformityCheck, strongFilter,
                          maxConditionsForDispersionEstimate) 
  }
)


setGeneric( name = "splicingReport",
            def = function (bdu, 
                            jdu,
                            maxBinFDR = 0.2, 
                            maxJunctionFDR = 0.2
            ) standardGeneric("splicingReport") )


setMethod(
  f = "splicingReport",
  signature = "ASpliDU",
  definition = function (
    bdu,
    jdu, 
    maxBinFDR = 0.2, 
    maxJunctionFDR = 0.2
  ) {
    
    .splicingReport( bdu, jdu, maxBinFDR = 0.2, maxJunctionFDR = 0.2 ) 
  }
)

setGeneric( name = "integrateSignals",
            def = function (sr=NULL,
                            asd = NULL, 
                            bin.fdr=0.05,
                            unif=0.1,
                            dPIN=0.05,
                            dPIR=0.05,
                            j.fdr=0.05,
                            j.particip=0.1,
                            usepvalBJS=FALSE,
                            bjs.fdr=0.1,
                            otherSources = NULL
            ) standardGeneric("integrateSignals") )


setMethod(
  f = "integrateSignals",
  signature = "ASpliSplicingReport",
  definition = function (
    sr=NULL,
    asd = NULL,
    bin.fdr=0.05,
    unif=0.1,
    dPIN=0.05,
    dPIR=0.05,
    j.fdr=0.05,
    j.particip=0.1,
    usepvalBJS=FALSE,
    bjs.fdr=0.1,
    otherSources = NULL
  ) {
    .integrateSignals(sr, asd, bin.fdr, unif, dPIN, dPIR, j.fdr, j.particip, usepvalBJS, bjs.fdr, otherSources) 
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
    
    #paths <- list()
    # Creates output folder structure
    #if ( containsGenesAndBins( du ) ) {
    #  genesFile     <- file.path( output.dir, "genes","gene.de.tab" )       
    #  exonsFile     <- file.path( output.dir, "exons", "exon.du.tab")       
    #  intronsFile   <- file.path( output.dir, "introns", "intron.du.tab")
    #  paths <- append( paths, list( genesFile, exonsFile, intronsFile  ))
    #}
    
    #if ( containsJunctions( du ) ) {
    #  junctionsFile <- file.path( output.dir, "junctions", "junction.du.tab")
    #  paths <- append( paths , list( junctionsFile ) )
    #}
    
    file.exists( output.dir ) || dir.create( output.dir )
    output.dir <- paste(output.dir, paste(names(du@contrast)[du@contrast != 0], collapse="-"), sep="/")
    file.exists( output.dir ) || dir.create( output.dir )


    #for (filename in paths ) {
    #  dir.create( dirname( filename ) )
    #}
    
    if ( containsGenesAndBins( du ) ) {
      # Export Genes  
      write.table( genesDE( du ), paste(normalizePath(output.dir), "gene.de.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names = NA )
      
      # Export Exons
      exonBins <- binsDU(du)[binsDU(du)$feature == "E",]
      exonBins <- exonBins[exonBins$event !="IR",]
      write.table( exonBins, paste(normalizePath(output.dir), "exon.du.tab", sep="/"), sep="\t", quote=FALSE, col.names=NA)
      
      
      # Export Introns 
      intronBins <- rbind( 
        binsDU(du)[binsDU(du)$feature == "I" ,], 
        binsDU(du)[binsDU(du)$feature == "Io",],
        binsDU(du)[binsDU(du)$event   == "IR",])
      write.table( intronBins, paste(normalizePath(output.dir), "intron.du.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names = NA )
    }
    # Export Junctions
    if ( containsJunctions( du ) ) {
      write.table( junctionsDU( du ), paste(normalizePath(output.dir), "junction.du.tab", sep="/"), sep = "\t", quote = FALSE, 
                   col.names=NA )
    }
  }
)

setGeneric( name = "writeJDU", 
            def = function ( jdu, output.dir="jdu"  ) standardGeneric( "writeJDU" ) )

setMethod(
  f = "writeJDU",
  signature = "ASpliJDU",
  definition = function( jdu, output.dir="jdu" ) {
    
    file.exists( output.dir ) || dir.create( output.dir )
    output.dir <- paste(output.dir, paste(names(jdu@contrast)[jdu@contrast != 0], collapse="-"), sep="/")
    file.exists( output.dir ) || dir.create( output.dir )
    
    for(slotName in slotNames(jdu)){
      # Export Genes  
      write.table( slot(jdu, slotName), paste(output.dir, paste0(slotName, ".txt"), sep="/"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE )
    }
    
  }
)

setGeneric( name = "writeSplicingReport", 
            def = function ( sr, output.dir="sr"  ) standardGeneric( "writeSplicingReport" ) )

setMethod(
  f = "writeSplicingReport",
  signature = "ASpliSplicingReport",
  definition = function( sr, output.dir="sr" ) {
    
    file.exists( output.dir ) || dir.create( output.dir )
    #output.dir <- paste(output.dir, paste(names(jdu@contrast)[jdu@contrast != 0], collapse="-"), sep="/")
    #file.exists( output.dir ) || dir.create( output.dir )
    
    for(slotName in slotNames(sr)){
      # Export Genes  
      if(class(slot(sr, slotName)) == "data.frame"){
        write.table( slot(sr, slotName), paste(output.dir, paste0(slotName, ".txt"), sep="/"), sep = "\t", quote = FALSE, col.names = T, row.names = F )
      }
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


setGeneric( name = "exportSplicingReports", 
            def = function ( sr, output.dir="sr"  ) standardGeneric( "exportSplicingReports" ) )

setMethod(
  f = "exportSplicingReports",
  signature = "ASpliSplicingReport",
  definition = function( sr, output.dir="sr" ) {
    output.dir <- paste0(normalizePath(output.dir), "/", paste0(names(sr@contrast)[sr@contrast != 0], collapse="-"))
    file.exists( output.dir ) || dir.create( output.dir , recursive = T)
    
    
    for(s in slotNames(sr)){
      if(class(slot(sr, s)) == "data.frame"){
        b <- slot(sr, s)
        if(s == "binbased"){
            b$feature <- as.factor(b$feature)
            b$bin.event <- as.factor(b$bin.event)
            #b[, c(10:12, 15:ncol(b))] <- apply(b[, c(10:12, 15:ncol(b))], 2, function(s){return(signif(as.numeric(s), digits = 4))})
        }
	columnas_numericas <- which(sapply(b, class) == "numeric")
	b[, columnas_numericas] <- apply(b[, columnas_numericas], 2, function(s){return(signif(as.numeric(s), digits = 4))})
        titulo <- paste0('ASpli: ', s, ". Contrasts: ", paste(names(sr@contrast)[sr@contrast != 0], collapse = " - "))
        y <- datatable(b,
                       escape = TRUE,
                       filter ="top",
                       extensions = c('Buttons', 'KeyTable'), 
                       options = list(dom = 'lfrtBip',
                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print', I('colvis')),
                                      columnDefs = list(
                                        #  list(visible = FALSE, targets = c(0, 2, 3)),
                                        list(orderable = FALSE, className =
                                               'details-control', targets = 1)
                                      ),
                                      keys = TRUE
                       ),   caption = htmltools::tags$caption(
                         style = 'caption-side: top; text-align: left;',
                         htmltools::h1(titulo)
                       ))    
        ffile <- paste0(normalizePath(output.dir), "/", s, "Report.html")
        suppressWarnings(saveWidget(y, file = ffile, title = paste0(names(sr@contrast)[sr@contrast != 0], collapse="-")))
        browseURL(ffile)
      }
    }    
  }
)

setGeneric( name = "exportIntegratedSignals", 
            def = function ( is, output.dir="is", 
                             sr, counts, features, 
                             mergedBams, 
                             jCompletelyIncluded = FALSE, zoomRegion = 1.5, 
                             useLog = FALSE, tcex = 1, ntop = NULL
                             ) standardGeneric( "exportIntegratedSignals" ) )

setMethod(
  f = "exportIntegratedSignals",
  signature = "data.table",
  definition = function( is, output.dir="is", sr, counts, features, mergedBams, jCompletelyIncluded = FALSE, zoomRegion = 1.5, useLog = FALSE, tcex = 1, ntop = NULL) {
    
    if(!all(colnames(is) %in% c('region','locus','b','bjs','ja','jl','bin','feature','bin.event','J3','binreg','locus_overlap','b.fdr','b.logfc','bjs.fdr','bjs.logfc','bjs.nonuniformity','bjs.inclussion','a.fdr','a.logfc','a.nonuniformity','a.participation','a.dparticipation','l.fdr','l.logfc','l.participation','l.dparticipation'))){
      stop("is must be a data.table generated by integrateSignals method. Not all integratedSignal columns are present")
    }
    
    if(class(sr) != "ASpliSplicingReport"){
      stop("sr must be an ASpliSplicingReport object")
    }

    if(class(counts) != "ASpliCounts"){
      stop("counts must be an ASpliCounts object")
    }

    if(class(features) != "ASpliFeatures"){
      stop("features must be an ASpliFeatures object")
    }
    mergedBams <- mergedBams[mergedBams$condition %in% names(sr@contrast)[sr@contrast != 0], ]
    if(nrow(mergedBams) == 0){
      stop("Merged bams dont match with contrasts")  
    }
    output.dir <- paste0(normalizePath(output.dir), "/", paste0(names(sr@contrast)[sr@contrast != 0], collapse="-"))
    file.exists( output.dir ) || dir.create( output.dir, recursive = T )
    file.exists( paste0(output.dir, "/img") ) || dir.create( paste0(output.dir, "/img") )

    is[,b:=as.factor(b)]
    is[,bjs:=as.factor(bjs)]
    is[,ja:=as.factor(ja)]
    is[,jl:=as.factor(jl)]
    is$feature[is.na(is$feature)] <- "-" 
    is[,feature:=as.factor(feature)]
    is[,bin.event:=as.factor(bin.event)]
    is[,b.fdr:=signif(as.numeric(b.fdr), 4)]
    is[,b.logfc:=signif(as.numeric(b.logfc), 4)]
    is[,bjs.fdr:=signif(as.numeric(bjs.fdr), 4)]
    is[,bjs.logfc:=signif(as.numeric(bjs.logfc), 4)]    
    is[,bjs.nonuniformity:=signif(as.numeric(bjs.nonuniformity), 4)]
    is[,bjs.inclussion:=signif(as.numeric(bjs.inclussion), 4)]    
    is[,a.fdr:=signif(as.numeric(a.fdr), 4)]
    is[,a.logfc:=signif(as.numeric(a.logfc), 4)]
    is[,a.nonuniformity:=signif(as.numeric(a.nonuniformity), 4)]
    is[,a.participation:=signif(as.numeric(a.participation), 4)]
    is[,a.dparticipation:=signif(as.numeric(a.dparticipation), 4)]
    is[,l.fdr:=signif(as.numeric(l.fdr), 4)]
    is[,l.logfc:=signif(as.numeric(l.logfc), 4)]
    is[,l.participation:=signif(as.numeric(l.participation), 4)]
    is[,l.dparticipation:=signif(as.numeric(l.dparticipation), 4)]

    message("Generating graphs...")
    if(!is.numeric(ntop)){
        ntop <- length(is$region)
    }else{
      if(ntop < 1){
        ntop <- length(is$region)
      }
    }
    #for(i in 1:length(dev.list))dev.off()
    for(i in 1:ntop){
      r <- is$region[i]
      if(i %% 10 == 0){
        message(paste0(signif(i/ntop, 2)*100, "% completed"))
      }
      #tryCatch({
      print(r)
        png(width = 1400, height=700, filename = paste0(normalizePath(output.dir), "/img/", r, "_gene.png"))
        .plotSplicingPattern(r, is, counts, features, mergedBams, sr, genePlot = TRUE, jCompletelyIncluded, zoomRegion, useLog, tcex)
        dev.off()
        png(width = 1400, height=700, filename = paste0(normalizePath(output.dir), "/img/", r, ".png"))
        .plotSplicingPattern(r, is, counts, features, mergedBams, sr, genePlot = FALSE, jCompletelyIncluded, zoomRegion, useLog, tcex)
        dev.off()
      # }, warning = function(warning_condition) {
      #     message(warning_condition)   
      #     dev.off()
      # }, error = function(error_condition) {
      #     message(error_condition)
      #     dev.off()
      # }, finally={
      #   
      # })
    }
    titulo <- paste0('ASpli: integrated signals. Contrasts: ', paste(names(sr@contrast)[sr@contrast != 0], collapse = " - "))
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'View'),
          th(rowspan = 2, 'Region'),
          th(rowspan = 2, 'Locus'),
          th(rowspan = 2, 'Locus overlap'),
          th(rowspan = 2, 'Bin Evidence'),
          th(rowspan = 2, 'Bin SJ Evidence'),
          th(rowspan = 2, 'Anchor Evidence'),
          th(rowspan = 2, 'Locale Evidence'),
          th(rowspan = 2, 'Bin'),
          th(rowspan = 2, 'Feature'),
          th(rowspan = 2, 'Event'),
          th(colspan = 2, 'Bins', bgcolor="#DCDCDC"),
          th(colspan = 4, 'Bin Supporting Junctions', bgcolor="#C0C0C0"),
          th(colspan = 5, 'Anchor Junctions', bgcolor="#DCDCDC"),
          th(colspan = 4, 'Locale Junctions', bgcolor="#C0C0C0")
        ),
        tr(
          lapply(c("logFC", "FDR", "logFC", "FDR", "Non Uniformity", "Inclussion", "logFC", "FDR", "Non uniformity", "Participation", "dParticipation", "logFC", "FDR", "Participation", "dParticipation"), th)
        )
      )
    ))
    y <- datatable(cbind('&oplus;', is[1:ntop,c("region", "locus", "locus_overlap", "b", "bjs", "ja", "jl", "bin", "feature", "bin.event", "b.logfc", "b.fdr", "bjs.logfc", "bjs.fdr", "bjs.nonuniformity", "bjs.inclussion", "a.logfc", "a.fdr", "a.nonuniformity", "a.participation", "a.dparticipation", "l.logfc", "l.fdr", "l.participation", "l.dparticipation")]),
              rownames = FALSE,
              escape = -1,
              filter ="top",
              fillContainer = F,
              extensions = c('Buttons', 'KeyTable'), 
              options = list(dom = 'lfrtBip',
                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print', I('colvis')),
                             columnDefs = list(
                               #  list(visible = FALSE, targets = c(0, 2, 3)),
                               list(orderable = FALSE, className =
                                      'details-control', targets = 0)
                             ),
                             keys = TRUE
              ),   caption = htmltools::tags$caption(
                style = 'caption-side: top; text-align: left;',
                htmltools::h1(titulo)
              ), container = sketch,
              callback = JS("
            table.order([10, 'asc']).draw();
            table.column(0).nodes().to$().css({cursor: 'pointer'});
            var format = function(d) {
             return '<div><h4>Gene view</h4></br><img src=\"img/' + d[1] + '_gene.png\" height=\"700\"></img></div><div></br></div><div><h4>Region view</h4></br><img src=\"img/' + d[1] + '.png\" height=\"700\"></img></div>';
            };
            table.on('click', 'td.details-control', function() {
               var td = $(this), row = table.row(td.closest('tr'));
               if (row.child.isShown()) {
                  row.child.hide();
                  td.html('&oplus;');
              } else {
                  row.child(format(row.data())).show();
                  td.html('&CircleMinus;');
              }
           });")
    )    
    
    ffile <- paste0(normalizePath(output.dir), "/integratedSignals.html")
    suppressWarnings(saveWidget(y, file = ffile, title = paste(names(sr@contrast)[sr@contrast != 0], collapse = " - ")))
    browseURL(ffile)
    
  }
)

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

