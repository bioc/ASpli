# Make plot for a set of bins
.plotBins <- function( counts, as, bin, factorsAndValues, targets, main ,
    colors ,  panelTitleColors, panelTitleCex, innerMargins , 
    outerMargins, useBarplots, barWidth, barSpacer, las.x, useHCColors,
    legendAtSide, outfolder, outfileType, deviceOpt) {
  
  for ( cBin in bin) {
    print(cBin)
    filename <- if ( ! is.null( outfolder )) {
          file.path( outfolder , .makeValidFileName( paste0(cBin,'.pr.',outfileType ) ) )
        } else {
          outfolder
        }
    
    .plotSingleBin( counts, as, cBin, factorsAndValues, targets, main, colors, 
        panelTitleColors, panelTitleCex, innerMargins, outerMargins, 
        useBarplots, barWidth, barSpacer, las.x, useHCColors, legendAtSide,
        filename, outfileType, deviceOpt )
  }
}

# ---------------------------------------------------------------------------- #
# .getHCColor return a color that should have high contrast with a another given
# color and also with the white background 
.getHCColor <- function ( color ) {
  rgbCol <- col2rgb ( color , alpha = TRUE)
  brightness <- mean( rgbCol [ 1:3,1 ] )
  if ( brightness < 95 ) {
    if ( rgbCol [ 1,1 ] > rgbCol [ 2,1 ] & rgbCol [ 1,1 ] > rgbCol [ 3,1 ] ) {
      return( rgb( 0.3,0.5,1,1 ) )
    }
    if ( rgbCol [ 2,1 ] > rgbCol [ 1,1 ] & rgbCol [ 2,1 ] > rgbCol [ 3,1 ] ) {
      return( rgb( 1, 0.1,0.1,1 ) )
    }
    return( rgb( 0.3,1, 0.5,1 ) )
  } else {
    return( rgb( 0, 0, 0 ,1 ) )
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# This function get a vector data for a bin property, with many values as
# samples. Values are averaged by condition a then returned as a list of
# vectors grouped by the main factor values.
.prepareData <- function( samplesProfileData, gridPanels, targets, factorsAndValues ) {
  
  mainFactorIndex <- length( factorsAndValues )
  targets$samplesNames <- rownames( targets)
  
  data <- lapply ( 1:nrow( gridPanels ) , function ( i ) {
        
        samplesData <- merge( targets, gridPanels[ i, , drop = FALSE ], 
            by = colnames( gridPanels[ i, , drop = FALSE ] ) )
        
        mainFactorValues <- factorsAndValues[[ mainFactorIndex ]]
        
        plotData <- sapply ( mainFactorValues, 
            function( x ) { 
              samples <- samplesData[ 
                  samplesData[ , names( factorsAndValues[ mainFactorIndex ] ) ] == x  , 
                  'samplesNames' ]
              rowMeans ( samplesProfileData[ rownames( targets ) %in% samples ] )
            } )
        
        plotData[ is.na( plotData )] <- 0
        return( plotData )
      } )
  return(data)
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .makeLegends creates a top-left legend over a single plot.
.makeLegends <- function( main , useHCColors, col, panelTitleColor, panelTitleCex ) {
  legend( 
      x = "topleft", 
      legend = main,
      adj = c( 0 , 0 ),
      yjust = 0.5,
      bty = "n",
      cex = panelTitleCex,
      text.col = if ( useHCColors ) .getHCColor( col ) else { panelTitleColor } )
}
# ---------------------------------------------------------------------------- #


# -------------------------------------------------------------------------- #
# Function to make a single xy plot 
.plotLines <- function ( data, dataCons, gridPanels, nPoints, xTicksLabels, 
    factorsAndValues, mainFactorIndex, minValue, maxValue, main, col, 
    legendAtSide, las.x , panelTitleCex ) {
  
  plot( x = 1:length( dataCons ) , 
      y = dataCons,
      pch = 20,
      type = 'p',
      col = 'gray',
      ylim = c( minValue, maxValue),
      xlab='',
      xaxt = "n",
      ylab = if ( legendAtSide ) main else '',
      yaxt = 'n',
      cex.lab = panelTitleCex )
  
  secFactorTable <- table( gridPanels[ , ncol(gridPanels)] )
  st = 0
  for ( i in 1:nrow( secFactorTable )) {
    from <-  (st + 1)
    to <- st + secFactorTable[i] * nPoints
    at <- c(from:to)
    axis( side = 1, 
        at = at, 
        labels= xTicksLabels[from:to], 
        las=las.x , 
        cex.axis=0.7, tck=0,
        mgp=c(3,0.3,0),
        col= rgb(0.9,0.9,0.9) )
    st <- to
  }
  
  st = 0
  for( i in 1:nrow( gridPanels )) { 
    lines( x = st + 1:length( factorsAndValues[[ mainFactorIndex ]] ) , 
        y = data[[i]],
        col = col)
    st <- st + length( factorsAndValues[[mainFactorIndex]] )
  }
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to make a single bar plot 
.plotBars <- function( dataCons, nPoints, barWidth, barSpacer, maxValue, main, 
    useHCColors, panelTitleColor, xTicksLabels, col, las.x, legendAtSide, panelTitleCex) {
  spacers <- c( (1 - barWidth) / barWidth, barSpacer)
  dataCons <- matrix( dataCons, ncol = nPoints, byrow = TRUE ) 
  barplot ( height = t( dataCons ),
      col = col,
      ylim = c( 0, maxValue),
      xlab = '',
      ylab = if ( legendAtSide ) main else '',
      yaxt = 'n',
      width = barWidth,
      space = spacers,
      border = '#DDDDDDFF',
      beside = TRUE,
      las = 2,
      cex.lab = panelTitleCex)
  
  ticks <- matrix( 1, 
      ncol = ncol(dataCons),
      nrow = nrow(dataCons))
  spacerTick <- barWidth * barSpacer - ( 1 - barWidth )
  ticks <- cbind( ticks, spacerTick )
  
  ticks <- as.vector( t( ticks ) )
  
  ac <- 0
  tickSum <- c()
  for ( i in 1 : (length( ticks ) -1) ) {
    tickSum <- append( tickSum , ac + ticks[i])
    ac <- tickSum[i]
  }
  
  tickSum <- tickSum - ( 1 - ( spacerTick + (1-barWidth) ) ) + barWidth / 2
  
  
  tlblmat <- matrix ( xTicksLabels, byrow=TRUE, ncol = ncol(dataCons) ) 
  tlblmat <- cbind( tlblmat, '' )
  xTicksLabels <- as.vector( t(tlblmat) )[1:( length( tickSum ) )]
  
  axis( 
      side = 1, 
      at = tickSum , 
      labels = xTicksLabels, 
      las = las.x , 
      cex.axis = 0.7, 
      tck=0,
      mgp = c(3,0.3,0),
      col = rgb( 0.9,0.9,0.9, 0 ),
      padj = 0.5 )
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Create the labels for factor corresponding to the x axis of plots. 
.makeXticksLabels <- function( justOneFactor, gridPanels, factorsAndValues, 
    mainFactorIndex, data ) {
  
  if ( ! justOneFactor ) {
    
    xTicksLabels <- rep( apply( gridPanels, 1, 
            function( x ) paste0(x,collapse=  ".") ), sapply( data, length  ) )
    xTicksLabels <- apply( gridPanels, 1, function( x ) paste0(x,collapse=  ".") )
    xTicksLabels <- expand.grid( factorsAndValues[[ mainFactorIndex]], xTicksLabels )
    xTicksLabels <- paste0( xTicksLabels[,2] ,".",  xTicksLabels[,1] )
    
  } else {
    xTicksLabels <- factorsAndValues[[mainFactorIndex]]
  }
  return( xTicksLabels )
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Main plotting function
.makePlot <- function ( data, useBarPlot = FALSE, main, col, factorsAndValues, 
    justOneFactor, gridPanels, nPoints, legendAtSide, las.x, panelTitleCex,
    panelTitleColor, useHCColors, barWidth, barSpacer ) {
  
  mainFactorIndex <- length( factorsAndValues )

  # Get the labels of ticks in the x axis.
  xTicksLabels <- .makeXticksLabels( justOneFactor , gridPanels, 
      factorsAndValues, mainFactorIndex, data )
  
  # Get the y axis minimum and maximum values. Used to scale the plots correctly.
  maxValue <- max ( sapply ( data, max ))
  minValue <- min ( sapply ( data, min ))
  
  # Make a unique vector will all data
  dataCons <- Reduce( function( a, b) { a <- append( a, b )}, data )

  # Remove box around plots
  par( bty = 'n')
  
  # Make the main plot
  if ( ! useBarPlot ) {
    .plotLines (data, dataCons , gridPanels , nPoints, xTicksLabels, 
        factorsAndValues, mainFactorIndex, minValue, maxValue, main, 
        col = col, legendAtSide = legendAtSide , las.x, panelTitleCex )
  } else {
    .plotBars(dataCons, nPoints, barWidth, barSpacer, maxValue, main, 
        useHCColors, panelTitleColor, xTicksLabels, col, las.x, legendAtSide, 
        panelTitleCex )
  }
  
  # Add legend inside plot is required
  if ( ! legendAtSide ) {
    .makeLegends( main, useHCColors, col, panelTitleColor, panelTitleCex )
  }
  
  # Draw y axis ticks
  axis( 2, las=1, cex.axis=0.7 )
  
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# .plotbins draw all plots for a given bin
.plotSingleBin <- function( 
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
    outfile         = NULL,
    outfileType     = 'png',
    deviceOpt       = NULL
    ) {
  
  # -------------------------------------------------------------------------- #
  # Manipulate factors
  # If just one factor, add another fictional factor with only one value to
  # make computation easier.
  justOneFactor <- FALSE
  if ( length ( factorsAndValues) == 1 ) {
    factorsAndValues <- append( factorsAndValues, list( single = 1 ) )
    colnames( joint(as) ) [ 
        match( unique( targets[,names(factorsAndValues)[[1]]] ), colnames( joint(as) ))
    ] <- paste( unique( targets[,names(factorsAndValues)[[1]]] ), '1',sep='_')
    targets <- cbind( targets, single = 1 )
    justOneFactor <- TRUE
  }
  # Create condition names
  targets <- .condenseTargetsConditions(targets)
  # factorsAndValues is reversed, this is required later in expand.grid.
  factorsAndValues <- factorsAndValues[ rev( names( factorsAndValues ) ) ]
  # Get the index of the first factor before factor list is reversed.
  mainFactorIndex  <- length( factorsAndValues )
  # Get the number of points of main factor
  nPoints          <- length( factorsAndValues[[ mainFactorIndex ]] )
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Get the combination of all factors other than the main 
  gridPanels <- expand.grid( factorsAndValues[-mainFactorIndex] )
  
  # Subset from all factor combinations those that are present in the targets
  condInTargets <-  unique( apply( targets[ , names( gridPanels ) , drop= FALSE], 1,
          function(x) paste0(x,collapse = "_") ) )
  
  condInGrid <- apply( gridPanels, 1,
      function(x) paste0(x,collapse = "_") )
  
  gridPanels <- gridPanels[ condInGrid %in% condInTargets , , drop=FALSE ]
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Extract data from bin
  binCounts <- .extractCountColumns( countsb( counts )[ bin,] , targets )
  geneLocus <- countsb( counts )[ bin,'locus']
  
  geneCounts <- .extractCountColumns( countsg( counts )[ geneLocus,] , targets )
  
  psir <- joint( as )[ bin, as.character( unique( targets$condition ) ) ]
  psir [ is.na(psir) ] <- 0
  
  J123 <- joint( as )[ bin,  colnames( joint( as ) ) %in%  rownames(targets)]
  ac <- 0
  J1   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J2   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J3   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # creates the graphic device
  outputIsAValidFile <- ( ! is.null( outfile ) ) && 
                        ( file.access( dirname( outfile ) , 2 ) == 0 ) 
                    
  if ( outputIsAValidFile ) {

    devicesNames <- c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')
    devices      <- list( png, bmp, jpeg, tiff, pdf)

    deviceIndex <- match( tolower( outfileType ), devicesNames )
    
    if ( is.na(deviceIndex) ) { 
      dev.new()
      message( paste('Format',outfileType,'is not recognized.',
              'Select png, bmp, jpeg, tiff or pdf') )

    } else {
      device <- devices[[ deviceIndex ]]
      undefDeviceOpt <-  is.null( deviceOpt )
      
      if ( (! is.na( deviceIndex )) && deviceIndex %in%  c(1:4) ) {
        fileNameOpt <- list( filename = outfile)
        defaultDeviceOpt <- list( res = 200, height = 1500, width = 600, 
            pointsize = 12, units = 'px')
      }
      
      if ( (! is.na( deviceIndex )) && deviceIndex ==  5 ) {
        fileNameOpt <- list( file = outfile)
        defaultDeviceOpt <- list( height = 7.5, width = 3, 
            pointsize = 12)        
      }
      
      deviceOpt <- if ( undefDeviceOpt ) 
            append( defaultDeviceOpt, fileNameOpt ) 
          else 
            append( deviceOpt, fileNameOpt )
      do.call( device, deviceOpt )    
      
    }
    
  } else {
    dev.new()
    if ( ! is.null ( outfile ) ) {
      message( paste( "File:",outfile,"cannot be created." ) )
    }
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Draw the plot

  # Define colors. Repeat colors if required
  colors <- rep_len( colors, 6 )
  panelTitleColors <- rep_len( panelTitleColors, 6)

  # Sets plot layout  
  par( mfrow= c( 6 , 1 ) )
  
  # Set inner and outer margins
  par( oma = outerMargins,
       mar = innerMargins )
   
  # If useBarplot is not defined then set it as true if main factor has two or
  # less points to show
  useBarplots <- ( ( ! is.null ( useBarplots ) )  &&  useBarplots ) |
                 ( is.null( useBarplots ) ) && ( ( nPoints <= 2 ) ) 
  
  # Add Bin count data
  data <- .prepareData( binCounts, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Bin counts", col = colors[1], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[1], useHCColors, barWidth, barSpacer )
  
  # Add Gene count data
  data <- .prepareData( geneCounts, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Gene counts", col = colors[2], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[2], useHCColors, barWidth, barSpacer )
  
  # Add PSI/PIR count data
  # There is only one psir value for each condition ( instead for each sample as
  # other profiles being plotted. The value for condition is repeated for each
  # sample in that condition and then processed as the other data.
  psir <- psir[, rep( colnames(psir), as.vector( table(targets$condition)[ unique(targets$condition) ] ))]
  colnames(psir) <- rownames( targets )
  data <- .prepareData( psir , gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "PSI / PIR", col = colors[3], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[3], useHCColors, barWidth, barSpacer )
  
  # Add J1 Exc count data
  data <- .prepareData( J1, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Incl. junc 1", col = colors[4], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[4], useHCColors, barWidth, barSpacer )
  
  # Add J2 Exc count data
  data <- .prepareData( J2, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Incl. junc 2", col = colors[5], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[5], useHCColors, barWidth, barSpacer )
  
  # Add J3 Exc count data
  data <- .prepareData( J3, gridPanels, targets, factorsAndValues )
  .makePlot( data, useBarplots, main = "Exclu. junc", col = colors[6], 
      factorsAndValues, justOneFactor, gridPanels, nPoints, legendAtSide,
      las.x, panelTitleCex, panelTitleColors[6], useHCColors, barWidth, barSpacer )
  
  # Draw the main Title of the plot
  if ( is.null(main) ) { main = bin }
  title( main = bin, outer=TRUE)
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Close graphic device if required
  if ( outputIsAValidFile ) {
    output <- capture.output( dev.off() )
  } 
  # -------------------------------------------------------------------------- #
  
}


#-------------------------------------------------
#Plot splicing pattern
#-------------------------------------------------
# colnames(reportes$col_16_pcp_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_16_pcp_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$col_23_col_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_23_col_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$pcp_23_pcp_16@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$pcp_23_pcp_16@anchorbased)[7] <- "junction.nonuniformity"
# colnames(reportes$col_23_pcp_23@binbased)[21]   <- "junction.nonuniformity"
# colnames(reportes$col_23_pcp_23@anchorbased)[7] <- "junction.nonuniformity"

#region=NULL;exones=NULL;genePlot=TRUE;zoomRegion=1.5;chrMap=NULL;hCov=0.7;hJun=0.3;useLog=FALSE
#bamFiles <- c("col0_16.star.bam","pcp_16.star.bam")
#mergedBAMs <- data.frame(bam=bamFiles,condition=c("col_16","pcp_16"))
#region <- iss$region[1]
#iss <- integrateSignals(sr, asd)
#chrMap <- as.character(1:5)
#names(chrMap) <- paste0("Chr",1:5)
#.plotSplicingPattern(region, iss, counts, f, mergedBAMs, sr, exones, chrMap = chrMap)
.plotSplicingPattern<-function(region=NULL,iss,counts,f,mergedBAMs,sr,exones=NULL,genePlot=TRUE,zoomRegion=1.5,chrMap=NULL,hCov=0.7,hJun=0.3,useLog=FALSE){
  
  greg <- region
  nConditions <- nrow(mergedBAMs)
  
  if(is.null(region)) warning("Out.\n")
  iiss  <- iss[iss$region %in% greg,]
  
  #encuentra x: GRanges para graficar
  if(!genePlot){
    
    #nombre de cromosoma en aspli
    aspli.chr <- strsplit2(iiss$region,":")[1]
    
    #nombre de cromosomas en features
    features.chr<-levels(seqnames(featuresb(f))@values)
    
    if(is.na(match(aspli.chr,features.chr))){
      if(is.null(chrMap)){
        warning(paste("No se pudo mapear nombres de cromosomas.",
                      "\n aspli.chr=",aspli.chr,
                      "\n features.chr=",paste(features.chr,collapse="/")))
      }else{
        chr <- features.chr[match(chrMap[aspli.chr],features.chr)]
      }
    }else{
      chr <- aspli.chr
    }
    
    #si hay una J3 != NA la uso para definir el rango
    if(!is.na(iiss$J3)){
      #si J3 se movio en un cluster, uso el cluster para definir el rango
      #  if(iiss$jl==1){
      #    #TODO
      #  }else{
      roi   <- as.numeric(strsplit2(as.character(iiss$J3),".",fixed=TRUE)[2:3])
      #  }
    }else{
      roi   <- as.numeric(strsplit2(strsplit2(iiss$region,":")[2],"-"))
    }
    delta <- roi[2]-roi[1]
    
    zroi  <- c(roi[1]-delta*zoomRegion/2 , roi[2]+delta*zoomRegion/2)
    zdelta<-zroi[2]-zroi[1]
    
    gr    <- as(paste0(chr,":",zroi[1],"-",zroi[2]),"GRanges")
    hits  <- findOverlaps(gr,featuresb(f))
    bins     <-featuresb(f)[subjectHits(hits)]
    
    
  }else{
    geneName<-iiss[,"locus"]
    
    #identifico coordenadas de biones ebinsonicos
    #iE <- which(mcols(featuresb(f))$locus%in%geneName & mcols(featuresb(f))$feature=="E"
    iE <- which(mcols(featuresb(f))$locus%in%geneName)
    if(length(iE)==0){
      warning("Something terribly wrong happened...No annotation data
              for",geneName,"was found!\n")
      return()
    }
    bins <- featuresb(f)[iE,]
    bins <- bins[mcols(bins)$feature!="Io",]
    
    zroi <- roi <- c(start(bins)[1],end(bins)[length(bins)])
  }
  
  if(is.null(exones)){
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),0.1)
    hh <- hh/sum(hh)
  }else{
    hh <- c(rep(c(hCov,hJun)/nConditions,nConditions),(hCov+hJun)/nConditions)
    hh <- hh/sum(hh)
  }
  layout(matrix((2*nConditions+1):1,(2*nConditions+1),1),height=hh)
  par(mar=c(1, 4.1, 1, 2.1))
  
  
  rownames(mergedBAMs)<-mergedBAMs[,1]
  
  
  #Transcriptos
  if(!is.null(exones)){
    transcriptGene<-strsplit2(names(exones),".",fixed=TRUE)[,1]
    
    ig<-which(transcriptGene%in%geneName)
    tex <- exones[ig]
    tex <- tex[order(names(tex))]
    
    limExones<-unlist(lapply(tex,function(x){
      return(c(start(x),end(x)))
    }))
    
    plot(0,typ="n",xlim=range(limExones),ylim=c(0,length(tex)+3),axes=FALSE,xlab="",ylab="")
    
    if(length(tex)>1){
      for(itex in 1:length(tex)){
        bbins <- tex[[itex]]
        y<-length(tex)-(itex-1)
        lines(c(start(bbins[1]),end(bbins[length(bbins)])),c(y,y),col="gray")
        rect(start(bbins),y-.3,end(bbins),y+0.3,col="white")
        
      }
    }
    ycollapsed <- length(tex)+2
  }else{
    ycollapsed <- 0
  }
  
  
  #bines del genoma anotado
  iE   <- mcols(bins)$feature=="E"
  iroi <- (start(bins)>=roi[1] & start(bins)<roi[2])
  
  nbines <- length(bins)
  delta <- end(bins[nbines])-start(bins[1])
  
  if(is.null(exones)){
    plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(-.5,.5),axes=FALSE,xlab="",ylab="")
  }
  lines(c(start(bins[1]),end(bins[nbines])),c(ycollapsed,ycollapsed),col="gray")
  rect(start(bins)[iE],ycollapsed-.45,end(bins)[iE],ycollapsed+.45,col="white")

  
  #bines diferenciales
  if(FALSE){
    ss <- iss[iss$region==greg,]
    if(ss$b==1){
      if(ss$feature=="I"){
        lines(c(start(bins[ss$bin]),end(bins[ss$bin])),c(0,0),col="orange",lwd=3)
      }else{
        rect(start(bins[ss$bin]),-1,end(bins[ss$bin]),1,col="orange",border=NA)
      }
    }
  }else{
    #que bines diferenciales hay en la region?
    iE<-which(names(bins)%in%iss$bin)
    if(length(iE)>0){
      for(iie in seq_along(iE)){
        if(mcols(bins)$feature[iE[iie]]=="I"){
          lines(c(start(bins[iE[iie]]),end(bins[iE[iie]])),c(ycollapsed,ycollapsed),col="orange",lwd=3)
        }else{
          rect(start(bins[iE[iie]]),ycollapsed-.45,end(bins[iE[iie]]),ycollapsed+.45,col="orange",border=NA)
        }
      }
    }
  }
  
  
  
  
  
  # junturas
  # dibujo todas las junturas de la region de interes zroi,
  # deberia trabajar con jdu...que tiene todas las analizadas, pero por ahora uso anchorage(sr)
  # para diseniar el ploteado
  aspli.chr <- strsplit2(iiss$region,":")[1]
  gr    <- as(paste0(chrMap[aspli.chr],":",zroi[1],"-",zroi[2]),"GRanges")
  
  js <- unique(anchorbased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  ijs <-which(aux[,1]==chrMap[aspli.chr] &
                as.numeric(aux[,2])>=zroi[1] &
                as.numeric(aux[,2])<=zroi[2])
  nj<-0
  jcoords<-c()
  if(length(ijs)>0){
    jcoords          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    rownames(jcoords)<-js[ijs]
    nj               <- nrow(jcoords)
    
    if(FALSE){
      iasr <- match(rownames(jcoords),anchorbased(sr)$junction)
      anchorbased(sr)[iasr,]
    }
  }
  
  
  #analizo si alguna locale no aparece...la agrego
  js <- unique(localebased(sr)$junction)
  aux <- strsplit2(js,".",fixed=TRUE)
  ijs <-which(aux[,1]==chrMap[aspli.chr] &
                as.numeric(aux[,2])>=zroi[1] &
                as.numeric(aux[,2])<=zroi[2])
  if(length(ijs)>0){
    jcoords2          <- matrix(as.numeric(aux[ijs,]),ncol=3)
    rownames(jcoords2)<-js[ijs]
    nj2               <- nrow(jcoords2)
    idiff<-setdiff(rownames(jcoords2),rownames(jcoords))
    if(length(idiff)>0){
      rnames <- c(rownames(jcoords),idiff)
      jcoords<- rbind(jcoords,jcoords2[idiff,])
      rownames(jcoords)<-rnames
      nj     <- nj + length(idiff)
      
      #TODO jclusters
      #ilsr <- match(rownames(jcoords),localebased(sr)$junction)
      #jclus <- unique(localebased(sr)[ilsr,"junction.cluster"])
      
      
    }
  }
  
  if(nj>0){
    abline(v=unique(c(jcoords[,2],jcoords[,3])),col="gray",lty=3)
    
    #paneles
    #cuento junturas por condicion
    x<-countsj(counts)[rownames(jcoords),]
    mjsum <-c()
    for(ifila in 1:nrow(x)){
      jsum<-c()
      for(ccond in counts@condition.order){
        icols<-grep(ccond,colnames(x))
        jsum <- c(jsum,sum(x[ifila,icols]))
      }
      mjsum<-rbind(mjsum,jsum)
    }
    colnames(mjsum)<-counts@condition.order
    rownames(mjsum)<-rownames(x)
  }
  
  
  for(icond in 1:nConditions){
    #ppath <- "/home/ariel/Projects/RNAseq/Marcelo/PCP/03_STAR/mergedBAMS/"
    ad <- system(paste0("samtools depth -r ",
                        paste0(chrMap[aspli.chr],":",zroi[1],"-",zroi[2]," "),
                        bamFiles[icond]), intern = T)
    ad <- matrix(as.numeric(strsplit2(ad,"\t")),ncol=3)
    yylim <- range(ad[,3])
    
    
    
    #panel junturas
    plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=c(1,nj+1),axes=FALSE,xlab="",ylab="")
    if(nj>0){
      jcounts <- mjsum[,mergedBAMs[icond,2]]
      ww <- jcounts/max(jcounts)*3
      
      abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightgray",lty=3)
      for(ij in 1:nj){
        lines(jcoords[ij,2:3],rep(ij,2),lwd=ww[ij])
        points(jcoords[ij,2:3],rep(ij,2),pch=18)
        text(mean(jcoords[ij,2:3]),ij,jcounts[ij],pos=3)
        
      }
    }
    
    #coverage
    if(useLog){
      plot(0.01,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=yylim,axes=FALSE,xlab="",ylab="",log="y")
    }else{
      plot(0,typ="n",xlim=c(start(bins[1]),end(bins[nbines])),ylim=yylim,axes=FALSE,xlab="",ylab="")
    }
    polygon(c(ad[1,2],ad[,2],ad[nrow(ad),2]), c(0.01,ad[,3],0.01)
            ,border=NA,col=topo.colors(nConditions,0.5)[icond],fillOddEven = TRUE)
    lines(c(ad[1,2],ad[,2],ad[nrow(ad),2]),  c(0.01,ad[,3],0.01),
          col=topo.colors(nConditions,0.7)[icond])
    text(ad[1,2],0.01,paste0("[0-",max(ad[,3]),"]"),cex=1.2,adj=c(-.25,-.5))
    
    if(nj>0)abline(v=unique(c(jcoords[,2],jcoords[,3])),col="lightgray",lty=3)
    
    
  }
  mtext(paste(iiss$locus,iiss$region),line=-.5)
  
  
}


#plotSplicingPattern(region,iss,counts,f=features,mergedBAMs,genePlot=TRUE,chrMap=chrMap)
#plotSplicingPattern(region,iss,counts,f=features,mergedBAMs,exones,genePlot=TRUE,chrMap=chrMap)


#for(i in 1:length(iss$region)){
#  plotSplicingPattern(region=iss$region[i],iss,counts,f=features,mergedBAMs,exones,genePlot=TRUE,chrMap=chrMap,useLog=FALSE)
#  readline("Press <ENTER>")
#}







