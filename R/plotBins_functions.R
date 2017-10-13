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










