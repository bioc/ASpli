# This file contains a collections of functions used ubiquitously throughout 
# the code.

# Replace all illegal chars in filenames to underscore chars.
# The filename provided should not have any folder name.
.makeValidFileName <- function( filename ) {
  filename <- gsub( '[<>:"/\\|?*]', '_', filename )
  return( filename )
}

# Subset a dataframe to contain only the columns that matchs the samples in the
# targets. That columns are the ones with the count data.
.extractCountColumns <- function ( aDataframe, targets ) {
  result <- aDataframe[ , match( row.names(targets), colnames( aDataframe ) ) , F]
  colnames( result ) <- as.character( row.names(targets) )
  return( result )
}

# Subset a dataframe to contain only the columns that do not matchs the samples
# in the targets. That columns are the ones that do not contain count data.
.extractDataColumns <- function ( aDataframe, targets ) {
  result <- aDataframe[ , - match( row.names(targets) , colnames( aDataframe ) ) ]
  return( result )
}

# create the names of the conditions of a targets by their factors
.condenseTargetsConditions <- function ( targets ) {
  if( ! "condition" %in% colnames( targets ) ) {
    targets <- data.frame( 
        targets, 
        condition = apply( targets[ , -1 , drop = FALSE] ,1 ,paste,collapse="_"),
        stringsAsFactors = FALSE)
    
  }
  return( targets )
}

# This function sums counts of a data frame by condition.
# The conditions are given in the targets data.frame.
# the dataframe to be summed must have the same number of columns as samples 
# in the targets, and they must have the same order.
.sumByCond <- function( countDf, targets ) {
  countDf[ is.na( countDf )] <- 0
  uniqueConditions <- unique( targets$condition )
  nConditions <- length( uniqueConditions )
  result <- matrix( 
      data = 0, 
      nrow = nrow( countDf) , 
      ncol = nConditions )
  
  for( i in 1:nConditions ) {
    result[ , i ] <- rowSums( countDf[ , targets$condition == uniqueConditions[i], drop = FALSE ] )
  }
  colnames( result ) <- uniqueConditions
  return ( result )
}
