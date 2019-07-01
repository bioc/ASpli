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
.condenseTargetsConditions <- function ( targets, collapse = "_" ) {
  if( ! "condition" %in% colnames( targets ) ) {
    targets <- data.frame( 
        targets, 
        condition = apply( targets[ , -1 , drop = FALSE] ,1 ,paste, collapse = collapse),
        stringsAsFactors = FALSE)
    
  }
  return( targets )
}

# create the names of samples
.generateSamplesNames <- function ( targets, collapse = "_" ) {
  
  #Auxiliary functions
  is.sequential <- function(x){
    all(diff(x) == diff(x)[1])
  }
  
  my.make.unique <- function(s, sep = "_"){
    tab <- unique(s)
    tab <- setNames(rep(1, length(tab)), tab)
    sapply(s, function(ss){
      sss <- paste(ss, tab[ss], sep = sep)
      tab[ss] <<- tab[ss] + 1
      return(sss)
    })
  }
  
  #Do we need to generate the names or do they already exist? If they dont exists, they are a sequence of numbers from 1 to nrow(targets)
  r <- suppressWarnings(as.numeric(rownames(targets)))
  if(all(!is.na(r))){
    if(is.sequential(r)){
      if( ! "condition" %in% colnames( targets ) ) {
        targets <- .condenseTargetsConditions( targets, collapse )
      }      
      rownames(targets) <- my.make.unique(targets$condition, sep = collapse)      
    }
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

####################################################
##This functions embed one data.table inside another
####################################################
.datatable2 <- function(x, vars = NULL, opts = NULL, ...) {
  
  names_x <- names(x)
  if (is.null(vars)) stop("'vars' must be specified!")
  pos <- match(vars, names_x)
  #if (any(map_chr(x[, pos], typeof) == "list"))
  #  stop("list columns are not supported in datatable2()")
  
  pos <- pos[pos <= ncol(x)] + 1
  rownames(x) <- NULL
  if (nrow(x) > 0) x <- cbind(' ' = '&oplus;', x)
  
  # options
  opts <- c(
    opts, 
    list(
      columnDefs = list(
        list(visible = FALSE, targets = c(0, pos)),
        list(orderable = FALSE, className = 'details-control', targets = 1),
        list(className = 'dt-left', targets = 1:3),
        list(className = 'dt-right', targets = 4:ncol(x))
      )
    ))
  
  datatable(
    x, 
    ...,
    escape = -2,
    options = opts,
    callback = JS(.callback2(x = x, pos = c(0, pos)))
  )
}
####################################################
##This functions embed one data.table inside another
####################################################
.callback2 <- function(x, pos = NULL) {
  
  part1 <- "table.column(1).nodes().to$().css({cursor: 'pointer'});"
  
  part2 <- .child_row_table2(x, pos = pos)
  
  part3 <- 
    "
  table.on('click', 'td.details-control', function() {
  var td = $(this), row = table.row(td.closest('tr'));
  if (row.child.isShown()) {
  row.child.hide();
  td.html('&oplus;');
  } else {
  row.child(format(row.data())).show();
  td.html('&ominus;');
  }
  });"
  
  paste(part1, part2, part3)
} 
####################################################
##This functions embed one data.table inside another
####################################################
.child_row_table2 <- function(x, pos = NULL) {
  
  names_x <- paste0(names(x), ":")
  text <- "
  var format = function(d) {
  text = '<div><table >' + 
  "

  for (i in seq_along(pos)) {
    text <- paste(text, glue::glue(
        "'<tr>' +
        '<td>' + '{names_x[pos[i]]}' + '</td>' +
        '<td>' + d[{pos[i]}] + '</td>' +
        '</tr>' + " ))
  }

  paste0(text,
    "'</table></div>'
    return text;};"
  )
}
