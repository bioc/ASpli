.junctionDUreportExt <- function(
    asd, 
    minAvgCounts              = 5, 
    contrast                  = NULL,
    filterWithContrasted      = FALSE,
    runUniformityTest         = FALSE,
    mergedBams                = NULL,
    maxPValForUniformityCheck = 0.2,
    strongFilter              = FALSE,
    maxConditionsForDispersionEstimate = 24 
  ){
  
  targets <- asd@targets
  
  # Generate conditions combining experimental factors
  if(!"condition" %in% colnames(targets)){
    targets             <- .condenseTargetsConditions(targets)
  }

  #If no contrast provided takes last two and warns
  if(is.null(contrast)){
    contrast <- c(rep(0, times=length(unique(targets$condition))-2), -1, 1)
    warning(paste0("Null contrast povided, using following constrast: ", paste(contrast, collapse=",")))
  }
  
  # Create result object                   
  jdu <- new( Class="ASpliJDU" )
  jdu@contrast <- setNames(contrast, getConditions(targets))
  
  ##############
  #junctionsPJU#
  ##############
  message("Running junctionsPJU test")
  data                  <- junctionsPJU(asd)
  start_J1              <- grep("StartHit", colnames(data)) + 1
  start_J2              <- grep("EndHit", colnames(data)) + 1
  start_J3              <- 9
  end_J3                <- start_J3 + nrow(targets) - 1
  
  junctions_of_interest <- .filterJunctionBySampleWithContrast(data[,start_J3:end_J3], targets=targets, threshold = minAvgCounts, filterWithContrasted, contrast )
  
  J1                    <- as.character(data$StartHit[rownames(data) %in% rownames(junctions_of_interest)])
  J2                    <- as.character(data$EndHit[rownames(data) %in% rownames(junctions_of_interest)])
  J3                    <- rownames(junctions_of_interest)
  clusters              <- .makeClusters(J1, J2, J3, strongFilter)
  
  countData             <- .makeCountDataWithClusters(data[names(clusters$membership),start_J3:end_J3], clusters)
  
  ltsp                  <- .binsJDUWithDiffSplice(countData, targets, contrast, maxConditionsForDispersionEstimate = maxConditionsForDispersionEstimate)
  
  jPSI                  <- ltsp[["junction"]]
  mean.counts           <- rowMeans(countData[rownames(jPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]])
  jPSI$log.mean         <- log2(mean.counts)
  
  mean.counts.per.condition     <- sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(countData[rownames(jPSI), rownames(targets)[targets$condition %in% i]]))})
  participation                 <- aggregate(mean.counts.per.condition ~ jPSI$cluster, FUN = function(r){return(rowSums(t(r)))})
  rownames(participation)       <- participation$cluster
  participation                 <- participation[jPSI$cluster, -1]
  participation                 <- mean.counts.per.condition/participation
  maxparticipation              <- apply(participation, 1, max)
  minparticipation              <- apply(participation, 1, min)
  deltapariticipation           <- maxparticipation - minparticipation 
  participation                 <- maxparticipation
  
  jPSI                  <- jPSI[, c("cluster", "log.mean", "logFC", "P.Value", "FDR")]
  jPSI$annotated        <- ifelse(data[rownames(jPSI), "junction"] != "noHit", "Yes", "No")
  jPSI$participation    <- participation
  jPSI$dParticipation   <- deltapariticipation
  colnames(jPSI)        <- c("cluster", "log.mean", "logFC", "pvalue", "FDR", "annotated", "participation", "dParticipation")
  
  
  jPSI                  <- cbind(jPSI, counts=mean.counts.per.condition)
  
  localej(jdu)          <- jPSI

  ltsp                  <- ltsp[["cluster"]]
  colnames(ltsp)        <- c("size", "cluster.LR", "pvalue", "FDR")
  localec(jdu)          <- ltsp
            
  
  ##############
  #junctionsPIR# 
  ##############
  message("Running junctionsPIR test")
  data                  <- junctionsPIR(asd)
  start_J1              <- 3
  start_J2              <- 3+nrow(targets)
  start_J3              <- 3+2*nrow(targets)
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast, strongFilter)
  countData             <- .makeCountData(Js$J3, Js$J1, Js$J2)
  
  #We reduce data so dispersion estimates can be computed in a razonable ammount of time
  reduxData             <- .makeReduxData(countData, targets, contrast, maxConditionsForDispersionEstimate)
  ltsp                  <- .binsJDUWithDiffSplice(reduxData$countData, reduxData$targets, reduxData$contrast)#, test = "gene")

  jPIR                  <- ltsp[["junction"]]
  #jPIR                  <- jPIR[-(grep("[.][1-2]$", rownames(jPIR))), ] #Saco las junturas que no son J3
  #hist(jPIR$P.Value)
  #jPIR$P.Value          <- p.adjust(jPIR$P.Value, "fdr") #Vuelvo a calcular los fdr solo sobre las junturas J3
  mean.counts           <- rowMeans(countData[rownames(jPIR), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]])
  jPIR$log.mean         <- log2(mean.counts)
  

  #Corremos Uniformity solamente sobre los elementos con menor pvalue
  data_unif             <- data[rownames(jPIR), getConditions(targets)[contrast != 0]]
  data_unif$pvalue      <- jPIR$P.Value
  
  if(runUniformityTest){
    message("Testing uniformity in junctionsPIR")
    jPIR$NonUniformity       <- .testUniformity(data_unif, mergedBams, maxPValForUniformityCheck, targets, contrast)
  }else{
    jPIR$NonUniformity       <- rep(NA, nrow(jPIR)) 
  }

  #Sacamos "cluster" de anchorj 
  jPIR      <- cbind(jPIR, countsJ1 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J1[paste0(rownames(jPIR), ".1"), rownames(targets)[targets$condition %in% i]]))}))
  jPIR      <- cbind(jPIR, countsJ2 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J2[paste0(rownames(jPIR), ".2"), rownames(targets)[targets$condition %in% i]]))}))
  jPIR      <- cbind(jPIR, countsJ3 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J3[rownames(jPIR), rownames(targets)[targets$condition %in% i]]))}))
  participation      <- jPIR[, grep("countsJ3", colnames(jPIR))]/(jPIR[, grep("countsJ1", colnames(jPIR))] + jPIR[, grep("countsJ2", colnames(jPIR))] + jPIR[, grep("countsJ3", colnames(jPIR))])
  jPIR$participation <- apply(participation, 1, max)
  jPIR$dParticipation <- jPIR$participation-apply(participation, 1, min)

  
  jPIR$annotated        <- ifelse(!is.na(data[rownames(jPIR), "hitIntron"]), "Yes", "No")
  jPIR                  <- jPIR[, c("log.mean", "logFC", "P.Value", "FDR", "NonUniformity", "participation", "dParticipation", "annotated", colnames(jPIR)[grep("counts", colnames(jPIR))])]

  colnames(jPIR)        <- c("log.mean", "logFC", "pvalue", "FDR", "NonUniformity", "participation", "dParticipation", "annotated", colnames(jPIR)[grep("counts", colnames(jPIR))])
           
  anchorj(jdu)          <- jPIR
  ltsp                  <- ltsp[["cluster"]][,!colnames(ltsp[["cluster"]])%in%"size"]
  colnames(ltsp)        <- c("cluster.LR", "pvalue", "FDR")
  anchorc(jdu)          <- ltsp
  
  #######
  #irPIR#
  #######
  message("Running irPIR test")
  data                  <- irPIR(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast, strongFilter)
  countData             <- .makeCountData(Js$J3, Js$J1,  Js$J2)
  
  #We reduce data so dispersion estimates can be computed in a razonable ammount of time
  reduxData             <- .makeReduxData(countData, targets, contrast, maxConditionsForDispersionEstimate)
  ltsp                   <- .binsDUWithDiffSplice(reduxData$countData, reduxData$targets, reduxData$contrast)
  tsp <- ltsp[["junction"]]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jirPIR                <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jirPIR$log.mean       <- log2(rowMeans(countData[rownames(jirPIR), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))


  data_unif             <- data[rownames(jirPIR), getConditions(targets)[contrast != 0]]
  rownames(data_unif)   <- data[rownames(jirPIR), "J3"]
  data_unif$pvalue      <- pmin(jirPIR$pvalue,tsp$pvalue)  #corro unif test si bin o juntura son suficientemente significativos
  
  if(runUniformityTest){
    message("Testing uniformity in irPIR")
    jirPIR$NonUniformity     <- .testUniformity(data_unif, mergedBams, maxPValForUniformityCheck, targets, contrast)
  }else{
    jirPIR$NonUniformity     <- rep(NA, nrow(jirPIR))
  }

  
  jirPIR$J3             <- data[rownames(jirPIR), "J3"] 
  jirPIR                <- jirPIR[, c("J3", "logFC", "log.mean", "pvalue", "bin.fdr", "NonUniformity")]
  
  dpir                  <- irPIR(asd)[rownames(jirPIR),getConditions(targets)]
  jirPIR$dPIR           <- apply(dpir,1,function(x){sum(x*contrast)}) 
  colnames(jirPIR)      <- c("J3", "logFC", "log.mean", "pvalue", "FDR", "NonUniformity","dPIR")
  jirPIR$multiplicity   <- "No"
  jirPIR$multiplicity[grep(";", jirPIR$J3)] <- "Yes"
  
  #Sacamos "cluster" de anchorj 
  jirPIR        <- cbind(jirPIR, countsJ1 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J1[paste0(rownames(jirPIR), ".1"), rownames(targets)[targets$condition %in% i]]))}))
  jirPIR        <- cbind(jirPIR, countsJ2 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J2[paste0(rownames(jirPIR), ".2"), rownames(targets)[targets$condition %in% i]]))}))
  jirPIR        <- cbind(jirPIR, countsJ3 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J3[rownames(jirPIR), rownames(targets)[targets$condition %in% i]]))}))
  
  #jirPIR$dPIR   <- apply(jirPIR,1,function(x){sum(x*contrast)}) 
  #participation        <- jirPIR$countsJ3/(jirPIR$countsJ1 + jirPIR$countsJ2 + jirPIR$countsJ3)
  #jirPIR$participation <- apply(participation, 1, max)
  
  jir(jdu)              <- jirPIR
  
  ########
  #ES PSI#
  ########
  message("Running esPIN test")
  data                  <- esPIN(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast, strongFilter)
  
  countData             <- .makeCountData(Js$J3, Js$J1, Js$J2)

  #We reduce data so dispersion estimates can be computed in a razonable ammount of time
  reduxData             <- .makeReduxData(countData, targets, contrast, maxConditionsForDispersionEstimate)
  ltsp                   <- .binsDUWithDiffSplice(reduxData$countData, reduxData$targets, reduxData$contrast)#, test = "gene")
  tsp <- ltsp[["junction"]]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jesPSI                <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jesPSI$log.mean       <- log2(rowMeans(countData[rownames(jesPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))
  jesPSI$event          <- data[rownames(jesPSI), "event"]
  jesPSI$J3             <- data[rownames(jesPSI), "J3"]
  jesPSI                <- jesPSI[, c("event", "J3", "logFC", "log.mean", "pvalue", "bin.fdr")]
  
  dpsi                  <- esPIN(asd)[rownames(jesPSI),getConditions(targets)]
  jesPSI$dPSI           <- apply(dpsi,1,function(x){sum(x*contrast)}) 
  colnames(jesPSI)      <- c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPIN")
  jesPSI$multiplicity        <- "No"
  jesPSI$multiplicity[grep(";", jesPSI$J3)] <- "Yes"
  
  jesPSI      <- cbind(jesPSI, countsJ1 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J1[paste0(rownames(jesPSI), ".1"), rownames(targets)[targets$condition %in% i]]))}))
  jesPSI      <- cbind(jesPSI, countsJ2 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J2[paste0(rownames(jesPSI), ".2"), rownames(targets)[targets$condition %in% i]]))}))
  jesPSI      <- cbind(jesPSI, countsJ3 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J3[rownames(jesPSI), rownames(targets)[targets$condition %in% i]]))}))
  
  jes(jdu)              <- jesPSI
  
    
  #########
  #ALT PSI#
  #########
  message("Running altPIN test")
  data                  <- altPIN(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast, strongFilter, alt = T)
  countData             <- .makeCountData(Js$J3, Js$J1 + Js$J2)
  
  ltsp                   <- .binsDUWithDiffSplice(countData, targets, contrast)#, test = "gene")
  tsp <- ltsp[["junction"]]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jaltPSI               <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jaltPSI$log.mean      <- log2(rowMeans(countData[rownames(jaltPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))  
  jaltPSI$event         <- data[rownames(jaltPSI), "event"]
  jaltPSI$J3            <- data[rownames(jaltPSI), "J3"]
  jaltPSI               <- jaltPSI[, c("event", "J3", "logFC", "log.mean", "pvalue", "bin.fdr")]
  
  
  dpsi                  <- altPIN(asd)[rownames(jaltPSI),getConditions(targets)]
  jaltPSI$dPSI          <- apply(dpsi,1,function(x){sum(x*contrast)}) 
  colnames(jaltPSI)     <- c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPIN")
  jaltPSI$multiplicity        <- "No"
  jaltPSI$multiplicity[grep(";", jaltPSI$J3)] <- "Yes"

  
  jaltPSI      <- cbind(jaltPSI, countsJ1 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J1[paste0(rownames(jaltPSI), ".1"), rownames(targets)[targets$condition %in% i]]))}))
  jaltPSI      <- cbind(jaltPSI, countsJ2 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J2[paste0(rownames(jaltPSI), ".2"), rownames(targets)[targets$condition %in% i]]))}))
  jaltPSI      <- cbind(jaltPSI, countsJ3 = sapply(getConditions(targets)[contrast != 0], function(i){return(rowMeans(Js$J3[rownames(jaltPSI), rownames(targets)[targets$condition %in% i]]))}))

  jalt(jdu)             <- jaltPSI
  
  return(jdu)
}

.filterJunctionBySampleWithContrast <- function( 
  df0,
  targets,
  threshold,
  filterWithContrasted = FALSE,
  contrast = NULL 
) {
  
  # Simple function to compute the minimum of values of two vectors, by position
  # Example:
  # vecMin( c(1,2,3), c(2,0,0) ) -> c(1,0,0)
  vecMin <- function ( a, b ) {
    at <- a < b
    bt <- a >= b
    av <- rep( NA, length( a ) )
    bv <- rep( NA, length( b ) )
    av[ at ] <- a[ at ]
    av[ !at ] <- 0
    bv[ bt ] <- b[ bt ]
    bv[ !bt ] <- 0
    return( av + bv )
  }
  
  
  if(filterWithContrasted){
    if(is.null(contrast)) warning("Constrast should be provided if filterWithContrast=TRUE")
    uniqueConditions <- getConditions(targets)[contrast!=0]
    auxtargets       <- targets[targets$condition%in%uniqueConditions,]
    cropped          <- .extractCountColumns( df0, auxtargets ) 
  }else{
    uniqueConditions <- getConditions(targets)
    auxtargets       <- targets
  }
  cropped          <- .extractCountColumns( df0, auxtargets )
  
  
  # Creates an matrix with Inf ( the neutral element for Min function)
  filter <- matrix( Inf ,
                    ncol = length( uniqueConditions ) ,
                    nrow = nrow( cropped ) )
  
  # Iterates over conditions and over samples of each condition
  # uses filter matrix to compute partial min operations
  for ( i in 1:length( uniqueConditions ) ) {
    byCond <- cropped[ , auxtargets$condition == uniqueConditions[ i ] , drop = FALSE]
    for ( j in 1:ncol( byCond ) ) {
      filter[ , i ] <- vecMin( filter[ , i ] , byCond[ , j ]  )
    }
  }
  
  # Filter the initial dataframe
  filter <- rowSums( filter > threshold ) > 0
  return( df0[ filter,  ] )
  
}

.makeClusters <- function(
  J1, 
  J2, 
  J3, 
  bStrongFilter = FALSE
){
  junctions_of_interest <- !is.na(J3) & J3 != "-"
  J1       <- J1[junctions_of_interest]
  J2       <- J2[junctions_of_interest]
  J3       <- J3[junctions_of_interest]
  
  #JJ3      <- J3[junctions_of_interest]
  JJ3      <- strsplit(J3, ";")
  main_junctions <- unlist(lapply(JJ3, function(s){return(s[1])}))
  nJJ3     <- unlist(lapply(JJ3, length))-1
  JJ3      <- unlist(lapply(JJ3[nJJ3 > 0], function(s){return(s[-1])}))
  resMed   <- cbind(rep(main_junctions, times=nJJ3), JJ3)
  
  junctions_of_interest <- J1 != "-" & !is.na(J1)
  J1       <- J1[junctions_of_interest]
  J1       <- strsplit(J1, ";")
  nJ1      <- unlist(lapply(J1, length))
  J1       <- unlist(J1)
  resStart <- cbind(rep(main_junctions[junctions_of_interest], times=nJ1), J1)
  
  junctions_of_interest <- J2 != "-" & !is.na(J2)
  J2       <- J2[junctions_of_interest]
  J2       <- strsplit(J2, ";")
  nJ2      <- unlist(lapply(J2, length))
  J2       <- unlist(J2)
  resEnd   <- cbind(rep(main_junctions[junctions_of_interest], times=nJ2), J2)
  
  #saco nodos que no pasaron el filtro (y que estaban en el string list starthit o endhit)
  if(bStrongFilter){
    resStart <- resStart[which(resStart[,1] %in% J3 & resStart[,2] %in% J3),]
    resEnd   <- resEnd[which(resEnd[,1] %in% J3 & resEnd[,2] %in% J3),]
  }
  
  g     <- graph_from_edgelist(rbind(resEnd, resStart, resMed),directed = FALSE)
  gclus <- clusters(g)
  
  return(gclus)
  
}

.binsJDUWithDiffSplice <- function( 
  countData,
  targets,
  contrast,
  ignoreExternal = FALSE,
  ignoreIo = TRUE,
  ignoreI = FALSE,
  maxConditionsForDispersionEstimate = 4
) {
  
  # Filter bins
  countData = countData[ ! ignoreExternal | countData$event != "external" ,] 
  countData = countData[ ! ignoreIo | countData$feature != "Io" ,] 
  countData = countData[ ! ignoreI | countData$feature != "I" ,] 
  
  # Define group and contrast
  group <- targets$condition
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts( group )
  
  # make DU analysis
  y <- DGEList( counts = .extractCountColumns( countData, targets ),
                group = factor( group, levels = getConditions(targets), ordered = TRUE) ,
                genes = .extractDataColumns(countData, targets) ) #Must use regular targets in order to effectively remove the target columns
  
  # TODO: Este filtro es muy resctrictivos
  #  keep <- rowSums( cpm( y ) > 1) >= 2
  #  y <- y[ keep, , keep.lib.sizes = FALSE ]
  y <- calcNormFactors( y )
  
  
  # model.matrix sort columns alphabetically if formula has characters instead 
  # of factors. Therefore to preserve order group is converted to ordered factors 
  groupFactor <- factor( group, unique( group ), ordered = TRUE )
  
  design <- model.matrix( ~0 + groupFactor, data = y$samples )
  
  y   <- estimateDisp( y, design )
  fit <- glmFit( y, design )
  ds  <- diffSpliceDGE( fit, contrast = contrast, geneid = "locus", 
                        exonid = NULL, verbose = FALSE )
  tsp <- topSpliceDGE( ds, test = "exon", FDR = 2, number = Inf )
  
  tspg <- topSpliceDGE( ds, test = "gene", FDR = 2, number = Inf )
  
  # make column names equal to the results of DUReport method
  #colnames( tsp )[ match( 'FDR', colnames( tsp )) ] <- 'bin.fdr'
  #colnames( tsp )[ match( 'P.Value', colnames( tsp )) ] <- 'pvalue'
  #tsp$exon.LR <- NULL
  
  tsp <- tsp[,c("locus","logFC","P.Value","FDR")]
  colnames(tsp)[1]<-"cluster"
  
  rownames(tspg) <- tspg$locus
  tspg <- tspg[,c("NExons","gene.LR","P.Value","FDR")]
  colnames(tspg)[1:2] <- c("size","cluster.LR")
  
  return( list(junction=tsp,cluster=tspg) )
  
} 

.makeCountDataWithClusters <- function(
  countData, 
  clusters
){
  countData                <- countData[names(clusters$membership), ]
  countData                <- cbind(length = "", countData)
  countData                <- cbind(end = "", countData)
  countData                <- cbind(start = "", countData)
  countData                <- cbind(gene_coordinates = "", countData)
  countData                <- cbind(symbol =  "", countData)
  countData                <- cbind(locus_overlap = "", countData)
  countData                <- cbind(locus = clusters$membership, countData)
  countData                <- cbind(event = "", countData)
  countData                <- cbind(feature = "", countData)
  countData                <- countData[order(rownames(countData)), ]
  return(countData)
}

.makeJunctions <- function(
  data,
  targets,
  start_J1,
  start_J2,
  start_J3,
  minAvgCounts,
  filterWithContrasted = FALSE,
  contrast = NULL,
  strongFilter = TRUE,
  alt = FALSE
){
  
  #Pasamos todos los NA a 0
  data[, start_J1:(start_J1 + nrow(targets) - 1)][is.na(data[, start_J1:(start_J1 + nrow(targets) - 1)])] <- 0
  data[, start_J2:(start_J2 + nrow(targets) - 1)][is.na(data[, start_J2:(start_J2 + nrow(targets) - 1)])] <- 0
  data[, start_J3:(start_J3 + nrow(targets) - 1)][is.na(data[, start_J3:(start_J3 + nrow(targets) - 1)])] <- 0
  
  J1 <- data[, start_J1:(start_J1 + nrow(targets) - 1)]
  J2 <- data[, start_J2:(start_J2 + nrow(targets) - 1)]
  J3 <- data[, start_J3:(start_J3 + nrow(targets) - 1)]  
  
  #reliables <- rep(FALSE, times=nrow(J3))
  #for(condition in unique(targets$condition)[contrast != 0]){
  #  reliables <- reliables | (rowMeans(J1[, grep(condition, colnames(J1))]) > minAvgCounts &
  #                              rowMeans(J2[, grep(condition, colnames(J2))]) > minAvgCounts &
  #                              rowMeans(J3[, grep(condition, colnames(J3))]) > minAvgCounts)
  #}
  reliables <- list()
  reliables[["J1"]] <- rownames(.filterJunctionBySampleWithContrast( J1, targets, minAvgCounts, filterWithContrasted = filterWithContrasted, contrast ))
  reliables[["J2"]] <- rownames(.filterJunctionBySampleWithContrast( J2, targets, minAvgCounts, filterWithContrasted = filterWithContrasted, contrast ))
  reliables[["J3"]] <- rownames(.filterJunctionBySampleWithContrast( J3, targets, minAvgCounts, filterWithContrasted = filterWithContrasted, contrast ))
  #No queremos filtrar todas las junturas de alt, solo que tenga o J1 o J2
  if(!alt){
    reliables <- Reduce(intersect, reliables)
  }else{
    reliables <- c(intersect(reliables[["J1"]], reliables[["J3"]]), 
                   intersect(reliables[["J2"]], reliables[["J3"]]))
  }
  #for(condition in unique(targets$condition)[contrast != 0]){
  #  reliables <- reliables | rowMeans(J3[, grep(condition, colnames(J3))]) > minAvgCounts
  #}  
  #reliables <- reliables & apply(data, 1, function(s){return(!any(is.na(s)))})
  rownames(J1) <- paste0(rownames(J3), ".1") 
  rownames(J2) <- paste0(rownames(J3), ".2") 
  J1 <- J1[paste0(reliables, ".1"), ]
  J2 <- J2[paste0(reliables, ".2"), ]
  J3 <- J3[reliables, ]
  return(list(J1 = J1, J2 = J2, J3 = J3))
}

.makeCountData <- function(
  J3,
  Jaux1,
  Jaux2 = NULL
){
  reps                 <- 2
  countData            <- rbind(J3, Jaux1)
  if(!is.null(Jaux2)){
    countData            <- rbind(countData, Jaux2)
    reps         <- 3
  }
  countData            <- cbind(length = "", countData)
  countData            <- cbind(end = "", countData)
  countData            <- cbind(start = "", countData)
  countData            <- cbind(gene_coordinates = "", countData)
  countData            <- cbind(symbol =  "", countData)
  countData            <- cbind(locus_overlap = "", countData)
  countData            <- cbind(locus = rep(rownames(J3), times=reps), countData)
  countData            <- cbind(event = "", countData)
  countData            <- cbind(feature = "", countData)
  countData            <- countData[order(rownames(countData)), ]
  return(countData)
}

.testUniformity <- function(
  data, 
  mergedBams,
  maxPValForUniformityCheck,
  targets,
  contrast
){
  
  if(is.null(mergedBams)){
    warning("Can't run uniformity test without merged bams.")
    return(rep(NA, length=nrow(data))) 
  }
  if(class(mergedBams)=="data.frame"){
    mergedBams <- as.character(mergedBams[,1])
  }
  if(!is.character(mergedBams)){
    warning("Can't run uniformity test. Problems with merged bams specifications.")
    return(rep(NA, length=nrow(data))) 
  }
  
  if(!all(file.exists(c(mergedBams, paste0(mergedBams, ".bai"))))){
    warning("Can't run uniformity test. Couldn't find merged bams or their index.")
    return(rep(NA, length=nrow(data))) 
  }

  ii                   <- rownames(data)[data$pvalue < maxPValForUniformityCheck]
  Uniformity           <- rep(NA, length=nrow(data))
  names(Uniformity)    <- rownames(data)
  
  # +     pb <- txtProgressBar(style=3)
  # +     setTxtProgressBar(pb, i)
  # +     close(pb
  i = 0
  pb <- txtProgressBar(min=1,max=length(ii),style=3)
  Uniformity[ii] <- sapply(ii, function(p){
    #if(i %% 100 == 0) print(paste(i/length(elementos), "%"))
    i <<- i + 1
    setTxtProgressBar(pb,i)
    # if(i %% 10 == 0){
    #   message(paste0("Uniformity test: ", round(i/length(ii)*100), "% completed"))
    # }
    
    #reader <- bamReader(mergedBams[2], idx=TRUE)
    #reader <- readers[[which.max(data[p, getConditions(targets)])]]
    
    pp      <- strsplit2(p, "[.]")[1, ]
    ad <- system(paste0("samtools depth -r ", paste0(pp[1], ":", (as.numeric(pp[2])-10), "-",  (as.numeric(pp[3])+10)), " ", mergedBams[which.max(data[p, getConditions(targets)[contrast != 0]])]), intern = T)
    if(length(ad) > 40){
      ad <- as.numeric(strsplit2(ad, "\t")[, 3])
      #if(refid[pp[1]] < 5){
      #for(i in 1:100){
      #reader  <- bamReader(mergedBams[which.max(data[p, getConditions(targets)])], idx=TRUE)
      #brange  <- bamRange(reader,c(refid[pp[1]],(as.numeric(pp[2])-10), (as.numeric(pp[3])+10)), complex = FALSE)
      #brange  <- bamRange(reader,c(0, 10, 10000))
      #ad      <- alignDepth(brange)
      #bamClose(reader)
      #}
      
      b      <- sd(ad[11:(length(ad)-12)])
      a      <- mean(c(ad[1:10], ad[(length(ad)-11):length(ad)]))
      return(b/a)
    }else{
      return(NA)
    }
    #bamClose(reader)
  })
  close(pb)
  
  #for(i in 1:length(readers)){
  #  bamClose(readers[[i]])
  #}
  
  return(Uniformity)
  
}

#Para que el test corra en un tiempo razonable en datasets muy grandes, evaluamos en un subconjunto
.makeReduxData <- function(countData, targets, contrast, maxConditionsForDispersionEstimate){
  group                 <- targets$condition 
  ugroup                <- unique(group)
  names(contrast)       <- ugroup
  contrast_groups       <- ugroup[contrast != 0]
  other_groups          <- ugroup[contrast == 0]
  other_groups          <- sample(other_groups, min(length(other_groups), abs(maxConditionsForDispersionEstimate-length(contrast_groups))))
  final_groups          <- c(contrast_groups, other_groups)
  redux_targets         <- targets[targets$condition %in% final_groups, ]
  contrast              <- contrast[final_groups]
  targets               <- targets[targets$condition %in% redux_targets$condition, ]
  
  countData             <- cbind(.extractDataColumns(countData, targets), .extractCountColumns( countData, redux_targets ))
  reduxData             <- list(countData = countData, targets = targets, contrast = contrast)
  return(reduxData)
}

