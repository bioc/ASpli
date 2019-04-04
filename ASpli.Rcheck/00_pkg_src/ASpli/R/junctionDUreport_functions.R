.junctionDUreportExt <- function(
    asd, 
    targets, 
    minAvgCounts              = 5, 
    contrast                  = NULL,
    filterWithContrasted      = FALSE,
    runUniformityTest         = FALSE,
    mergedBams                = NULL,
    maxPValForUniformityCheck = 0.2
  ){
  
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
  #junctionsPSI#
  ##############
  message("Runing junctionsPSI test")
  data                  <- junctionsPSI(asd)
  start_J1              <- grep("StartHit", colnames(data)) + 1
  start_J2              <- grep("EndHit", colnames(data)) + 1
  start_J3              <- 9
  end_J3                <- start_J3 + nrow(targets) - 1
  
  junctions_of_interest <- .filterJunctionBySampleWithContrast(data[,start_J3:end_J3], targets=targets, threshold = minAvgCounts, filterWithContrasted, contrast )
  
  J1                    <- as.character(data$StartHit[rownames(data) %in% rownames(junctions_of_interest)])
  J2                    <- as.character(data$EndHit[rownames(data) %in% rownames(junctions_of_interest)])
  J3                    <- rownames(junctions_of_interest)
  clusters              <- .makeClusters(J1, J2, J3)
  
  countData             <- .makeCountDataWithClusters(data[names(clusters$membership),start_J3:end_J3], clusters)
  
  ltsp                  <- .binsJDUWithDiffSplice(countData, targets, contrast)
  
  jPSI                  <- ltsp[["junction"]]
  jPSI$log.mean         <- log2(rowMeans(countData[rownames(jPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))
  
  jPSI                  <- jPSI[, c("cluster", "log.mean", "logFC", "P.Value", "FDR")]
  colnames(jPSI)        <- c("cluster", "log.mean", "logFC", "pvalue", "FDR")
  localej(jdu)          <- jPSI

  ltsp                  <- ltsp[["cluster"]]
  colnames(ltsp)        <- c("size", "cluster.LR", "pvalue", "FDR")
  localec(jdu)          <- ltsp
            
  
  ##############
  #junctionsPIR# 
  ##############
  message("Runing junctionsPIR test")
  data                  <- junctionsPIR(asd)
  start_J1              <- 3
  start_J2              <- 3+nrow(targets)
  start_J3              <- 3+2*nrow(targets)
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast)
  countData             <- .makeCountData(Js$J3, Js$J1, Js$J2)
  
  ltsp                  <- .binsJDUWithDiffSplice(countData, targets, contrast)

  jPIR                  <- ltsp[["junction"]]
  jPIR                  <- jPIR[-(grep("[.][1-2]$", rownames(jPIR))), ] #Saco las junturas que no son J3
  jPIR$P.Value          <- p.adjust(jPIR$P.Value, "fdr") #Vuelvo a calcular los fdr solo sobre las junturas J3
  jPIR$log.mean         <- log2(rowMeans(countData[rownames(jPIR), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))
  
  #Corremos Uniformity solamente sobre los elementos con menor pvalue
  data_unif             <- data[rownames(jPIR), getConditions(targets)]
  data_unif$pvalue      <- jPIR$P.Value
  
  if(runUniformityTest){
    message("Testing uniformity in junctionsPIR")
    jPIR$Uniformity       <- .testUniformity(data_unif, mergedBams, maxPValForUniformityCheck, targets)
  }else{
    jPIR$Uniformity       <- rep(NA, nrow(jPIR)) 
  }

  #Sacamos "cluster" de anchorj 
  jPIR                  <- jPIR[, c("log.mean", "logFC", "P.Value", "FDR", "Uniformity")]
  colnames(jPIR)        <- c("log.mean", "logFC", "pvalue", "FDR", "Uniformity")
  anchorj(jdu)          <- jPIR
  ltsp                  <- ltsp[["cluster"]][,!colnames(ltsp[["cluster"]])%in%"size"]
  colnames(ltsp)        <- c("cluster.LR", "pvalue", "FDR")
  anchorc(jdu)          <- ltsp
  
  #######
  #irPIR#
  #######
  message("Runing irPIR test")
  data                  <- irPIR(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast)
  countData             <- .makeCountData(Js$J3, Js$J1,  Js$J2)
  
  tsp                   <- .binsDUWithDiffSplice(countData, targets, contrast)
  tsp                   <- tsp[-(grep("[.][1-2]$", rownames(tsp))), ]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jirPIR                <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jirPIR$log.mean       <- log2(rowMeans(countData[rownames(jirPIR), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))


  data_unif             <- data[rownames(jirPIR), getConditions(targets)]
  rownames(data_unif)   <- data[rownames(jirPIR), "J3"]
  data_unif$pvalue      <- jirPIR$pvalue
  
  if(runUniformityTest){
    message("Testing uniformity in irPIR")
    jirPIR$Uniformity     <- .testUniformity(data_unif, mergedBams, maxPValForUniformityCheck, targets)
  }else{
    jirPIR$Uniformity     <- rep(NA, nrow(jirPIR))
  }

  jirPIR$J3             <- data[rownames(jirPIR), "J3"] 
  jirPIR                <- jirPIR[, c("J3", "logFC", "log.mean", "pvalue", "bin.fdr", "Uniformity")]
  
  dpir                  <- irPIR(asd)[rownames(jirPIR),getConditions(targets)]
  jirPIR$dPIR           <- apply(dpir,1,function(x){sum(x*contrast)}) 
  colnames(jirPIR)      <- c("J3", "logFC", "log.mean", "pvalue", "FDR", "Uniformity","dPIR")
  jir(jdu)              <- jirPIR
  
  ########
  #ES PSI#
  ########
  message("Runing esPSI test")
  data                  <- esPSI(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast)
  
  countData             <- .makeCountData(Js$J3, Js$J1, Js$J2)
  
  tsp                   <- .binsDUWithDiffSplice(countData, targets, contrast)
  tsp                   <- tsp[-(grep("[.][1-2]$", rownames(tsp))), ]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jesPSI                <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jesPSI$log.mean       <- log2(rowMeans(countData[rownames(jesPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))
  jesPSI$event          <- data[rownames(jesPSI), "event"]
  jesPSI$J3             <- data[rownames(jesPSI), "J3"]
  jesPSI                <- jesPSI[, c("event", "J3", "logFC", "log.mean", "pvalue", "bin.fdr")]
  
  dpsi                  <- esPSI(asd)[rownames(jesPSI),getConditions(targets)]
  jesPSI$dPSI           <- apply(dpsi,1,function(x){sum(x*contrast)}) 
  colnames(jesPSI)      <- c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPSI")
  jes(jdu)              <- jesPSI
  
    
  #########
  #ALT PSI#
  #########
  message("Runing altPSI test")
  data                  <- altPSI(asd)
  start_J1              <- grep("J1", colnames(data)) + 1
  start_J2              <- grep("J2", colnames(data)) + 1
  start_J3              <- grep("J3", colnames(data)) + 1
  
  data                  <- data[!is.na(data$J3), ]
  
  Js                    <- .makeJunctions(data, targets, start_J1, start_J2, start_J3, minAvgCounts, filterWithContrasted, contrast)
  countData             <- .makeCountData(Js$J3, Js$J1, Js$J2)
  
  tsp                   <- .binsDUWithDiffSplice(countData, targets, contrast)
  tsp                   <- tsp[-(grep("[.][1-2]$", rownames(tsp))), ]
  tsp$bin.fdr           <- p.adjust(tsp$pvalue, "fdr")
  
  jaltPSI               <- tsp[, c("logFC", "pvalue", "bin.fdr")]
  jaltPSI$log.mean      <- log2(rowMeans(countData[rownames(jaltPSI), rownames(targets)[targets$condition %in% getConditions(targets)[contrast != 0]]]))  
  jaltPSI$event         <- data[rownames(jaltPSI), "event"]
  jaltPSI$J3            <- data[rownames(jaltPSI), "J3"]
  jaltPSI               <- jaltPSI[, c("event", "J3", "logFC", "log.mean", "pvalue", "bin.fdr")]
  
  
  dpsi                  <- altPSI(asd)[rownames(jaltPSI),getConditions(targets)]
  jaltPSI$dPSI          <- apply(dpsi,1,function(x){sum(x*contrast)}) 
  colnames(jaltPSI)     <- c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPSI")
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
  ignoreI = FALSE
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
                genes = .extractDataColumns(countData, targets) )       
  
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
  contrast = NULL
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
  reliables <- .filterJunctionBySampleWithContrast( J3, targets, minAvgCounts, filterWithContrasted = filterWithContrasted, contrast )
  reliables <- rownames(reliables)
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
  targets
){
  
  if(is.null(mergedBams)){
    warning("Can't run uniformity test without merged bams.")
    return(rep(NA, length=nrow(data))) 
  }
  if(!all(file.exists(mergedBams))){
    warning("Couldn't find merged bams. Can't run uniformity test.")
    return(rep(NA, length=nrow(data))) 
  }

  ii                    <- rownames(data)[data$pvalue < maxPValForUniformityCheck]
  Uniformity           <- rep(NA, length=nrow(data))
  names(Uniformity)    <- rownames(data)
  
  i = 0
  Uniformity[ii] <- sapply(ii, function(p){
    #if(i %% 100 == 0) print(paste(i/length(elementos), "%"))
    i <<- i + 1
    print(i/length(ii))
    #reader <- bamReader(mergedBams[2], idx=TRUE)
    
    #reader <- readers[[which.max(data[p, getConditions(targets)])]]
    pp      <- strsplit2(p, "[.]")[1, ]
    ad <- system(paste0("samtools depth -r ", paste0(pp[1], ":", (as.numeric(pp[2])-10), "-",  (as.numeric(pp[3])+10)), " ", mergedBams[which.max(data[p, getConditions(targets)])]), intern = T)
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
    #}else{
    #  return(NA)
    #}
    #bamClose(reader)
  })
  
  #for(i in 1:length(readers)){
  #  bamClose(readers[[i]])
  #}
  
  return(Uniformity)
  
}

.mergeReports <- function(bdu, jdu){
  
  #Sanity check
  if(class(bdu) != "ASpliDU"){
    stop("bdu must be an ASpliDU object, try running binDUreport first") 
  }
  if(class(jdu) != "ASpliJDU"){
    stop("jdu must be an ASpliJDU object, try running junctionDUreport first") 
  }
  
  mr <- new( Class="ASpliMergedReports" )
  
  #Mergin bin based junctions
  jir                    <- jir(jdu)
  jir$event              <- NA
  jir$dPSI               <- NA
  jir                    <- jir[, c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPSI", "Uniformity")]
  colnames(jir)[c(6, 9)] <- c("junction.fdr", "uniformity")
  #jir$bin                <- rownames(jir)
  
  jes                    <- jes(jdu)
  jes$uniformity         <- NA
  jes$dPIR               <- NA
  jes                    <- jes[, c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPSI", "uniformity")]
  colnames(jes)[6]       <- "junction.fdr"
  #jes$bin                <- rownames(jes)
  
  jalt                   <- jalt(jdu)
  jalt$uniformity        <- NA
  jalt$dPIR              <- NA
  jalt                   <- jalt[, c("event", "J3", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPSI", "uniformity")]
  colnames(jalt)[6]      <- "junction.fdr"
  #jalt$bin               <- rownames(jalt)
  
  j                      <- data.table(rbind(jir, jes, jalt), keep.rownames = T)
  bins                   <- data.table(binsDU(bdu), keep.rownames = T)
  aux                    <- merge(bins, j, by="rn", all=T)
  colnames(aux)          <- c("bin", "feature", "bin.event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start",
                              "end", "length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction.event", "junction", "junction.logFC",
                              "junction.log.mean", "junction.pvalue", "junction.fdr", "junction.dPIR", "junction.dPSI", "junction.uniformity")
  binbased(mr)           <- aux
  
  ########################################
  localej   <- localej(jdu)
  junctions <- strsplit2(rownames(localej), "[.]")
  seqnames  <- junctions[, 1]
  start     <- as.numeric(junctions[, 2])
  end       <- as.numeric(junctions[, 3])
  
  grjunctions <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames <- strsplit2(binsDU(bdu)$gene_coordinates, "[:]")[, 1]
  start    <- binsDU(bdu)$start
  end      <- binsDU(bdu)$end
  
  grbines      <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap      <- data.frame(findOverlaps(grbines, grjunctions))
  
  bins <- binsDU(bdu)
  colnames(bins)[10] <- "bin.logFC"
  
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(localej(jdu))[overlap$subjectHits], localej(jdu)[overlap$subjectHits, ])
  
  bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  colnames(bins)[1]      <- "bin"
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(localej(jdu)[!rownames(localej(jdu)) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[1]  <- "J3"
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", "junction.cluster", "junction.log.mean",
                                         "junction.logFC", "junction.pvalue", "junction.fdr")
  
  junctionbased_junctions$junction.cluster <- as.character(junctionbased_junctions$junction.cluster)
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  localec <- data.table(localec(jdu), keep.rownames = T)
  colnames(localec) <- c("rn", "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr")
  fulldt <- merge(junctionbased_junctions, localec, by.x = "junction.cluster", by.y = "rn", all = T)
  
  fulldt <- fulldt[, 
                   c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", 
                     "bin.start", "bin.end",
                     "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", 
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", "junction.cluster", 
                     "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr")
                   ]
  
  localebased(mr) <- fulldt[order(fulldt$bin.fdr, fulldt$junction.fdr, fulldt$cluster.fdr), ]
  
  
  ########################################
  anchorj  <- anchorj(jdu)
  junctions <- strsplit2(rownames(anchorj), "[.]")
  seqnames <- junctions[, 1]
  start    <- as.numeric(junctions[, 2])
  end      <- as.numeric(junctions[, 3])
  
  grjunctions <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames <- strsplit2(bdu@bins$gene_coordinates, "[:]")[, 1]
  start    <- binsDU(bdu)$start
  end      <- binsDU(bdu)$end
  
  grbines      <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap      <- data.frame(findOverlaps(grbines, grjunctions))
  
  bins <- binsDU(bdu)
  colnames(bins)[10] <- "bin.logFC"
  
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(anchorj(jdu))[overlap$subjectHits], anchorj(jdu)[overlap$subjectHits, ])
  
  bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  colnames(bins)[1]      <- "bin"
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(anchorj(jdu)[!rownames(anchorj(jdu)) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[1]  <- "J3"
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", "junction.log.mean",
                                         "junction.logFC", "junction.pvalue", "junction.fdr", "junction.uniformity")
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  anchorc <- data.table(anchorc(jdu), keep.rownames = T)
  colnames(anchorc) <- c("rn", "cluster.LR", "cluster.pvalue", "cluster.fdr")
  fulldt <- merge(junctionbased_junctions, anchorc, by.x = "junction", by.y = "rn", all = T)
  
  fulldt <- fulldt[, 
                   c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", 
                     "bin.start", "bin.end",
                     "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", 
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", "junction.uniformity", 
                     "cluster.LR", "cluster.pvalue", "cluster.fdr")
                   ]
  
  anchorbased(mr) <- fulldt[order(fulldt$bin.fdr, fulldt$junction.fdr, fulldt$cluster.fdr), ]  
  return(mr)
}
