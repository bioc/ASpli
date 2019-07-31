#################################################################################
#Splicing Report
#
#
#
#This function integrates differential usage information from different sources.
#Receives bin usage and junction usage and merges them in three steps.
#First merges all bin based junction information (jir, jes, jalt).
#Then merges locale information with coverage information in the regions.
#Finally does the same with anchors.
#################################################################################
.splicingReport <- function(bdu, jdu){
  
  #Sanity check
  if(class(bdu) != "ASpliDU"){
    stop("bdu must be an ASpliDU object, try running gbDUreport first") 
  }
  if(class(jdu) != "ASpliJDU"){
    stop("jdu must be an ASpliJDU object, try running jDUreport first") 
  }
  
  mr <- new( Class="ASpliSplicingReport" )
  mr@contrast <- jdu@contrast
  
  #Mergin bin based junctions
  jir                    <- jir(jdu)
  jir$event              <- NA
  jir$dPSI               <- NA
  jir                    <- jir[, c("event", "J3", "multiplicity", "logFC", "log.mean", "LR", "pvalue", "FDR", "dPIR", "dPSI", "NonUniformity", colnames(jir)[grep("counts", colnames(jir))])]
  colnames(jir)[colnames(jir) %in% c("FDR", "NonUniformity")] <- c("junction.fdr", "nonuniformity")
  #jir$bin                <- rownames(jir)
  
  jes                    <- jes(jdu)
  jes$nonuniformity      <- NA
  jes$dPIR               <- NA
  jes                    <- jes[, c("event", "J3", "multiplicity", "logFC", "log.mean", "LR", "pvalue", "FDR", "dPIR", "dPSI", "nonuniformity", colnames(jes)[grep("counts", colnames(jes))])]
  colnames(jes)[colnames(jes) %in% c("FDR")] <- c("junction.fdr")
  #jes$bin                <- rownames(jes)
  
  jalt                   <- jalt(jdu)
  jalt$nonuniformity     <- NA
  jalt$dPIR              <- NA
  jalt                   <- jalt[, c("event", "J3", "multiplicity", "logFC", "log.mean", "LR", "pvalue", "FDR", "dPIR", "dPSI", "nonuniformity", colnames(jalt)[grep("counts", colnames(jalt))])]
  colnames(jalt)[colnames(jalt) %in% c("FDR")] <- c("junction.fdr")
  #jalt$bin               <- rownames(jalt)
  
  j                      <- data.table(rbind(jir, jes, jalt), keep.rownames = T)
  bins                   <- data.table(binsDU(bdu), keep.rownames = T)
  aux                    <- data.frame(merge(bins, j, by="rn", all=T))
  colnames(aux)          <- c("bin", "feature", "bin.event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start",
                              "end", "length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction.event", "J3", "J3.multiplicity", "junction.logFC",
                              "junction.log.mean", "junction.LR", "junction.pvalue", "junction.fdr", "junction.dPIR", "junction.dPSI", "junction.nonuniformity", colnames(aux)[grep("counts", colnames(aux))])
  
  aux                    <- aux[, !colnames(aux) %in% c("symbol", "junction.event")]
  
  #aux                    <- aux[union(which(aux$bin.fdr < maxBinFDR), which(aux$junction.fdr < maxJunctionFDR)), ]
  aux                    <- aux[order(aux$bin), ]
  rownames(aux)          <- NULL
  binbased(mr)           <- aux
  
  ########################################
  localej   <- localej(jdu)
  junctions <- strsplit2(rownames(localej), "[.]")
  seqnames  <- junctions[, 1]
  start     <- as.numeric(junctions[, 2]) + 1
  end       <- as.numeric(junctions[, 3]) - 1
  
  grjunctions <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames <- strsplit2(binsDU(bdu)$gene_coordinates, "[:]")[, 1]
  start    <- binsDU(bdu)$start
  end      <- binsDU(bdu)$end
  
  grbines      <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap      <- data.frame(findOverlaps(grbines, grjunctions))
  
  bins <- binsDU(bdu)
  colnames(bins)[colnames(bins) %in% c("logFC")] <- "bin.logFC"
  
  #Unique identifier por locale pvalue
  colnames(localej)[colnames(localej) %in% c("pvalue")] <- "junction.pvalue"
  
  #Merge overlaping bins and locales
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(localej)[overlap$subjectHits], localej[overlap$subjectHits, ])
  
  
  #bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #colnames(bins)[1]      <- "bin"
  
  #junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(localej[!rownames(localej) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[colnames(junctions) %in% c("rn")]  <- c("J3")
  
  #Adds non overlaping locales
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr","junction", "junction.cluster", "junction.log.mean",
                                         "junction.logFC", "junction.pvalue", "junction.fdr", "junction.annotated", "junction.participation", "junction.dparticipation", colnames(junctionbased_junctions)[grep("counts", colnames(junctionbased_junctions))])
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  #Change type of junction.cluster so we can merge with cluster information
  junctionbased_junctions$junction.cluster <- as.character(junctionbased_junctions$junction.cluster)
  

  localec <- data.table(localec(jdu), keep.rownames = T)

  #Tries to find cluster locus
  junctions  <- strsplit2(localec$range, "[.]")
  seqnames   <- junctions[, 1]
  start      <- as.numeric(junctions[, 2]) + 1
  end        <- as.numeric(junctions[, 3]) - 1
  
  grclusters <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames   <- strsplit2(genesDE(bdu)$gene_coordinates, "[:]")[, 1]
  start      <- genesDE(bdu)$start
  end        <- genesDE(bdu)$end
  
  grgenes    <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap    <- data.frame(findOverlaps(grclusters, grgenes))
  localec$locus <- NA
  localec$locus[overlap$queryHits] <- as.character(genesDE(bdu)$symbol[overlap$subjectHits])
  colnames(localec) <- c("rn", "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr", "cluster.range", "cluster.participation", "cluster.dparticipation", "cluster.locus")
  
  fulldt <- data.frame(merge(junctionbased_junctions, localec, by.x = "junction.cluster", by.y = "rn", all = T))
  
  fulldt <- fulldt[, 
                   c("junction", "junction.annotated",
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", 
                     colnames(fulldt)[grep("counts", colnames(fulldt))],
                     "junction.cluster", "junction.participation", "junction.dparticipation",
                     "cluster.locus", "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr", "cluster.range", "cluster.participation", "cluster.dparticipation", 
                     "bin", "bin.pvalue", "bin.fdr")
                   ]
  rownames(fulldt) <- NULL
  #fulldt           <- fulldt[union(which(fulldt$bin.fdr < maxBinFDR), which(fulldt$cluster.fdr < maxJunctionFDR)), ] #Traemos todos
  index.orden      <- strsplit2(fulldt$junction, "[.]")
  index.orden      <- order(GRanges(index.orden[, 1], IRanges(as.numeric(index.orden[, 2]), as.numeric(index.orden[, 3]))))
  fulldt           <- fulldt[index.orden, ] 
  rownames(fulldt) <- NULL
  localebased(mr)  <- fulldt
  
  ########################################
  anchorj  <- anchorj(jdu)
  junctions <- strsplit2(rownames(anchorj), "[.]")
  seqnames <- junctions[, 1]
  start    <- as.numeric(junctions[, 2]) + 1
  end      <- as.numeric(junctions[, 3]) - 1
  
  grjunctions <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames <- strsplit2(bdu@bins$gene_coordinates, "[:]")[, 1]
  start    <- binsDU(bdu)$start
  end      <- binsDU(bdu)$end
  
  grbines      <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap      <- data.frame(findOverlaps(grbines, grjunctions))
  
  bins <- binsDU(bdu)
  colnames(bins)[colnames(bins) %in% c("logFC")] <- "bin.logFC"

  colnames(anchorj)[colnames(anchorj) %in% c("pvalue")] <- "junction.pvalue"  
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(anchorj)[overlap$subjectHits], anchorj[overlap$subjectHits, ])
  
  #bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #colnames(bins)[1]      <- "bin"
  
  #junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(anchorj[!rownames(anchorj) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[colnames(junctions) %in% c("rn")]  <- c("J3")
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", "junction.log.mean",
                                         "junction.logFC", "junction.LR", "junction.pvalue", "junction.fdr", "J1.pvalue", "J2.pvalue",
                                         "junction.nonuniformity", "junction.dPIR", "junction.annotated",
                                         colnames(junctionbased_junctions)[grep("counts", colnames(junctionbased_junctions))])
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  anchorc <- data.table(anchorc(jdu), keep.rownames = T)  #Tries to find cluster locus
  junctions  <- strsplit2(anchorc$rn, "[.]")
  seqnames   <- junctions[, 1]
  start      <- as.numeric(junctions[, 2]) + 1
  end        <- as.numeric(junctions[, 3]) - 1
 
  grclusters <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  seqnames   <- strsplit2(genesDE(bdu)$gene_coordinates, "[:]")[, 1]
  start      <- genesDE(bdu)$start
  end        <- genesDE(bdu)$end
  
  grgenes    <- GRanges(seqnames, IRanges(start, end), strand="*")
  
  overlap    <- data.frame(findOverlaps(grclusters, grgenes))
  anchorc$locus <- NA
  anchorc$locus[overlap$queryHits] <- as.character(genesDE(bdu)$symbol[overlap$subjectHits])
  
  colnames(anchorc) <- c("rn", "cluster.LR", "cluster.pvalue", "cluster.fdr", "cluster.locus")
  fulldt <- data.frame(merge(junctionbased_junctions, anchorc, by.x = "junction", by.y = "rn", all = T))
  
  fulldt <- fulldt[, 
                   c("junction", "junction.annotated",
                     "junction.log.mean",
                     "junction.logFC", "junction.LR", "junction.pvalue", "junction.fdr", "J1.pvalue", "J2.pvalue",
                     "junction.nonuniformity", "junction.dPIR",
                     colnames(fulldt)[grep("counts", colnames(fulldt))],
                     "cluster.locus", "cluster.pvalue", "cluster.fdr",
                     "bin", "bin.pvalue", "bin.fdr")
                  ]
  rownames(fulldt) <- NULL
  #fulldt          <- fulldt[union(which(fulldt$bin.fdr < maxBinFDR), which(fulldt$cluster.fdr < maxJunctionFDR)), ]
  index.orden     <- strsplit2(fulldt$junction, "[.]")
  index.orden     <- order(GRanges(index.orden[, 1], IRanges(as.numeric(index.orden[, 2]), as.numeric(index.orden[, 3]))))
  #countsJ1 <- fulldt[, colnames(fulldt)[grep("countsJ1", colnames(fulldt))]]
  #countsJ2 <- fulldt[, colnames(fulldt)[grep("countsJ2", colnames(fulldt))]]
  #countsJ3 <- fulldt[, colnames(fulldt)[grep("countsJ3", colnames(fulldt))]]
  #pir             <- (countsJ1 + countsJ2)/(countsJ1 + countsJ2 + countsJ3)
  #colnames(pir)   <- strsplit2(colnames(pir), "[.]")[, 2]
  #dpir            <- rowSums(t(t(pir)*jdu@contrast[colnames(pir)]))
  #fulldt$dPIR     <- dpir
  fulldt          <- fulldt[index.orden, ] 
  rownames(fulldt)<- NULL
  anchorbased(mr) <- fulldt
  return(mr)
}

#Genera un data.table de regiones genicas con seniales diferenciales de diferente naturaleza
#                     region b bjs ja jl
# 1:       Chr1:45560-45645 1   1  0  0
# 2:     Chr1:387673-388267 1   1  1  1
# 3:     Chr1:406793-406902 1   1  0  0
# ...
#
# b: bin coverage signal / bjs: bin junction support signal / ja: junction anchor / jl: junction locale
#
#
# para la senial de junturas soporte de bin coverage se usan criterios basados en dPIN/dPI y 
# adicionalmente es posible usar, o no, el valor de fdr asociado a la juntura J3
# Aca lo que hacemos es ver el overlap entre regiones que tienen seniales diferentes. Para el caso del overlap entre 
# b o bjs y cualquier otro, usamos cualquier overlap, mientras que para el overlap entre ja y jl tienen que ser las mismas regiones
#
# Ademas, filtramos las ja y jl por junction.fdr unicamente. De esa forma, una region que tiene una juntura con uso diferencial, 
# por ejemplo, jl y que se solapa con un bin en al menos 3 pares de bases, aparece con soporte en bjs y en jl, mientras que si se solapa
# con b, solamente de coverage, aparece con b y jl. Las junturas que aparecen reportadas al final son las que aparecen en el bin en caso
# de tratarse de una region meramente "binica" o son la region en caso de venir de ja o jl.
#bin.fdr=0.1;nonunif=0.1;usenonunif=FALSE;dPSI=0.1;dPIR=0.1;j.fdr=0.1;j.particip=0.1;usepvalBJS=FALSE;bjs.fdr=0.1; otherSources = NULL; overlapType = "any"
.integrateSignals<-function(sr = NULL, asd = NULL, bin.fdr=0.1,nonunif=0.1,usenonunif=FALSE,dPSI=0.1,dPIR=0.1,j.fdr=0.1, j.particip=0.1,usepvalBJS=FALSE, bjs.fdr=0.1, otherSources = NULL, overlapType = "any"){
  
  if(class(sr) != "ASpliSplicingReport"){
    stop("sr must be an ASpliSplicingReport object")   
  }
  
  if(class(asd) != "ASpliAS"){
    stop("asd must be an ASpliAS object") 
  }
  
  #
  # Pongo todo lo que quiero comparar en GRanges
  #
  #bines significativos y uniformes 1.130149.130372 
  b  <- sr@binbased
  b  <- b[!is.na(b$start), ]
  if(usenonunif){
   b  <- b[ replace_na(b$bin.fdr < bin.fdr, FALSE) & (is.na(b$junction.dPIR) | replace_na(b$junction.nonuniformity < nonunif, FALSE)),]
  }else{
   b  <- b[ replace_na(b$bin.fdr < bin.fdr, FALSE),]
  }

  if(nrow(b) > 0){
    start   = as.numeric(b$start)#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(b$end)#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = strsplit2(b$gene_coordinates, ":")[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("binbased", times=length(seqnames))
    binbased <- GRanges(seqnames, IRanges(start, end), strand, tipo)
  }else{
    binbased <- GRanges()
  }
  
  #soporte de juntura para bines
  b  <- sr@binbased
  b  <- b[!is.na(b$start), ]
 
  if(!usenonunif){
    bunif <- TRUE
  }else{
    bunif <- replace_na(b$junction.nonuniformity < nonunif, FALSE)
  }
  if(usepvalBJS){
    b  <- b[ b$junction.fdr < bjs.fdr &
               (replace_na(abs(b$junction.dPSI) > dPSI, FALSE) | 
                 (replace_na(abs(b$junction.dPIR) > dPIR, FALSE) & bunif)), ]
    
  }else{  
    b  <- b[replace_na(abs(b$junction.dPSI) > dPSI, FALSE) | 
              (replace_na(abs(b$junction.dPIR) > dPIR, FALSE) & bunif), ]
  }
  if(nrow(b) > 0){
    start   = as.numeric(b$start)#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(b$end)#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = strsplit2(b$gene_coordinates, ":")[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("binbased", times=length(seqnames))
    binsupport <- GRanges(seqnames, IRanges(start, end), strand, tipo)
  }else{
    binsupport <- GRanges()
  }
  
  #Locales
  b <- sr@localebased
  b <- b[replace_na(b$cluster.fdr < j.fdr, FALSE) & 
           replace_na(b$cluster.participation > j.particip, FALSE), ] #Reemplazamos la participacion de la juntura por la participacion del cluster al que pertenece que es el dato importante
  if(nrow(b) > 0){
    aux <- strsplit2(b$cluster.range, "[.]") #Reemplazamos la juntura por el rango del cluster al que pertenece que es el dato importante
    start   = as.numeric(aux[, 2])#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(aux[, 3])#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = aux[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("localebased", times=length(aux[, 1]))
    localebased <- unique(GRanges(seqnames, IRanges(start, end), strand, tipo)) #Puede haber duplicados porque machean con varios bines
  }else{
    localebased <- GRanges()
  }
  
  b <- sr@anchorbased
  if(!usenonunif | all(is.na(b$junction.nonuniformity))){
    bunif <- TRUE
  }else{
    bunif <- replace_na(b$junction.nonuniformity < nonunif, FALSE)
  }
  b <- b[replace_na(b$junction.fdr < j.fdr, FALSE) & 
           replace_na(b$junction.dPIR > dPIR, FALSE) & 
           bunif, ]
  if(nrow(b) > 0){
    aux <- strsplit2(b$junction, "[.]")
    start   = as.numeric(aux[, 2]) + 1#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(aux[, 3]) - 1#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = aux[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("localebased", times=length(aux[, 1]))
    anchorbased <- unique(GRanges(seqnames, IRanges(start, end), strand, tipo)) #Puede haber duplicados porque machean con varios bines
  }else{
    anchorbased <- GRanges()
  }
  
  
  #
  # Overlap estimation
  #
  laux  <- list(b=binbased,bjs=binsupport,jl=localebased,ja=anchorbased)
  if(class(otherSources) == "GRanges") laux$otherSources = otherSources
  lover <- list()
  for(i in seq_along(laux)){
    for(j in (i+1):length(laux)){
      if(j>length(laux)) break
      saux <- paste0("overlaps_",paste0(names(laux)[c(i,j)],collapse="_"))

      ttype<-overlapType
      # if(names(laux)[i]%in%c("b","bjs","otherSources") | names(laux)[j]%in%c("b","bjs","otherSources")){
      #   ttype<-overlapType
      # }else{
      #   ttype<-"equal"
      # }

      suppressWarnings(taux <- as.data.table(findOverlaps(laux[[i]],laux[[j]],type=ttype, minoverlap = 3)))
      if(nrow(taux) > 0){
        taux[,queryHits:=paste(names(laux)[i],queryHits,sep=".")]
        taux[,subjectHits:=paste(names(laux)[j],subjectHits,sep=".")]
        colnames(taux) <- names(laux)[c(i,j)]
        lover[[saux]]           <-  taux
      }
    }
  }  
  
  #Genero overlaps_aux: con region y 0/1 si hay evidencia
  overlaps     <- rbindlist(lover, use.names = F)
  overlaps_aux <- data.table(region=character(), b=numeric(), bjs=numeric(), ja=numeric(), jl=numeric())
  if(class(otherSources) == "GRanges") overlaps_aux$otherSources = numeric()
  for(i in unique(c(overlaps$b, overlaps$bjs))){
    j <- c(i, overlaps$bjs[overlaps$b == i], overlaps$b[overlaps$bjs == i])
    d <- data.table(region = i,
                    b=as.numeric(length(grep("b.", j,fixed=TRUE)) > 0),
                    bjs=as.numeric(length(grep("bjs.", j,fixed=TRUE)) > 0),
                    ja=as.numeric(length(grep("ja.", j,fixed=TRUE)) > 0),
                    jl=as.numeric(length(grep("jl.", j,fixed=TRUE)) > 0))
    if(class(otherSources) == "GRanges") d$otherSources = as.numeric(length(grep("otherSources.", j,fixed=TRUE)) > 0)
    overlaps_aux <- rbindlist(list(overlaps_aux, d))
  }
  
  #junction anchorBased que no tienen overlap con otra clase de eventos
  if(length(anchorbased) > 0){
    nooverlap <- setdiff(paste0("ja.", (1:length(anchorbased))), overlaps_aux$region[grep("ja.", overlaps_aux$region,fixed=TRUE)])
    if(length(nooverlap) > 0){
      d <- data.table(region = nooverlap,
                      b=0,bjs=0,ja=1,jl=0)
      if(class(otherSources) == "GRanges") d$otherSources = 0    
      overlaps_aux <- rbindlist(list(overlaps_aux, d))
    }  
  }
  
  #junction localeBased que no tienen overlap con otra clase de eventos
  if(length(localebased) > 0){
    nooverlap <- setdiff(paste0("jl.", (1:length(localebased))), overlaps_aux$region[grep("jl.", overlaps_aux$region,fixed=TRUE)])
    if(length(nooverlap) > 0){
      d <- data.table(region = nooverlap,
                      b=0,bjs=0,ja=0,jl=1)
      if(class(otherSources) == "GRanges") d$otherSources = 0
      overlaps_aux <- rbindlist(list(overlaps_aux, d))
    }
  }
  
  #bin based que no tienen overlap con otra clase de eventos
  if(length(binbased) > 0){
    nooverlap <- setdiff(paste0("b.", (1:length(binbased))), overlaps_aux$region[grep("b.", overlaps_aux$region,fixed=TRUE)])
    if(length(nooverlap) > 0){
      d <- data.table(region = nooverlap,
                      b=1,bjs=0,ja=0,jl=0)
      if(class(otherSources) == "GRanges") d$otherSources = 0
      overlaps_aux <- rbindlist(list(overlaps_aux, d))
    }
  }
  
  #junction bin support que no tienen overlap con otra clase de eventos
  if(length(binsupport) > 0){
    nooverlap <- setdiff(paste0("bjs.", (1:length(binsupport))), overlaps_aux$region[grep("bjs.", overlaps_aux$region,fixed=TRUE)])
    if(length(nooverlap) > 0){
      d <- data.table(region = nooverlap,
                      b=0,bjs=1,ja=0,jl=0)
      if(class(otherSources) == "GRanges") d$otherSources = 0
      overlaps_aux <- rbindlist(list(overlaps_aux, d))
    }
  }
  
  #Other sources que no tienen overlap con otra clase de eventos
  if(class(otherSources) == "GRanges"){
    if(length(otherSources) > 0){
      nooverlap <- setdiff(paste0("otherSources.", (1:length(otherSources))), overlaps_aux$region[grep("otherSources.", overlaps_aux$region,fixed=TRUE)])
      if(length(nooverlap) > 0){
        d <- data.table(region = nooverlap,
                        b=0,bjs=0,ja=0,jl=0, otherSources = 1)
        overlaps_aux <- rbindlist(list(overlaps_aux, d))
      }
    }
  }
  #Hasta aca arme la tabla de 1/0 evidenciaL overlaps_aux
  
  
  # i <- grep("ja.", overlaps_aux$region,fixed=TRUE)
  # if(length(i)){
  #   L <- as.data.frame(anchorbased[as.numeric(strsplit2(overlaps_aux$region[i], "ja.")[, 2])])
  #   overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
  # }
  
  # i <- grep("jl.", overlaps_aux$region,fixed=TRUE)
  # L <- as.data.frame(localebased[as.numeric(strsplit2(overlaps_aux$region[i], "jl.")[, 2])])
  # overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start+1, "-", L$end-1)
  
  #
  #Cambio id de region por Rango genomico
  #
  i <- grep("b.", overlaps_aux$region,fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(binbased[as.numeric(strsplit2(overlaps_aux$region[i], "b.")[, 2])])
    overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
  }
  
  i <- grep("bjs.", overlaps_aux$region, fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(binsupport[as.numeric(strsplit2(overlaps_aux$region[i], "bjs.")[, 2])])
    overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
  }
  
  if(class(otherSources) == "GRanges"){
    i <- grep("otherSources.", overlaps_aux$region,fixed=TRUE)
    if(length(i)){
      L <- as.data.frame(otherSources[as.numeric(strsplit2(overlaps_aux$region[i], "otherSources.")[, 2])])
      overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
    }
  }
  
  i <- grep("jl.", overlaps_aux$region,fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(localebased[as.numeric(strsplit2(overlaps_aux$region[i], "jl.")[, 2])])
    overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
  }  
  
  #corregi rango de junturas para buscar bounderies de intrones, pero se reporta
  #el rango de la juntura original  
  # se consigna con extremos corregidos para matchear borde de intron del anchor con bin-intronico con senial,
  # si es un rango nuevo, se consigna el rango de la juntura original
  i <- grep("ja.", overlaps_aux$region,fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(anchorbased[as.numeric(strsplit2(overlaps_aux$region[i], "ja.")[, 2])])
    aux_region <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
    j <- which(aux_region %in% setdiff(aux_region, overlaps_aux$region))
    overlaps_aux$region[i] <- aux_region
    overlaps_aux$region[i[j]] <- paste0("Chr", L$seqnames[j], ":", L$start[j] - 1, "-", L$end[j] + 1)
  }  
  
  #comienzo a mergear rangos  
  if(nrow(overlaps_aux) > 0){
    overlaps_aux <- unique(overlaps_aux)  
    overlaps_aux <- overlaps_aux[order(overlaps_aux$region), ]
      
    #matcheo rango con bin
    roi  <- overlaps_aux[,region]
    rroi <- unlist(lapply(strsplit(roi,":",fixed=TRUE),function(x){return(x[2])}))
    #binbased(sr)
    ii<-match(rroi,paste0(binbased(sr)$start,"-",binbased(sr)$end))
    # table(is.na(ii))
    #table(binbased(sr)[ii,2:3])
    
    aa <- binbased(sr)[ii,c(1:3,7:8,13)]
    aa <- cbind(aa[,-c(4,5)],binreg=paste(aa$start,aa$end,sep="-"))
    
    aa  <-cbind(overlaps_aux,aa)
    roi <- gsub("[Chr]", "", roi)
    roi <- gsub("[:]", ".", roi)
    roi <- gsub("[-]", ".", roi)
    aa$J3 <- as.character(aa$J3)
    aa$J3[is.na(aa$J3)] <- roi[is.na(aa$J3)] #Hacemos esto para machear con los anchor. Recordar restar uno al start y sumar uno al end para machear correctamente
    
    #Se hizo esto para recalcular el locus de los locales que ahora vienen por cluster
    
    if(FALSE){
      if(class(asd) == "ASpliAS"){
        if(nrow(asd@junctionsPJU) > 0){
          aa$locus <- strsplit2(aa$bin, ":")[, 1]
          regiones <- gsub("[Chr]", "", aa$region)
          regiones <- gsub("[:]", ".", regiones)
          regiones <- gsub("[-]", ".", regiones)
          regiones <- sapply(regiones, function(r){
            i <- rownames(asd@junctionsPJU) == r
            if(sum(i) > 0){
              return(as.character(asd@junctionsPJU$symbol[i]))
            }else{
              return("") 
            }
          })
          aa$locus[regiones != ""] <- regiones[regiones != ""]
        }
      }
    }else{
      if(nrow(asd@junctionsPJU) > 0){
        regiones <- gsub("[Chr]", "", aa$region)
        regiones <- gsub("[:]", ".", regiones)
        regiones <- gsub("[-]", ".", regiones)
        regiones <- strsplit2(regiones, "[.]")
        regiones <- GRanges(regiones[, 1], ranges = IRanges(start = as.numeric(regiones[, 2]), end = as.numeric(regiones[, 3])))
        genes    <- gsub("[:]", ".", asd@junctionsPJU$gene_coordinates)
        genes    <- gsub("[-]", ".", genes)
        genes    <- strsplit2(genes, "[.]")
        rownames(genes) <- asd@junctionsPJU$symbol
        genes    <- unique(genes)
        grgenes  <- GRanges(genes[, 1], ranges = IRanges(start = as.numeric(genes[, 2]), end = as.numeric(genes[, 3])))
        overlap  <- data.frame(findOverlaps(regiones, grgenes, type="any"))
        aa$locus[overlap$queryHits] <- rownames(genes)[overlap$subjectHits]
      }else{
        aa$locus <- NA
      }
    }

    aa <- aa[, c(1, 11, 8, 2:7, 9:10)]
    
    aa$locus_overlap <- "-"
    locus_overlap <- binbased(sr)$locus_overlap[binbased(sr)$locus %in% aa$locus[!is.na(aa$locus)]]
    names(locus_overlap) <- binbased(sr)$locus[binbased(sr)$locus %in% aa$locus[!is.na(aa$locus)]]
    aa$locus_overlap[!is.na(aa$locus)] <- locus_overlap[aa$locus[!is.na(aa$locus)]]
    aa$locus_overlap[is.na(aa$locus_overlap)] <- "-"
       
    aa$b.fdr          <- NA
    aa$b.logfc        <- NA

    aa$bjs.lr            <- NA
    aa$bjs.fdr           <- NA
    aa$bjs.logfc         <- NA
    aa$bjs.nonuniformity <- NA
    aa$bjs.inclussion    <- NA
    
    aa$a.lr            <- NA
    aa$a.fdr           <- NA
    aa$a.logfc         <- NA
    aa$a.nonuniformity <- NA
    aa$a.dpir          <- NA

    aa$l.lr            <- NA
    aa$l.fdr           <- NA
    #aa$l.logfc         <- NA
    aa$l.participation <- NA
    aa$l.dparticipation <- NA

    #Generamos el rango para los bines para clasificar los locales
    bines       <- data.frame(bin = sr@binbased$bin, feature = sr@binbased$feature, gene_coordinates = sr@binbased$gene_coordinates, start = sr@binbased$start, end = sr@binbased$end, stringsAsFactors = F)
    bines       <- bines[complete.cases(bines), ]
    chromosomes <- strsplit2(bines$gene_coordinates, "[:]")[, 1]
    rango_bines <- GRanges(seqnames = chromosomes, ranges = IRanges(start = bines$start, end = bines$end))
    
    elementos_a_eliminar <- c() #Puede que podamos colapsar algunos elementos, por lo que los duplicados tienen que eliminarse
    
    for(b in 1:nrow(aa)){
      if(aa$b[b] != 0){
        i <- which(binbased(sr)$bin == aa$bin[b])
        if(length(i) > 0){
          aa$b.fdr[b] <- binbased(sr)$bin.fdr[i[1]]
          aa$b.logfc[b] <- binbased(sr)$bin.logFC[i[1]]
        }else{
          aa$b[b] <- "*"
        }
      }
      if(aa$bjs[b] != 0){
        i <- which(binbased(sr)$J3 == aa$J3[b])
        if(length(i) > 0){
          aa$bjs.lr[b] <- binbased(sr)$junction.LR[i[1]]
          aa$bjs.fdr[b] <- binbased(sr)$junction.fdr[i[1]]
          aa$bjs.logfc[b] <- binbased(sr)$junction.logFC[i[1]]
          aa$bjs.nonuniformity[b] <- binbased(sr)$junction.nonuniformity[i[1]]
          aa$bjs.inclussion[b] <- replace_na(binbased(sr)$junction.dPSI[i[1]], binbased(sr)$junction.dPIR[i[1]])  
        }else{
          aa$bjs[b] <- "*"
        }
      }
      if(aa$ja[b] != 0){
        i <- which(anchorbased(sr)$junction == aa$J3[b])
        if(length(i) > 0){
          aa$a.lr[b] <- anchorbased(sr)$junction.LR[i[1]]
          aa$a.fdr[b] <- anchorbased(sr)$junction.fdr[i[1]]
          aa$a.logfc[b] <- anchorbased(sr)$junction.logFC[i[1]]
          aa$a.nonuniformity[b] <- anchorbased(sr)$junction.nonuniformity[i[1]]
          aa$a.dpir[b] <- anchorbased(sr)$junction.dPIR[i[1]]
          #Si no tiene feature se lo tratamos de asignar
          if(is.na(aa$feature[b])){
            i <- which(binbased(sr)$J3 == aa$J3[b])
            if(length(i) > 0){
              if(length(grep("Io", binbased(sr)$bin)) > 0){
                aa$bin.event[b] <- "IoR" #Es un evento conocido, puede venir de un Io  
                aa$bin[b] <- binbased(sr)$bin[i[[1]]]
              }
            }else{            
              if(aa$b[b] == "*"){ #Es un evento complejo, csp
                aa$bin.event[b] <- "CSP" 
              }else if(aa$b[b] == "0"){ #Es un evento nuevo que no se asocia a ningun bin, new alternative splicing pattern
                aa$bin.event[b] <- "Novel ASP" 
              }
            }
          }
        }else{
          aa$ja[b] <- "*"
        }
      }
      if(aa$jl[b] != 0){
        if(FALSE){
          if(!is.na(localebased(sr)$J3[b])){
            i <- which(localebased(sr)$junction == aa$J3[b])
            if(length(i) > 0){
              aa$l.lr[b] <- localebased(sr)$cluster.LR[i[1]]
              aa$l.fdr[b] <- localebased(sr)$junction.fdr[i[1]]
              aa$l.logfc[b] <- localebased(sr)$junction.logFC[i[1]]
              aa$l.participation[b] <- localebased(sr)$junction.participation[i[1]]
              aa$l.dparticipation[b] <- localebased(sr)$junction.dparticipation[i[1]]
            }
          } 
        }else{
          #Aca busco el cluster por region en lugar de buscar por juntura especifica
          roi <- gsub("[Chr]", "", aa$region[b])
          roi <- gsub("[:]", ".", roi)
          roi <- gsub("[-]", ".", roi)
          i <- which(localebased(sr)$cluster.range == roi)
          if(length(i) > 0){
            aa$l.lr[b] <- localebased(sr)$cluster.LR[i[1]]
            aa$l.fdr[b] <- localebased(sr)$cluster.fdr[i[1]]
            aa$l.participation[b] <- localebased(sr)$cluster.participation[i[1]]
            aa$l.dparticipation[b] <- localebased(sr)$cluster.dparticipation[i[1]]
          }else{ #Vemos si machea con la juntura quizas
            i <- which(localebased(sr)$cluster.range == aa$J3[b])
            if(length(i) > 0){
              aa$l.lr[b] <- localebased(sr)$cluster.LR[i[1]]
              aa$l.fdr[b] <- localebased(sr)$cluster.fdr[i[1]]
              aa$l.participation[b] <- localebased(sr)$cluster.participation[i[1]]
              aa$l.dparticipation[b] <- localebased(sr)$cluster.dparticipation[i[1]]
            }else{
              aa$jl[b] <- "*"
            }
          }
          
          if(aa$jl[b] == 1){
            if(is.na(aa$feature[b])){ #Tratamos de darle sentido a la juntura
              junturas <- strsplit2(unique(localebased(sr)$junction[i]), "[.]")
              #Vemos si hay un unico pivot.
              punta_1 <- length(table(junturas[, 2]))
              punta_2 <- length(table(junturas[, 3]))
              
              #Vemos si es alt o complejo              
              if(punta_1 == 1 | punta_2 == 1){
                
                if(punta_1 == 1){
                  #Buscamos el rango de las junturas
                  rango_min    <- min(as.numeric(junturas[, 3]))
                  rango_max    <- max(as.numeric(junturas[, 3]))
                  rango        <- GRanges(junturas[1, 1], IRanges(rango_min, rango_max))
                  rango_intron <- GRanges(junturas[1, 1], IRanges(rango_min, rango_max - 1))
                  
                }else if (punta_2 == 1){ 
                  
                  #Buscamos el rango de las junturas
                  rango_min    <- min(as.numeric(junturas[, 2]))
                  rango_max    <- max(as.numeric(junturas[, 2]))
                  rango        <- GRanges(junturas[1, 1], IRanges(rango_min, rango_max))
                  rango_intron <- GRanges(junturas[1, 1], IRanges(rango_min + 1, rango_max))
 
                }
                
                #Vemos que tipo de bines atraviezan las junturas. Si son todos del mismo tipo, es alt
                hits        <- data.frame(findOverlaps(rango, rango_bines))
                if(length(table(bines$feature[hits$subjectHits])) == 1){
                  aa$bin.event[b] <- "Alt 5'/3'"
                }else{
                  #Vemos si se trata de un intron alt
                  hits        <- data.frame(findOverlaps(rango_intron, rango_bines))
                  if(length(table(bines$feature[hits$subjectHits])) == 1){
                    aa$bin.event[b] <- "Alt 5'/3'"
                  }else{ #No parece ser alt, es complejo entonces
                    aa$bin.event[b] <- "CSP"
                  }
                }
              }else{
                #Vemos si es ES
                if(nrow(junturas) == 3){
                  i <- which(localebased(sr)$cluster.range == aa$J3[b])
                  if(length(i) > 0){ 
                    aa$bin.event[b] <- "ES"
                    #Vemos si podemos encontrar el bin asociado a este ES
                    j <- which(aa$region == aa$J3[b])
                    if(length(j) > 0){ #Lo encontramos, sacamos directamente esta fila porque coincide con la que ya existe y tiene mas datos.
                      elementos_a_eliminar <- c(elementos_a_eliminar, b)
                      #aa$b.fdr[b]   <- binbased(sr)$bin.fdr[j[1]]
                      #aa$b.logfc[b] <- binbased(sr)$bin.logFC[j[1]]
		                  #aa$b[b]       <- 1
                    }              
                  }else{
                    aa$bin.event[b] <- "CSP"
                  }
                }else{ #Si no es nada de lo anterior entonces es complejo
                  aa$bin.event[b] <- "CSP"
                }
              }
              #Vemos si es un novel event
              if(any(asd@junctionsPJU$junction[rownames(asd@junctionsPJU) %in% localebased(sr)$junction[i]] == "noHit")){
                aa$bin.event[b] <- paste("Novel", aa$bin.event[b])
              }
            }
          }
        }
      }  
      
      #Tratamos de darle sentido a otros eventos - de bines
      if((aa$b[b] == 1 | aa$bjs[b] == 1 | aa$jl[b] == 1) & aa$bin.event[b] == "-"){
        if(aa$feature[b] == "E"){
          aa$bin.event[b] <- "ASCE"
        }else if(aa$feature[b] == "I"){ 
          if(aa$b[b] != "*" & aa$bjs[b] != "*" & aa$ja[b] != "*" & aa$jl[b] != "*"){
            aa$bin.event[b] <- "IR"
          }else{
            aa$bin.event[b] <- "CSP"
          }
        }
      }      
    }  
    
    #Descartamos elementos duplicados por ES
    if(!is.null(elementos_a_eliminar)){
      aa <- aa[-elementos_a_eliminar, ]
    }
    
    aa <- aa[order(aa$b.fdr), ]
    r <- strsplit2(unique(aa$region), "Chr")[, 2]
    chr <- strsplit2(r, ":")[, 1]
    chr <- sapply(chr, function(s){return(is.na(suppressWarnings(as.numeric(s))))})
    aa$region[chr] <- r[chr] 
    
    #Finds events with bjs and locale and tries to merge them
    bjs_jl <- aa$b == 0 & aa$bjs == 1 & aa$ja == 0 & aa$jl == 1
    if(any(bjs_jl)){
      regions <- strsplit2(strsplit2(aa$region[bjs_jl], ":")[, 2], "-")
      regions <- cbind(regions, locus=aa$locus[bjs_jl])
      regions <- cbind(regions, bin=aa$bin[bjs_jl])
      regions <- cbind(regions, region_original=aa$region[bjs_jl])
      regions <- as.data.table(regions)
      regions$V1 <- as.numeric(regions$V1)
      regions$V2 <- as.numeric(regions$V2)
      a_descartar <- c()
      for(i in unique(regions$locus)){
        a_ver <- regions[regions$locus == i, ]
        posible_a_descartar <- which(is.na(a_ver$bin))
        if(length(posible_a_descartar) > 0){
          for(j in posible_a_descartar){
            if(any(a_ver$V1 == a_ver$V1[j] + 1 | a_ver$V1 == a_ver$V1[j] - 1 | a_ver$V2 == a_ver$V2[j] - 1 | a_ver$V2 == a_ver$V2[j] - 1)){
              a_descartar <- c(a_descartar, which(aa$region == a_ver$region_original[j] & is.na(aa$bin) & aa$b == 0 & aa$bjs == 1 & aa$ja == 0 & aa$jl == 1))
            }
          }
        }
      }
      if(length(a_descartar)>0){
        aa <- aa[-a_descartar, ]
      }
    }
    
  }else{
    aa <- data.table(region = character(), locus = character(), b = numeric(), bjs = numeric(), ja = numeric(),
                     jl = numeric(), bin = character(), feature = character(), bin.event = character(),
                     J3 = character(), binreg = character(), locus_overlap = numeric(), b.fdr = numeric(),
                     b.logfc = numeric(), bjs.lr = numeric(), bjs.fdr = numeric(), bjs.logfc = numeric(), bjs.nonuniformity = numeric(),
                     bjs.inclussion = numeric(), a.lr = numeric(), a.fdr = numeric(), a.logfc = numeric(), a.nonuniformity = numeric(),
                     a.participation = numeric(),
                     a.dparticipation = numeric(), l.lr = numeric(), l.fdr = numeric(), l.logfc = numeric(), l.participation = numeric(),l.dparticipation = numeric())
  }
  

  return(aa)
}
