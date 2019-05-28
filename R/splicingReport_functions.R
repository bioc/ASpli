
.splicingReport <- function(bdu, jdu, maxBinFDR, maxJunctionFDR){
  
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
  jir$dPIN               <- NA
  jir                    <- jir[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "NonUniformity", colnames(jir)[grep("counts", colnames(jir))])]
  colnames(jir)[c(7, 10)] <- c("junction.fdr", "nonuniformity")
  #jir$bin                <- rownames(jir)
  
  jes                    <- jes(jdu)
  jes$nonuniformity      <- NA
  jes$dPIR               <- NA
  jes                    <- jes[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "nonuniformity", colnames(jes)[grep("counts", colnames(jes))])]
  colnames(jes)[7]       <- "junction.fdr"
  #jes$bin                <- rownames(jes)
  
  jalt                   <- jalt(jdu)
  jalt$nonuniformity     <- NA
  jalt$dPIR              <- NA
  jalt                   <- jalt[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "nonuniformity", colnames(jalt)[grep("counts", colnames(jalt))])]
  colnames(jalt)[7]      <- "junction.fdr"
  #jalt$bin               <- rownames(jalt)
  
  j                      <- data.table(rbind(jir, jes, jalt), keep.rownames = T)
  bins                   <- data.table(binsDU(bdu), keep.rownames = T)
  aux                    <- data.frame(merge(bins, j, by="rn", all=T))
  colnames(aux)          <- c("bin", "feature", "bin.event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start",
                              "end", "length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction.event", "J3", "J3.multiplicity", "junction.logFC",
                              "junction.log.mean", "junction.pvalue", "junction.fdr", "junction.dPIR", "junction.dPIN", "junction.nonuniformity", colnames(aux)[grep("counts", colnames(aux))])
  
  aux                    <- aux[, !colnames(aux) %in% c("symbol", "junction.event")]
  
  aux                    <- aux[union(which(aux$bin.fdr < maxBinFDR), which(aux$junction.fdr < maxJunctionFDR)), ]
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
  colnames(bins)[10] <- "bin.logFC"
  
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(localej(jdu))[overlap$subjectHits], localej(jdu)[overlap$subjectHits, ])
  colnames(junctionbased_junctions)[18] <- "junction.pvalue"
  
  #bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #colnames(bins)[1]      <- "bin"
  
  #junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(localej(jdu)[!rownames(localej(jdu)) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[c(1, 5)]  <- c("J3", "junction.pvalue")
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr","junction", "junction.cluster", "junction.log.mean",
                                         "junction.logFC", "junction.pvalue", "junction.fdr", "junction.annotated", "junction.participation", colnames(junctionbased_junctions)[grep("counts", colnames(junctionbased_junctions))])
  
  junctionbased_junctions$junction.cluster <- as.character(junctionbased_junctions$junction.cluster)
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  localec <- data.table(localec(jdu), keep.rownames = T)
  colnames(localec) <- c("rn", "cluster.size", "cluster.LR", "cluster.pvalue", "cluster.fdr")
  fulldt <- data.frame(merge(junctionbased_junctions, localec, by.x = "junction.cluster", by.y = "rn", all = T))
  
  fulldt <- fulldt[, 
                   c("junction", "junction.annotated",
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", 
                     colnames(fulldt)[grep("counts", colnames(fulldt))],
                     "junction.cluster", "junction.participation",
                     "cluster.size", "cluster.pvalue", "cluster.fdr", 
                     "bin", "bin.pvalue", "bin.fdr")
                   ]
  rownames(fulldt) <- NULL
  fulldt          <- fulldt[union(which(fulldt$bin.fdr < maxBinFDR), which(fulldt$junction.fdr < maxJunctionFDR)), ]
  index.orden     <- strsplit2(fulldt$junction, "[.]")
  index.orden     <- order(GRanges(index.orden[, 1], IRanges(as.numeric(index.orden[, 2]), as.numeric(index.orden[, 3]))))
  rownames(fulldt)<- NULL
  localebased(mr) <- fulldt[index.orden, ]    
  
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
  colnames(bins)[10] <- "bin.logFC"
  
  junctionbased_junctions <- data.table(bin=rownames(bins)[overlap$queryHits], bins[overlap$queryHits, ], 
                                        J3 = rownames(anchorj(jdu))[overlap$subjectHits], anchorj(jdu)[overlap$subjectHits, ])
  colnames(junctionbased_junctions)[17] <- "junction.pvalue"
  #bins                   <- data.table(bins[!rownames(bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #colnames(bins)[1]      <- "bin"
  
  #junctionbased_junctions <- rbindlist(list(junctionbased_junctions, bins), use.names = T, fill = T)
  
  junctions               <- data.table(anchorj(jdu)[!rownames(anchorj(jdu)) %in% junctionbased_junctions$J3, ], keep.rownames = T)
  colnames(junctions)[c(1, 4)]  <- c("J3", "junction.pvalue")
  
  junctionbased_junctions <- rbindlist(list(junctionbased_junctions, junctions), use.names = T, fill = T)
  
  colnames(junctionbased_junctions) <- c("bin", "feature", "event", "locus", "locus_overlap", "symbol", "gene_coordinates", "bin.start", "bin.end",
                                         "bin.length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction", "junction.log.mean",
                                         "junction.logFC", "junction.pvalue", "junction.fdr", "junction.nonuniformity", "junction.participation", "junction.annotated",
                                         colnames(junctionbased_junctions)[grep("counts", colnames(junctionbased_junctions))])
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  anchorc <- data.table(anchorc(jdu), keep.rownames = T)
  colnames(anchorc) <- c("rn", "cluster.LR", "cluster.pvalue", "cluster.fdr")
  fulldt <- data.frame(merge(junctionbased_junctions, anchorc, by.x = "junction", by.y = "rn", all = T))
  
  fulldt <- fulldt[, 
                   c("junction", "junction.annotated",
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", "junction.nonuniformity", "junction.participation",
                     colnames(fulldt)[grep("counts", colnames(fulldt))],
                     "cluster.pvalue", "cluster.fdr",
                     "bin", "bin.pvalue", "bin.fdr")
                   ]
  rownames(fulldt) <- NULL
  fulldt          <- fulldt[union(which(fulldt$bin.fdr < maxBinFDR), which(fulldt$junction.fdr < maxJunctionFDR)), ]
  index.orden     <- strsplit2(fulldt$junction, "[.]")
  index.orden     <- order(GRanges(index.orden[, 1], IRanges(as.numeric(index.orden[, 2]), as.numeric(index.orden[, 3]))))
  countsJ1 <- fulldt[, colnames(fulldt)[grep("countsJ1", colnames(fulldt))]]
  countsJ2 <- fulldt[, colnames(fulldt)[grep("countsJ2", colnames(fulldt))]]
  countsJ3 <- fulldt[, colnames(fulldt)[grep("countsJ3", colnames(fulldt))]]
  pir             <- (countsJ1 + countsJ2)/(countsJ1 + countsJ2 + countsJ3)
  colnames(pir)   <- strsplit2(colnames(pir), "[.]")[, 2]
  dpir            <- rowSums(t(t(pir)*jdu@contrast[colnames(pir)]))
  fulldt$dPIR     <- dpir
  rownames(fulldt)<- NULL
  anchorbased(mr) <- fulldt[index.orden, ]  
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

.integrateSignals<-function(sr = NULL, asd = NULL, bin.fdr=0.05,unif=0.1,dPIN=0.05,dPIR=0.05,j.fdr=0.05,j.particip=0.1,usepvalBJS=FALSE,bjs.fdr=0.1, otherSources = NULL){
  
  if(class(sr) != "ASpliSplicingReport"){
    stop("sr must be an ASpliSplicingReport object") 
  }
  
  if(class(asd) != "ASpliAS"){
    stop("asd must be an ASpliAS object") 
  }
  
  #
  # Pongo todo lo que quiero compara en GRanges
  #
  #bines significativos y uniformes
  b  <- sr@binbased
  b  <- b[!is.na(b$start), ]
  b  <- b[ replace_na(b$bin.fdr < bin.fdr, FALSE) & (is.na(b$junction.dPIR) | replace_na(b$junction.nonuniformity < unif, FALSE)),]
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
  
  if(usepvalBJS){
    b  <- b[ b$junction.fdr < bjs.fdr &
               (replace_na(abs(b$junction.dPIN) > dPIN, FALSE) | 
                  (replace_na(abs(b$junction.dPIR) > dPIR, FALSE) & replace_na(b$junction.nonuniformity < unif, FALSE))), ]
    
  }else{  
    b  <- b[replace_na(abs(b$junction.dPIN) > dPIN, FALSE) | 
              (replace_na(abs(b$junction.dPIR) > dPIR, FALSE) & replace_na(b$junction.nonuniformity < unif, FALSE)), ]
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
  
  b <- sr@localebased
  b <- b[replace_na(b$junction.fdr < j.fdr, FALSE) & 
           replace_na(b$junction.participation > j.particip, FALSE), ]
  if(nrow(b) > 0){
    aux <- strsplit2(b$junction, "[.]")
    start   = as.numeric(aux[, 2])#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(aux[, 3])#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = aux[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("localebased", times=length(aux[, 1]))
    localebased <- GRanges(seqnames, IRanges(start, end), strand, tipo)
  }else{
    localebased <- GRanges()
  }
  
  b <- sr@anchorbased
  b <- b[replace_na(b$junction.fdr < j.fdr, FALSE) & 
           replace_na(b$junction.participation > j.particip, FALSE) & 
           replace_na(b$junction.nonuniformity < unif, FALSE), ]
  if(nrow(b) > 0){
    aux <- strsplit2(b$junction, "[.]")
    start   = as.numeric(aux[, 2])#aggregate(start ~ cluster, data=b, FUN=min)
    end     = as.numeric(aux[, 3])#aggregate(end ~ cluster, data=b, FUN=max)
    seqnames = aux[, 1]
    strand  = rep("*", times=length(seqnames))
    tipo    = rep("localebased", times=length(aux[, 1]))
    anchorbased <- GRanges(seqnames, IRanges(start, end), strand, tipo)
  }else{
    anchorbased <- GRanges()
  }
  
  
  #
  # Overlaps
  #
  laux  <- list(b=binbased,bjs=binsupport,jl=localebased,ja=anchorbased)
  if(class(otherSources) == "GRanges") laux$otherSources = otherSources
  lover <- list()
  for(i in seq_along(laux)){
    for(j in (i+1):length(laux)){
      if(j>length(laux)) break
      saux <- paste0("overlaps_",paste0(names(laux)[c(i,j)],collapse="_"))
      #if(names(laux)[i]%in%c("b","bjs") | names(laux)[j]%in%c("b","bjs")){
      if(names(laux)[i]%in%c("b","bjs","otherSources") | names(laux)[j]%in%c("b","bjs","otherSources")){
        ttype<-"any"
      }else{
        ttype<-"equal"
      }
      suppressWarnings(taux           <- as.data.table(findOverlaps(laux[[i]],laux[[j]],type=ttype, minoverlap = 3)))
      if(nrow(taux) > 0){
        taux[,queryHits:=paste(names(laux)[i],queryHits,sep=".")]
        taux[,subjectHits:=paste(names(laux)[j],subjectHits,sep=".")]
        colnames(taux) <- names(laux)[c(i,j)]
        lover[[saux]]           <-  taux
      }
    }
  }  
  
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
  
  # i <- grep("ja.", overlaps_aux$region,fixed=TRUE)
  # if(length(i)){
  #   L <- as.data.frame(anchorbased[as.numeric(strsplit2(overlaps_aux$region[i], "ja.")[, 2])])
  #   overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start, "-", L$end)
  # }
  
  # i <- grep("jl.", overlaps_aux$region,fixed=TRUE)
  # L <- as.data.frame(localebased[as.numeric(strsplit2(overlaps_aux$region[i], "jl.")[, 2])])
  # overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start+1, "-", L$end-1)
  
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
  
  
  #corregi rango de junturas para buscar bounderies de intrones, pero se reporta
  #el rango de la juntura original
  i <- grep("jl.", overlaps_aux$region,fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(localebased[as.numeric(strsplit2(overlaps_aux$region[i], "jl.")[, 2])])
    aux_region <- paste0("Chr", L$seqnames, ":", L$start + 1, "-", L$end - 1)
    j <- which(aux_region %in% setdiff(aux_region, overlaps_aux$region))
    overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start + 1, "-", L$end - 1)
    overlaps_aux$region[i[j]] <- paste0("Chr", L$seqnames[j], ":", L$start[j], "-", L$end[j])
  }  
  
  i <- grep("ja.", overlaps_aux$region,fixed=TRUE)
  if(length(i)){
    L <- as.data.frame(localebased[as.numeric(strsplit2(overlaps_aux$region[i], "ja.")[, 2])])
    aux_region <- paste0("Chr", L$seqnames, ":", L$start + 1, "-", L$end - 1)
    j <- which(aux_region %in% setdiff(aux_region, overlaps_aux$region))
    overlaps_aux$region[i] <- paste0("Chr", L$seqnames, ":", L$start + 1, "-", L$end - 1)
    overlaps_aux$region[i[j]] <- paste0("Chr", L$seqnames[j], ":", L$start[j], "-", L$end[j])
  }  
  
  
  
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
    
    aa <-cbind(overlaps_aux,aa)
    roi <- gsub("[Chr]", "", roi)
    roi <- gsub("[:]", ".", roi)
    roi <- gsub("[-]", ".", roi)
    aa$J3 <- as.character(aa$J3)
    aa$J3[is.na(aa$J3)] <- roi[is.na(aa$J3)]
    
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
        aa$locus[!is.na(regiones)] <- regiones[!is.na(regiones)]
        aa <- aa[, c(1, 11, 2:10)]
      }
    }
    aa$locus_overlap <- "-"
    locus_overlap <- binbased(sr)$locus_overlap[binbased(sr)$locus %in% aa$locus[aa$locus != ""]]
    names(locus_overlap) <- binbased(sr)$locus[binbased(sr)$locus %in% aa$locus[aa$locus != ""]]
    for(i in 1:length(locus_overlap)){
      j <- aa$locus == names(locus_overlap)[i]
      aa$locus_overlap[j] <- locus_overlap[i]
    }
        
    aa$bfdr          <- 1
    aa$blogfc        <- 0

    aa$bjsfdr        <- 1
    aa$bjslogfc      <- 0
    aa$bjsnonuniformity <- Inf
    aa$bjsinclussion <- 0
    
    aa$afdr           <- 1
    aa$alogfc         <- 0
    aa$anonuniformity  <- Inf
    aa$aparticipation <- 0
    
    aa$lfdr           <- 1
    aa$llogfc         <- 0
    aa$lparticipation <- 0
    
        
    for(b in 1:nrow(aa)){
      if(aa$b[b] != 0){
        i <- which(binbased(sr)$bin == aa$bin[b])
        if(length(i) > 0){
          aa$bfdr[b] <- binbased(sr)$bin.fdr[i[1]]
          aa$blogfc[b] <- binbased(sr)$bin.logFC[i[1]]
        }
      }
      if(aa$bjs[b] != 0){
        i <- which(binbased(sr)$J3 == aa$J3[b])
        if(length(i) > 0){
          aa$bjsfdr[b] <- binbased(sr)$junction.fdr[i[1]]
          aa$bjslogfc[b] <- binbased(sr)$junction.logFC[i[1]]
          aa$bjsnonuniformity[b] <- binbased(sr)$junction.nonuniformity[i[1]]
          aa$bjsinclussion[b] <- replace_na(binbased(sr)$junction.dPIN[i[1]], binbased(sr)$junction.dPIR[i[1]])  
        }
      }      
      if(aa$ja[b] != 0){
        i <- which(anchorbased(sr)$junction == aa$J3[b])
        if(length(i) > 0){
          aa$afdr[b] <- anchorbased(sr)$junction.fdr[i[1]]
          aa$alogfc[b] <- anchorbased(sr)$junction.logFC[i[1]]
          aa$anonuniformity[b] <- anchorbased(sr)$junction.nonuniformity[i[1]]
          aa$aparticipation[b] <- anchorbased(sr)$junction.participation[i[1]]
        }
      }
      if(aa$jl[b] != 0){
        i <- which(localebased(sr)$junction == aa$J3[b])
        if(length(i) > 0){
          aa$lfdr[b] <- localebased(sr)$junction.fdr[i[1]]
          aa$llogfc[b] <- localebased(sr)$junction.logFC[i[1]]
          aa$lparticipation[b] <- localebased(sr)$junction.participation[i[1]]
        }
      }          
    }  
    aa <- aa[order(aa$bfdr), ]
    r <- strsplit2(unique(aa$region), "Chr")[, 2]
    chr <- strsplit2(r, ":")[, 1]
    chr <- sapply(chr, function(s){return(is.na(suppressWarnings(as.numeric(s))))})
    aa$region[chr] <- r[chr] 
  }else{
    aa <- data.table(region = character(), locus = character(), b = numeric(), bjs = numeric(), ja = numeric(),
                     jl = numeric(), bin = character(), feature = character(), bin.event = character(),
                     J3 = character(), binreg = character())
  }
  

  return(aa)
}
