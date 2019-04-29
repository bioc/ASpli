
.splicingReport <- function(bdu, jdu, maxBinFDR, maxJunctionFDR){
  
  #Sanity check
  if(class(bdu) != "ASpliDU"){
    stop("bdu must be an ASpliDU object, try running gbDUreport first") 
  }
  if(class(jdu) != "ASpliJDU"){
    stop("jdu must be an ASpliJDU object, try running jDUreport first") 
  }
  
  mr <- new( Class="ASpliSplicingReport" )
  
  #Mergin bin based junctions
  jir                    <- jir(jdu)
  jir$event              <- NA
  jir$dPIN               <- NA
  jir                    <- jir[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "Uniformity", colnames(jir)[grep("counts", colnames(jir))])]
  colnames(jir)[c(7, 10)] <- c("junction.fdr", "uniformity")
  #jir$bin                <- rownames(jir)
  
  jes                    <- jes(jdu)
  jes$uniformity         <- NA
  jes$dPIR               <- NA
  jes                    <- jes[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "uniformity", colnames(jes)[grep("counts", colnames(jes))])]
  colnames(jes)[7]       <- "junction.fdr"
  #jes$bin                <- rownames(jes)
  
  jalt                   <- jalt(jdu)
  jalt$uniformity        <- NA
  jalt$dPIR              <- NA
  jalt                   <- jalt[, c("event", "J3", "multiplicity", "logFC", "log.mean", "pvalue", "FDR", "dPIR", "dPIN", "uniformity", colnames(jalt)[grep("counts", colnames(jalt))])]
  colnames(jalt)[7]      <- "junction.fdr"
  #jalt$bin               <- rownames(jalt)
  
  j                      <- data.table(rbind(jir, jes, jalt), keep.rownames = T)
  bins                   <- data.table(binsDU(bdu), keep.rownames = T)
  aux                    <- data.frame(merge(bins, j, by="rn", all=T))
  colnames(aux)          <- c("bin", "feature", "bin.event", "locus", "locus_overlap", "symbol", "gene_coordinates", "start",
                              "end", "length", "bin.logFC", "bin.pvalue", "bin.fdr", "junction.event", "J3", "J3.multiplicity", "junction.logFC",
                              "junction.log.mean", "junction.pvalue", "junction.fdr", "junction.dPIR", "junction.dPIN", "junction.uniformity", colnames(aux)[grep("counts", colnames(aux))])
  
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
                                         "junction.logFC", "junction.pvalue", "junction.fdr", "junction.uniformity", "junction.participation", "junction.annotated",
                                         colnames(junctionbased_junctions)[grep("counts", colnames(junctionbased_junctions))])
  
  #bins <- data.table(bdu@bins[rownames(bdu@bins) %in% junctionbased_junctions$bin, ], keep.rownames = T)
  #completo <- rbindlist(junctionbased_junctions, bins, use.names = T, fill = T)
  
  anchorc <- data.table(anchorc(jdu), keep.rownames = T)
  colnames(anchorc) <- c("rn", "cluster.LR", "cluster.pvalue", "cluster.fdr")
  fulldt <- data.frame(merge(junctionbased_junctions, anchorc, by.x = "junction", by.y = "rn", all = T))
  
  fulldt <- fulldt[, 
                   c("junction", "junction.annotated",
                     "junction.log.mean",
                     "junction.logFC", "junction.pvalue", "junction.fdr", "junction.uniformity", "junction.participation",
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