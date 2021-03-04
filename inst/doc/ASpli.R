### R code from vignette source 'ASpli.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()
options(continue=" ")


###################################################
### code chunk number 2: quickStart1
###################################################
library(ASpli)
library(GenomicFeatures)

# gtf preprocessing ----
gtfFileName <- aspliExampleGTF()
genomeTxDb  <- makeTxDbFromGFF( gtfFileName )

# feature extraction ----
features    <- binGenome( genomeTxDb )

#bams and target file ----
BAMFiles <- aspliExampleBamList()
targets  <- data.frame(row.names = paste0('Sample',c(1:6)),
                       bam = BAMFiles[1:6],
                       f1  = c( 'control','control','control','treatment','treatment','treatment'),
                       stringsAsFactors = FALSE)
mBAMs <- data.frame( bam = sub("_[012]","",targets$bam[c(1,4)]),
                     condition = c("control","treatment"))



###################################################
### code chunk number 3: quickStart2
###################################################
gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 100, maxISize = 50000)
gbcounts


###################################################
### code chunk number 4: quickStart2b
###################################################
asd <- jCounts(counts=gbcounts, features=features, minReadLength=100)
asd


###################################################
### code chunk number 5: quickStart3
###################################################
gb  <- gbDUreport(gbcounts, contrast = c(-1,1))
gb


###################################################
### code chunk number 6: quickStart3b
###################################################
jdur <- jDUreport(asd, contrast=c(-1,1))
jdur


###################################################
### code chunk number 7: quickStart4
###################################################
sr <- splicingReport(gb, jdur, counts=gbcounts)


###################################################
### code chunk number 8: quickStart4b
###################################################
is <- integrateSignals(sr,asd)


###################################################
### code chunk number 9: quickStart5
###################################################
exportIntegratedSignals(is,sr=sr,
                          output.dir = "aspliExample",
                          counts=gbcounts,features=features,asd=asd,
                          mergedBams = mBAMs)


###################################################
### code chunk number 10: installation (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## 
## # The following line initializes usage of Bioc devel branch and should not be necessary after
## # the next official release scheduled for October 2020
## if(Sys.Date()<"2020-11-01") BiocManager::install(version='devel')
## BiocManager::install("ASpli")
## library(ASpli)


###################################################
### code chunk number 11: makeTx (eval = FALSE)
###################################################
## library(GenomicFeatures)
## TxDb <- makeTxDbFromGFF(
##   file="genes.gtf",
##   format="gtf")


###################################################
### code chunk number 12: saveTx (eval = FALSE)
###################################################
##  saveDb(TxDB,file="gene.sqlite")


###################################################
### code chunk number 13: binGenome (eval = FALSE)
###################################################
## annFile       <- aspliExampleGTF()
## aTxDb         <- makeTxDbFromGFF(annFile)
## features      <- binGenome( aTxDb ) 
## geneCoord     <- featuresg( features )
## binCoord      <- featuresb( features )
## junctionCoord <- featuresj( features )


###################################################
### code chunk number 14: binGenome2 (eval = FALSE)
###################################################
## symbols       <- data.frame( row.names = genes( aTxDb ), 
##                              symbol = paste( 'This is symbol of gene:',
##                                              genes( aTxDb ) ) )
## features      <- binGenome( aTxDb, geneSymbols = symbols ) 


###################################################
### code chunk number 15: targetsDF2 (eval = FALSE)
###################################################
## BAMFiles <- c("path_to_bams/CT_time1_rep1.BAM", "path_to_bams/CT_time1_rep2.BAM",
##               "path_to_bams/CT_time2_rep1.BAM", "path_to_bams/CT_time2_rep2.BAM",
##               "path_to_bams/TR_time1_rep1.BAM", "path_to_bams/TR_time1_rep2.BAM",
##               "path_to_bams/TR_time2_rep1.BAM", "path_to_bams/TR_time2_rep2.BAM")
## (targets <- data.frame( bam = BAMFiles,
##                         genotype = c( 'CT', 'CT', 'CT',  'CT', 
##                                       'TR', 'TR', 'TR', 'TR' ),
##                         time     = c( 't1', 't1', 't2', 't2', 
##                                       't1', 't1', 't2', 't2' ),
##                          stringsAsFactors = FALSE ))


###################################################
### code chunk number 16: targetsDF2b (eval = FALSE)
###################################################
## getConditions( targets )


###################################################
### code chunk number 17: gbCounts (eval = FALSE)
###################################################
## gbcounts  <- gbCounts( features = features, 
##                            targets = targets, 
##                            minReadLength = 100, maxISize = 50000,
##                            libType="SE", 
##                            strandMode=0)


###################################################
### code chunk number 18: gbCountAccessors (eval = FALSE)
###################################################
## GeneCounts <- countsg(counts)
## GeneRd <- rdsg(counts)
## BinCounts <- countsb(counts)
## BinRd <- rdsb(counts)
## JunctionCounts <- countsj(counts)


###################################################
### code chunk number 19: gbCountsWrite (eval = FALSE)
###################################################
## writeCounts(counts=counts, output.dir = "example")
## writeRds(counts=counts, output.dir = "example")


###################################################
### code chunk number 20: gbCountAccessors2 (eval = FALSE)
###################################################
## e1iCounts <- countse1i(counts)
## ie2Counts <- countsie2(counts)


###################################################
### code chunk number 21: asd (eval = FALSE)
###################################################
## asd         <-  jCounts(counts = gbcounts, 
##                      features = features, 
##                      minReadLength = 100,
##                      libType="SE", 
##                      strandMode=0)


###################################################
### code chunk number 22: asdAccesor (eval = FALSE)
###################################################
## irPIR  <- irPIR( asd )
## altPSI <- altPSI( asd )
## esPSI  <- esPSI( asd )
## allBins      <- joint( asd )
## 
## junctionsPJU <- junctionsPIR( asd )
## junctionsPIR <- junctionsPIR( asd )


###################################################
### code chunk number 23: asdAccesor2 (eval = FALSE)
###################################################
## writeAS(as=asd, output.dir="example")


###################################################
### code chunk number 24: gbDUreport (eval = FALSE)
###################################################
##  gb <- gbDUreport( counts, 
##              minGenReads = 10, 
##              minBinReads = 5,
##              minRds = 0.05, 
##              contrast = NULL, 
##              ignoreExternal = TRUE, 
##              ignoreIo = TRUE, 
##              ignoreI = FALSE,
##              filterWithContrasted = TRUE,
##              verbose = TRUE,
##              formula = NULL,
##              coef = NULL)


###################################################
### code chunk number 25: asdAccesor (eval = FALSE)
###################################################
## geneX  <- genesDE( gb )
## binDU  <- binsDU( gb )


###################################################
### code chunk number 26: gbCountsWrite (eval = FALSE)
###################################################
## writeDU(gb, output.dir = "example")


###################################################
### code chunk number 27: jDUreport (eval = FALSE)
###################################################
## jdu <- jDUreport(asd, 
##             minAvgCounts                       = 5, 
##             contrast                           = NULL,
##             filterWithContrasted               = TRUE,
##             runUniformityTest                  = FALSE,
##             mergedBAMs                         = NULL,
##             maxPValForUniformityCheck          = 0.2,
##             strongFilter                       = TRUE,
##             maxConditionsForDispersionEstimate = 24,
##             formula                            = NULL,
##             coef                               = NULL,
##             maxFDRForParticipation             = 0.05,
##             useSubset                          = FALSE)


###################################################
### code chunk number 28: jDUreportAccesors (eval = FALSE)
###################################################
##  localej( jdu )
##  localec( jdu )
##  
##  anchorj( jdu )
##  anchorc( jdu )
##  
##  jir( jdu )
##  jes( jdu )
##  jalt( jdu )


###################################################
### code chunk number 29: splicingReportA (eval = FALSE)
###################################################
##  binbased( sr )
##  localebased( sr )
##  anchorbased( sr )


###################################################
### code chunk number 30: splicingReportA (eval = FALSE)
###################################################
## writeSplicingReport( sr, output.dir = "test")


###################################################
### code chunk number 31: integrateSignals (eval = FALSE)
###################################################
## is    <- integrateSignals(sr, asd,
##                  bin.FC = 3, bin.fdr = 0.05, bin.inclussion = 0.2,
##                  nonunif = 1, usenonunif = FALSE, 
##                  bjs.inclussion = 10, bjs.fdr = 0.01, 
##                  a.inclussion   = 0.3, a.fdr = 0.01, 
##                  l.inclussion   = 0.3, l.fdr = 0.01,
##                  otherSources = NULL, overlapType = "any")


###################################################
### code chunk number 32: integrateSignalsA (eval = FALSE)
###################################################
##  signals( is )
##  filters( is )


###################################################
### code chunk number 33: exportSplicingReports (eval = FALSE)
###################################################
## exportSplicingReports( sr, 
##                        output.dir="sr",
##                        openInBrowser = FALSE, 
##                        maxBinFDR = 0.2, 
##                        maxJunctionFDR = 0.2 )


###################################################
### code chunk number 34: exportIntegratedSignals (eval = FALSE)
###################################################
## exportIntegratedSignals( is, output.dir="is", 
##                          sr, counts, features, asd,
##                          mergedBams, 
##                          jCompletelyIncluded = FALSE, zoomRegion = 1.5, 
##                          useLog = FALSE, tcex = 1, ntop = NULL, openInBrowser = F, 
##                          makeGraphs = T, bforce=FALSE
##                         )


###################################################
### code chunk number 35: Ex02.a (eval = FALSE)
###################################################
## #library( ASpli )
## library( GenomicFeatures )
## gtfFileName <- aspliExampleGTF()
## genomeTxDb <- makeTxDbFromGFF( gtfFileName )
## features <- binGenome( genomeTxDb )


###################################################
### code chunk number 36: Ex02.e (eval = FALSE)
###################################################
## BAMFiles <- aspliExampleBamList()
## 
## targets <- data.frame( 
##   row.names = paste0('Sample',c(1:12)),
##   bam = BAMFiles,
##   f1 = c( 'A','A','A','A','A','A',
##           'B','B','B','B','B','B'),
##   f2 = c( 'C','C','C','D','D','D',
##           'C','C','C','D','D','D'),
##   stringsAsFactors = FALSE)


###################################################
### code chunk number 37: Ex02.g (eval = FALSE)
###################################################
## getConditions(targets)


###################################################
### code chunk number 38: Ex02.f (eval = FALSE)
###################################################
## mBAMs <- data.frame(bam      = sub("_[02]","",targets$bam[c(1,4,7,10)]),
##                     condition= c("A_C","A_D","B_C","B_D"))


###################################################
### code chunk number 39: Ex02.i (eval = FALSE)
###################################################
## gbcounts  <- gbCounts( features = features, 
##                            targets = targets, 
##                            minReadLength = 100, maxISize = 50000,
##                            libType="SE", 
##                            strandMode=0)


###################################################
### code chunk number 40: Ex02.j (eval = FALSE)
###################################################
## asd<- jCounts(counts = gbcounts, 
##                      features = features, 
##                      minReadLength = 100,
##                      libType="SE", 
##                      strandMode=0,
##                       threshold     = 5,
##                       minAnchor     = 10)


###################################################
### code chunk number 41: Ex02.k (eval = FALSE)
###################################################
## gb      <- gbDUreport(gbcounts,contrast = c( 1, -1, -1, 1 ) )
## 


###################################################
### code chunk number 42: Ex02.k2 (eval = FALSE)
###################################################
## genesDE(gb)[1:5,]
## binsDU(gb)[1:5,]


###################################################
### code chunk number 43: Ex02.l (eval = FALSE)
###################################################
## jdur    <- jDUreport(asd, contrast =  c( 1, -1, -1, 1 ))


###################################################
### code chunk number 44: Ex02.kk (eval = FALSE)
###################################################
## localec(jdur)[1:5,]


###################################################
### code chunk number 45: Ex02.kkk (eval = FALSE)
###################################################
## localej(jdur)[1:5,1:8]


###################################################
### code chunk number 46: Ex02.m (eval = FALSE)
###################################################
## sr      <- splicingReport(gb, jdur, counts = gbcounts)
## is      <- integrateSignals(sr,asd)


###################################################
### code chunk number 47: Ex03 (eval = FALSE)
###################################################
## exportIntegratedSignals( is, output.dir="aspliExample", 
##                          sr, gbcounts, features, asd, mBAMs)
## 


###################################################
### code chunk number 48: Ex04 (eval = FALSE)
###################################################
##  form <- formula(~f1+f2+f1:f2)


###################################################
### code chunk number 49: Ex05 (eval = FALSE)
###################################################
##  model.matrix(form,targets)


###################################################
### code chunk number 50: Ex02.k (eval = FALSE)
###################################################
## gb      <- gbDUreport(counts, formula = form , coef = 4)
## jdur    <- jDUreport(asd, formula = form, coef = 4 ,
##                      runUniformityTest = TRUE,
##                      mergedBams = mBAMs)


###################################################
### code chunk number 51: Ex02.kp (eval = FALSE)
###################################################
## targetPaired <- targets[c(1,4,7,10),]
## 
## gbcounts  <- gbCounts( features = features, 
##                            targets = targets, 
##                            minReadLength = 100, maxISize = 50000,
##                            libType="SE", 
##                            strandMode=0)
## asd    <- jCounts(counts = gbcounts, 
##                      features = features, 
##                      minReadLength = 100,
##                      libType="SE", 
##                      strandMode=0)
## form   <- formula(~f1+f2)
## gb     <- gbDUreport(gbcounts, formula = form)
## #jdur   <- jDUreport(asd    , formula = form)
## #sr <- splicingReport(gb, jdur, counts = counts)
## #is <- integrateSignals(sr,asd,bjs.fdr = 0.1, bjs.inclussion = 0.1, l.inclussion=0.001,l.fdr = 1)
## # exportSplicingReports(sr,output.dir="paired")
## #exportIntegratedSignals(is,output.dir="paired",sr,counts,features,asd,mBAMs,tcex=2)


