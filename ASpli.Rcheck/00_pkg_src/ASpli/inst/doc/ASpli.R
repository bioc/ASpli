### R code from vignette source 'ASpli.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()
options(continue=" ")


###################################################
### code chunk number 2: installation (eval = FALSE)
###################################################
## source("https://www.bioconductor.org/biocLite.R")
## biocLite("ASpli")


###################################################
### code chunk number 3: loadASpli (eval = FALSE)
###################################################
## library(ASpli)


###################################################
### code chunk number 4: makeTx (eval = FALSE)
###################################################
## library(GenomicFeatures)
## TxDb <- makeTxDbFromGFF(
##   file="genes.gtf",
##   format="gtf")


###################################################
### code chunk number 5: binGenome (eval = FALSE)
###################################################
## library(GenomicFeatures)
## annFile       <- aspliExampleGTF()
## aTxDb         <- makeTxDbFromGFF(annFile)
## # extract features from annotation
## features      <- binGenome( aTxDb ) 
## # Accesors of ASpliFeatures object 
## geneCoord     <- featuresg( features )
## binCoord      <- featuresb( features )
## junctionCoord <- featuresj( features )
## # Extract metadata annotation, such as bins names, event and feature type
## binMetadata   <- mcols( binCoord )


###################################################
### code chunk number 6: binGenome (eval = FALSE)
###################################################
## # extract features and map symbols to genes
## symbols       <- data.frame( row.names = genes( aTxDb ), 
##                              symbol = paste( 'This is symbol of gene:',
##                                              genes( aTxDb ) ) )
## features      <- binGenome( aTxDb, geneSymbols = symbols ) 
## features


###################################################
### code chunk number 7: targetsDF
###################################################
bamFiles <- c( "Ct1.bam", "Ct2.bam", "Mut1.bam","Mut2.bam" )
targets <- data.frame( row.names =  c("CT_rep1","CT_rep2", "Mut_rep1", "Mut_rep2"),
                       bam = bamFiles,
                       genotype = c("CT","CT", "Mut", "Mut") , 
                       stringsAsFactors = FALSE )
targets


###################################################
### code chunk number 8: targetsDF2
###################################################
bamFiles <- c("Ct_time1_rep1.bam", "Ct_time1_rep2.bam",
              "Ct_time2_rep1.bam", "Ct_time2_rep2.bam",
              "Mut_time1_rep1.bam","Mut_time1_rep2.bam",
              "Mut_time2_rep1.bam","Mut_time2_rep2.bam")
targets_2 <- data.frame( row.names = c( 'CT_t1_r1',  'CT_t1_r2',
                                        'CT_t2_r1',  'CT_t2_r2',
                                        'Mut_t1_r1', 'Mut_t1_r2',
                                        'Mut_t2_r1', 'Mut_t2_r2' ),
                         bam = bamFiles,
                         genotype = c( 'CT',  'CT',  'CT',  'CT', 
                                       'MUT', 'MUT', 'MUT', 'MUT' ),
                         time     = c( 't1', 't1', 't2', 't2', 
                                       't1', 't1', 't2', 't2' ),
                         stringsAsFactors = FALSE )
targets_2


###################################################
### code chunk number 9: targetsDF2 (eval = FALSE)
###################################################
## # Show the name of each unique condition in the console
## getConditions( targets_2 )


###################################################
### code chunk number 10: loadBam (eval = FALSE)
###################################################
## targets <- aspliTargetsExample()
## bam <- loadBAM(targets)


###################################################
### code chunk number 11: readCounts (eval = FALSE)
###################################################
## counts <- readCounts ( 
##     features, 
##     bam, 
##     targets, 
##     cores = 1, 
##     readLength = 100L, 
##     maxISize = 50000,
##     minAnchor = 10 )


###################################################
### code chunk number 12: readCountAccessors (eval = FALSE)
###################################################
## # Accessing count and read density data
## GeneCounts <- countsg(counts)
## GeneRd <- rdsg(counts)
## 
## BinCounts <- countsb(counts)
## BinRd <- rdsb(counts)
## 
## JunctionCounts <- countsj(counts)
## 
## # Export count data to text files 
## writeCounts(counts=counts, output.dir = "example")
## writeRds(counts=counts, output.dir = "example")


###################################################
### code chunk number 13: readCountAccessors2 (eval = FALSE)
###################################################
## # Accessing counts for intron flanking regions
## e1iCounts <- countse1i(counts)
## ie2Counts <- countsie2(counts)
## 


###################################################
### code chunk number 14: AsDiscover (eval = FALSE)
###################################################
## as <- AsDiscover( counts, targets, features, bam, readLength=100L, 
##                   threshold = 5)


###################################################
### code chunk number 15: writeAS (eval = FALSE)
###################################################
## # Access result from an ASpliAS object 
## irPIR  <- irPIR( as )
## altPSI <- altPSI( as )
## esPSI  <- esPSI( as )
## junctionsPIR <- junctionsPIR( as )
## junctionsPSI <- junctionsPSI( as )
## allBins      <- joint( as )
## # Export results to files
## writeAS(as=as, output.dir="example")


###################################################
### code chunk number 16: export (eval = FALSE)
###################################################
## contrast <- c( -1, 1 ) 


###################################################
### code chunk number 17: write (eval = FALSE)
###################################################
## du <- DUreport( counts,
##                 targets,
##                 minGenReads = 10, 
##                 minBinReads = 5,
##                 minRds = 0.05, 
##                 offset = FALSE, 
##                 offsetAggregateMode =  c( "geneMode", "binMode" )[2],
##                 offsetUseFitGeneX = TRUE,
##                 contrast = NULL, 
##                 forceGLM = FALSE, 
##                 ignoreExternal = TRUE, 
##                 ignoreIo = FALSE, 
##                 ignoreI = FALSE,
##                 filterWithConstrasted = FALSE,
##                 verbose = FALSE )


###################################################
### code chunk number 18: export (eval = FALSE)
###################################################
## du <- DUreportBinSplice( counts, 
##                          targets, 
##                          minGenReads = 10, 
##                          minBinReads = 5,
##                          minRds = 0.05, 
##                          contrast = NULL, 
##                          forceGLM = FALSE, 
##                          ignoreExternal = TRUE, 
##                          ignoreIo = TRUE, 
##                          ignoreI = FALSE,
##                          filterWithContrasted = FALSE )


###################################################
### code chunk number 19: export (eval = FALSE)
###################################################
## du <- junctionDUReport(  counts, 
##                          targets, 
##                          appendTo = NULL, 
##                          minGenReads = 10,
##                          minRds = 0.05,
##                          threshold = 5,
##                          offset   = FALSE,
##                          offsetUseFitGeneX = TRUE,
##                          contrast = NULL,
##                          forceGLM = FALSE )


###################################################
### code chunk number 20: write (eval = FALSE)
###################################################
## writeCounts(counts, "example_counts")
## writeDU(du, output.dir="example_du");
## writeAS(as=as, output.dir="example_as");
## writeAll(counts=counts, du=du, as=as, output.dir="example_all")


###################################################
### code chunk number 21: targetsPlot (eval = FALSE)
###################################################
## library(GenomicFeatures)
## bamfiles <- system.file( 'extdata', 
##     c('A_C_0.bam', 'A_C_1.bam', 'A_C_2.bam',
##       'A_D_0.bam', 'A_D_1.bam', 'A_D_2.bam',
##       'B_C_0.bam', 'B_C_1.bam', 'B_C_2.bam',
##       'B_D_0.bam', 'B_D_1.bam', 'B_D_2.bam' ),
##     package = "ASpli") 
## 
## targets <- data.frame(
##   row.names = c( 'A_C_0', 'A_C_1', 'A_C_2',
##                  'A_D_0', 'A_D_1', 'A_D_2',
##                  'B_C_0', 'B_C_1', 'B_C_2',
##                  'B_D_0', 'B_D_1', 'B_D_2' ),
##   bam = bamfiles,
##   f1  = rep( c('A','B'), each=6 ),
##   f2  = rep( c('C','D'), 2, each=3 ),
##   stringsAsFactors = FALSE )
## 
## genomeTxDb <- makeTxDbFromGFF( system.file( 'extdata', 'genes.mini.gtf',
## package='ASpli'))
## features <- binGenome( genomeTxDb )
## 
## # Draw a single plot on a window
## # xIsBin arguments is used to specify if names given are bins or genes.
## # When multiple names are given, all must correspond to bins or correspond to
## # genes. Mixed names are not allowed.
## 
## plotGenomicRegions( features, 'GENE02:E002' , genomeTxDb, targets, 
##     xIsBin = TRUE, verbose = TRUE)
## 
## # Draw a single plot on a window, with a custom layout
## plotGenomicRegions( features, 'GENE02:E002' , genomeTxDb, targets, xIsBin =
##    TRUE, layout = matrix( c('A_C','A_D','B_C','B_D') ,nrow = 1 ), 
##    verbose = TRUE, color = '#AA8866')
## 
## # Draw many plots on files.
## # Argument outfileType define the file format of the generated files.
## # Argument useTransparency specified that transparency can be used to draw the
## # plot. No all graphic devices support transparency.
## # Argument annotationHeight specifies the proportional height of the plot used
## # to draw the annotation.
## # Argument deviceOpt gives additional arguments (width and height in this
## # example) to the underlying graphic device (pdf in this example)
## # One of the plots generated here is shown in a figure below.
## plotGenomicRegions( features, c( 'GENE02:E002', 'GENE01:E002', 'GENE02:E003' ), 
##     genomeTxDb, targets, xIsBin = TRUE, outfolder = 'grplots', verbose = TRUE,
##     color = '#334466', outfileType='pdf', useTransparency = TRUE,
##     annotationHeight = 0.12, deviceOpt = c( width = 8, height = 6) )
##     
## # Draw a single plot on a window, with premerged bams
## mergedBams <- data.frame(
##   row.names = c( 'A_C','A_D','B_C','B_D' ),
##   bam = c( 'A_C.bam',  # Warning, these files do not exists.
##            'A_D.bam',  # Bam files should be indexed.
##            'B_C.bam',
##            'B_D.bam' ))
## plotGenomicRegions( features, 'GENE02:E002' , genomeTxDb, targets,
##    xIsBin = TRUE, verbose = TRUE, preMergedBAMs = mergedBams )    


###################################################
### code chunk number 22: plotbins (eval = FALSE)
###################################################
## # Defines an experiment with one experimental factor (genotype, with two
## # values: wild-type (WT) and mutant (MT) ), and three replicate samples for each
## # condition.
## targets <- data.frame(
##   row.names = c( 'WT1', 'WT2', 'WT3', 
##                  'MT1', 'MT2', 'MT3' ),
##   bam = c( 'WT1.bam', 'WT2.bam', 'WT3.bam', 
##            'MT1.bam', 'MT2.bam', 'MT3.bam' ),
##   genotype = c( 'WT','WT','WT','MT','MT','MT' ),
##   stringsAsFactors = FALSE )
## # Specifies what factors and values will be plotted, in this example all of
## # them.
## fv = list( genotype = c( 'WT', 'MT' ) ) 
## plotbins(
##         counts, 
##         as,
##         'GENE02:E002', 
##         factorsAndValues = fv, 
##         targets )


###################################################
### code chunk number 23: plotbins (eval = FALSE)
###################################################
## # Defines an experiment with two experimental factor (treat, with four
## # values: A, B, C and D; and time, with twelve values from T1 to T12.), and
## # two replicate samples for each condition. In the definition, many values are
## # not shown to reduce the space of the example, and are replaced by '...' . 
## targets <- data.frame(
##   row.names = c( 'A_T1.r1', 'A_T1.r2', 
##                  'A_T2.r1', 'A_T2.r2',
##                  ... , 
##                  'D_T12.r1', 'D_T12.r2' ),
##   bam = c( 'A_T1.r1.bam', 'A_T1.r2.bam', 
##            'A_T2.r1.bam', 'A_T2.r2.bam',
##            'D_T12.r1.bam', 'D_T12.r2.bam' ),
##   treat = c( 'A', 'A',
##              'A', 'A',
##              ... ,
##              'D', 'D' ),
##   time = c( 'T1', 'T1',
##             'T2', 'T2',
##             ... ,
##             'T12', 'T12' ),
##   stringsAsFactors = FALSE )
## # Draw the plots.
## plotbins(
##         counts, 
##         as,
##         'GENE02:E002', 
##         factorsAndValues = fv, 
##         targets )
## 
## # Specifies what factors and values will be plotted. In this example there are
## # two factor, the first one in 'fv' list is the main factor, values will
## # be grouped by this factor. 
## fv = list( time = c( 'T1', 'T2', ... , 'T12' ) ,
##            treat = c( 'A', 'B', 'C', 'D' ) ) 
## # Draw the plots
## plotbins(
##         counts, 
##         as,
##         'GENE02:E002', 
##         factorsAndValues = fv, 
##         targets )


###################################################
### code chunk number 24: containsDUfunctions (eval = FALSE)
###################################################
## du <- aspliDUexample1()
## containsJunctions( du )
## containsGenesAndBins( du )


###################################################
### code chunk number 25: filterDU1 (eval = FALSE)
###################################################
## # Filter genes that changes at least a 50 percent their expression with a 
## # FDR value of 0.05 or less, and bins changes at least a 50 percent their expression with a 
## # FDR value of 0.05 
## du <- aspliDUexample1()
## duFilt <- filterDU(du, 'genes', '0.05', log(1.5,2))
## duFilt <- filterDU(duFilt, 'bins', '0.1', log(1.5,2))


###################################################
### code chunk number 26: filterDU2 (eval = FALSE)
###################################################
## # Filter genes that reduces their expression to the half
## du <- aspliDUexample1()
## duFilt <- filterDU( du, 'genes', '0.05', -1, absLogFC = FALSE,
##                     logFCgreater = FALSE  ) 


###################################################
### code chunk number 27: subset (eval = FALSE)
###################################################
## # Subset counts
## counts <- aspliCountsExample()
## countsSmall <- subset(counts, targets, select = c("A_C", "A_D") )
## # Subset AS results. Note that PIR/PSI metrics are not recalculated.
## as <- aspliASexample()
## asSmall <- subset(as, targets, select = c("A_C_0", "A_C_1") )
## # Subset targets
## targets      <- aspliTargetsExample()
## targetsSmall <- subsetTargets( targets, c("A_C", "A_C") )
## # Subset BAMs
## bams     <- aspliBamsExample()
## bamSmall <- subsetBams( bams, targets, c("A_C", "A_C") )


###################################################
### code chunk number 28: mergeBinDUAS (eval = FALSE)
###################################################
## du <- aspliDUexample1()
## as <- aspliASexample()
## mergeBinDUAS(du,as, targets)


###################################################
### code chunk number 29: Ex01.a (eval = FALSE)
###################################################
## 
## library( GenomicFeatures )
## # gtfFileName contains the full path to the example gtf file in your system.
## gtfFileName <- aspliExampleGTF()
## 
## # Create a TxDb object using makeTxDbFromGFF function from GenomicFeatures 
## # package
## genomeTxDb <- makeTxDbFromGFF( gtfFileName )
## 
## # Extract the genomic features
## features <- binGenome( genomeTxDb )
## 


###################################################
### code chunk number 30: Ex01.b (eval = FALSE)
###################################################
## 
## # bamFiles contains the full path of the bam files in your system
## bamFiles <- aspliExampleBamList( )
## 
## # Define the targets
## targets <- data.frame( 
##   row.names = paste0('Sample',c(1:12)),
##   bam = bamFiles,
##   f1 = c( 'A','A','A','A','A','A',
##           'B','B','B','B','B','B'),
##   f2 = c( 'C','C','C','D','D','D',
##           'C','C','C','D','D','D'),
##   stringsAsFactors = FALSE
## )
## 
## # Examinate the condition names
## getConditions(targets)
## 
## # Load the bam files
## bams = loadBAM(targets)


###################################################
### code chunk number 31: Ex01.c (eval = FALSE)
###################################################
## counts <- readCounts( features, bams, targets, readLength = 100, 
##                       maxISize = 5000 )


###################################################
### code chunk number 32: Ex01.d (eval = FALSE)
###################################################
## # DU/DE analysis for genes and bins.
## du <- DUreport( counts, targets, contrast = c( 1, -1, -1, 1 ) )
## # DU analysis for junctions
## du <- junctionDUReport( counts, targets, appendTo = du, 
##                         contrast = c( 1, -1, -1, 1 ) )


###################################################
### code chunk number 33: Ex01.e (eval = FALSE)
###################################################
## # AS analysis 
## as <- AsDiscover( counts, targets, features, bams, readLength = 100)


###################################################
### code chunk number 34: Ex01.f (eval = FALSE)
###################################################
## # Select top tags from DU
## topTagsBins <- filterDU( du, what = 'bins', fdr = 0.05, logFC = log(1.5,2) )
## # DU results are merged with AS results
## topTagsBins <- mergeBinDUAS( du, as, targets, contrast = c( 1, -1, -1, 1 ) )
## # We filter those top tags that have also junction evidence that supports the
## # differential usage, a change of at least 0.1 between observed and expected 
## # in the PSI or PIR metric is required.
## topTagsBins <- topTagsBins[ abs(topTagsBins[,'delta' ]) > 0.10,]


###################################################
### code chunk number 35: Ex01.f (eval = FALSE)
###################################################
## plotGenomicRegions( features,'GENE10:E002', genomeTxDb, targets )


###################################################
### code chunk number 36: Ex02.a (eval = FALSE)
###################################################
## library(RNAseqData.HNRNPC.bam.chr14)


###################################################
### code chunk number 37: Ex02.b (eval = FALSE)
###################################################
## chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
## genome <- loadDb(chr14)


###################################################
### code chunk number 38: Ex02.c (eval = FALSE)
###################################################
## features <- binGenome(genome) 


###################################################
### code chunk number 39: Ex02.d (eval = FALSE)
###################################################
## targets <- data.frame( 
##               bam = RNAseqData.HNRNPC.bam.chr14_BAMFILES,
##               treat = c( "CT", "CT", "CT", "CT", 
##                          "KD", "KD", "KD", "KD" ),
##               stringsAsFactors = FALSE )


###################################################
### code chunk number 40: Ex02.e (eval = FALSE)
###################################################
## bam <- loadBAM(targets)


###################################################
### code chunk number 41: Ex02.f (eval = FALSE)
###################################################
## counts <- readCounts( features, bam, targets, readLength=100L, maxISize=50000 )


###################################################
### code chunk number 42: Ex02.g (eval = FALSE)
###################################################
## du <- DUreport(counts, targets)


###################################################
### code chunk number 43: Ex02.g (eval = FALSE)
###################################################
## as <- AsDiscover (counts, targets, features, bam, readLength=100L, threshold = 5, cores = 1 )


###################################################
### code chunk number 44: Ex02.h (eval = FALSE)
###################################################
## duFiltered <- filterDU( du, 'bins', fdr = 0.05, logFC = log(1.5,2) ) 
## merged     <- mergeBinDUAS( duFiltered, as, targets )


###################################################
### code chunk number 45: Ex02.i (eval = FALSE)
###################################################
## # Show the complete gene, with the bin hightlighted
## plotGenomicRegions( features, 'ZNF410:E013', genome, targets, verbose=TRUE,
##                     colors = matrix(c('#223344','#223344'), ncol=2 ), outfileType='pdf', outfolder='grplots')
## # Show a zoom in in the bin
## plotGenomicRegions( features, 'ZNF410:E013', genome, targets, verbose=TRUE,
##                     colors = matrix(c('#223344','#223344'), ncol=2 ),
##                     zoomOnBins = 0.05)


