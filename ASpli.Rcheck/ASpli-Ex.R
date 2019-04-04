pkgname <- "ASpli"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ASpli')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ASpli-package")
### * ASpli-package

flush(stderr()); flush(stdout())

### Name: ASpli-package
### Title: Analysis of alternative splicing using RNAseq
### Aliases: ASpli-package ASpli
### Keywords: alternative splicing, RNA-seq, junctions

### ** Examples
 
  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes, bins and junctions 
  du       <- DUreport.norm( counts, targets )

  # Calculate PSI / PIR for bins and junction.
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )



cleanEx()
nameEx("AsDiscover-methods")
### * AsDiscover-methods

flush(stderr()); flush(stdout())

### Name: AsDiscover
### Title: Report PSI and PIR using experimental junctions
### Aliases: AsDiscover

### ** Examples


  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes, bins and junctions 
  du       <- DUreport.norm( counts, targets )

  # Calculate PSI / PIR for bins and junction.
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )
  
  writeAS( as = as, output.dir = "only_as" )




cleanEx()
nameEx("DUreport.norm")
### * DUreport.norm

flush(stderr()); flush(stdout())

### Name: DUreport.norm
### Title: Differential gene expression and differential bin usage
###   estimation
### Aliases: DUreport.norm

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  du       <- DUreport.norm( counts, targets )

  # Export results  
  writeDU( du = du, output.dir = "only_du" )



cleanEx()
nameEx("DUreport.offset")
### * DUreport.offset

flush(stderr()); flush(stdout())

### Name: DUreport.offset
### Title: Differential gene expression and differential bin usage
###   estimation
### Aliases: DUreport.offset

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  du       <- DUreport.offset( counts, targets )

  # Export results  
  writeDU( du = du, output.dir = "only_du" )



cleanEx()
nameEx("acc_AS")
### * acc_AS

flush(stderr()); flush(stdout())

### Name: AS accessors
### Title: Accessors for ASpliAS object
### Aliases: altPSI esPSI irPIR joint junctionsPIR junctionsPSI altPSI<-
###   esPSI<- irPIR<- joint<- junctionsPIR<- junctionsPSI<-

### ** Examples


# Accessing data tables from an ASpliAS object

as <- aspliASexample()

ap <- altPSI(as)
ep <- esPSI(as)
ip <- irPIR(as)
j  <- joint(as)
jpi <- junctionsPIR(as)
jps <- junctionsPSI(as)

# Setting data tables to an ASpliAS object

as2 <- new( 'ASpliAS' )

altPSI( as2 )  <- ap
esPSI( as2 )  <- ep
irPIR( as2 )  <- ip
joint( as2 )  <- j
junctionsPIR( as2 )  <- jpi
junctionsPSI( as2 )  <- jps




cleanEx()
nameEx("acc_DU")
### * acc_DU

flush(stderr()); flush(stdout())

### Name: DU accessors
### Title: Accessors for ASpliDU object
### Aliases: genesDE genesDE<- binsDU binsDU<- junctionsDU junctionsDU<-

### ** Examples

  
  # Get data tables from an ASpliDU object
  
  du <- aspliDUexample1()
  
  gde <- genesDE( du )
  bdu <- binsDU( du )
  jdu <- junctionsDU( du )
  
  # Set data tables to an ASpliDU object 
  
  genesDE( du )     <- gde
  binsDU( du )      <- bdu
  junctionsDU( du ) <- jdu




cleanEx()
nameEx("acc_counts")
### * acc_counts

flush(stderr()); flush(stdout())

### Name: Counts accesors
### Title: Accessors for ASpliCounts object
### Aliases: countsb countse1i countsg countsie2 countsj rdsg rdsb
###   countsb<- countse1i<- countsg<- countsie2<- countsj<- rdsg<- rdsb<-

### ** Examples


# Get data tables from an ASpliCounts object

counts <- aspliCountsExample() 

cb1  <- countsb(counts)
ce1i <- countse1i(counts)
cg   <- countsg(counts)
cie2 <- countsie2(counts)
cj   <- countsj(counts)
rg   <- rdsg(counts)
rb   <- rdsb(counts)

# Set data tables to an ASpliCounts object

countsb(counts)   <- cb1
countse1i(counts) <- ce1i
countsg(counts)   <- cg
countsie2(counts) <- cie2
countsj(counts)   <- cj
rdsg(counts)      <- rg
rdsb(counts)      <- rb




cleanEx()
nameEx("acc_features")
### * acc_features

flush(stderr()); flush(stdout())

### Name: Features accesors
### Title: Accessors for ASpliFeatures object
### Aliases: featuresb featuresg featuresj featuresb<- featuresg<-
###   featuresj<-

### ** Examples

  # Get data from an ASpliFeatures object
  
  features <- aspliFeaturesExample()
  
  fg <- featuresg( features )
  fb <- featuresb( features )
  fj <- featuresj( features )
  
  # Set data to an ASpliFeatures object
  
  featuresg( features ) <- fg 
  featuresb( features ) <- fb 
  featuresj( features ) <- fj 




cleanEx()
nameEx("binDUreport-method")
### * binDUreport-method

flush(stderr()); flush(stdout())

### Name: binDUreport
### Title: Differential gene expression and differential bin usage
###   estimation
### Aliases: binDUreport

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  du       <- binDUreport( counts, targets )

  # Export results  
  writeDU( du = du, output.dir = "only_du" )



cleanEx()
nameEx("binGenome")
### * binGenome

flush(stderr()); flush(stdout())

### Name: binGenome
### Title: Feature coordinates extraction
### Aliases: binGenome

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )

  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Extract gene, bin and junctions features
  GeneCoord <- featuresg(features)
  BinCoord <- featuresb(features)
  JunctionCoord <- featuresj(features)
  



cleanEx()
nameEx("containsJunctions-function")
### * containsJunctions-function

flush(stderr()); flush(stdout())

### Name:  Examine ASpliDU objects 
### Title: Examine ASpliDU objects
### Aliases: containsJunctions containsGenesAndBins

### ** Examples


  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Create an ASpliDU object.
  # 
  du       <- DUreport.norm( counts, targets )
  
  # Verify if du contains results for genes and bins.
  containsGenesAndBins( du )

  # Verify if du contains results for genes and bins.
  containsJunctions( du )
  
  



cleanEx()
nameEx("example_data")
### * example_data

flush(stderr()); flush(stdout())

### Name: Example data
### Title: Example Aspli objects
### Aliases: aspliASexample aspliBamsExample aspliCountsExample
###   aspliDUexample1 aspliDUexample2 aspliExampleBamList aspliExampleGTF
###   aspliFeaturesExample aspliJunctionDUexample aspliTargetsExample

### ** Examples

  as <- aspliASexample()
  bams <- aspliBamsExample()
  counts <- aspliCountsExample()
  du1 <- aspliDUexample1()
  du2 <- aspliDUexample2()
  bamfiles <- aspliExampleBamList()
  gtffile <- aspliExampleGTF()
  features <- aspliFeaturesExample()
  jdu <- aspliJunctionDUexample()
  targets <-aspliTargetsExample()



cleanEx()
nameEx("filterDU-method")
### * filterDU-method

flush(stderr()); flush(stdout())

### Name: filterDU
### Title: Filtering ASpliDU objects
### Aliases: filterDU

### ** Examples


 # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of junctions only 
  du       <- DUreport.norm( counts, targets )
  
  # Filter by FDR
  duFiltered1 <- filterDU( du, what=c('genes','bins'), 
     fdr = 0.01 )


  # Filter by logFC, only those that were up-regulated 
  duFiltered2 <- filterDU( du, what=c('genes','bins'), 
     logFC = log( 1.5, 2 ),  absLogFC = FALSE )




cleanEx()
nameEx("getConditions-function")
### * getConditions-function

flush(stderr()); flush(stdout())

### Name: getConditions
### Title: Retrieve condition names from a targets data frame.
### Aliases: getConditions

### ** Examples

  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files.
  # Return value is c('C', 'D') in this example. 
  conditions <- getConditions(targets)


  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", "A_C_3.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam", "A_D_3.bam"  )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:8)),
               bam = file.path( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','C','D','D','D','D'),
               factor2 = c( 'E','E','F','F','E','E','F','F')  )
  
  # Load reads from bam files.
  # Return value is c("C_E", "C_F", "D_E", "D_F") in this example. 
  conditions <- getConditions(targets)




cleanEx()
nameEx("junctionDUreport")
### * junctionDUreport

flush(stderr()); flush(stdout())

### Name: junctionDUreport
### Title: Differential junction usage estimation
### Aliases: junctionDUreport

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  
  #as <- new("ASpliAS")
  #targets <- 1
  #a <- junctionDUreportExt(as, targets)
  



cleanEx()
nameEx("loadBAM")
### * loadBAM

flush(stderr()); flush(stdout())

### Name: loadBAM
### Title: Load BAM files
### Aliases: loadBAM

### ** Examples


  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
               
  # Load reads from bam files
  bams <- loadBAM( targets )
  



cleanEx()
nameEx("mergeBinDUAS-method")
### * mergeBinDUAS-method

flush(stderr()); flush(stdout())

### Name: mergeBinDUAS
### Title: Differential usage of bins and PSI/PIR.
### Aliases: mergeBinDUAS

### ** Examples


  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam",
                     "B_C_0.bam", "B_C_1.bam", "B_C_2.bam", 
                     "B_D_0.bam", "B_D_1.bam", "B_D_2.bam" )
                     
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:12)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'A','A','A','A','A','A','B','B','B','B','B','B'),
               factor2 = c( 'C','C','C','D','D','D','C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  du <- DUreport.norm( counts, targets , contrast = c(1,-1,-1,1))
  
  # Calculate PSI / PIR for bins and junction.
  as <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )

  mas <- mergeBinDUAS( du, as, targets, contrast =  c(1,-1,-1,1) )                     
  



cleanEx()
nameEx("mergeReports")
### * mergeReports

flush(stderr()); flush(stdout())

### Name: mergeReports
### Title: Differential junction usage estimation. COMPLETAR
### Aliases: mergeReports

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored.
  
  #as <- new("ASpliAS")
  #targets <- 1
  #a <- junctionDUreportExt(as, targets)
  



cleanEx()
nameEx("plotBins")
### * plotBins

flush(stderr()); flush(stdout())

### Name: plotBins
### Title: Draw plots of gene counts, bin counts, PSI/PIR value, inclusion
###   and exclusion junctions for selected bins.
### Aliases: plotBins

### ** Examples


  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam",
                     "B_C_0.bam", "B_C_1.bam", "B_C_2.bam", 
                     "B_D_0.bam", "B_D_1.bam", "B_D_2.bam" )
                     
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:12)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'A','A','A','A','A','A','B','B','B','B','B','B'),
               factor2 = c( 'C','C','C','D','D','D','C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes, bins and junctions 
  du       <- DUreport.norm( counts, targets , contrast = c(1,-1,-1,1))

  # Calculate PSI / PIR for bins and junction.
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )
  
  # Plot bin data. Factor2 is the main factor for graphic representation in
  # this example as it is the first in factorsAndValues argument.
  # This makes a bar plot comparing four conditions, grouped by factor1.
  plotBins( counts, as, 'GENE03:E002', 
    factorsAndValues = list( 
      factor2 = c('C','D'), 
      factor1 = c('A','B') ),
    las.x = 1,
    legendAtSide = TRUE,
    useHCColors = TRUE,   
    targets = targets,
    barWidth = 0.95,
    innerMargins = c( 2.1, 4.1, 1.1, 1.1 ) )
    
    
  # Redefine targets  
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:12)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'A','A','B','B','C','C','D','D','E','E','F','F') )
  
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )
  
  # This makes a line plot for six conditions, grouped by factor1.                       
  plotBins( counts, as, 'GENE03:E002', 
    factorsAndValues = list( 
      factor1 = c('A','B','C','D','E','F') ),
    las.x = 1,
    legendAtSide = FALSE,
    targets = targets,
    innerMargins = c( 2.1, 4.1, 1.1, 1.1 ) )                        
  



cleanEx()
nameEx("plotGenomicRegions")
### * plotGenomicRegions

flush(stderr()); flush(stdout())

### Name: plotGenomicRegions
### Title: Create genomic regions coverage plots
### Aliases: plotGenomicRegions

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D'),
               stringsAsFactors = FALSE )
  
  # Plot a single bin to a window
  plotGenomicRegions( 
    features, 
    'GENE01:E002', 
    genomeTxDb, 
    targets, 
    sashimi = FALSE,
    colors = '#AA4444', 
    annotationHeight = 0.1, 
    tempFolder = 'tmp', 
    verbose = TRUE , 
    avoidReMergeBams = FALSE, 
    useTransparency = FALSE ) 
    
  # plot two bins to pdf files.
  plotGenomicRegions( 
    features, c( 'GENE01:E002', 'GENE02:E002' ), 
    genomeTxDb, 
    targets, 
    layout = matrix( c( 'C', 'D'), ncol = 1),
    colors = matrix( c( '#663243', '#363273'), ncol = 1),
    plotTitles = matrix( c( 'C condition', 'D condition'), ncol = 1),
    sashimi = FALSE,
    mainFontSize = 12,
    annotationHeight = 0.1, 
    tempFolder = 'tmp', 
    verbose = TRUE , 
    avoidReMergeBams = FALSE, 
    useTransparency = TRUE,
    outfolder = '.',
    outfileType = 'pdf',
    deviceOpt = list( height = 6, width = 5, paper = 'a4r' ) ) 



cleanEx()
nameEx("rds")
### * rds

flush(stderr()); flush(stdout())

### Name: rds
### Title: Read density of gene and bins
### Aliases: rds

### ** Examples

    # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
                          
  # Calculates read densities
  counts <- rds( counts, targets )



cleanEx()
nameEx("readCounts")
### * readCounts

flush(stderr()); flush(stdout())

### Name: readCounts
### Title: Summarize read overlaps
### Aliases: readCounts

### ** Examples

  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  # Export data
  writeCounts( counts, output.dir = "only_counts" )
  



cleanEx()
nameEx("subset-methods")
### * subset-methods

flush(stderr()); flush(stdout())

### Name: Subset ASpli objects
### Title: Subset ASpli objects
### Aliases: subsetBams subsetTargets subset

### ** Examples


  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Create ASpliAS object                        
  as       <- AsDiscover( counts, targets, features, bams, readLength = 100, 
                          threshold = 5, cores = 1 )
                          
  # Define selection
  select <- c('Sample_1', 'Sample_2', 'Sample_4', 'Sample_5')                       
  
  # Subset target 
  targets2 <- subsetTargets( targets, select )
  
  # Subset bams 
  bams2 <- subsetBams( bams, targets, select )
  
  # Subset ASpliCounts object 
  counts2 <- subset( counts, targets, select )
  
  # Subset ASpliAS object 
  as2 <- subset( as, targets, select )
  



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
