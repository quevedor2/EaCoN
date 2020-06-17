#### Libraries ####
## Instal ASCAT and FACETS
install.pkg <- FALSE
if(install.pkg){
  devtools::install_github("quevedor2/EaCoN")
  devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
  devtools::install_github("mskcc/facets")
  
  ## Instal Bioconductor packages
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("affxparser", "Biostrings", "aroma.light", "BSgenome", "copynumber", "GenomicRanges", "limma", "rhdf5", "sequenza"))
  
  ## Install the most recent STABLE version (@master)
  #devtools::install_github("gustaveroussy/EaCoN")
  
  # The affy.CN.norm package provides pre-computed GC% and wave-effect 
  # (re)normalization datasets for all compatible Affymetrix designs, 
  # for both NA33/NA35 (hg19) and NA36 (hg38) human genome builds.
  install.packages("https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz", repos = NULL, type = "source")
  
  # Contains the embedded APT tool
  devtools::install_github("gustaveroussy/apt.snp6.1.20.0")
  
  # Annotations for SNP6.0
  install.packages("https://partage.gustaveroussy.fr/pydio_public/152397?dl=true&file=/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz", repos = NULL, type = "source")
  
  # Rawcopy subset for BAF normalization
  install.packages( "https://partage.gustaveroussy.fr/pydio_public/e6fe22?dl=true&file=/rcnorm_0.1.5.tar.gz", repos = NULL, type = "source")
}

##############
#### Main ####

#devtools::install_github("quevedor2/EaCoN", ref = 'tads')
library(optparse)
library(EaCoN)
library(parallel)
library(foreach)
library(Biobase)
library(foreach)

option_list <- list(
  make_option(c("-i", "--idx"), type="integer", default=NULL,
              help="Index of sample group [default= %default]", metavar="integer"),
  make_option(c("-g", "--grpsize"), type="integer", default=10,
              help="Size of the groups to run in one job [default= %default]", metavar="integer"),
  make_option(c("-d", "--dataset"), type="character", default='CCLE',
              help="Dataset to use, either 'GDSC' or 'GNE' or 'CCLE'"),
  make_option(c("-p", "--pdir"), type="character", default='/mnt/work1/users/pughlab/projects/cancer_cell_lines',
              help="Parent directory path"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


.splitSamples <- function(sample.ids, range.idx, grp.size=10){
  if(((grp.size * range.idx) - grp.size) > length(sample.ids)){
    NULL 
  } else {
    if(grp.size > 1){
      range.s <- seq(1, length(sample.ids), by=grp.size)
      range.e <- c((range.s + (grp.size-1)), length(sample.ids))
      
      range.se <- c(range.s[range.idx], range.e[range.idx])
      print(paste0(range.idx, " : ", paste(range.se, collapse="-")))
      sample.ids[range.se[1]:range.se[2]]
    } else {
      sample.ids[range.idx]
    }
  }
}

.removeRedundantFiles <- function(pattern1, pattern2, unlink.path=NULL, 
                                  astep="L2R", segmenter='ASCAT', overwrite=FALSE,
                                  unlink.stat=TRUE){
  #p1.files <- list.files(pattern=pattern1, recursive = TRUE, full.names = TRUE)
  dirs <- function(f, astep, segmenter='ASCAT'){
    switch(astep,
         "L2R"={
           read.dir <- file.path(f)
           out.dir <- file.path(f, toupper(segmenter), astep)
         },
         "ASCN"={
           read.dir <- file.path(f, toupper(segmenter), 'L2R')
           out.dir <- file.path(f, toupper(segmenter), astep)
         })
    return(list("read"=read.dir, "out"=out.dir))
  }
  p1.files <- as.character(sapply(list.files(), function(f){
    list.files(path = dirs(f, astep, segmenter)$read, pattern=pattern1, full.names = TRUE)
  }))
  p2.files <- as.character(sapply(list.files(), function(f){
    list.files(path = dirs(f, astep, segmenter)$out, pattern=pattern2, recursive=FALSE, full=T)
  }))
  zero1.idx <- sapply(p1.files, function(i) i == "character(0)") 
  zero2.idx <- sapply(p2.files, function(i) i == "character(0)") 
  if(any(zero1.idx)) p1.files <- p1.files[-which(zero1.idx)]
  if(any(zero2.idx)) p2.files <- p2.files[-which(zero2.idx)]
  
  if(overwrite){
    rm.ids <- sapply(strsplit(p2.files, "/"), function(x) x[[1]])
    rm.idx <- unlist(sapply(rm.ids, function(x) grep(paste0("\\/", x, "\\/"), p1.files)))
    if(length(rm.idx) > 0){
      p.files <- p1.files[-rm.idx]
    } else {
      p.files <- p1.files
    }
  } else {
    p2.ids <- sapply(strsplit(p2.files, "/"), function(x) x[[1]])
    p1.ids <- sapply(strsplit(p1.files, "/"), function(x) x[[1]])
    p.ids <- setdiff(p1.ids, p2.ids)
    p.files <- as.character(setNames(p1.files, p1.ids)[p.ids])
  }
  
  
  if(!is.null(unlink.path)){
    print(paste0("Unlinking files in path: <Sample>/", unlink.path))
    files.to.unlink <- sapply(strsplit(p.files, "/"), function(x) x[[1]])
    sapply(files.to.unlink, function(ftu){
      print(file.path(ftu, unlink.path))
      if(unlink.stat) unlink(x = file.path(ftu, unlink.path), recursive = T)
    })
  }
  p.files
}


##################
#### 0) Setup ####
segmenter <- 'ASCAT'
cl <- 10 # Number of cores
run_parallel <- FALSE 
dataset <- opt$dataset  #'GDSC', 'CCLE', 'GNE', '1000G'
pdir <- file.path(opt$pdir, dataset)
CEL.dir <- file.path(pdir, "data")

message(paste0("Searching for CEL$ files in: ", file.path(pdir, "data")))
sample.paths <- list.files(CEL.dir, pattern="CEL$", recursive = T, 
                         ignore.case = T, full.names = T)
if(any(grepl("symlinks", c("a", "b", "symlinks")))) sample.paths <- sample.paths[-grep("symlinks", sample.paths)]

sample.ids <- gsub(".cel$", "", basename(sample.paths), ignore.case = TRUE)
data.dir <- dirname(sample.paths)
regm <- regexpr(".cel", basename(sample.paths), ignore.case=T)
cel.suffix <- regmatches(x = basename(sample.paths), m = regm)


dir.create(file.path(pdir, "eacon"), recursive = TRUE)
setwd(file.path(pdir, "eacon"))

#################################
#### 1) Signal Normalization ####
# For Affy6 analyses, this pipeline uses the 'rawcopy' method of
# signal normalization. 
#   APT:    v1.20.0
#   Genome build:   na35.r1
# Using the apt.snp6.1.20.0::apt.snp6.process() function, 
# this creates an intermediatry CNCHP file using the apt-copynumber-workflow
# in the APT toolkit (https://www.affymetrix.com/support/developer/powertools/changelog/apt-copynumber-workflow.html)
# This is then converted into a OSCHP format for standardization purposes using
# the apt2-dset-util APT tool (https://media.affymetrix.com/support/developer/powertools/changelog/apt2-dset-util.html).
# This function employs the  standard APT tools which are embedded 
# in the apt.snp6.1.20.0 R package, where the annotations are found in the
# GenomeWideSNP.6.na35.r1_0.1.0 R package.
#
# The OSCHP is read in to the rcnorm package where after a little metadata
# addition, the BAF is normalized using rcnorm::rcnorm.snp(). The L2R
# undergoes wave and GC-normalization using EaCoN:::renorm.go() function, 
# which uses the pre-computed GC% and wave-effect (re)normalization datasets 
# stored in the affy.CN.norm.data_0.1.2 package. Post-normalized L2Rs are then
# median centered.
#
# Finally, an ASCAT-like object is built to standardize the analyses for all 
# subsequent steps between different data types.

# sample.ids <- setdiff(sample.ids, list.files())
if(dataset=='GNE'){
  illumina.dir <- file.path(pdir, 'illumina', 'GNE_Matrices')
  EaCoN::Build.OMNI25(illumina.dir=illumina.dir, parent.dir=pdir)
} else {
  if(run_parallel){
    mclapply(sample.ids, function(sample){
      setwd(file.path(pdir, "eacon"))
      print(paste0(sample, ": (", match(sample, sample.ids), "/", 
                   length(sample.ids), ")"))
      tryCatch({
        idx <- grep(sample, basename(sample.paths))
        EaCoN:::SNP6.Process(CEL = file.path(data.dir[idx], paste0(sample, cel.suffix[idx])), 
                             samplename = sample)
      }, error=function(e){NULL})
    }, mc.cores = cl)
  } else {
    for(sample in sample.ids){
      print(paste0(sample, ": (", match(sample, sample.ids), "/", 
                   length(sample.ids), ")"))
      tryCatch({
        setwd(file.path(pdir, "eacon"))
        idx <- grep(sample, basename(sample.paths))
        EaCoN:::SNP6.Process(CEL = file.path(data.dir[idx], paste0(sample, cel.suffix[idx])), 
                             samplename = sample)
      }, error=function(e){NULL})
    }
  }
  
}


##############################
#### 2) L2R Segmentation: ####
# Takes in the _processed.RDS file from the apt/rcnorm pre-processing
# Uses the copynumber::winsorize to winsorize L2R at smooth.k MADs
# Alelle-specific piecewise fit segmentation using ASCAT::ascat.aspcf
# Recenters based on the 'recenter' parameter:
#   - l2r.centeredpeak (Default): Most centred peak in L2R density
#   - l2r.mainpeak: Heighest peak of the L2R density
#   - l2r.median: Median of L2R distribution
# Winsorization of segmented L2R using a MAD k=5
# Rescues small CN events using the changepoint::cpt.mean PELT method with
# a penalization factor of 40
#
for(segmenter in c("ASCAT")){
  message("Running L2R segmentation...")
  
  ## Select non-processed files

  RDS.files <- .removeRedundantFiles(pattern1="_processed.RDS$", 
                                     pattern2=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"),
                                     unlink.path = file.path(toupper(segmenter), "L2R"),
                                     astep="L2R", segmenter=toupper(segmenter))
  RDS.files <- .splitSamples(RDS.files, opt$idx, opt$grpsize)
  
  EaCoN::Segment.ff.Batch(RDS.file = RDS.files,  segmenter = segmenter, nthread=3)
}

########################
#### 3) ASCN Calls: ####
# Takes in the *.SEG.ASCAT.RDS file from the ASCAT aspcf analysis
# Scans a gamma range of (gammaRange[1], gammaRange[2], by=0.05) (0.35, 0.95, 0.05)
#  Gamma parameter: (https://www.crick.ac.uk/research/labs/peter-van-loo/software)
#    The parameter represents the drop in LogR for a change from two copies to 
#    one copy in 100% of cells. Gamma theoretically should equal 1, but due to 
#    array background signal and bespoke array normalisation procedures, in 
#    practice it is often significantly lower. Its default setting of 0.55 works 
#    for many but not all SNP arrays
# Runs ASCAT, ASCAT::ascat.runAscat, using the gamma parameter
# Basic visualization and reports after this function
# 
for(segmenter in c("ASCAT")){
  message("Obtaining allele-specific copy number (ASCN) calls...")
  ## CN Estimation:
  # Provides ASCN calls from ASCAT
  l2r.rds <- .removeRedundantFiles(pattern1=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"), 
                                   pattern2="gammaEval.txt$",
                                   unlink.path = file.path(toupper(segmenter), "ASCN"), 
                                   astep="ASCN", segmenter=toupper(segmenter)) #, unlink.stat=FALSE
  l2r.rds <- .splitSamples(l2r.rds, opt$idx, opt$grpsize)
  print(l2r.rds)
  
  fit.val <- EaCoN::ASCN.ff.Batch(RDS.files = l2r.rds, nthread=4, force=TRUE)
}


#### Output builder: ####
for(segmenter in c("ASCAT")){
  message("Selecting the best gamma fit for ASCAT, annotating, and generating output files...")
  if(toupper(segmenter)=='ASCAT'){
    
    fit.vals.path <- unlist(sapply(list.files(), function(x){
      if(file.exists(file.path(x, segmenter, "ASCN", 
                               paste0(x, ".gammaEval.txt")))){
        file.path(".", x, segmenter, "ASCN", paste0(x, ".gammaEval.txt"))
      }
    }))

    if(any(grepl("bkup", fit.vals.path))) fit.vals.path <- fit.vals.path[-grep("bkup", fit.vals.path)]
    all.fits <- lapply(fit.vals.path, function(fvp){
      list(fit=read.table(fvp, sep="\t", header=T, stringsAsFactors = F, 
                          check.names = F, fill=F),
           sample=strsplit(fvp, "/")[[1]][2])
    })
    names(all.fits) <- sapply(all.fits, function(x) x$sample)
    
    meta.l <- switch(dataset,
                     GNE={
                       data(gne.meta)
                       meta <- gne.meta[,c('INVENTORY_SAMPLE_NAME', 'TCGA_code')]
                       colnames(meta) <- c('sample', 'TCGA_code')
                       list('meta.tcga'=meta, 'meta'=gne.meta)
                     },
                     CCLE={
                       data(CCLE_meta)
                       meta <- ccle.meta[,c('SNP arrays', 'tcga_code')]
                       colnames(meta) <- c('sample', 'TCGA_code')
                       list('meta.tcga'=meta, 'meta'=ccle.meta)
                     },
                     GDSC={
                       data(GDSC_meta)
                       meta <- gdsc.meta[,c('Sample Name', 'Cancer Type (matching TCGA label)')]
                       colnames(meta) <- c('sample', 'TCGA_code')
                       list('meta.tcga'=meta, 'meta'=gdsc.meta)
                     },
                     '1000G'={
                       meta <- data.frame("sample"=names(all.fits),
                                          "TCGA_code"=rep("Normal", length(all.fits)))
                       list('meta.tcga'=meta, 'meta'=meta)
                     })
    
    data(pancanPloidy.noSegs)
    pancan.ploidy <- pancan.obj.segless$breaks$ploidy
    # plot(pancan.ploidy[,c('breaks', 'BRCA')], type='l')
    
    max.process <- opt$grpsize # 15
    split.range <- seq(1, length(all.fits), by=max.process)
    split.range <- data.frame("start"=split.range,
                              "end"=c(split.range[-1]-1, length(all.fits)))
    r <- apply(split.range, 1, function(r){
      print(paste(r, collapse="-"))
      start_time <- Sys.time()
      
      gr.cnv <- annotateRDS.Batch(all.fits[r['start']:r['end']], 
                                  toupper(segmenter), nthread=4,
                                  gamma.method='score', gamma.meta=meta.l$meta.tcga,
                                  pancan.ploidy=pancan.ploidy, 
                                  feature.set=c('bins'),
                                  bin.size=5000)
      
      # cbio.path=file.path("out", "cBio")
      # buildCbioOut(gr.cnv, cbio.path="./out/cBio")
      
      ## Build standard bin and gene PSets
      pset.path=file.path("out", "PSet")
      col.ids <- colnames(gr.cnv[[1]]$bins$genes)[c(-1)]
      col.ids <- c("seg.mean", col.ids[-c(grep("seg.mean", col.ids), length(col.ids))])
      col.ids <- c("Log2Ratio", "L2Rraw")
      
      buildPSetOut(gr.cnv, anno.name=dataset, pset.path, 
                   seg.id='Log2Ratio',
                   meta=meta.l$meta, 
                   out.idx=c(r['start'], r['end']),
                   cols=col.ids)
      
      end_time <- Sys.time()
      t1 <- end_time - start_time
      gc()
    })
    
  }
}


## Outdated? Maybe?
qsub.split <- FALSE
if(qsub.split){
  if(file.exists(file.path("..", "scripts", "samples.RData"))){
    EaCoN:::tmsg("Getting sample IDs...")
    load(file.path("..", "scripts", "samples.RData"))
    range.s <- seq(1, length(sample.ids), by=10)
    range.e <- c((range.s + 9), length(sample.ids))
    range.idx <- opt$idx
    
    range.se <- c(range.s[range.idx], range.e[range.idx])
    print(paste0(range.idx, " : ", paste(range.se, collapse="-")))
    
    sample.ids <- sample.ids[range.se[1]:range.se[2]]
    #sample.ids <- sample.ids[opt$idx]
  } else {
    # fit.vals.path <- list.files(pattern="gammaEval.txt", recursive=T, full=T)
    # rm.idx <- unlist(sapply(sample.ids, function(x) any(grepl(x, fit.vals.path))))
    
    rm.idx <- sapply(sample.ids, function(s){
      if(dir.exists(s)){
        if(any(grepl("temp", list.files(s)))){
          print(paste("Incomplete sample [temp present]:", s))
          unlink(file.path(s), recursive = TRUE)
          temp.idx <- FALSE
        } else if(any(grepl("RDS$", list.files(s)))) {
          temp.idx <- TRUE
        } else {
          print(paste("Incomplete sample [RDS missing]:", s))
          unlink(file.path(s), recursive = TRUE)
          FALSE
        }
      } else {
        temp.idx <- FALSE
      }
    })
    
    #names(rm.idx) <- sample.ids
    if(sum(rm.idx) > 0) sample.ids <- sample.ids[-which(rm.idx)]
    #sample.ids <- sample.ids[sample(x = c(1:length(sample.ids)), size = length(sample.ids), replace = F)]
    
    EaCoN:::tmsg("Saving samples...")
    save(sample.ids, file=file.path("..", "scripts", "samples.RData"))
  }
}
