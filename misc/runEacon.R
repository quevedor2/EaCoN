#### Libraries ####
## Instal ASCAT and FACETS
install.pkg <- FALSE
if(install.pkg){
  devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
  devtools::install_github("mskcc/facets")
  
  ## Instal Bioconductor packages
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("affxparser", "Biostrings", "aroma.light", "BSgenome", "copynumber", "GenomicRanges", "limma", "rhdf5", "sequenza"))
  
  ## Install the most recent STABLE version (@master)
  #devtools::install_github("gustaveroussy/EaCoN")
  
  install.packages("https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz", repos = NULL, type = "source")
  devtools::install_github("gustaveroussy/apt.snp6.1.20.0")
  install.packages("https://partage.gustaveroussy.fr/pydio_public/152397?dl=true&file=/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz", repos = NULL, type = "source")
  install.packages( "https://partage.gustaveroussy.fr/pydio_public/e6fe22?dl=true&file=/rcnorm_0.1.5.tar.gz", repos = NULL, type = "source")
  
}

##############
#### Main ####

#devtools::install_github("quevedor2/EaCoN", ref = 'tads')
library(optparse)
library(EaCoN)
library(parallel)
library(Biobase)

option_list <- list(
  make_option(c("-i", "--idx"), type="integer", default=NULL,
              help="Index of sample group [default= %default]", metavar="integer"),
  make_option(c("-g", "--grpsize"), type="integer", default=10,
              help="Size of the groups to run in one job [default= %default]", metavar="integer"),
  make_option(c("-d", "--dataset"), type="character", default='GDSC',
              help="Dataset to use, either 'GDSC' or 'GNE' or 'CCLE'"),
  make_option(c("-p", "--pdir"), type="character", default='/mnt/work1/users/pughlab/projects/cancer_cell_lines',
              help="Parent directory path"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


.splitSamples <- function(sample.ids, range.idx, grp.size=10){
  if((grp.size * range.idx) > length(sample.ids)){
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



segmenter <- 'ASCAT'
dataset <- opt$dataset  #'GDSC', 'CCLE', 'gCSI'
pdir <- file.path(opt$pdir, dataset)

## Normalization
# Outputs a ./YT_4941/YT_4941_GenomeWideSNP_6_hg19_processed.RDS file
message(paste0("Searching for CEL$ files in: ", file.path(pdir, "data")))
CEL.dir <- file.path(pdir, "data")
sample.paths <- list.files(CEL.dir, pattern="CEL$", recursive = T, 
                         ignore.case = T, full.names = T)
sample.ids <- gsub(".cel", "", basename(sample.paths), ignore.case = TRUE)
regm <- regexpr(".cel", basename(sample.paths), ignore.case=T)
cel.suffix <- regmatches(x = basename(sample.paths), m = regm)


dir.create(file.path(pdir, "eacon"), recursive = TRUE)
setwd(file.path(pdir, "eacon"))



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

if(dataset=='GNE'){
  illumina.dir <- file.path(pdir, 'illumina', 'GNE_Matrices')
  EaCoN::Build.OMNI25(illumina.dir=illumina.dir, parent.dir=pdir)
} else {
  mclapply(sample.ids, function(sample, sample.paths){
    idx <- grep(sample, basename(sample.paths))
    SNP6.Process(CEL = file.path(CEL.dir, paste0(sample, cel.suffix[idx])), 
                 samplename = sample)
  }, sample.paths=sample.paths, mc.cores = 2)
}



#### L2R Segmentation: ####
# Takes in the _processed.RDS file 
for(segmenter in c("ASCAT")){
  message("Running L2R segmentation...")
  
  ## Select non-processed files

  RDS.files <- .removeRedundantFiles(pattern1="_processed.RDS$", 
                                     pattern2=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"),
                                     unlink.path = file.path(toupper(segmenter), "L2R"),
                                     astep="L2R", segmenter=toupper(segmenter))
  RDS.files <- .splitSamples(RDS.files, opt$idx, opt$grpsize)
  
  EaCoN:::Segment.ff.Batch(RDS.file = RDS.files,  segmenter = segmenter, nthread=2)
}

#### ASCN Calls: ####
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
  
  fit.val <- EaCoN:::ASCN.ff.Batch(RDS.files = l2r.rds, nthread=2, force=TRUE)
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
                                  toupper(segmenter), nthread=3,
                                  gamma.method='score', gamma.meta=meta.l$meta.tcga,
                                  pancan.ploidy=pancan.ploidy, 
                                  feature.set=c('bins'),
                                  bin.size=5000)
      
      cbio.path=file.path("out", "cBio")
      # buildCbioOut(gr.cnv, cbio.path="./out/cBio")
      
      ## Build standard bin and gene PSets
      pset.path=file.path("out", "PSet")
      buildPSetOut(gr.cnv, dataset, pset.path, meta=meta.l$meta, out.idx=c(r['start'], r['end']))
      
      end_time <- Sys.time()
      t1 <- end_time - start_time
      print(paste0("Time to create jaccard matrix: ", round(t1,2), "s"))
      
      gc()
      r
    })
    
  }
}
