#### Libraries ####
## Instal ASCAT and FACETS
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

###################
#### Variables ####
library(EaCoN)
library(foreach)
sample <- "YT_4941"
segmenter <- 'ASCAT'
pdir <- '/mnt/work1/users/home2/quever/snp6tmp'

## ASCN.ASCAT:

tmsg <- function (text = NULL) {
  message(paste0(" [", Sys.info()[["nodename"]], ":", Sys.getpid(), 
                 "] ", text))
}
write.table.fast <- function (x, file = NULL, header = TRUE, sep = "\t", fileEncoding = "", 
                              row.names = FALSE, ...) {
  if (header) 
    write.table(x = x[NULL, ], file = file, sep = "\t", quote = FALSE, 
                row.names = FALSE, fileEncoding = fileEncoding)
  if (!row.names) 
    rownames(x) <- NULL
  trychk <- try(iotools::write.csv.raw(x = x, file = file, 
                                       sep = sep, col.names = FALSE, fileEncoding = fileEncoding, 
                                       append = header, ...))
  if (!is.null(trychk)) {
    print("Fast write failed, using canonical write.table ...")
    write.table(x = x, file = file, sep = sep, row.names = row.names, 
                quote = FALSE)
  }
  gc()
}
EaCoN.Rorschard.plot <- function (data = NULL, cnpTotal = NULL) {
  k.sqrt <- ceiling(sqrt(length(data$data$chrs)))
  par(mar = c(1, 1, 1, 1), mfrow = c(k.sqrt, k.sqrt))
  for (k in 1:length(data$data$ch)) {
    graphics::plot(data$data$Tumor_BAF[[1]][data$germline$germlinegenotypes], 
                   data$data$Tumor_LogR[[1]][data$germline$germlinegenotypes], 
                   pch = ".", cex = 2, xlim = c(0, 1), ylim = c(-2, 
                                                                2), col = "grey95")
    points(data$data$Tumor_BAF[[1]][!data$germline$germlinegenotypes], 
           data$data$Tumor_LogR[[1]][!data$germline$germlinegenotypes], 
           pch = ".", cex = 2, col = "grey50")
    kin <- data$data$ch[[k]][!(data$data$ch[[k]] %in% which(data$germline$germlinegenotypes))]
    r.col <- if (is.null(cnpTotal)) 
      3
    else cnpTotal[kin] + 1
    points(data$data$Tumor_BAF[[1]][kin], data$data$Tumor_LogR[[1]][kin], 
           pch = ".", cex = 4, col = r.col)
    try(text(x = 0.5, y = 2, labels = data$data$chrs[k], 
             pos = 1, cex = 2))
  }
}


RDS.file <- list.files(file.path(pdir, sample, segmenter, "L2R"), pattern=".RDS$")
data <- readRDS(file.path(pdir, sample, segmenter, 'L2R', RDS.file))
out.dir <- file.path(pdir, sample, segmenter, 'ASCN')

gammaRange = c(.35,.7)
nsubthread = 1
cluster.type = "PSOCK"
force = FALSE


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
  make_option(c("-d", "--dataset"), type="character", default=NULL,
              help="Dataset to use, either 'GDSC' or CCLE'"))
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

.removeRedundantFiles <- function(pattern1, pattern2, unlink.path=NULL){
  p1.files <- list.files(pattern=pattern1, recursive = TRUE, full.names = TRUE)
  p2.files <- list.files(pattern=pattern2, recursive=T, full=T)
  rm.ids <- sapply(strsplit(p2.files, "/"), function(x) x[[2]])
  rm.idx <- unlist(sapply(rm.ids, function(x) grep(paste0("\\/", x, "\\/"), p1.files)))
  if(length(rm.idx) > 0){
    p.files <- p1.files[-rm.idx]
  } else {
    p.files <- p1.files
  }
  
  if(!is.null(unlink.path)){
    print(paste0("Unlinking files in path: <Sample>/", unlink.path))
    files.to.unlink <- sapply(strsplit(p.files, "/"), function(x) x[[2]])
    sapply(files.to.unlink, function(ftu){
      unlink(x = file.path(ftu, unlink.path), recursive = T)
    })
  }
  p.files
}



segmenter <- 'ASCAT'
dataset <- opt$dataset  #'GDSC'
pdir <- file.path('/mnt/work1/users/pughlab/projects/cancer_cell_lines', dataset)

## Normalization
# Outputs a ./YT_4941/YT_4941_GenomeWideSNP_6_hg19_processed.RDS file
CEL.dir <- file.path(pdir, "data")
sample.paths <- list.files(CEL.dir, pattern="CEL$", recursive = T, 
                         ignore.case = T, full.names = T)
sample.ids <- gsub(".cel", "", basename(sample.paths), ignore.case = TRUE)
regm <- regexpr(".cel", basename(sample.paths), ignore.case=T)
cel.suffix <- regmatches(x = basename(sample.paths), m = regm)




dir.create(file.path(pdir, "eacon"), recursive = TRUE)
setwd(file.path(pdir, "eacon"))



qsub.split <- TRUE
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

mclapply(sample.ids, function(sample, sample.paths){
  idx <- grep(sample, basename(sample.paths))
  SNP6.Process(CEL = file.path(CEL.dir, paste0(sample, cel.suffix[idx])), 
               samplename = sample)
}, sample.paths=sample.paths, mc.cores = 2)


#### L2R Segmentation: ####
# Takes in the _processed.RDS file 
for(segmenter in c("ASCAT")){
  print(segmenter)
  
  ## Select non-processed files
  RDS.files <- .removeRedundantFiles(pattern1="_processed.RDS$", 
                                     pattern2=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"),
                                     unlink.path = file.path(toupper(segmenter), "L2R"))
  RDS.files <- .splitSamples(RDS.files, opt$idx, opt$grpsize)
  
  Segment.ff.Batch(RDS.file = RDS.files,  segmenter = segmenter, nthread=5)
}

#### ASCN Calls: ####
for(segmenter in c("ASCAT")){
  print(segmenter)
  ## CN Estimation:
  # Provides ASCN calls from ASCAT
  l2r.rds <- .removeRedundantFiles(pattern1=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"), 
                                   pattern2="gammaEval.txt$",
                                   unlink.path = file.path(toupper(segmenter), "ASCN"))
  l2r.rds <- .splitSamples(l2r.rds, opt$idx, opt$grpsize)
  print(l2r.rds)
  
  fit.val <- ASCN.ff.Batch(RDS.files = l2r.rds, nthread=2)
}


#### Output builder: ####
for(segmenter in c("ASCAT")){
  print(segmenter)
  if(toupper(segmenter)=='ASCAT'){
    
    fit.vals.path <- unlist(sapply(list.files(), function(x){
      if(file.exists(file.path(x, segmenter, "ASCN", 
                               paste0(x, ".gammaEval.txt")))){
        file.path(".", x, segmenter, "ASCN", paste0(x, ".gammaEval.txt"))
      }
    }))
    # fit.vals.path <- list.files(pattern="gammaEval.txt", recursive=T, full=T)
    
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
                                  pancan.ploidy=pancan.ploidy, feature.set=c('bins', 'tads'))
      
      cbio.path=file.path("out", "cBio")
      buildCbioOut(gr.cnv, cbio.path="./out/cBio")
      
      ## Build standard bin and gene PSets
      pset.path=file.path("out", "PSet")
      buildPSetOut(gr.cnv, dataset, pset.path, meta=meta.l$meta, out.idx=c(r['start'], r['end']))
      
      ## Build TAD and CRE PSets
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
