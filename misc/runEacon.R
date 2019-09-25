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

#devtools::install_github("quevedor2/EaCoN", ref = 'preprocess-builder')
library(optparse)
library(EaCoN)
library(parallel)

# option_list <- list(make_option(c("-i", "--idx"), type="integer", default=NULL,
#               help="Index of sample group [default= %default]", metavar="integer"))
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)

pdir <- '/mnt/work1/users/pughlab/projects/CCLE'

## Normalization
# Outputs a ./YT_4941/YT_4941_GenomeWideSNP_6_hg19_processed.RDS file
CEL.dir <- file.path(pdir, "data", "symlinks")
sample.ids <- list.files(CEL.dir, pattern="CEL$")
sample.ids <- gsub(".cel", "", sample.ids, ignore.case = TRUE)


dir.create(file.path(pdir, "eacon"), recursive = TRUE)
setwd(file.path(pdir, "eacon"))

qsub.split <- TRUE
if(qsub.split){
  if(file.exists(file.path("..", "scripts", "samples.RData"))){
    EaCoN:::tmsg("Getting sample IDs...")
    load(file.path("..", "scripts", "samples.RData"))
    # range.s <- seq(1, length(sample.ids), by=10)
    # range.e <- c((range.s + 9), length(sample.ids))
    # range.idx <- opt$idx
    # 
    # range.se <- c(range.s[range.idx], range.e[range.idx])
    # print(paste0(range.idx, " : ", paste(range.se, collapse="-")))
    # 
    # sample.ids <- sample.ids[range.se[1]:range.se[2]]
    sample.ids <- sample.ids[opt$idx]
  } else {
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
    
    names(rm.idx) <- sample.ids
    sample.ids <- sample.ids[-which(rm.idx)]
    sample.ids <- sample.ids[sample(x = c(1:length(sample.ids)), size = length(sample.ids), replace = F)]
    
    EaCoN:::tmsg("Saving samples...")
    save(sample.ids, file=file.path("..", "scripts", "samples.RData"))
  }
}

mclapply(sample.ids, function(sample){
  SNP6.Process(CEL = file.path(CEL.dir, paste0(sample, ".CEL")), 
               samplename = sample)
}, mc.cores = 10)


## Segmentation:
# Takes in the _processed.RDS file 
for(segmenter in c("ASCAT")){
  print(segmenter)
  RDS.files <- list.files(pattern="_processed", recursive = TRUE, full.names = TRUE)
  Segment.ff.Batch(RDS.file = RDS.files,  segmenter = segmenter, nthread=5)
  
  ## CN Estimation:
  # Provides ASCN calls from ASCAT
  l2r.rds <- list.files(pattern=paste0("\\.SEG\\.", toupper(segmenter), ".*\\.RDS$"), recursive=T, full=T)
  if(any(grepl('bkup', l2r.rds))) l2r.rds <- l2r.rds[-grep("bkup", l2r.rds)]
  fit.val <- ASCN.ff.Batch(RDS.files = l2r.rds, nthread=5)
}

for(segmenter in c("ASCAT")){
  print(segmenter)
  if(toupper(segmenter)=='ASCAT'){
    fit.vals.path <- list.files(pattern="gammaEval.txt", recursive=T, full=T)
    fit.vals.path <- fit.vals.path[-grep("bkup", fit.vals.path)]
    all.fits <- lapply(fit.vals.path, function(fvp){
      list(fit=read.table(fvp, sep="\t", header=T, stringsAsFactors = F, 
                          check.names = F, fill=F),
           sample=strsplit(fvp, "/")[[1]][2])
    })
    names(all.fits) <-sapply(all.fits, function(x) x)
    
    gr.cnv <- annotateRDS.Batch(all.fits, toupper(segmenter), nthread=3)
    
    cbio.path=file.path("out", "cBio")
    buildCbioOut(gr.cnv, cbio.path="./out/cBio", overwrite=sample)
    
    pset.path=file.path("out", "PSet")
    buildPSetOut(gr.cnv, "CGP", pset.path, meta=cell.line.anno)
  }
}
FIEFS_p_NCLE_DNA_Affy5_GenomeWideSNP_6_C06_411046


x <- list("A"=1:10,
          "B"=5:15,
          "C"=10:20)
boxplot(x)
sapply(seq_along(x), function(x.pos){
  require(scales)
  spread <- 0.1
  y.vals <- x[[x.pos]]
  x.vals <- x.pos + runif(n=length(y.vals), min = -1*spread, max = spread)
  points(x=x.vals, y=y.vals, pch=16, col=alpha("black", 0.7))
})
