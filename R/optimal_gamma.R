#### setup ####
library(EaCoN)
pdir <- "/mnt/work1/users/pughlab/projects/CCLE/eacon"
setwd(pdir)
sample <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_A05_256010'
segmenter <- 'ASCAT'
gamma <- 0.7
#my.data <- EaCoN:::loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
my.data <- loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
pancan.dir <- '/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/raw'
pancan.obj <- createPancanSegRef(pancan.dir, out.dir='/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/cleaned')

#### Functions ####
createPancanSegRef <- function(ref.dir, out.dir=NULL, verbose=T){
  setwd(ref.dir)
  seg.f <- list.files(pattern="whitelisted.seg$")
  anno.f <- list.files(pattern="annotations.tsv$")
  ploidy.f <- list.files(pattern="JSedit.fixed.txt$")
  pancan.files <- list("seg"=seg.f, "anno"=anno.f, "ploidy"=ploidy.f)
  
  pancan.data <- lapply(pancan.files, read.table, sep="\t", header=T, 
                        stringsAsFactors = F, check.names = F, fill=T,
                        quote='', comment.char='')
  
  ## Remove problematic samples
  fail.idx <- grep("failed", pancan.data[['anno']]$`cancer type`)
  blank.idx <- grep("^$", pancan.data[['anno']]$`cancer type`, perl=TRUE)
  dnu.idx <- grep("True", pancan.data[['anno']]$`Do_not_use`, perl=TRUE)
  dnu.blank.idx <- grep("^$", pancan.data[['anno']]$`Do_not_use`, perl=TRUE)
  pancan.data[['anno']] <- pancan.data[['anno']][-c(fail.idx, blank.idx,
                                                    dnu.idx, dnu.blank.idx),]
  
  ## Remove samples with no ploidy call
  blank.idx <- grep("^$", pancan.data[['ploidy']]$`call status`)
  pancan.data[['ploidy']] <- pancan.data[['ploidy']][-blank.idx,]
  
  ## Merge Ploidy with Annotations
  anno.ploidy <- Reduce(function(x,y) merge(x,y, by.x='sample', by.y='aliquot_barcode', all=F), 
                        pancan.data[c('ploidy', 'anno')])
  
  ## Identify the .Seg data for Annotated-ploidy samples 
  seg <- split(pancan.data[['seg']], f=pancan.data[['seg']]$Sample)
  seg <- seg[match(anno.ploidy$sample, names(seg))]
  na.idx <- which(is.na(names(seg)))
  seg <- seg[-na.idx]
  anno.ploidy.seg <- anno.ploidy[-na.idx,]
  
  if(verbose){
    ap.cnt <- nrow(anno.ploidy)
    aps.cnt <- nrow(anno.ploidy.seg)
    a.cnt <- nrow(pancan.data[['anno']])
    p.cnt <- nrow(pancan.data[['ploidy']])
    s.cnt <- length(unique(pancan.data[['seg']]$Sample))
    
    cat('> Overlap between Annotations and Ploidy:\n')
    cat(paste0("\t[", ap.cnt, " / ", a.cnt, "] samples have an associated ploidy/purity value \n"))
    cat(paste0("\t[", ap.cnt, " / ", p.cnt, "] samples with purity/ploidy values have an associated annotation \n"))
    
    cat('> Overlap between Annotations, Ploidy and Seg-data:\n')
    cat(paste0("\t[", aps.cnt, " / ", ap.cnt, "] annotated-ploidy samples also have seg-data \n"))
    cat(paste0("\t[", aps.cnt, " / ", s.cnt, "] samples with seg-data have an associated annotation and ploidy value \n"))
  }
  
  pancan.obj <- list('AP-meta'=anno.ploidy,
                     'APS-meta'=anno.ploidy.seg,
                     'seg'=seg)
  if(!is.null(out.dir)){
    if(!dir.exists(out.dir)){
      dir.create(path = out.dir, recursive = T)
      warning(paste("Outputing to directory: ", out.dir))
    }
    saveRDS(pancan.obj, file = file.path(out.dir, "pancanPloidy.RDS"))
  }
  return(pancan.obj)
}

scoreSegsABC <- function(seg1.gr, seg2.gr){
  require(GenomicRanges)
  .expandRle <- function(rle.obj){
    as.character(rep(rle.obj@values, rle.obj@lengths))
  }
  .ovSeg <- function(gr0, gr1, col, id){
    ov.idx <- findOverlaps(gr1, gr0)
    gr2 <- gr0[subjectHits(ov.idx),]
    mcols(gr2)[,id] <- mcols(gr1[queryHits(ov.idx), ])[,col]
    gr2
  }
  .getABC <- function(gr, col.ids){
    abc <- (mcols(gr)[,col.ids[1]] - mcols(gr)[,col.ids[2]]) * width(gr)
    sum(abs(abc))
  }
  
  comb.gr <- sort(GRanges(seqnames=c(.expandRle(seqnames(seg1.gr)), 
                                     .expandRle(seqnames(seg2.gr))),
                          IRanges(start=c(start(seg1.gr), start(seg2.gr)), 
                                  end=c(end(seg1.gr), end(seg2.gr)))))
  suppressWarnings(seqlevelsStyle(comb.gr) <- 'UCSC')
  
  comb.gr <- .ovSeg(comb.gr, seg1.gr, col='seg.mean', id='seg')
  comb.gr <- .ovSeg(comb.gr, seg2.gr, col='Segment_Mean', id='pan.seg')
  
  abc <- .getABC(comb.gr, c('seg', 'pan.seg'))
  abc
}

#### Load ####
seg <- my.data$segments_raw
seg$width <- round(with(seg, endpos - startpos)/100000,0) + 1
seg$nABraw <- with(seg, nAraw + nBraw)
seg$seg.mean <- round(log2(seg$nABraw / 2),3)
seg$seg.mean <- seg$seg.mean - median(rep(seg$seg.mean, seg$width))
seg$seg.mean[seg$seg.mean < -2] <- -2

## Make GRanges objects:
#seg1
seg1.gr <- makeGRangesFromDataFrame(seg, keep.extra.columns = T,
                                    start.field = c('startpos'),end.field = c('endpos'))
seqlevelsStyle(seg1.gr) <- 'UCSC'

#seg2
pancan.seg <- as.data.frame(do.call(rbind, pancan.obj[['seg']]))
pancan.seg$Chromosome <- gsub("23$", "X", pancan.seg$Chromosome) %>%
  gsub("24$", "X", .)
pancan.grl <- makeGRangesListFromDataFrame(pancan.seg, split.field='Sample', keep.extra.columns=T)
seqlevelsStyle(pancan.grl) <- 'UCSC'

## Calculate ABC
nthread <- 10
pancan.abc <- mclapply(pancan.grl, function(seg2) scoreSegsABC(seg1.gr=seg1.gr, seg2.gr=seg2), mc.cores = nthread)




## Visualize periodicity
par(mfrow=c(2,1))
with(seg, hist(rep(seg.mean, width), breaks = 50))
with(seg, hist(rep(nABraw, width), breaks = 50))

seg.spl <- split(seg, f=seg$chr)
seg.spl <- seg.spl[paste0("chr", c(1:22, "X",))]
par(mfrow=c(length(seg.spl),1), mar=c(0.5, 5.1, 0.5, 4.1))
x <- sapply(seg.spl, function(seg.i, ...){
  require(pracma)
  require(scales)
  h.data <- with(seg.i, hist(rep(seg.mean, width), breaks=seq(-2, 3, by=0.02),
                             xlim=c(-2, 3), main='', xlab='', col="grey"))
  
  pk <- findpeaks(h.data$counts, nups = 0, ndowns = 0, zero = "0", ...)
  if(!is.null(pk)){
    br <- h.data$breaks[-1]
    apply(pk, 1, function(p){rect(xleft = br[p[3]], ybottom = 0, xright = br[p[4]], ytop = p[1],
                                  border=NA, col=alpha("blue", 0.5))})
    abline(v = br[pk[,2]], col="red")
  }
}, minpeakheight = 50, minpeakdistance = 10)
