#### setup ####
library(EaCoN)
library(dplyr)
library(GenomicRanges)
pdir <- "/mnt/work1/users/pughlab/projects/CCLE/eacon"
setwd(pdir)
sample <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_A05_256010'
segmenter <- 'ASCAT'
gamma <- 0.7
#my.data <- EaCoN:::loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
my.data <- loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
pancan.dir <- '/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/raw'
pancan.obj <- createPancanSegRef(pancan.dir, out.dir='/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/cleaned',
                                 write.seg=T)

#### Functions ####
loadBestFitRDS <- function(gamma, ...){
  RDS.file <- file.path(sample, toupper(segmenter), 'ASCN', paste0('gamma', format(gamma, nsmall=2)),
                        paste(sample, 'ASCN', toupper(segmenter), 'RDS', sep="."))
  rds <- readRDS(RDS.file)
  return(rds)
}

createPancanSegRef <- function(ref.dir, out.dir=NULL, verbose=T, write.seg=FALSE){
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
  
  ## Output stuff
  if(!is.null(out.dir)){
    if(!dir.exists(out.dir)){
      dir.create(path = out.dir, recursive = T)
      warning(paste("Outputing to directory: ", out.dir))
    }
    saveRDS(pancan.obj, file = file.path(out.dir, "pancanPloidy.RDS"))
  }
  
  ## Writes cancer-specific seg files
  if(write.seg & !is.null(out.dir)){
    aps.cancer.type <- split(anno.ploidy.seg, f=anno.ploidy.seg$`cancer type`)
    
    sapply(names(aps.cancer.type), function(cancer.type){
      cancer.seg <- seg[aps.cancer.type[[cancer.type]]$sample]
      cancer.seg <- do.call(rbind, cancer.seg)
      write.table(cancer.seg, file=file.path(out.dir, paste0(cancer.type, "_pancan.seg")),
                  col.names=T, row.names=F, quote = F, sep = '\t')
    })
  }
  
  return(pancan.obj)
}

makeRefGRanges <- function(pancan.gr, gr1){
  total.gr <- unlist(GRangesList(pancan.gr, granges(gr1)))
  pancan.chr.grl <- split(total.gr, f=seqnames(total.gr))
  comb.grl <- lapply(pancan.chr.grl, function(p.gr){
    pos <- unique(sort(c(start(p.gr), end(p.gr))))
    GRanges(seqnames=rep(unique(as.character(seqnames(p.gr))), length(pos)-1),
            ranges = IRanges(start=pos[-length(pos)], 
                             end=pos[-1]))
  })
  comb.gr <- unlist(as(comb.grl, 'GRangesList'))
  return(comb.gr)
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

combineGr <- function(gr){
  gr.chr <- split(gr, seqnames(gr))
  gr.c <- as(lapply(gr.chr, function(gr0){
    pos <- unique(sort(c(start(gr0), end(gr0))))
    GRanges(seqnames = rep(unique(as.character(seqnames(gr0))), length(pos)-1),
            IRanges(start=pos[-length(pos)],
                    end=pos[-1]))
  }), 'GRangesList')
  gr.c <- unlist(gr.c)
  gr.c
}

populateGr <- function(grl, ref.gr, col.id){
  for(each.id in names(grl)){
    ov.idx <- findOverlaps(grl[[each.id]], ref.gr)
    mcols(ref.gr)[,each.id] <- NA
    mcols(ref.gr)[subjectHits(ov.idx),each.id] <- mcols(grl[[each.id]])[queryHits(ov.idx),col.id]
  }
  ref.gr
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

## Make the combined reference GRanges obj
pancan.gr <- sort(makeGRangesFromDataFrame(pancan.seg))
seqlevelsStyle(pancan.gr) <- 'UCSC'
ref.gr <- makeRefGRanges(pancan.gr, seg1.gr)

## Create seg.mean matrix

pancan.abc <- mclapply(pancan.grl, function(seg, col) {
  ov.idx <- findOverlaps(seg, ref.gr)
  mcols(ref.gr)[,id] <- mcols(seg[queryHits(ov.idx),])[,'Segment_Mean']
})

## Calculate ABC
nthread <- 10
pancan.abc <- mclapply(pancan.grl, function(seg2) scoreSegsABC(seg1.gr=seg1.gr, seg2.gr=seg2), mc.cores = nthread)
pancan.abc <- as.matrix(unlist(pancan.abc))

## Group ABC based on tumor type
aps.spl <- split(pancan.obj[['APS-meta']], f=pancan.obj[['APS-meta']]$`cancer type`)
abc.spl <- lapply(aps.spl, function(i) pancan.abc[i$sample,])
boxplot(abc.spl[order(sapply(abc.spl, mean))], horizontal=T, las=1)


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
