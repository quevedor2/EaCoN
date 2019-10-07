library(DNAcopy)
library(GenomicRanges)
library(scales)
library(plyr)
library(optparse)

option_list <- list(make_option(c("-c", "--chr"), type="character", default=NULL,
                                help="Chromosome (e.g. chr1)"),
                    make_option(c("-o", "--outdir"), type="character", default="./",
                                help="Out directory"),
                    make_option(c("-p", "--plot"), action="store_true", default=TRUE,
                                help="Generate plots [default]"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
                    
PDIR <- "/mnt/work1/users/pughlab/references/ENCODE/cancer_cell_lines/HiC/raw/TAD"
setwd(PDIR)

#### Functions ####
.addSD <- function(x, y, sd){
  segments(x0 = c(x, x), y0 = c(y,y), 
           x1 = c(x, x), y1 = c(y+sd, y-sd))
}

.addSeg <- function(x, col){
  for(i in seq(1:dim(x$output)[1])) {
    #guideline at CNV median
    segments(x$output[i,'loc.start'],
             as.numeric(x$output[i,'seg.mean']),
             x$output[i,'loc.end'],
             as.numeric(x$output[i,'seg.mean']),
             col=col,
             lwd=2)       
  }          
}

plotGr0 <- function(gr0, max.y,  ...){
  plot(gr0$idx, sort(start(gr0)), col=alpha("firebrick1", 0.4), 
       ylim=c(0, max.y), pch=16, ...)
  st.cbs <- getCbsSeg(chr='chr21', val=sort(start(gr0)),
                      pos = gr0$idx, name="start", 
                      undo.splits='prune', min.width=4)
  .addSeg(st.cbs, 'firebrick4')
  
  points(gr0$idx, sort(end(gr0)), col=alpha("dodgerblue", 0.4), pch=16)
  end.cbs <- getCbsSeg(chr='chr21', val=sort(end(gr0)),
                       pos = gr0$idx, name="end",
                       undo.splits='prune', min.width=4)
  .addSeg(end.cbs, 'dodgerblue4')
}

.jaccard <- function(g0, g1){
  AnB <- intersect(g0, g1)
  AoB <- sapply(seq_along(g1), function(idx){
    width(reduce(GRanges(c(g0, g1[idx]))))
  })
  
  round((width(AnB) / AoB),3)
}

getCbsSeg <- function(chr, val, pos, name, aval=0.01, 
                      nperm=1000, weights=NULL, ...){
  chr.id <- unlist(regmatches(chr, regexec(pattern = "[0-9X]+", chr)))
  if(!chr.id %in% 'X'){
    chr.num <- as.numeric(chr.id)
  } else {
    chr.num <- chr.id
  }
  
  CNA.object <- CNA(as.numeric(val),
                    rep(chr.num, length(pos)),
                    pos, data.type="logratio", sampleid=name)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  if(!is.null(weights)){
    if(any(is.na(weights))) weights[is.na(weights)] <- 0.0001
    if(any(weights == 0)) weights[weights == 0] <- 0.0001
    segment.smoothed.CNA.object <- segment(CNA.object, verbose=1, 
                                           alpha=aval, nperm=nperm,
                                           weights=weights, ...)
  } else {
    segment.smoothed.CNA.object <- segment(CNA.object, verbose=1, 
                                           alpha=aval, nperm=nperm, ...)
  }
  return(segment.smoothed.CNA.object)
}


#### Setup Data ####
## Read in and sort all the TAD data
all.tads <- lapply(list.files(pattern="bed$"), read.table, header=F,
                   stringsAsFactors=F, check.names=F, fill=F,
                   col.names=c('chr', 'start', 'end', 'tad', 'score'))
names(all.tads) <- gsub("_TAD.*", "", list.files(pattern='bed$'))
tads.gr <- lapply(names(all.tads), function(id){
  all.tads[[id]]$ID <- id
  makeGRangesFromDataFrame(all.tads[[id]], keep.extra.columns=T)
})
tads.gr <- sort(as(tads.gr, 'GRangesList'))

# Collapse and split by chromosome
gr.t <- unlist(tads.gr)
chr.gr <- split(gr.t, seqnames(gr.t))



#### Jaccard Finding TADs ####

chr <- opt$chr


print(paste0(" > ", chr))
gr0 <- sort(chr.gr[[chr]])
max.y <- max(end(gr0))
gr0$idx <- seq_along(gr0)
jacc.cutoff <- 0.95
consensus.grl <- list()
if(opt$plot) pdf(file.path(opt$outdir, paste0("tadPlot_", chr, ".pdf")))
if(opt$plot) plotGr0(gr0, max.y, main=paste0("Initial: ", chr))

while(length(gr0) > 0 & jacc.cutoff > 0.5){
  c.gr0 <- split(gr0, gr0$ID)
  
  # Isolate a single sample to use a reference
  cl.id <- sample(names(c.gr0), 1)
  print(paste(cl.id, "-", length(gr0), "-", jacc.cutoff))
  s.gr0 <- c.gr0[[cl.id]]
  
  # Compare the single sample TADs to all other TADs
  # Assess the reproducibility of the TAD using a jaccard index
  print("Generating Jaccard matrix...")
  start_time <- Sys.time()
  jacc.mat <- mclapply(c.gr0, function(g){
    ov.idx <- as.data.frame(findOverlaps(s.gr0, g))
    spl.ov.idx <- split(ov.idx, ov.idx$queryHits)
    
    max.j <- sapply(spl.ov.idx, function(ov){
      max.jacc <- max(.jaccard(g0=s.gr0[ov$queryHits[1]],
                               g1=g[ov$subjectHits]))
      return(max.jacc)
    })
    
    return(as.data.frame(t(max.j)))
  }, mc.cores = 4)
  end_time <- Sys.time()
  t1 <- end_time - start_time
  print(paste0("Time to create jaccard matrix: ", round(t1,2), "s"))
  jacc.mat <- rbind.fill(jacc.mat)
  jacc.mat[is.na(jacc.mat)] <- 0
  avg.jacc <- round(colMeans(jacc.mat, na.rm=T),3)
  
  
  # Flag high confidence regions and collapse into consensus TADs
  print("Identifying concordant TADs...")
  consensus.ref <- s.gr0[which(avg.jacc > jacc.cutoff)]
  # Aggregate TADs for the confident regions
  ov.idx <- as.data.frame(findOverlaps(consensus.ref, gr0))
  spl.ov.idx <- split(ov.idx, ov.idx$queryHits)
  grl <- lapply(spl.ov.idx, function(ov){
    GRanges(seqnames = unique(as.character(seqnames(gr0))),
            ranges = IRanges(start = median(start(gr0[ov$subjectHits])),
                             end=median(end(gr0[ov$subjectHits]))))
  })
  consensus.gr <- unlist(as(grl, "GRangesList"))
  consensus.grl[[as.character(jacc.cutoff)]] <- consensus.gr
  
  # Remove TADs that overlap the consensus TADs
  print("plotting...")
  ov.idx <- findOverlaps(consensus.gr, gr0)
  if(length(subjectHits(ov.idx)) > 0) gr0 <- gr0[-subjectHits(ov.idx)]
  if(opt$plot) try(plotGr0(gr0, max.y, main=paste0("jaccard: ", jacc.cutoff)))
  jacc.cutoff <- (jacc.cutoff - 0.02)
}

if(opt$plot) dev.off()
print("saving...")
c.gr <- unlist(as(consensus.grl, "GRangesList"))
c.gr$jaccard <- gsub("\\.[0-9]*$", "", names(c.gr))
save(c.gr, file=file.path(opt$outdir, paste0(chr, "_TAD.rda")))
print("Saved!")


if(1==0){
  chr <- 'chr21'
  load(paste0(chr, "_TAD.rda"))
  df <- as.data.frame(c.gr)
  df$ID <- rep('consensus', nrow(df))
  df <- df[,c('ID', 'seqnames', 'start', 'end', 'width', 'jaccard')]
  colnames(df) <- c('ID', 'chrom', 'loc.start', 'loc.end',
                    'num.mark', 'seg.mean')
  write.table(df, file = paste0(chr, "_TAD.seg"), sep="\t",
              row.names = F, col.names = T, quote = F)
}



if(1==0){

  #### Identify optimal "collapse" search span ####
  ## Create 'Saturation curves' to find optimal search span size
  col <- setNames(rainbow(n = 23), paste0("chr", c(1:22, "X")))
  plot(0, type='n', xlim=c(0, 1100000), ylim=c(0, 16),
       ylab='median', xlab='span')
  search.size <- seq(10000, 1010000, by=10000)
  meds <- sapply(paste0("chr", c(1:22, "X")), function(chr){
    st0 <- start(chr.gr[[chr]])
    end0 <- end(chr.gr[[chr]])
    cnt.saturation <- t(sapply(search.size, function(ss){
      num.breaks <- max(end(chr.gr[[chr]])) / ss
      cnts <- hist(st0, breaks=num.breaks, plot = F)$counts
      pos.cnts <- cnts[which(cnts != 0)]
      c("span"=ss, 'median'=mean(pos.cnts), 'mad'=sd(pos.cnts))
    }))
    cnt.saturation <- as.data.frame(cnt.saturation)
    with(cnt.saturation, points(span, median, pch=16, col=alpha(col[chr], 0.5)))
    #with(cnt.saturation, .addSD(span, median, mad))
    loess.fit <- loess(formula = median ~ span, data = cnt.saturation, span = 0.1)$fitted
    lines(cnt.saturation$span, loess.fit, col=alpha(col[chr], 0.5))
    text(x = 1000000, y=max(cnt.saturation$median), labels = chr, 
         col=alpha(col[chr], 0.5), adj=-1)
    return(cnt.saturation$median)
  })
  
  ## Plot inflection points
  rownames(meds) <- search.size
  meds <- round(meds, 2)
  diff.meds <- apply(meds, 2, diff)
  row.x <- as.numeric(rownames(diff.meds))
  plot(0, type='n', xlim=c(0, max(row.x)), ylim=c(0, 10))
  sapply(colnames(diff.meds), function(chr){
    rect(xleft = row.x-1000, ybottom = 0, 
         xright = row.x+1000, ytop = diff.meds[,chr], 
         col=alpha(col[chr], 0.5))
    row.idx <- which(diff.meds[,chr] > 0)
    pos.diff <- diff.meds[row.idx,chr]
    text(x = row.x[row.idx], y=pos.diff, 
         labels=rep(chr, length(row.idx)), srt=90, cex=0.7, adj=-1)
  })
  idx <- diff.meds[,'chr22'] > 2 &  diff.meds[,'chr22'] < 6
  span <- row.x[which(idx)]
  
  #### PRACMA peak finding TADs ####
  lapply(chr.gr, function(gr0){
    num.bins <- ceiling((max(end(gr0)) / span*2.5) + 1)
    h.data <- hist(start(gr0), plot=T, breaks=num.bins, col=alpha("red", 0.6), border = NA)
    h.data <- hist(end(gr0), add=T, breaks=num.bins, col=alpha("blue", 0.6), border = NA)
    
    pk <- findpeaks(h.data$counts, nups = 0, ndowns = 0, zero = "0", 
                    minpeakheight = 4, minpeakdistance = 1)
    if(!is.null(pk)){
      br <- h.data$breaks[-1]
      apply(pk, 1, function(p){rect(xleft = br[p[3]], ybottom = 0, xright = br[p[4]], ytop = p[1],
                                    border=NA, col=alpha("blue", 0.5))})
      abline(v = br[pk[,2]], col="red")
    }
    return(br[pk[,2]])
  }, minpeakheight = 50, minpeakdistance = 10)
  
  
}