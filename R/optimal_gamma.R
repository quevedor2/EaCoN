.junk <- function(){
  #### setup ####
  library(EaCoN)
  library(dplyr)
  library(GenomicRanges)
  pdir <- "/mnt/work1/users/pughlab/projects/CCLE/eacon"
  gistic.dir <- '/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/gistic/output'
  setwd(pdir)
  
  all.samples <- list.files(pattern="gammaEval.txt$", include.dirs = T, recursive = T)
  all.samples <- sapply(strsplit(all.samples, "/"), function(x) x[1])
  #sample <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_A05_256010' # VMRCRCZ_KIDNEY	VMRC-RCZ	renal_cell_carcinoma
  #sample <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_B05_256034' # NCIH1975_LUNG non_small_cell_carcinoma
  #sample <- 'ARLES_p_NCLE_DNAAffy2_S_GenomeWideSNP_6_C07_256062' # OVCAR8_OVARY ovary	carcinoma
  
  segmenter <- 'ASCAT'
  gamma <- 0.7
  #my.data <- EaCoN:::loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
  my.data <- loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
  pancan.dir <- '/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/raw'
  pancan.obj <- createPancanSegRef(pancan.dir, out.dir='/mnt/work1/users/pughlab/references/TCGA/TCGA_Pancan_ploidyseg/cleaned',
                                   write.seg=T)
}

#### Functions ####
loadBestFitRDS <- function(gamma, segmenter, sample){
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
  
  ## Create a data structure of ploidy probabilities per TCGA onco-codes
  cancer.type.ploidy <- split(anno.ploidy.seg, f=anno.ploidy.seg$`cancer type`)
  col.of.interest <- c('purity', 'ploidy', 'Cancer DNA fraction', 'Subclonal genome fraction')
  ## Get the counts per breakpoint
  ctp.breaks <- lapply(cancer.type.ploidy, function(ctp0){
    tcga.code <- unique(ctp0$`cancer type`)
    coi.list <- lapply(col.of.interest, function(coi){
      breaks <- switch(coi,
                       ploidy=seq(0.05, 10.05, by=0.1),
                       seq(-0.005, 1.005, by=0.01))
      gap <- diff(breaks)[1] / 2
      coi.breaks <- hist(ctp0[,coi], breaks=breaks, plot = F)
      data.frame("breaks"= coi.breaks$breaks[-length(coi.breaks$breaks)] + gap,
                 "P"=round(coi.breaks$counts / sum(coi.breaks$counts),3))
    })
    names(coi.list) <- col.of.interest
    coi.list
  })
  ## Reduce into a single matrix
  ctp.break.mats <- lapply(col.of.interest, function(coi){
    coi.mat <- Reduce(function(x,y) merge(x,y,by='breaks'),
                      lapply(ctp.breaks, function(cb) cb[[coi]]))
    colnames(coi.mat) <- c('breaks', names(cancer.type.ploidy))
    coi.mat$AVG <- rowMeans(coi.mat[,-1])
    plot(coi.mat[,c('breaks', 'AVG')])
    coi.mat
  })
  names(ctp.break.mats) <- col.of.interest
  
  
  
  pancan.obj <- list('AP-meta'=anno.ploidy,
                     'APS-meta'=anno.ploidy.seg,
                     'seg'=seg,
                     'breaks'=ctp.break.mats)
  
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
  if(class(gr) == 'list'){
    gr <- unlist(as(gr, 'GRangesList'))
  } else if(class(gr) == 'GRangesList'){
    gr <- unlist(gr)
  }
  
  gr.chr <- split(gr, seqnames(gr))
  gr.c <- as(lapply(gr.chr, function(gr0){
    pos <- unique(sort(c(start(gr0), end(gr0)+1)))
    GRanges(seqnames = rep(unique(as.character(seqnames(gr0))), length(pos)-1),
            IRanges(start=pos[-length(pos)],
                    end=pos[-1] - 1))
  }), 'GRangesList')
  gr.c <- unlist(gr.c)
  gr.c
}

populateGr <- function(grl, ref.gr, col.id=NULL){
  if(is.null(names(grl))){
    names(grl) <- letters[1:length(grl)]
  }
  
  for(each.id in names(grl)){
    ov.idx <- findOverlaps(grl[[each.id]], ref.gr)
    if(is.null(col.id)){
      mcol.mat <- as.data.frame(matrix(as.numeric(),
                         nrow = length(ref.gr),
                         ncol = ncol(mcols(grl[[each.id]]))))
      
      mcol.mat.tmp <- as.data.frame(mcols(grl[[each.id]])[queryHits(ov.idx),,drop=F])
      mcol.mat[subjectHits(ov.idx),] <- mcol.mat.tmp
      colnames(mcol.mat) <- colnames(mcol.mat.tmp)

      mcols(ref.gr) <- cbind(mcols(ref.gr), mcol.mat)
    } else {
      mcols(ref.gr)[,each.id] <- NA
      mcols(ref.gr)[subjectHits(ov.idx),each.id] <- mcols(grl[[each.id]])[queryHits(ov.idx),col.id]
    }
    
  }
  ref.gr
}

reduceGr <- function(gr, sig=1){
  gr <- sort(gr)
  mcols(gr) <- round(as.matrix(mcols(gr)), sig)
  
  runs <- apply(mcols(gr),1,paste, collapse="_")
  rle.runs <- rle(runs)
  
  idx <- cumsum(rle.runs$lengths)
  merge.df <- data.frame("start"=c(1, idx[-length(idx)]+1),
                          "end"=idx)
  
  nonmerge.idx <- which(rle.runs$length == 1)
  nonmerge.gr <- gr[merge.df[nonmerge.idx,1]]
  
  merge.idx <- which(rle.runs$length > 1)
  m.gr <- unlist(reduce(as(apply(merge.df[merge.idx,], 1, function(x){
    gr[x[1]:x[2],]
  }), "GRangesList")))
  m.gr <- sort(unlist(as(list(nonmerge.gr, m.gr), "GRangesList")))
  
  mcols(m.gr) <- mcols(gr)[merge.df$start,]
  m.gr
}

sampleGisticProfiles <- function(gr, bootstrap=1000, seed=1234){
  amplitude.idx <- grep("A_", colnames(mcols(gr)))
  frequency.idx <- grep("F_", colnames(mcols(gr)))
  
  freq.mat <- as.matrix(mcols(gr)[,frequency.idx])
  amp.mat <- as.matrix(mcols(gr)[,amplitude.idx])
  
  set.seed(seed)
  ## Select an element from each row using the probability matrix (row by row)
  # Random uniform sampling for each row [0-1]
  random <- matrix(round(runif(nrow(freq.mat) * bootstrap),3), ncol = bootstrap)
  
  # Calculate cumulative sum for each row 
  cumul.w <- freq.mat %*% upper.tri(diag(ncol(freq.mat)), diag = TRUE) / rowSums(freq.mat)
  cumul.w[,3] <- 1
  # Identify which range the uniform sampling falls in
  s.bootstrap <- apply(random, 2, function(x){
    rowSums(x > cumul.w) + 1L
  })
  
  amp.s.mat <- apply(s.bootstrap, 2, function(s) amp.mat[cbind(seq_along(s), s)])
  freq.s.mat <- apply(s.bootstrap, 2, function(s) freq.mat[cbind(seq_along(s), s)])
  
  list("A"=amp.s.mat,
       "F"=freq.s.mat)
}

fillMissingGistic <- function(gr, col.check='A_Amp'){
  na.rows <- which(is.na(mcols(gr)[,col.check]))
  na.cols <- grep("[FA]_", colnames(mcols(gr)), invert = F)
  mcols(gr)[na.rows,na.cols] <- 0
  mcols(gr)[na.rows, 'F_Neut'] <- 1
  
  ## Remove any NA rows for the given seg
  seg.col <- grep("[FA]_", colnames(mcols(gr)), invert = T)
  if(any(is.na(mcols(gr)[,seg.col]))){
    gr <- gr[-which(is.na(mcols(gr)[,seg.col]))]
  }
  gr
}

calcABC <- function(gr, bootstrap.gr){
  seg.col <- grep("[FA]_", colnames(mcols(gr)), invert = T)
  abc <- sweep(bootstrap.gr[['A']], 1, as.matrix(mcols(gr)[,seg.col]), "-")
  abc <- sweep(abc, 1, width(gr), "*")
  freq <- colSums(log2(bootstrap.gr[['F']]+0.01), na.rm=T)
  
  ## Calculate the best fit model
  amplitude.idx <- grep("A_", colnames(mcols(gr)))
  frequency.idx <- grep("F_", colnames(mcols(gr)))
  seg.mean.idx <- grep("seg.mean", colnames(mcols(gr)))
  
  freq.mat <- as.matrix(mcols(gr[,c(frequency.idx)]))
  amp.mat <- as.matrix(mcols(gr[,c(amplitude.idx)]))
  seg.mat <- as.matrix(mcols(gr[,c(seg.mean.idx)]))
  
  abc.seg.mat <- abs(sweep(amp.mat, 1, seg.mat, '-'))
  abc.seg.mat <- sweep(abc.seg.mat, 1, width(gr), "*")
  min.idx <- apply(abc.seg.mat, 1, which.min)
  min.abc <- data.frame("A"=sum(abc.seg.mat[cbind(seq_along(min.idx), min.idx)], na.rm=T), 
                        "F"=sum(log2(freq.mat[cbind(seq_along(min.idx), min.idx)]+0.01), na.rm = T))
  
  list("bootstrap"=data.frame("A"=colSums(abs(abc)),
                              "F"=freq),
       "best.fit"=min.abc)
}


#### Load in GISTIC Cancer Type ####
.test <- function(){
  if(!file.exists(file.path(gistic.dir, "gistic_scores.RData"))){
    tumor.types <- list.files(gistic.dir)
    tumor.types <- tumor.types[-grep("LUAD", tumor.types)]
    all.gistics <- lapply(tumor.types, function(tumor.type){
      print(tumor.type)
      gistic.path <- file.path(gistic.dir, tumor.type, "scores.gistic")
      gistic <- read.table(gistic.path, header=T, stringsAsFactors = F, check.names = F, sep="\t")
      gr <- makeGRangesFromDataFrame(gistic, keep.extra.columns = T)
      seqlevelsStyle(gr) <- 'UCSC'
      grl <- split(gr, gr$Type)
      gr.c.raw <- combineGr(gr)
      gr.c.amp <- populateGr(grl, gr.c.raw, col.id = 'average amplitude')
      gr.c.amp$Neut <- 0
      gr.c.freq <- populateGr(grl, gr.c.raw, col.id = 'frequency')
      gr.c <- populateGr(list(gr.c.amp, gr.c.freq), gr.c.raw)
      colnames(mcols(gr.c)) <- paste0(c(rep("A_", 3), rep("F_",2)), 
                                      colnames(mcols(gr.c)))
      
      gr.c <- reduceGr(gr.c, sig=1)
      gr.c$F_Neut <- apply(mcols(gr.c), 1, function(x) 1-sum(x[c('F_Del', 'F_Amp')]))
      gr.c
    })
    names(all.gistics) <- tumor.types
    
    save(all.gistics, file=file.path(gistic.dir, "gistic_scores.RData"))
  } else {
    load(file=file.path(gistic.dir, "gistic_scores.RData"))
  }
  
  #### Load in Seg File####
  sample.abcs <- lapply(all.samples[1:10], function(sample){
    print(paste0("Processing ", sample, "..."))
    my.data <- loadBestFitRDS(gamma=gamma,sample=sample, segmenter=segmenter)
    
    seg <- my.data$segments_raw
    seg$width <- round(with(seg, endpos - startpos)/100000,0) + 1
    seg$nABraw <- with(seg, nAraw + nBraw)
    seg$seg.mean <- round(log2(seg$nABraw / 2),3)
    seg$seg.mean <- seg$seg.mean - median(rep(seg$seg.mean, seg$width))
    seg$seg.mean[seg$seg.mean < -2] <- -2
    seg.gr <- makeGRangesFromDataFrame(seg, start.field = 'startpos', 
                                       end.field = 'endpos', keep.extra.columns = T)
    mcols(seg.gr) <- mcols(seg.gr)[,'seg.mean',drop=F]
    
    #### Compare ABCs ####
    all.abcs <- lapply(all.gistics, function(gr.c, seg.gr){
      print("...")
      gr.c.raw <- combineGr(list(seg.gr, gr.c))
      gr.gistic.seg <- populateGr(list(gr.c, seg.gr), gr.c.raw)
      gr.gistic.seg <- reduceGr(gr.gistic.seg, sig=1)
      gr.gistic.seg$A_Del <- -1 * gr.gistic.seg$A_Del
      gr.gistic.seg <- fillMissingGistic(gr.gistic.seg)
      
      
      #### Exhaustively compare seg to GISTIC Cancer type ####
      gr.c$A_Del <- -1 * gr.c$A_Del
      sampled.gistic <- sampleGisticProfiles(gr.gistic.seg, bootstrap=1000, seed=1234)
      
      #gr.gistic.seg <- mapSegToRef(seg.gr, gr.c)
      bootstrap.abc <- calcABC(gr.gistic.seg, sampled.gistic)
      
      bootstrap.abc
    }, seg.gr=seg.gr)
    
    
    mean.abc <- as.data.frame(t(sapply(all.abcs, function(x) colMeans(x[[1]]))))
    min.abc <- as.data.frame(t(sapply(all.abcs, function(x) unlist(x[[2]]))))
    z.abc <- t(sapply(seq_along(all.abcs), function(i, mean.abc, min.abc){
      .l2 <- function(x, y, c1, c2){
        sqrt((x-c1)^2 + (y-c2)^2)
      }
      d <- .l2(x=all.abcs[[i]][[1]][,1], y=all.abcs[[i]][[1]][,2],  
               mean.abc[i,1], mean.abc[i,2])
      sd <- sum(d)/nrow(all.abcs[[i]][[1]])
      
      min.d <- .l2(min.abc[i,1], min.abc[i,2], mean.abc[i,1], mean.abc[i,2])
      c("sd"=sd,
        "d"=min.d,
        "z"=(min.d / sd))
    }, mean.abc=mean.abc, min.abc=min.abc))
    z.abc <- as.data.frame(z.abc)
    rownames(z.abc) <- names(all.abcs)
    
    return(list("mean"=mean.abc, "min"=min.abc, "z"=z.abc))
  })
  names(sample.abcs) <- all.samples[1:10]
  
  lapply(sample.abcs, function(x) {
    mean.abc <- x[['mean']]
    min.abc <- x[['min']]
    z.abc <- x[['z']]
    
    ranks <- rank(mean.abc$A) + rank(abs(mean.abc$F)) + rank(z.abc$z) + rank(min.abc$A)
    names(ranks) <- rownames(mean.abc)
    head(sort(ranks),3)
    
    head(z.abc[order(z.abc$z),], 3)
  })
  
  
  
  
  
  
  
  
  
  mean.abc[order(mean.abc$A),]
  mean.abc[order(mean.abc$F, decreasing = T),]
  min.abc[order(min.abc$A),]
  z.abc[order(z.abc$z),]
  
  ranks <- rank(mean.abc$A) + rank(abs(mean.abc$F)) + rank(z.abc$z) + rank(min.abc$A)
  names(ranks) <- rownames(mean.abc)
  sort(ranks)
  
  
  #### Machine Learning on Pancan Seg ####
  sample.ids <- names(pancan.obj[['seg']])
  pancan.gr <- lapply(pancan.obj[['seg']], function(x){
    seg.df <- x[,c('Chromosome', 'Start', 'End', 'Segment_Mean')]
    seg.df$Segment_Mean <- round(seg.df$Segment_Mean, 2)
    makeGRangesFromDataFrame(seg.df, keep.extra.columns=T)
  })
  pancan.gr.c <- combineGr(pancan.gr)
  min.bin.size <- 500
  pancan.gr.c <- pancan.gr.c[-which(width(pancan.gr.c) < min.bin.size),]
  pancan.gr.seg <- populateGr(pancan.gr, pancan.gr.c)
  pancan.gr.seg
  
  
  
  #### OLD METHOD ####
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
  
}


