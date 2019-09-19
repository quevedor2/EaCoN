loadBestFitRDS <- function(gamma, ...){
  RDS.file <- file.path(sample, toupper(segmenter), 'ASCN', paste0('gamma', format(gamma, nsmall=2)),
                        paste(sample, 'ASCN', toupper(segmenter), 'RDS', sep="."))
  rds <- readRDS(RDS.file)
  return(rds)
}

cleanGR <- function(gr0){
  for(i in c(1:ncol(elementMetadata(gr0)))){
    icol <- elementMetadata(gr0)[,i]
    if(class(icol)=='numeric'){
      elementMetadata(gr0)[,i] <- round(icol, 3)
    }
  }
  gr0$TCN <- rowSums(as.matrix(elementMetadata(gr0)[,c('nMajor', 'nMinor')]))
  gr0$seg.mean <- round(log2(gr0$TCN / 2),3)
  gr0$seg.mean[gr0$seg.mean < log2(1/50)] <- round(log2(1/50), 2)
  gr0
}

getGenes <- function(genome.build="hg19"){
  switch(genome.build,
         hg18={ 
           suppressPackageStartupMessages(require("TxDb.Hsapiens.UCSC.hg18.knownGene"))
           package <- TxDb.Hsapiens.UCSC.hg18.knownGene 
         },
         hg19={ 
           suppressPackageStartupMessages(require("TxDb.Hsapiens.UCSC.hg19.knownGene"))
           package <- TxDb.Hsapiens.UCSC.hg19.knownGene 
         },
         hg38={ 
           suppressPackageStartupMessages(require("TxDb.Hsapiens.UCSC.hg38.knownGene"))
           package <- TxDb.Hsapiens.UCSC.hg38.knownGene 
         },
         stop("genome must be hg18, hg19, or hg38"))
  
  list(txdb=package,
       txdb.genes=genes(package))
}

getMapping <- function(in.col='ENTREZID', 
                       out.cols=c("SYMBOL", "ENSEMBL")){
  gene.map <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db, in.col), 
                     keytype="ENTREZID", columns=out.cols)
  gene.map
}

ASCAT.selectBestFit <- function(fit.val, method='GoF'){
  idx <- switch(method,
                "GoF"=which.max(fit.val$GoF),
                'psi'=which.min(fit.val$psi))
  return(fit.val$gamma[idx])
}

genWindowedBed <- function(bin.size=1000000, seq.style="UCSC"){
  chrs <- seqlengths(Hsapiens)[paste0("chr", c(1:22,"X", "Y"))]
  
  ## Construct intervals across the genome of a certain bin size
  start.points <- seq(1, 500000000, by=bin.size)
  grl <- lapply(names(chrs), function(chr.id){
    chr <- chrs[chr.id]
    ir <- IRanges(start=start.points[start.points < chr], width=bin.size)
    end(ir[length(ir),]) <- chr
    gr <- GRanges(seqnames = chr.id, ir)
    gr
  })
  
  ## Assemble all GRanges and set seq level style
  grl <- as(grl, "GRangesList")
  suppressWarnings(seqlevelsStyle(grl) <- seq.style)
  gr <- unlist(grl)
  gr
}

segmentCNVs <- function(cnv, bed, reduce='mean'){
  olaps = findOverlaps(cnv, bed)
  
  # Flag BED bins that map to multiple CNVs
  dup.idx <- which(duplicated(subjectHits(olaps), fromLast=TRUE))
  dup.idx <- c(dup.idx, which(duplicated(subjectHits(olaps), fromLast=FALSE)))
  dup.idx <- sort(dup.idx)
  
  # Use a summary metric (Default=mean) to reduce the CNV information that
  # spans multiple bed windows
  dup.df <- as.data.frame(olaps[dup.idx,])
  dup.spl <- split(dup.df, dup.df$subjectHits)
  dup.em <- lapply(dup.spl, function(i) {
    em.mat <- as.matrix(elementMetadata(cnv[i$queryHits,]))
    switch(reduce,
           mean=colMeans(em.mat),
           min=do.call(pmin, lapply(1:nrow(em.mat), function(i) em.mat[i,])),
           max=do.call(pmax, lapply(1:nrow(em.mat), function(i) em.mat[i,])),
           median=do.call(median, lapply(1:nrow(em.mat), function(i) em.mat[i,])))
    
  })
  dup.em <- do.call(rbind, dup.em)
  
  # Initialize a metadata matrix and populate it for the BED GRanges object
  em  <- matrix(nrow=length(bed), 
                ncol=ncol(elementMetadata(cnv)), 
                dimnames = list(NULL,colnames(elementMetadata(cnv))))
  dedup.olaps <- olaps[-dup.idx,]
  em[subjectHits(dedup.olaps),] <- as.matrix(elementMetadata(cnv)[queryHits(dedup.olaps),])
  em[unique(subjectHits(olaps)[dup.idx]),] <- dup.em 
  
  # Append metadata and return
  em <- as.data.frame(em)
  bed$ID <- em$ID <- paste0("bin_", c(1:nrow(em)))
  
  return(list(seg=bed, genes=em))
}

annotateCNVs <- function(cnv, txdb, anno=NULL,
                         cols=c("seg.mean", "nA", "nB")){
  stopifnot(is(cnv, "GRanges"), is(txdb, "TxDb"))
  
  ## Assign EntrezID to each segment 
  if(is.null(anno)) anno = genes(txdb)
  olaps = findOverlaps(cnv, anno)
  mcols(olaps)$gene_id = anno$gene_id[subjectHits(olaps)]  # Fixed the code here
  cnv_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  cnv$gene_id = IRanges::splitAsList(mcols(olaps)$gene_id, cnv_factor)
  
  ## Assign EntrezID to each segment 
  seg.entrez <- apply(as.data.frame(mcols(cnv)), 1, function(i){
    ids <- unlist(strsplit(x = as.character(unlist(i[['gene_id']])), split=","))
    segs <- do.call(rbind, replicate(length(ids), round(unlist(i[cols]),3), simplify = FALSE))
    
    as.data.frame(cbind(segs, 'ENTREZ'=ids))
  })
  seg.entrez <- do.call(rbind, seg.entrez)
  if(any(duplicated(seg.entrez$ENTREZ))) seg.entrez <- seg.entrez[-which(duplicated(seg.entrez$ENTREZ)),]
  
  
  ## Map ensembl and HUGO IDs to the ENTREZ ids
  seg.anno <- merge(seg.entrez, getMapping(),
                    by.x="ENTREZ", by.y="ENTREZID", all.x=TRUE)
  seg.anno <- seg.anno[-which(duplicated(seg.anno$ENTREZ)),]
  for(each.col in cols){
    seg.anno[,each.col] <- as.numeric(as.character(seg.anno[,each.col]))
  }
  
  list("seg"=cnv, "genes"=seg.anno)  
}

annotateRDS <- function(fit.val, sample, segmenter, build='hg19', 
                        bin.size=50000, ...){
  ## Assemble the ASCAT Seg file into a CNV GRanges object
  gamma <- ASCAT.selectBestFit(fit.val, method='GoF')
  my.data <- loadBestFitRDS(gamma)
  genes <- getGenes(build)
  
  tmsg(paste0("Annotating sample: ", sample, "..."))
  cnv <- makeGRangesFromDataFrame(my.data$segments_raw, keep.extra.columns=TRUE, 
                                  start.field='startpos', end.field='endpos')
  cnv <- cleanGR(cnv)
  
  ## Annotate the CNVs based on:
  cl.anno <- list()
  # Bins
  windowed.bed <- genWindowedBed(bin.size=bin.size)
  cl.anno[['bins']] <- segmentCNVs(cnv, windowed.bed, reduce='min')
  
  # Genes
  cols <- c('nMajor', 'nMinor', 'nAraw', 'nBraw', 'TCN', 'seg.mean')
  cl.anno[['genes']] <- suppressMessages(annotateCNVs(cnv, genes$txdb, 
                                                      anno=genes$txdb.genes, cols=cols))
  cl.anno$genes$genes <- cl.anno$genes$genes[-which(is.na(cl.anno$genes$genes$SYMBOL)),]
  
  # Raw seg
  cl.anno[['seg']] <- as.data.frame(cnv)
  
  return(cl.anno)
}



all.fits <- list()
all.fits[[sample]] <- list("fit"=fit.val,
                           "sample"=sample)
gr.cnv <- lapply(all.fits, function(x) annotateRDS(x$fit, x$sample, segmenter, 
                                                   build='hg19', bin.size=50000))
cbio.path=file.path("out", "cBio")

buildCbioOut(gr.cnv, cbio.path="./out/cBio", overwrite=sample)


buildCbioOut <- function(gr.cnv, cbio.path="./out/cBio", pattern="_CNA", 
                         cbio.cna.file=NULL, cbio.linear.file=NULL, cbio.seg.file=NULL,
                         amp.thresh=5, add.on.to.existing=TRUE, ...){
  .checkFile <- function(cbio.path, file.id, pat){
    idx <- grep(list.files(cbio.path), pattern=pat, perl=TRUE)[1]
    if(!is.na(idx)){
      cbio.file <- list.files(cbio.path)[idx]
      exists.stat <- TRUE
    } else {
      cbio.file <- file.id
      exists.stat <- FALSE
    }
    c('file'=cbio.file, 'exists'=exists.stat)
  }
  
  .adjustCnaMat <- function(cnv.mat, mat.type, amp.thresh=NULL, ord=NULL){
    tcn.mat <- cnv.mat[,-c(1,2),drop=FALSE]
    if(mat.type=='CNA'){
      ## Convert Total CN to the -2, -1, 0, 1, 2, standards
      tcn.mat[tcn.mat >= amp.thresh] <- amp.thresh
      tcn.mat <- tcn.mat - 2
      for(i in c(1:(amp.thresh-3))){ tcn.mat[tcn.mat == i] <- 1 }
      tcn.mat[tcn.mat == (amp.thresh-2)] <- 2
    } else if(mat.type=='linear') {
      tcn.mat <- round(tcn.mat, 2)
    }
    
    ## Recombine the CNV mat
    colnames(tcn.mat) <- names(gr.cnv)
    colnames(cnv.mat)[c(1,2)] <- c('Hugo_Symbol', 'Entrez_Gene_Id')
    cnv.mat <- cbind(cnv.mat[,c(1:2)], tcn.mat)
    
    ## Set the order
    if(is.null(ord)){
      cnv.mat <- cnv.mat[order(cnv.mat$Hugo_Symbol),]
    } else {
      cnv.mat <- cnv.mat[match(ord, cnv.mat$Hugo_Symbol),]
      na.idx <- apply(cnv.mat, 1, function(x) all(is.na(x)))
      if(any(na.idx)) cnv.mat <- cnv.mat[-which(na.idx),]
    }
    return(cnv.mat)
  }
  
  tmsg(paste0("Building a cBioportal Object..."))
  ## Locating existing cBioportal Objects
  suppressWarnings(dir.create(cbio.path, recursive = TRUE))
  if(is.null(cbio.cna.file)){
    cbio.cna.file <- .checkFile(cbio.path, 'data_CNA.txt',   pat=paste0("(?<!linear)", pattern))
  }
  if(is.null(cbio.linear.file)){
    cbio.linear.file <- .checkFile(cbio.path, 'data_linear_CNA.txt',   pat=paste0("linear", pattern))
  }
  if(is.null(cbio.seg.file)){
    cbio.seg.file <- .checkFile(cbio.path, 'data_cna_hg19.seg',   pat='seg$')
  }

  ## Create the data_CNA.txt file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#discrete-copy-number-data
  cna.mat <- Reduce(function(x,y) merge(x,y, by=c('SYMBOL', 'ENTREZ'), all.x=TRUE, all.y=TRUE),
                    lapply(gr.cnv, function(cnv){ cnv$genes$genes[,c('SYMBOL', 'ENTREZ', 'TCN')]}))
  cna.mat <- .adjustCnaMat(cna.mat, mat.type = 'CNA', amp.thresh = amp.thresh, ord = NULL)

  ## Create the data_linear_CNA.txt file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#continuous-copy-number-data
  linear.mat <- Reduce(function(x,y) merge(x,y, by=c('SYMBOL', 'ENTREZ'), all.x=TRUE),
                       lapply(gr.cnv, function(cnv){ cnv$genes$genes[,c('SYMBOL', 'ENTREZ', 'seg.mean')]}))
  linear.mat <- .adjustCnaMat(linear.mat, mat.type = 'linear', ord = cna.mat$Hugo_Symbol)
  
  ## Create the data_cna_hg19.seg data file: https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#segmented-data
  segs <- do.call("rbind", lapply(gr.cnv, function(x) x$seg))
  segs$ID <- gsub(".[0-9]*$", "", rownames(segs))
  segs <- segs[,c('ID', 'seqnames', 'start', 'end', 'width', 'seg.mean')]
  colnames(segs) <- c('ID', 'chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')
  
  ## If existing cBio objects exist, append to the existing data structure
  if(add.on.to.existing){
    if(as.logical(cbio.cna.file['exists'])){
      cna.mat <- .appendToCbioMat(cbio.path, cbio.cna.file, cna.mat, ...)
    }
    if(as.logical(cbio.linear.file['exists'])){
      linear.mat <- .appendToCbioMat(cbio.path, cbio.linear.file, linear.mat, ...)
    }
    if(as.logical(cbio.seg.file['exists'])){
      segs <- .appendToCbioSeg(cbio.path, cbio.seg.file, segs, ...)
    }
  }
  
  ## Write cBioportal Matrices
  .write <- function(...){
    write.table(..., sep="\t", col.names=TRUE, row.names=FALSE, quote=F)
  }
  .write(x=segs, file=file.path(cbio.path, cbio.seg.file['file']))
  .write(x=cna.mat, file=file.path(cbio.path, cbio.cna.file['file']))
  .write(x=linear.mat, file=file.path(cbio.path, cbio.linear.file['file']))
}

.appendToCbioSeg <- function(cbio.path, cbio.file, seg, overwrite=NULL){
  exist.seg <- read.table(file.path(cbio.path, cbio.file['file']), sep="\t", header=T,
                          stringsAsFactors = F, check.names = F, fill=F)
  exist.spl <- split(exist.seg, f=exist.seg$ID)
  if(!is.null(overwrite)) {
    ov.idx <- sapply(overwrite, function(id) grep(id, names(exist.spl)))
    tmsg(paste0("Overwriting samples: ", paste(names(exist.spl)[ov.idx], collapse=",")))
    exist.spl[ov.idx] <- NULL
  }
  
  new.spl <- split(seg, f=seg$ID)
  new.ids <- which(!names(new.spl) %in% names(exist.spl))
  
  if(length(new.ids) > 0){
    tmsg(paste0("New samples being added to cBio Seg file : ", 
                paste(names(new.spl)[new.ids], collapse=",")))
    seg <- do.call(rbind, append(new.spl[new.ids], exist.spl))
  } else {
    tmsg("No new samples to add.  If you want to replace an existing sample, please specify 
         using overwrite=c('SampleA', 'SampleB')")
  }
  return(seg)
}

.appendToCbioMat <- function(cbio.path, cbio.file, cnv.mat, overwrite=NULL){
  exist.cna <- read.table(file.path(cbio.path, cbio.file['file']), sep="\t", header=T,
                          stringsAsFactors = F, check.names = F, fill=F)
  if(!is.null(overwrite)) {
    ov.idx <- sapply(overwrite, function(id) grep(id, colnames(exist.cna)))
    tmsg(paste0("Overwriting samples: ", paste(colnames(exist.cna)[ov.idx], collapse=",")))
    exist.cna <- exist.cna[,-ov.idx]
  }
  new.cols <- which(!colnames(cnv.mat) %in% colnames(exist.cna))
  
  if(length(new.cols) > 0){
    tmsg(paste0("New samples being added to cBio CNA matrices: ", 
                paste(colnames(cnv.mat)[new.cols], collapse=",")))
    cnv.mat <- merge(exist.cna, cnv.mat[,c(1,2, new.cols)], 
                     by=c('Hugo_Symbol', 'Entrez_Gene_Id'), all.x=TRUE)
  } else {
    tmsg("No new samples to add.  If you want to replace an existing sample, please specify 
         using overwrite=c('SampleA', 'SampleB')")
  }
  return(cnv.mat)
}

buildPSetOut <- function(){

  
}


