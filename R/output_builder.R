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
  gr0
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


annotateRDS <- function(data, build='hg19', bin.size=50000, ...){
  gamma <- ASCAT.selectBestFit(fit.val, method='GoF')
  my.data <- loadBestFitRDS(gamma)
  genes <- getGenes(build)
  
  tmsg(paste0("Annotating sample: ", sample, "..."))
  cnv <- makeGRangesFromDataFrame(my.data$segments_raw, keep.extra.columns=TRUE, 
                                  start.field='startpos', end.field='endpos')
  cnv <- cleanGR(cnv)
  
  windowed.bed <- genWindowedBed(bin.size=bin.size)
  cl.anno <- segmentCNVs(cnv, windowed.bed, reduce='min')
  
  cols <- c('nMajor', 'nMinor', 'nAraw', 'nBraw', 'TCN', 'seg.mean')
  cl.anno <- suppressMessages(annotateCNVs(cnv, genes$txdb, 
                                           anno=genes$txdb.genes, cols=cols))
  cl.anno$genes <- cl.anno$genes[-which(is.na(cl.anno$genes$SYMBOL)),]
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

annotateSegments <- function(cn.data, genes, ...){
  suppressPackageStartupMessages(require(org.Hs.eg.db))
  suppressPackageStartupMessages(require(GenomicRanges))
  suppressPackageStartupMessages(require(AnnotationDbi))
  gr0 <- cn.data
  if(class(cn.data) != 'GRanges') gr0 <- makeGRangesFromDataFrame(cn.data, keep.extra.columns=TRUE, ...)
  #if(class(gr0) != 'GRanges') stop("Input data could not be converted to a GRanges object")
  seqlevelsStyle(gr0) <- "UCSC"  # current: NCBI
  if(is(genes, "TxDb")) anno = genes(genes) else anno=genes
  
  olaps <- findOverlaps(anno, gr0, type="any")
  idx <- factor(subjectHits(olaps), levels=seq_len(subjectLength(olaps)))
  gr0$gene_ids <- splitAsList(anno$gene_id[queryHits(olaps)], idx)
  gr0$gene_ids <- lapply(gr0$gene_ids, function(input.id) {
    if(length(input.id) > 0){ 
      tryCatch({
        suppressMessages(mapIds(org.Hs.eg.db,
                                keys=input.id,
                                column="SYMBOL",
                                keytype="ENTREZID",
                                multiVals="first"))
      }, error=function(e){NULL})
    } else { NA }
  })
  return(gr0)
}






buildCbioOut <- function(){
  
}

buildPSetOut <- function(){

  
}


