loadBestFitRDS <- function(gamma, ...){
  RDS.file <- file.path(sample, toupper(segmenter), 'ASCN', paste0('gamma', format(gamma, nsmall=2)),
                        paste(sample, 'ASCN', toupper(segmenter), 'RDS', sep="."))
  rds <- readRDS(RDS.file)
  return(rds)
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


annotateRDS <- function(data, build='hg19', ...){
  gamma <- ASCAT.selectBestFit(fit.val, method='GoF')
  my.data <- loadBestFitRDS(gamma)
  genes <- getGenes(build)
  
  tmsg(paste0("Annotating sample: ", sample, "..."))
  gr0 <- makeGRangesFromDataFrame(my.data$segments_raw, keep.extra.columns=TRUE, 
                                  start.field='startpos', end.field='endpos')
  
  seg <- annotateSegments(my.data$segments_raw, genes,
                          start.field='startpos', end.field='endpos')
  
  bin.size <- 50000
  windowed.bed <- genWindowedBed(bin.size=bin.size)
  
  if(map.to=='bin'){
    # Map CNV segments to a reference bed
    cl.anno <- segmentCNVs(gr0, windowed.bed, reduce='min')
  } else if(map.to == 'genes'){
    # Map CNV segments to genes
    cl.anno <- suppressMessages(annotateCnvs(cnv, txdb, 
                                             anno=txdb.genes,
                                             cols=cols))
    names(cl.anno) <- c('seg', 'genes')
  }
  
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
  
  genes0 <- genes(package)
  idx <- rep(seq_along(genes0), elementNROWS(genes0$gene_id))
  genes <- granges(genes0)[idx]
  genes$gene_id = unlist(genes0$gene_id)
  genes
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


