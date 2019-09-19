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

annotateRDS <- function(data, build='hg19', ...){
  gamma <- ASCAT.selectBestFit(fit.val, method='GoF')
  my.data <- loadBestFitRDS(gamma)
  genes <- getGenes(build)
  
  tmsg(paste0("Annotating sample: ", sample, "..."))
  seg <- annotateSegments(my.data$segments_raw, genes,
                          start.field='startpos', end.field='endpos')
  
  
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


