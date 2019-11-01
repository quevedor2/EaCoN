#### Objective ####
# Create the Pset for Signature and CN components
library(NMF)
library(Biobase)
library(EaCoN)

#### Functions ####
getColvOrder <- function(bmat, ref){
  cor.mat <- apply(ref,2,function(x){
    apply(bmat,2,cor,x,method="pearson")
  })
  ord <- setNames(rep(NA, ncol(cor.mat)), nm = colnames(cor.mat))
  
  while(!all(is.na(cor.mat))){
    unique.val <- F
    while(!unique.val){
      idx <- which(cor.mat == max(cor.mat, na.rm=T), arr.ind = T)[1,]
      if(idx[1] %in% ord){
        cor.mat[idx[1],idx[2]] <- NA
        unique.val <- F
      } else {
        unique.val <- T
      }
    }
    
    cor.mat[,idx[2]] <- NA
    ord[idx[2]] <- idx[1]
  }
  ord
}


#### Load in Data ####
pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
ccl.comp <- readRDS(file.path(pdir, "rds", "ccl_comp.RDS"))
ccl.nmf <- readRDS(file.path(pdir, "rds", "ccl_nmf.RDS"))
basis.mats <- lapply(ccl.nmf, basis)
coef.mats <- lapply(ccl.nmf, coef)

## Synchronize the signatures
sig.ord <- lapply(basis.mats, getColvOrder, ref=basis.mats[[1]])
coef.mats <- lapply(names(coef.mats), function(dataset){
  cm <- t(coef.mats[[dataset]][sig.ord[[dataset]],])
  colnames(cm) <- paste0("cnsig_", 1:ncol(cm))
  cm
})
names(coef.mats) <- names(basis.mats)

#### Create PSet ####
dataset <- 'CCLE'
for(dataset in c("CCLE", "GDSC")){
  print(dataset)
  ## Load in Metadata
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
  
  ## Assemble assayData
  eset.env <- new.env()
  exprs <- cbind(coef.mats[[dataset]], ccl.comp[[dataset]])
  assign("exprs", t(exprs), envir=eset.env)
  
  
  ## Assemble phenoData 
  if(verbose) print("Assembling phenoData...")
  meta <- EaCoN:::.overlapMetaWithExprs(exprs=eset.env$exprs, meta=meta.l$meta)
  cl.phenoData <- new("AnnotatedDataFrame", data=meta)
  
  ## Assemble featureData
  if(verbose) print("Assembling featureData [Genes]...")
  gene.fdata <- AnnotatedDataFrame(data=as.data.frame(matrix(nrow=nrow(eset.env$exprs),ncol=0)),
                                   varMetadata=data.frame(labelDescription=c()))
  rownames(gene.fdata) <- rownames(eset.env$exprs)
  
  comp.eset <- ExpressionSet(assayData=eset.env,
                             phenoData=cl.phenoData,
                             annotation=dataset,
                             featureData=gene.fdata)
  saveRDS(comp.eset, file=file.path(pdir, "rds", paste0(dataset, "_CN.comp.RDS")))
}
