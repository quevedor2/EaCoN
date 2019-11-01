combine_psets <- function(){
  require(Biobase)
  dataset <- 'GDSC'
  pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
  pset.dir <- file.path(pdir, dataset, "eacon", "out", "PSet")
  for(data.type in c('bins', 'tads')){
    print(data.type)
    all.files <- list.files(pset.dir, pattern=paste0(data.type, ".*RData"))
    psets <- lapply(all.files[1:10], function(file){
      load(file.path(pset.dir, file))
      eset <- switch(data.type,
                     gene=gene.eset,
                     bins=out.eset,
                     tads=out.eset)
      return(eset)
    })
    
    pset1 <- Reduce(function(x,y) combine(x,y), psets)
    rm(psets)
    
    psets <- lapply(all.files[11:length(all.files)], function(file){
      load(file.path(pset.dir, file))
      eset <- switch(data.type,
                     gene=gene.eset,
                     bins=out.eset,
                     tads=out.eset)
      return(eset)
    })
    
    pset2 <- Reduce(function(x,y) combine(x,y), psets)
    rm(psets)
    
    pset <- combine(pset1, pset2)
    saveRDS(pset, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".RDS")))
    rm(pset1, pset2)
  }
}
# 
# sample.names <- lapply(psets, sampleNames)
# features <- Reduce(function(x,y) combine(x,y), lapply(psets, featureData))
# x <- combine(psets[[1]], psets[[1]])
# phenoData
# 
# .appendToPset <- function(path.to, old.pset.f, new.pset, overwrite=T, verbose=T){
#   .load_obj <- function(f){
#     env <- new.env()
#     nm <- load(f, env)[1]
#     env[[nm]]
#   }
#   old.pset <- .load_obj(file.path(path.to, old.pset.f))
#   updated.pset <- old.pset
#   
#   #sampleNames(new.pset)[3] <- sampleNames(old.pset)[5]
#   dup.ids <- sampleNames(old.pset) %in% sampleNames(new.pset)
#   
#   ## Updating featureData obj
#   if(!all(rownames(fData(old.pset)) == rownames(fData(new.pset)))){
#     if(verbose) print("Updating featureData...")
#     featureData(updated.pset) <- combine(featureData(old.pset), featureData(new.pset))
#   } else {
#     if(verbose) print("Using existing featureData...")
#     featureData(updated.pset) <- featureData(old.pset)
#   }
#   
#   ## Updating phenoData obj
#   if(any(dup.ids) & overwrite){
#     ov.ids <- sampleNames(old.pset)[which(dup.ids)]
#     if(verbose) print(paste0("Removing duplicated sample: ", paste(ov.ids, collapse=",")))
#     
#     pd.tmp <- phenoData(old.pset)[-which(dup.ids),]
#     phenoData(updated.pset) <- combine(pd.tmp, phenoData(new.pset))
#   } else if(any(dup.ids) & !overwrite){
#     stop("Non-overwrite mode is not implemented yet")
#   } else {
#     if(verbose) print("Updating phenoData...")
#     phenoData(updated.pset) <- combine(phenoData(old.pset), phenoData(new.pset))
#   }
#   
#   ## Updating assayData obj
#   feature.df <- data.frame('features'=rownames(featureData(updated.pset)))
#   if(verbose) print("Updating assayData...")
#   env.mat <- ls(assayData(old.pset))
#   new.env <- lapply(env.mat, function(em, dup.ids){
#     # Cycle through each assayData matrix and update to match
#     # the order of featureData and phenoData
#     if(any(dup.ids)){
#       old.mat <- as.data.frame(assayData(old.pset)[[em]][,-which(dup.ids)])
#     } else {
#       old.mat <- as.data.frame(assayData(old.pset)[[em]])
#     }
#     
#     new.mat <- as.data.frame(assayData(new.pset)[[em]])
#     updated.mat <- as.data.frame(assayData(updated.pset)[[em]])
#     
#     old.mat$features <- rownames(featureData(old.pset))
#     new.mat$features <- rownames(featureData(new.pset))
#     
#     upd.mat <- Reduce(function(x,y) merge(x,y,by.x='features', all.x=T), 
#                       list(feature.df, old.mat, new.mat))
#     rownames(upd.mat) <- upd.mat$features
#     if(all(colnames(upd.mat)[-1] == rownames(pData(updated.pset)))){
#       as.matrix(upd.mat[,-1])
#     } else {
#       stop("Colnames in assayData does not match the phenoData")
#     }
#   }, dup.ids=dup.ids)
#   names(new.env) <- env.mat
#   assayData(updated.pset) <- list2env(new.env)
#   
#   return(updated.pset)
# }
