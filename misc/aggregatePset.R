library(Biobase)
dataset <- 'CCLE'
pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
pset.dir <- file.path(pdir, dataset, "eacon", "out", "PSet")
for(data.type in c('bins', 'tads')){
  print(data.type)
  all.files <- list.files(pset.dir, pattern=paste0(data.type, ".*RData"))
  splits <- seq(1, length(all.files), by=10)
  splits.df <- data.frame("start"=splits,
                          "end"=c(splits[-1]-1, length(all.files)))
  
  all.psets <- apply(splits.df, 1, function(idx){
    print(paste(idx, collapse="-"))
    
    psets <- lapply(all.files[idx['start']:idx['end']], function(file){
      load(file.path(pset.dir, file))
      eset <- switch(data.type,
                     gene=gene.eset,
                     bins=out.eset,
                     tads=out.eset)
      return(eset)
    })
    pset1 <- Reduce(function(x,y) combine(x,y), psets)
    rm(psets)
    return(pset1)
  })
  save(all.psets, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".rda")))
  gc()
  
  max.comb.size <- 5
  psets <- list(); i <- 1
  while(length(all.psets) > 0){
    if(length(all.psets) < max.comb.size) max.comb.size <- length(all.psets)
    psets[[i]] <- Reduce(function(x,y) combine(x,y), all.psets[1:max.comb.size])
    all.psets <- all.psets[-c(1:max.comb.size)]
    i <- i+1
    gc()
  }
  rm(all.psets)
  pset <- Reduce(function(x,y) combine(x,y), psets)
  rm(psets)
  
  saveRDS(pset, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".RDS")))
}