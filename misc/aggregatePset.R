library(Biobase)
dataset <- 'GNE'
pdir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines'
pset.dir <- file.path(pdir, dataset, "eacon", "out", "PSet")

###################
#### FUNCTIONS ####
combinePsets <- function(psets){
  ## Uses the BiocGenerics::combine function to merge psets
  while(length(psets) > 1){
    print(length(psets))
    pset1 <- lapply(seq(1, length(psets), by=2), function(pidx){
      if((pidx + 1) > length(psets)){
        psets[[pidx]]
      } else {
        combine(psets[[pidx]], psets[[pidx+1]])
      }
    })
    psets <- pset1
    rm(pset1); gc()
  }
  return(psets)
}

roundAssayData <- function(pset, sigdig=2){
  assay.d <- ls(assayData(pset))
  for(a in assay.d){
    assayData(pset)[[a]] <- round(assayData(pset)[[a]], sigdig)
  }
  return(pset)
}

combinePsetsLomem <- function(all.ps){
  assay.d <- ls(assayData(all.ps[[1]]))
  while(length(all.ps) > 1){
    print(paste0(length(all.ps), "..."))
    for(a in assay.d){
      assayData(all.ps[[1]])[[a]] <- cbind(assayData(all.ps[[1]])[[a]],
                                           assayData(all.ps[[2]])[[a]])
    }
    phenoData(all.ps[[1]])@data <- rbind(phenoData(all.ps[[1]])@data,
                                         phenoData(all.ps[[2]])@data)
    all.ps <- all.ps[-2]
  }
  return(all.ps[[1]])
}

###############################
#### Combine Bins and TADs ####
for(data.type in c('bins', 'tads')){
  print(data.type)
  all.files <- list.files(pset.dir, pattern=paste0(data.type, ".*[0-9].*RData"))
  idx  <- seq(1, length(all.files), by=5)
  se.df <- data.frame("start"=idx,
                      "end"=c(idx[-1]-1, length(all.files)))
  
  all.ps <- apply(se.df, 1, function(i, round.psets=FALSE){
    psets <- lapply(all.files[i['start']:i['end']], function(file){
      load(file.path(pset.dir, file))
      eset <- switch(data.type,
                     gene=gene.eset,
                     bins=out.eset,
                     tads=out.eset)
      return(eset)
    })
    psets <- combinePsetsLomem(psets)
    if(round.psets) psets <- roundAssayData(psets)
    return(psets)
  })
  pset <- combinePsetsLomem(all.ps) 

  rm(all.ps)
  save(pset, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".rda")))
  saveRDS(pset, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".RDS")))
}

#######################
#### Combine Genes ####
for(data.type in c('gene')){
  print(data.type)
  all.files <- list.files(pset.dir, pattern=paste0(data.type, ".*[0-9].*RData"))
  idx  <- seq(1, length(all.files), by=5)
  se.df <- data.frame("start"=idx,
                      "end"=c(idx[-1]-1, length(all.files)))
  
  all.ps <- apply(se.df, 1, function(idx, round.psets=FALSE){
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
  save(all.ps, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".rda")))
  gc()
  
  max.comb.size <- 5
  psets <- list(); i <- 1
  while(length(all.ps) > 0){
    if(length(all.ps) < max.comb.size) max.comb.size <- length(all.ps)
    psets[[i]] <- Reduce(function(x,y) combine(x,y), all.ps[1:max.comb.size])
    all.ps <- all.ps[-c(1:max.comb.size)]
    i <- i+1
    gc()
  }
  rm(all.ps)
  pset <- Reduce(function(x,y) combine(x,y), psets)
  rm(psets)
  
  saveRDS(pset, file=file.path(pset.dir, paste0(dataset, "_CN.", data.type, ".RDS")))
}




