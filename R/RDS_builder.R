Build.SNP6 <- function(ao.df, meta.b, affy.meta, CEL, sex.chr, samplename, cs){
  tmsg("Building normalized object ...")
  # my.ch <- sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) })
  my.ascat.obj <- list(
    data = list(
      Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, row.names = ao.df$ProbeSetName),
      Tumor_LogR = data.frame(sample = ao.df$L2R, row.names = ao.df$ProbeSetName),
      Tumor_BAF = data.frame(sample = ao.df$BAF, row.names = ao.df$ProbeSetName),
      Tumor_AD = data.frame(sample = ao.df[["Allele Difference"]], row.names = ao.df$ProbeSetName),
      Tumor_LogR_segmented = NULL,
      Tumor_BAF_segmented = NULL,
      Germline_LogR = NULL,
      Germline_BAF = NULL,
      SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, row.names = ao.df$ProbeSetName),
      ch = sapply(unique(ao.df$chr), function(x) { which(ao.df$chr == x) }),
      chr = sapply(unique(ao.df$chrgap), function(x) { which(ao.df$chrgap == x) }),
      chrs = unique(ao.df$chr),
      samples = samplename,
      gender = as.vector(meta.b$predicted.gender),
      sexchromosomes = sex.chr,
      failedarrays = NULL
    ), 
    germline = list(
      germlinegenotypes = matrix(as.logical(abs(ao.df$cluster2 - 2L)), ncol = 1),
      failedarrays = NULL
    )
  )
  
  colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
  my.ascat.obj$data$Tumor_BAF.unscaled = data.frame(sample = ao.df$BAF.unscaled, row.names = ao.df$ProbeSetName)
  colnames(my.ascat.obj$data$Tumor_BAF.unscaled) <- samplename
  my.ascat.obj$data$Tumor_BAF.unisomy = data.frame(sample = ao.df$uni, row.names = ao.df$ProbeSetName)
  colnames(my.ascat.obj$data$Tumor_BAF.unisomy) <- samplename
  rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
  genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
  rm(ao.df)
  gc()
  
  ## Adding meta
  my.ascat.obj$meta = list(
    basic = meta.b,
    affy = affy.meta
  )
  
  ## Adding CEL intensities
  my.ascat.obj$CEL = list(
    CEL1 = affxparser::readCel(filename = CEL)
  )
  my.ascat.obj$CEL$CEL1$intensities <- as.integer(my.ascat.obj$CEL$CEL1$intensities)
  
  return(my.ascat.obj)
}