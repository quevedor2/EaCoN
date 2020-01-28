#' Build.SNP6
#'
#' @param ao.df 
#' @param meta.b 
#' @param affy.meta 
#' @param CEL 
#' @param sex.chr 
#' @param samplename 
#' @param cs 
#'
#' @return
#' @export
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

#' Build.OMNI25
#'
#' @param illumina.dir 
#' @param parent.dir 
#'
#' @return
#' @export
#'
#' @examples Build.OMNI25(illumina.dir='/mnt/work1/users/pughlab/projects/cancer_cell_lines/GNE/illumina/GNE_Matrices',
#' eacon.dir='/mnt/work1/users/pughlab/projects/cancer_cell_lines/GNE')
Build.OMNI25 <- function(illumina.dir, parent.dir){
  #gne.mat.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/GNE/illumina/GNE_Matrices'
  #illumina.dir <- '/mnt/work1/users/pughlab/projects/cancer_cell_lines/GNE/illumina/GNE_Matrices'
  require(ASCAT)
  require(dplyr)
  
  # Log all Log2R and BAF files
  logr.files <- list.files(file.path(illumina.dir, 'LogR'), pattern="LogR")
  baf.files <- list.files(file.path(illumina.dir, 'BAF'), pattern="BAF")
  ids <- gsub("^.*_", "", logr.files) %>% gsub(".txt", "", .)
  
  # Split each sample into an individaul RDS
  lapply(ids, function(i){
    # Find sample in each baf/log2r matrix
    f.l2r.idx <- grep(paste0("_", i, ".txt"), logr.files)
    f.baf.idx <- grep(paste0("_", i, ".txt"), baf.files)
    
    ascat.obj.path <- file.path(illumina.dir, "ASCAT_obj", paste0("ascat_", i, ".rda"))
    if(!file.exists(ascat.obj.path)){
      tmsg("Creating ASCAT object from LogR and BAF...")
      ascat.bc = ASCAT::ascat.loadData(Tumor_LogR_file=file.path(illumina.dir, "LogR", logr.files[f.l2r.idx]),
                                       Tumor_BAF_file=file.path(illumina.dir, "BAF", baf.files[f.baf.idx]))
      gg <- ASCAT::ascat.predictGermlineGenotypes(ascat.bc, platform = "Illumina2.5M")
      dir.create(file.path(illumina.dir, "ASCAT_obj"), showWarnings = FALSE)
      save(ascat.bc, gg, file=ascat.obj.path)
    } else {
      tmsg("Reading in LogR/BAF ASCAT objects")
      load(ascat.obj.path)
    }
    
    samples <- gsub("\\.Log R.*", "", colnames(ascat.bc$Tumor_LogR))
    colnames(ascat.bc$Tumor_LogR) <- colnames(ascat.bc$Tumor_BAF) <- ascat.bc$samples <-samples
    colnames(gg$germlinegenotypes) <- samples
    
    ## Create a simple meta info sheet
    meta <- list(
      "basic"=list(
        samplename=samples,
        source='microarray',
        source.file='',
        type='Illumina2.5M',
        manufacturer='Infinium',
        species="Homo sapiens",
        genome='hg19',
        genome.pkg="BSgenome.Hsapiens.UCSC.hg19",
        predicted.gender='XY',
        wave.renorm='',
        gc.renorm='GC50,GC100'),
      "affy"=NULL
    )
    
    ## Loop through each sample and split into the ao.df data structure
    tmp <- lapply(samples, function(s){
      data.tmp <- ascat.bc
      gg.tmp <- gg
      meta.tmp <- meta
      
      data.tmp$Tumor_LogR = ascat.bc$Tumor_LogR[,s,drop=FALSE]
      data.tmp$Tumor_BAF = ascat.bc$Tumor_BAF[,s,drop=FALSE]
      data.tmp$samples = s
      data.tmp$gender = ascat.bc$gender[match(s, samples)]
      gg.tmp$germlinegenotypes = gg.tmp$germlinegenotypes[,s,drop=FALSE]
      meta.tmp$basic$samplename = s
      
      rds <- list("data"=data.tmp,
                  "meta"=meta.tmp,
                  "germline"=gg.tmp)
      sample.dir <- file.path(parent.dir, "eacon", s)
      dir.create(sample.dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(rds, file=file.path(sample.dir, paste0(s, "_processed.RDS")))
    })
  })
}

# segmenter = 'ASCAT'
# RDS.files = file.path(sample.dir, paste0(s, "_processed.RDS"))
# EaCoN:::Segment.ff.Batch(RDS.file = RDS.files,  segmenter = segmenter, nthread=5)
# 
# l2r.rds <- file.path(sample.dir, 'ASCAT', 'L2R', paste0(s, '.SEG.ASCAT.RDS'))
# file.path(sample.dir, paste0(s, "_processed.RDS"))
# fit.val <- EaCoN:::ASCN.ff.Batch(RDS.files = l2r.rds, nthread=2)