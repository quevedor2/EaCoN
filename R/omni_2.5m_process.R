## Performs CS CEL processing
Omni2.5m.Process <- function(CEL = NULL, samplename = NULL, 
                         l2r.level = "normal", method="rawcopy",
                         gc.renorm = TRUE, gc.rda = NULL, wave.renorm = TRUE, 
                         wave.rda = NULL, mingap = 1E+06, 
                         out.dir = getwd(), oschp.keep = FALSE, 
                         force.OS = NULL, apt.version = "1.20.0", 
                         apt.build = "na35.r1", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", 
                         return.data = FALSE, write.data = TRUE, 
                         plot = TRUE, force = FALSE) {
  
  # setwd("/home/job/WORKSPACE/EaCoN_tests/SNP6")
  # CEL <- "GSM820994.CEL.bz2"
  # samplename <- "BITES_TEST"
  # l2r.level <- "normal"
  # wave.renorm <- TRUE
  # wave.rda <- NULL
  # gc.renorm <- TRUE
  # gc.rda <- NULL
  # mingap <- 1E+06
  # out.dir <- getwd()
  # oschp.keep <- TRUE
  # force.OS <- NULL
  # apt.version <- "1.20.0"
  # apt.build <- "na35.r1"
  # genome.pkg <- "BSgenome.Hsapiens.UCSC.hg19"
  # return.data <- FALSE
  # write.data <- TRUE
  # plot <- TRUE
  # force <- FALSE
  # require(foreach)
  # source("~/git_gustaveroussy/EaCoN/R/mini_functions.R")
  # source("~/git_gustaveroussy/EaCoN/R/renorm_functions.R")
  
  
  ## Early checks
  if (is.null(CEL)) stop(tmsg("A CEL file is required !"), call. = FALSE)
  if (is.null(samplename)) stop(tmsg("A samplename is required !"), call. = FALSE)
  if (!file.exists(CEL)) stop(tmsg(paste0("Could not find CEL file ", CEL, " !")), call. = FALSE)
  if (gc.renorm) { if (!is.null(gc.rda)) { if (!file.exists(gc.rda)) stop(tmsg(paste0("Could not find gc.rda file ", gc.rda)), call. = FALSE) } }
  if (wave.renorm) { if (!is.null(wave.rda)) { if (!file.exists(wave.rda)) stop(tmsg(paste0("Could not find wave.rda file ", wave.rda)), call. = FALSE) } }
  if (is.null(genome.pkg)) stop(tmsg("A BSgenome package name is required !"), call. = FALSE)
  if (!genome.pkg %in% BSgenome::installed.genomes()) {
    if (genome.pkg %in% BSgenome::available.genomes()) {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), call. = FALSE)
    } else {
      stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), call. = FALSE)
    }
  }
  if (dir.exists(samplename)) { if (!force) stop(tmsg(paste0("A [", samplename, '] dir already exists !')), call. = FALSE) else unlink(samplename, recursive = TRUE, force = FALSE) }
  
  l2r.lev.conv <- list("normal" = "Log2Ratio", "weighted" = "SmoothSignal")
  if (!(l2r.level %in% names(l2r.lev.conv))) stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"), call. = FALSE)
  
  ## Handling compressed files
  CEL <- compressed_handler(CEL)
  
}
  