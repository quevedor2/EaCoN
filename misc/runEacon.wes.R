#### Libraries ####
## Instal ASCAT and FACETS
# devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
# devtools::install_github("mskcc/facets")
# biocLite(c("affxparser", "Biostrings", "aroma.light", "BSgenome", "copynumber", "GenomicRanges", "limma", "rhdf5", "sequenza"))
# devtools::install_github("gustaveroussy/EaCoN")
# devtools::install_github("quevedor2/EaCoN")

# devtools::install_github("gustaveroussy/apt.snp6.1.20.0")
# install.packages("https://partage.gustaveroussy.fr/pydio_public/083305?dl=true&file=/affy.CN.norm.data_0.1.2.tar.gz", repos = NULL, type = "source")
# install.packages("https://partage.gustaveroussy.fr/pydio_public/152397?dl=true&file=/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz", repos = NULL, type = "source")
# install.packages( "https://partage.gustaveroussy.fr/pydio_public/e6fe22?dl=true&file=/rcnorm_0.1.5.tar.gz", repos = NULL, type = "source")

##############
#### Main ####
#devtools::install_github("quevedor2/EaCoN", ref = 'tads')
library(optparse)
library(EaCoN)
library(parallel)
library(Biobase)

capture.ref.dir <- '/mnt/work1/users/pughlab/references/intervals'
capture <- 'V5utr'
genome <- 'hg19'
capture.bed <- switch(capture,
                      V5utr=file.path(capture.ref.dir, '/SureSelect_Human_All_Exon_V5+UTRs/S04380219_Covered.bed'),
                      V7=file.path(capture.ref.dir, '/SureSelect_Human_All_Exon_V7/hg38', 'S31285117_Covered.bed'))
genome.build <- switch(genome,
                       hg19='BSgenome.Hsapiens.UCSC.hg19',
                       hg38='BSgenome.Hsapiens.UCSC.hg38')

BINpack.Maker(bed.file = capture.bed, bin.size = 50, 
              genome.pkg = "BSgenome.Hsapiens.UCSC.hg38")

bam.dir <- '/mnt/work1/users/pughlab/data/INSPIRE_data_mordor/160418_D00331_0188_AC8AU8ANXX_cindy/processed_bam'
tumor <- 'INS-E-001-S-T.processed.bam'
normal <- 'INS-E-001-SB.processed.bam'
bin.rda <- '/mnt/work1/users/home2/quever/S04380219_Covered_hg19_b50.GC.rda'
WES.Bin(testBAM = file.path(bam.dir, tumor), refBAM = file.path(bam.dir, normal), 
        BINpack = bin.rda, samplename = 'INS-E-001')