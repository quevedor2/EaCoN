#' Pre-processed TCGA Pancancer Atlas data
#'
#' Metadata and matrices of TCGA samples that have ploidy/purity calls
#' on them as well, split by cancer type.
#' 
#' "The Cancer Genome Atlas (TCGA) Research Network has profiled and 
#' analyzed large numbers of human tumors to discover molecular aberrations 
#' at the DNA, RNA, protein and epigenetic levels. The resulting rich 
#' data provide a major opportunity to develop an integrated picture 
#' of commonalities, differences and emergent themes across tumor 
#' lineages. The Pan-Cancer Atlas initiative compares the 33 tumor 
#' types profiled by TCGA. Analysis of the molecular aberrations and 
#' their functional roles across tumor types will teach us how to extend 
#' therapies effective in one cancer type to others with a similar 
#' genomic profile."
#'
#' @docType data
#'
#' @usage data(pancan.obj.segless)
#'
#' @format A list with three elements:
#' \describe{
#'  \item{AP-meta}{A dataframe containing metadata of all annotated samples with ploidy}
#'  \item{APS-meta}{A dataframe containing metadata of all annotated samples with ploidy and seg files}
#'  \item{breaks}{A list of 4 matrices containing binned frequencies for: 'purity', 'ploidy', 'cancerDNA fraction', 'subclonal fraction'}
#' }
#' 
#' @keywords datasets
#'
#' @source \href{https://gdc.cancer.gov/about-data/publications/pancanatlas}{PanCancer Atlas}
#' ABSOLUTE-annotated seg file - TCGA_mastercalls.abs_segtabs.fixed.txt 
#' Copy Number - broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg
#' ABSOLUTE purity/ploidy file - TCGA_mastercalls.abs_tables_JSedit.fixed.txt 
#' Merged Sample Quality Annotations - merged_sample_quality_annotations.tsv
#' 
#' @examples
#' data(pancanPloidy.noSegs)
#' pancan.ploidy <- pancan.obj.segless$breaks$ploidy
#' pancan.ploidy[,c('breaks', 'ACC')]
"pancan.obj.segless"


#' CCLE Metadata
#'
#' @docType data
#'
#' @usage data(CCLE_meta)
#'
#' @format A dataframe containing CCLE metadata:
#' 
#' @keywords datasets
#'
#' @source \href{https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt}{CCLE metadata}
#' CCLE_sample_info_file_2012-10-18.txt
#' Cell_lines_annotations_20181226.txt
#' 
#' @examples
#' data(CCLE_meta)
"ccle.meta"