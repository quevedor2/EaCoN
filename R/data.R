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

#' GDSC Metadata
#'
#' @docType data
#'
#' @usage data(GDSC_meta)
#'
#' @format A dataframe containing GDSC metadata:
#' 
#' @keywords datasets
#'
#' @source \href{https://www.cancerrxgene.org/downloads/bulk_download}{GDSC metadata}
#' All cell lines screened 	Cell-line-annotation
#' Cell_Lines_Details.xlsx
#' GDSC_Cell_Lines_Details.csv
#' 
#' @examples
#' data(GDSC_meta)
"gdsc.meta"

#' Consensus TAD GenomicRanges
#'
#' @docType data
#'
#' @usage data(consensusTAD)
#'
#' @format A GenomicRanges object containing condensed TADs and their jaccard index cutoff
#' 
#' @keywords datasets
#'
#' @source \href{https://www.encodeproject.org/metadata/type%3DExperiment%26files.output_type%3Dnested%2Btopologically%2Bassociated%2Bdomains%26files.output_type%3Dtopologically%2Bassociated%2Bdomains%26files.file_type%3Dbed%2Bbed3%252B%26files.assembly%3Dhg19/metadata.tsv}{ENCODE data}
#' https://www.encodeproject.org/metadata/type%3DExperiment%26files.output_type%3Dnested%2Btopologically%2Bassociated%2Bdomains%26files.output_type%3Dtopologically%2Bassociated%2Bdomains%26files.file_type%3Dbed%2Bbed3%252B%26files.assembly%3Dhg19/metadata.tsv -X GET -H "Accept: text/tsv" -H "Content-Type: application/json" --data '{"elements": ["/experiments/ENCSR549MGQ/","/experiments/ENCSR312KHQ/","/experiments/ENCSR440CTR/","/experiments/ENCSR834DXR/","/experiments/ENCSR401TBQ/","/experiments/ENCSR862OGI/","/experiments/ENCSR489OCU/","/experiments/ENCSR998ZSP/","/experiments/ENCSR346DCU/","/experiments/ENCSR444WCZ/","/experiments/ENCSR213DHH/","/experiments/ENCSR343VKT/"]}'
#' https://www.encodeproject.org/files/ENCFF066NYU/@@download/ENCFF066NYU.bed.gz
#' https://www.encodeproject.org/files/ENCFF784LMI/@@download/ENCFF784LMI.bed.gz
#' https://www.encodeproject.org/files/ENCFF075QYD/@@download/ENCFF075QYD.bed.gz
#' https://www.encodeproject.org/files/ENCFF437EBV/@@download/ENCFF437EBV.bed.gz
#' https://www.encodeproject.org/files/ENCFF996RJN/@@download/ENCFF996RJN.bed.gz
#' https://www.encodeproject.org/files/ENCFF701HCM/@@download/ENCFF701HCM.bed.gz
#' https://www.encodeproject.org/files/ENCFF136AUY/@@download/ENCFF136AUY.bed.gz
#' https://www.encodeproject.org/files/ENCFF310FEU/@@download/ENCFF310FEU.bed.gz
#' https://www.encodeproject.org/files/ENCFF588KUZ/@@download/ENCFF588KUZ.bed.gz
#' https://www.encodeproject.org/files/ENCFF274VJU/@@download/ENCFF274VJU.bed.gz
#' https://www.encodeproject.org/files/ENCFF471EYL/@@download/ENCFF471EYL.bed.gz
#' https://www.encodeproject.org/files/ENCFF768YZZ/@@download/ENCFF768YZZ.bed.gz
#' https://www.encodeproject.org/files/ENCFF931RKD/@@download/ENCFF931RKD.bed.gz
#' https://www.encodeproject.org/files/ENCFF531PAV/@@download/ENCFF531PAV.bed.gz
#' https://www.encodeproject.org/files/ENCFF654GDV/@@download/ENCFF654GDV.bed.gz
#' https://www.encodeproject.org/files/ENCFF451MCF/@@download/ENCFF451MCF.bed.gz
#' https://www.encodeproject.org/files/ENCFF938WXQ/@@download/ENCFF938WXQ.bed.gz
#' https://www.encodeproject.org/files/ENCFF676WJO/@@download/ENCFF676WJO.bed.gz
#' https://www.encodeproject.org/files/ENCFF336WPU/@@download/ENCFF336WPU.bed.gz
#' https://www.encodeproject.org/files/ENCFF248DGB/@@download/ENCFF248DGB.bed.gz
#' https://www.encodeproject.org/files/ENCFF477VVJ/@@download/ENCFF477VVJ.bed.gz
#' https://www.encodeproject.org/files/ENCFF569KNS/@@download/ENCFF569KNS.bed.gz
#' Dataset was created using the createReferenceTad.R script
#' And condensed into a signal data structure using the collapsedTad.R script
#' 
#' @examples
#' data(consensusTAD)
"tad.gr"