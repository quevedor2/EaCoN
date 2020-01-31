##################
#### gne_meta ####
setwd("~/git/EaCoN/data-raw/")
gne.meta <- read.table('GNE_mapping.tsv', sep=",", check.names=F, 
                       header=T,  stringsAsFactors = F)
gne.meta$TCGA_code <- NA
usethis::use_data(gne.meta, overwrite = T)
