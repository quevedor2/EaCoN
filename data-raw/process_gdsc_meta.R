# Files downloaded from:
# https://www.cancerrxgene.org/downloads/bulk_download
# All cell lines screened 	Cell-line-annotation

setwd('/mnt/work1/users/bhklab/Projects/cell_line_clonality/total_cel_list/datasets/reference/GDSC')
gdsc.df <- read.table("GDSC_Cell_Lines_Details.csv", header=T,
                      sep=",", stringsAsFactors = F, check.names = F,
                      fill=F, quote = '"')

dataset <- 'GDSC'
pdir <- file.path('/mnt/work1/users/pughlab/projects', dataset)
CEL.dir <- file.path(pdir, "data")
sample.paths <- list.files(CEL.dir, pattern="CEL$", recursive = T, 
                           ignore.case = T, full.names = T)
sample.ids <- gsub(".cel", "", basename(sample.paths), ignore.case = TRUE)
regm <- regexpr(".cel", basename(sample.paths), ignore.case=T)
cel.suffix <- regmatches(x = basename(sample.paths), m = regm)
cel.files <- paste0(sample.ids, cel.suffix)

anno.df <- data.frame("sample"=sample.ids, "cel"=cel.files)
anno.df$sample <- gsub("ega-box-03_", "", anno.df$sample)
gdsc.meta <- merge(gdsc.df, anno.df, by.x='Sample Name', by.y='sample', all.y=T)

## Clean up TCGA labels
unable.idx <- gdsc.meta$`Cancer Type (matching TCGA label)` == 'UNABLE TO CLASSIFY'
blank.idx <- gdsc.meta$`Cancer Type (matching TCGA label)` == ''
na.idx <- which(as.logical(unable.idx + blank.idx))
gdsc.meta$`Cancer Type (matching TCGA label)`[na.idx] <- NA


save(gdsc.meta, file="GDSC_meta.rda")
write.table(gdsc.meta, file="GDSC_meta.tsv", sep="\t",
            col.names=T, row.names = F, quote = F)
