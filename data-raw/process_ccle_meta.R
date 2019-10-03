# CCLE data: https://portals.broadinstitute.org/ccle/data
# CCLE metadata: https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt

ccle.meta.files <- c('CCLE_sample_info_file_2012-10-18.txt', 'Cell_lines_annotations_20181226.txt')
ccle.meta <- lapply(ccle.meta.files, function(f){
  read.table(f, sep="\t", header=T, stringsAsFactors = F, check.names = F,
             fill=F, comment.char = '')
})
ccle.meta <- Reduce(function(x,y) merge(x,y,by.x='CCLE name', by.y='CCLE_ID', all=TRUE),
                    ccle.meta)
ids <- sapply(ccle.meta$`CCLE name`, function(x) strsplit(x, split="_")[[1]][1])
tissues <- sapply(ccle.meta$`CCLE name`, function(x) paste(strsplit(x, split="_")[[1]][-1], collapse="_"))

ccle.meta$Name <- ccle.meta$`Cell line primary name` <- ids
ccle.meta$site <- tolower(tissues)


write.table(ccle.meta, file="CCLE_meta.tsv", sep="\t",
            col.names=T, row.names = F, quote = F)
save(ccle.meta, file="CCLE_meta.RData")