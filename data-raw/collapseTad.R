library(GenomicRanges)
grl <- as(lapply(list.files(pattern="rda$"), function(i){
  load(i)
  c.gr
}), 'GRangesList')
names(grl) <- gsub("_TAD.*", "", list.files(pattern="rda$"))

gr <- unlist(grl)
seqlevelsStyle(gr) <- 'NCBI'
seqlevels(gr) <- c(1:22, "X", "Y")
gr <- sort(gr)
seqlevelsStyle(gr) <- 'UCSC'

seg <- as.data.frame(gr)
seg$ID <- rep("Consensus", nrow(seg))
seg <- seg[,c(7, 1, 2, 3, 4, 6)]
colnames(seg) <- c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
colour.line <- seg[1,]
colour.line[,c('loc.start', 'loc.end', 'seg.mean')] <- c(1, 2, -0.1)
seg <- rbind(colour.line, seg)

saveRDS(gr, file = "consensusTAD.RDS")
write.table(gr, file="consensusTAD.tsv", sep="\t",
            col.names = T, row.names = F, quote = F)
write.table(seg, file="consensusTAD.seg", sep="\t",
            col.names = T, row.names = F, quote = F)