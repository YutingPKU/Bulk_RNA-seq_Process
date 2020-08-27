################  normalize RNA-Seq raw data: GC content, gene length, sequencing depth
########## read count table and gene models

gclen <- fread("public/UCSC_GRCh37_gene_GC_lengths.txt", header = F)
names(gclen) <- c("id", "len", "gc")

load("results/allsamples.summarizeOverlaps.withcoldat.addH1.RData")
mat <- assay(all)
colnames(mat) <- unlist(as.character(colData(all)[,1]))
rownames(mat) <- unlist(as.character(rowData(all)[,1]))


#mat.filter <- mat[which(rowSums(mat) != 0), ]
loci <- match(rownames(mat), gclen$id)
mat.filter <- mat[which(!is.na(loci)),]
gclen.filter <- gclen[loci[which(!is.na(loci))], ]

cqn.subset = cqn(mat.filter, lengths = gclen.filter$len, x = gclen.filter$gc)
RPKM.cqn <- cqn.subset$y + cqn.subset$offset
colnames(RPKM.cqn) <- unlist(as.character(colData(all)[,1]))

#####filter by 
boxplot(RPKM.cqn)
write.table(RPKM.cqn, "results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.txt", sep = '\t', row.names = T, col.names = T, quote = F)
write.xlsx(RPKM.cqn, "results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.xlsx", rowNames =T, colNames = T )
test = read.xlsx("results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.xlsx", rowNames = T, colNames = T)
