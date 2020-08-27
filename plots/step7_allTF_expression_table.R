library("openxlsx")
#mat <- read.xlsx("results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.xlsx", 
             #    sheet = 1, colNames = T, rowNames = T)
mat <- read.xlsx("results/shareData/D190115_gene.fpkm_table.xlsx", rowNames = T, colNames = T)
#mat <- mat[,26:43]
tf <- fread("public/TF.name.txt", header = F)
loci <- match(tf$V1, rownames(mat))
loci <- loci[which(!is.na(loci))]
tfmat <- mat[loci,]
#tfmat <- 2^tfmat
write.xlsx(tfmat, "results/shareData/D190115.TF.gene.fpkm_table.xlsx", rowNames = T, colNames = T)

