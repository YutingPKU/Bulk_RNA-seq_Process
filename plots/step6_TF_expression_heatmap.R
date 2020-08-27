#### cluster analysis
library(openxlsx)
library(ggrepel)
mat <- read.xlsx("results/shareData/D181225_genes.fpkm_table.xlsx", rowNames = T,
                 colNames = T, sheet = 1)

smat <- mat[,c(48:53,59:64,66,81:83,67:68,74,89:90)]
#smat <- mat[,c(5:7,9:11,20:22,44,46)]

##### cluster based on TF
tff <- c("TAL","LYL1", "GATA1", "GATA2", "LMO2", "CMYB", "RUNX1", "SPL1", "GFL1","GFL1B",
         "EVL1","ETV2","ERG","FLI1","SOX17","SOX18","IKZF1","HOXA2","HOXA3","HOXA5","HOXA6",
         "HOXA7","HOXA9","HOXA10")
tf <- paste0("^", tff)
loci.tf <- grep(paste(tf,collapse="|"), rownames(smat))
rownames(mat)[loci.tf]
loci.tf <- loci.tf[-c(2:4,22:23,28)]
rownames(mat)[loci.tf]

tf.mat <- smat[loci.tf, ]

################# based on marker
marker <- c("KDR","CD34","CD31","CDH5","THY1","CXCR4","CD45","CD117","MPL","TEK","APLNR","PROCR","ACE")
marker <- paste0("^", marker)
loci.m <- grep(paste(marker,collapse="|"), rownames(smat))
rownames(smat)[loci.m]
loci.m <- loci.m[-c(3:5,12,15:20)]
rownames(smat)[loci.m]

mk.mat <- smat[loci.m,]
########## plot heatmap
pdf("results/plots/TF.surfacemarker.expression.heatmap.D181225.pdf", width = 10, height = 8)
pheatmap(tf.mat, main = "TF")
#pheatmap(tf.mat, cluster_cols = F, main = "TF")
pheatmap(mk.mat, main = "marker")
#pheatmap(mk.mat, cluster_cols = F, main = "marker")
dev.off()

########## save fpkm table
ws <- createWorkbook("matrix")

re <- pheatmap(tf.mat, cluster_cols = F, main = "TF")
addWorksheet(ws, "TF")
writeData(ws, sheet = "TF",tf.mat[re$tree_row$order,] , rowNames = T, colNames = T)

re <- pheatmap(mk.mat, cluster_cols = F, main = "marker")
addWorksheet(ws, "marker")
writeData(ws, sheet = "marker",mk.mat[re$tree_row$order,], rowNames = T, colNames = T )


saveWorkbook(ws, "results/select.heatmap.rawFPKM.D0525.matrix.xlsx", overwrite = T)

