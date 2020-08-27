####### select DEG in HSPC VS HPC, and show heatmap of expression in experiment groups

library("dplyr")
load("results/allsamples.withD180718.summarizeOverlaps.exonbygene.withcolD.RData")
#samples <- c("ABFCD34D22","B10CD34D22","BCDFCD34D22","s-D18CW","s-W100-1","CWTII45","CW37NN45","CWFTP45",
#             "CW45","TCW45","CWTB45")
samples <- c("ABFCD34D22","s-D18CW","T1","T6","T3","T8")
loci <- grep(paste(samples,collapse="|"), 
             colData(all)[,1], value=TRUE)
all = se[, loci]
colnames(all) = colData(all)[,1]
colData(all) <- cbind(colData(all), condition = c(rep("ABFCD34D22",3),"D18CW","T1","T3",
                                                  "T6","T8"))
colData(all)[,3] <- as.factor(colData(all)[,3])
## construct dds object
dds <- DESeqDataSet(all, design = ~ id)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]
nrow(dds)


## normalization
head(counts(dds),3)
vsd <- vst(dds, blind = T)
mat <- assay(vsd)
## DEG 
deg <- read.xlsx("results/xlsx/HPC_VS_HSPC/DEG_T3T8_VS_T1T6_UPINHSPC_FC2.xlsx", colNames = T, rowNames = T)
loci.deg <- match(rownames(deg), rownames(mat))
mat.deg <- mat[loci.deg,]
pheatmap(mat.deg[,c(1:3,7)], show_rownames = T)

test <- data.frame(mat.deg[,c(1:3,7)])
test$FC <- log2(test$s.D18CW/rowMeans(test[,1:3]))
test <- test[which(abs(test$FC) > 0.1),]
test <- test[order(test$FC),]
pheatmap(test[,-5], show_rownames = T, cluster_rows = T)






###### FPKM
exp <- read.xlsx("results/xlsx/quantile.normalized.FPKM.D180718.xlsx", colNames = T, rowNames = T)
deg <- read.xlsx("results/xlsx/HPC_VS_HSPC/DEG_T3T8_VS_T1T6_UPINHSPC_FC2.xlsx", colNames = T, rowNames = T)
#deg <- read.xlsx("results/tmp/select_DEG_T3T8.xlsx", colNames = T, rowNames = T)
loci.deg <- match(rownames(deg), rownames(exp))
samples <- c("ABFCD34D22","s-D18CW","T1","T8","T3","T6")
samples <- c("ABFCD34D22","s-D18CW","s-W100-1",
                         "T1","T3","T6","T8")
loci.sample <- grep(paste(samples,collapse="|"), colnames(exp), value=F)

mat <- exp[loci.deg[which(!is.na(loci.deg))], loci.sample]
mat <- mat[,c(4,6,5,7,8,1:3)]
#mat <- log2(mat+1)
mat$fct <- log2(rowMeans(mat[,3:4])/rowMeans(mat[,1:2]))
mat$fc <- log2(mat$`s-D18CW_0`/rowMeans(mat[,6:8]))
#mat <- mat[which(abs(mat$fc) >1),]
mat <- mat[order(  mat$fc),]
submat <- mat[which(abs(mat$fc)>1),1:8]
#pdf("results/plots/HPC_VS_HSPC/T1T6_VS_T3T8_DEG_FC2_withexperiments_heatmap_withgenename.pdf", width = 12, height = 28)
re <- pheatmap(log2(submat[,1:8]+1), cluster_rows = T, cluster_cols = F, show_rownames = T, show_colnames = F)
 pheatmap(log2(submat[match(retest, rownames(submat)),c(5:8,4,1:3)]+1), cluster_rows = F, cluster_cols = F, show_rownames = T)
#dev.off()


sd <-  apply(mat, 1, sd)
#mat <- mat[order(sd),]
sub.mat <- mat[which(sd>5),]
pheatmap(log2(sub.mat+1), show_rownames = F, cluster_cols = T)



