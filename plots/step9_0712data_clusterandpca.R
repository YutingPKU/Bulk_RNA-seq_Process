
omat <- read.xlsx("results/shareData/D190303_genes.fpkm_table.xslx", rowNames = T, colNames = T)
omat <- omat[which(rowSums(omat) > 10), ]
ymat <- read.xlsx("../tmpdata/genes.fpkm_table_180630.xlsx", rowNames = T, colNames = T)
colnames(omat)[46:63] <- substr(colnames(omat)[46:63],3, 30)
romat <- log2(omat+1)

######### compare with backup data
comm <- intersect(colnames(omat), colnames(ymat))
loci.o <- match(comm, colnames(omat))
loci.y <- match(comm, colnames(ymat))
comm.omat <- omat[, loci.o]
comm.ymat <- ymat[, loci.y]
test <- apply(comm.ymat, 2, as.numeric)

cor(comm.omat$T3, comm.ymat$T3)
cp.o <- cor(comm.omat, method = "spearman")
cp.y <- cor(test)




########### pca and cluster  by tophat fpkm
cp <- cor(romat, method = "pearson")
pdf("results/plots/allsamples_log2FPKM_spearman_correlation_D190303.pdf", width = 18, height = 14)
re <- pheatmap(cp)
dev.off()

ws <- createWorkbook("matrix")

addWorksheet(ws, "spearman correlation")
writeData(ws, sheet = "spearman correlation",cp[re$tree_row$order, re$tree_col$order] , rowNames = T, colNames = T)

saveWorkbook(ws, "results/xlsx/allsamples_FPKM_spearman_correlation_180712.xlsx", overwrite = T)



############## pac and cluster by vsd and rld

library("dplyr")
load("results/allsamples.backup.withD190128.withD190202.withD190303.summarizeOverlaps.exonbygene.withcolD.RData")
#loci = which(colData(se)[,3] == 3 | colData(se)[,3] == 4)
#loci <- match(c("s86T445jia-1","s86T445jia-2", "Thymus","T15","T22","CB34-7jia-1","CB34-7jia-2","H1-38-1","H1-38-2","BWF1","BWF2","BWF3"), rownames(colData(all)))
all = all[, loci]

## construct dds object
dds <- DESeqDataSet(all, design = ~ id)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 90
dds <- dds[keep,]
nrow(dds)


## normalization
head(counts(dds),3)
vsd <- vst(dds, blind = T)
head(assay(vsd),3)
rld <- rlog(dds, blind = T)
head(assay(rld),3)
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


## cluster and PCA
sampleDists.rld <- dist(t(assay(rld)))
sampleDistMatrix.rld <- as.matrix( sampleDists.rld )
rownames(sampleDistMatrix.rld) <- vsd$id
colnames(sampleDistMatrix.rld) <- NULL

sampleDists.vsd <- dist(t(assay(vsd)))
sampleDistMatrix.vsd <- as.matrix( sampleDists.vsd)
rownames(sampleDistMatrix.vsd) <- vsd$id
colnames(sampleDistMatrix.vsd) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("results/plots/allsamples_withD190303_VSDandRLDNormalized_cluster_distance.pdf", width = 16, height = 12)
re.vsd <- pheatmap(sampleDistMatrix.vsd,
                   clustering_distance_rows = sampleDists.vsd,
                   clustering_distance_cols = sampleDists.vsd,
                   col = colors, main = "vsd")
#re.rld <- pheatmap(sampleDistMatrix.rld,
#                   clustering_distance_rows = sampleDists.rld,
#                   clustering_distance_cols = sampleDists.rld,
#                   col = colors, main = "rld")
dev.off()








plotPCA(vsd, intgroup = c("id"))
pcaData <- plotPCA(vsd, intgroup = c( "type"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

cutvar <- function(mat, cutvalue){
  gvar <- apply(mat, 1, var)
  gcut <- quantile(gvar, cutvalue)
  loci <- which(gvar >= gcut)
  mcut <- mat[loci,]
  mcut <- as.matrix(mcut)
  colnames(mcut) = colnames(mat)
  return(mcut)
}
#mat <- cutvar(assay(vsd))
#loci = which(colData(all)[,3] == 3)
mat <- cutvar(assay(vsd),0.9)
res.pca <- prcomp(t(mat))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = colData(all)$id)) +
  geom_point(size =6) + geom_label_repel(aes(label = colData(all)$id),
                                         box.padding   = 0.35, 
                                         point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0(" PCA vsd (", nrow(mat), "genes) "))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),  
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28),
        legend.position = "none",
        plot.title = element_text(size=28, hjust = 0.5))

pdf("results/plots/allsamples_withD190303_vst_top10VarGenes_PCA.pdf", width = 16, height = 12)
g1
dev.off()

cs <- cor(mat, method = "pearson")
re <- pheatmap(cs, show_colnames = F)
pdf("results/plots/allsamples_withD190303_pearson_correlation.pdf", width = 16, height = 12)
re
dev.off()
write.xlsx(cs[re$tree_row$order, re$tree_col$order], file = "results/xlsx/allsample_withD180718_vsd_allgenes_spearman_correlation.xlsx", row.names = T, col.names = T)


loci <- match(c("H1-38-1","H1-1112-2","H1-1112-1","H1-37","H1139I","D2CW-1","D2C10","D2CW-2",
                "D2C5","CW37NN","T1","T3","T8","T6","s-CBG","s-CBZ1","s-CBZ2","Thymus","CB34-7jia-1","CB34-7jia-2"),
              colnames(mat))
boxplot(mat[338,loci], mat[338, -loci])
rat <- res.pca$rotation
rat <- rat[,1:2]
write.xlsx(rat, "results/shareData/D181225_genefpkm_PCAloading.xlsx", rowName = T, colName = T)
