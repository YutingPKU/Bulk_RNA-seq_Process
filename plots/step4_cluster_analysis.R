#### cluster analysis
library(openxlsx)
library(ggrepel)
library(dplyr)
#mat <- read.xlsx("results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.xlsx", rowNames = T, colNames = T)
#mat <- mat[,26:43]
#se = readRDS("results/allsamples.summarizeOverlaps.exonbygene.RData")
#all = se[,26:43]
#sampleData <- colData(all)
#add <- data.frame(cbind(id = "H1", type = "H1", batch = 5, replicate = 1))
#sampleData <- rbind(sampleData, add)

mat <- read.xlsx("results/xlsx/D180913_fpkm.table.xlsx", rowNames = T, colNames = T)
rmat <- log2(mat+1) 

#mat.filter <- mat[which(rowSums(mat) >= 0),]
## PCA
cutvar <- function(mat, cutvalue){
  gvar <- apply(mat, 1, var)
  gcut <- quantile(gvar, cutvalue)
  loci <- which(gvar >= gcut)
  mcut <- mat[loci,]
  mcut <- as.matrix(mcut)
  colnames(mcut) = colnames(mat)
  return(mcut)
}
mat.filter <- cutvar(rmat, 0.9)
res.pca <- prcomp(t(mat.filter))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)

sampleData = data.frame(colData(all))
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = colnames(mat))) +
  geom_point(size =6) + geom_label_repel(aes(label = colnames(mat)),
                                         box.padding   = 0.35, 
                                         point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0(" PCA (", nrow(mat.filter), "genes) "))+
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

pdf("results/CBbatch-RNA-seq-PCA-withname.pdf", width = 16, height = 10)
g1
dev.off()

## hcluster
cp = cor(mat.filter, method = "pearson")
c = cor(mat.filter, method = "spearman")
#annote = as.data.frame(cbind("region" = c("CP","GZ","CP","GZ"), "individual" = c("07456A","07456A","07456B","07456B")))
annote = as.data.frame(colData(all)[,2:4])
annote$type = as.factor(annote$type)
annote$replicate = as.factor(annote$replicate)
rownames(annote) <- colnames(cp)
pdf("results/CBbatch-RNA-Seq-pearson-correlation-cluster.pdf", width = 16, height = 12)
re = pheatmap(cp, main = "pearson")
res = pheatmap(c, main = "spearman")
#pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = 28, cluster_cols = T, annotation_col  = annote,
#         colorRampPalette(c("white","red"))(100), show_colnames = F)

dev.off()

ws <- createWorkbook("matrix")

addWorksheet(ws, "pearson")
writeData(ws, sheet = "pearson",cp[re$tree_row$order, re$tree_col$order] , rowNames = T, colNames = T)

addWorksheet(ws, "spearman")
writeData(ws, sheet = "spearman",c[res$tree_row$order, res$tree_col$order], rowNames = T, colNames = T )

saveWorkbook(ws, "results/CBbatch-RNASeq-FPKM-cluster.xlsx", overwrite = T)



pheatmap(as.matrix(cp[1:10,1:10]),fontsize_row = 20,cluster_rows = F, cluster_cols = F,show_colnames = T,
         display_numbers = T,number_format = "%2f", fontsize_number = 10,number_color = "black",cellwidth=50,cellheight=50,treeheight_row = 250,
         border_color = "white",fontsize_col = 20,color = colorRampPalette(c(rgb(72,111,167, max = 255), "white",rgb(219,110,49, max = 255)))(50))



## marker gene expression
markg = c()
loci = match(markg, mcols(rowRanges(all))[,2])
markgid <- as.character(mcols(rowRanges(all)[loci,])$id)
loci <- match(noquote(markgid), rownames(mat))
gmat = mat[loci,]
rownames(gmat) = markg
anno <- as.data.frame(colData(all)[, c("region","replicate")])
rownames(anno) <- colnames(gmat)
pheatmap(gmat, annotation_col = anno, cluster_rows = FALSE, fontsize = 18, main = " CP VS GZ")
gmat1 <- gmat[,1:6]
gmat1 <- gmat1 - rowMeans(gmat1)
pheatmap(gmat1[,1:6], annotation_col = anno, cluster_rows = FALSE, fontsize = 18, main = "11002B CP VS GZ")
