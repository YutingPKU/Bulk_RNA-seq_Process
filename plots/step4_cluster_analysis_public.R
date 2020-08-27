#### cluster analysis
library(openxlsx)
library(ggrepel)
library(dplyr)

mat <- read.xlsx("public/GSE74246_RNAseq/GSE74246_RNAseq_All_Counts.xlsx", rowNames = T, colNames = T)
mat <- mat[,1:49]
id <- colnames(mat)
ind <- substr(id, 1,4)
cellt <- substr(id,6, 20)
colD <- data.frame(cbind(id = id, individual = ind, celltype = cellt))
de <- DESeqDataSetFromMatrix(countData = mat, colData = colD, design = ~celltype)
saveRDS(de, file = "public/GSE74246_RNAseq/GSE74246_RNAseq_All_Counts_onlyhealth.RDS")




se = readRDS("public/GSE74246_RNAseq/GSE74246_RNAseq_All_Counts_onlyhealth.RDS")
dds = se
sampleData <- colData(dds)

####### vst 
## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds)


## normalization
head(counts(dds),3)
vsd <- vst(dds, blind = T)

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
mat <- assay(vsd)
mat.filter <- cutvar(mat, 0.9)
res.pca <- prcomp(t(mat.filter))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)

sampleData = data.frame(colData(all))
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$celltype)) +
  geom_point(size =6) + geom_label_repel(aes(label = sampleData$id),
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
cp = cor(mat, method = "pearson")
c = cor(mat.filter, method = "spearman")
#annote = as.data.frame(cbind("region" = c("CP","GZ","CP","GZ"), "individual" = c("07456A","07456A","07456B","07456B")))
annote = data.frame(ID = factor(sampleData$celltype))
rownames(annote) <- colnames(cp)
re = pheatmap(cp, main = "pearson")
res = pheatmap(c, main = "spearman")
aka3 = list(ID = c(HSC =rgb(0,68,27,maxColorValue = 255) , MPP= rgb(70,160,64,maxColorValue = 255),
                   LMPP = rgb(1,175,153,maxColorValue = 255), CMP =rgb(255,193,161,maxColorValue = 255),
                   GMP = rgb(255,163,0,maxColorValue = 255),MEP = rgb(246,49,62,maxColorValue = 255),
                   Mono = rgb(211,28,68,maxColorValue = 255), CD4Tcell = rgb(1,129,201,maxColorValue = 255),
                   CD8Tcell = rgb(0,21,136,maxColorValue = 255), NKcell = rgb(73,12,101,maxColorValue = 255),
                   Bcell = rgb(186,127,208,maxColorValue = 255), CLP=rgb(152,217,233,maxColorValue = 255),
                   Ery = rgb(143,19,54,maxColorValue = 255)))
pheatmap(as.matrix(cp),  cluster_rows = T, fontsize = 18, cluster_cols = T, annotation_row  = annote,
         annotation_colors =aka3[1] ,
         colorRamps::blue2green2red(100), show_colnames = F, show_rownames = F, border_color = NA)

dev.off()




############# combine with our data
public <- readRDS("public/GSE74246_RNAseq/GSE74246_RNAseq_All_Counts_onlyhealth.RDS")
load("results/allsamples.withD180718.summarizeOverlaps.exonbygene.withcolD.RData")
our <- all
rm(all)

comm.gene <- intersect(rownames(our), rownames(public))
comm.our <- assay(our)[match(comm.gene, rownames(our)),]
comm.pub <- counts(public)[match(comm.gene, rownames(public)),]
comm.all <- cbind(comm.our, comm.pub)

colD <- data.frame(cbind(id = c(as.character(colData(our)[,1]), colData(public)[,1]), 
                         source = c(rep("our", 70), rep("public",49)), 
                         type = c(as.character(colData(our)[,1]), as.character(colData(public)[,3]))))
de <- DESeqDataSetFromMatrix(countData = comm.all, colData = colD, design = ~type)


nrow(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
nrow(dds)


## normalization
head(counts(dds),3)
vsd <- vst(dds, blind = T)

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
mat <- assay(vsd)
mat.filter <- cutvar(mat, 0.9)
res.pca <- prcomp(t(mat.filter))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)

sampleData = data.frame(colData(dds))
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = sampleData$id)) +
  geom_point(size =6) + geom_label_repel(aes(label = sampleData$id),
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


## hcluster
cp = cor(mat, method = "pearson")
c = cor(mat.filter, method = "pearson")


re = pheatmap(cp, main = "pearson")
res = pheatmap(c, main = "pearson on Var genes")

ws <- createWorkbook("matrix")

addWorksheet(ws, "pearson")
writeData(ws, sheet = "pearson",cp[re$tree_row$order, re$tree_col$order] , rowNames = T, colNames = T)

addWorksheet(ws, "pearson on Var genes")
writeData(ws, sheet = "pearson on Var genes",c[res$tree_row$order, res$tree_col$order], rowNames = T, colNames = T )

saveWorkbook(ws, "public/GSE74246_RNAseq/compare_withOurData_cluster_correlation.xlsx", overwrite = T)



