library("dplyr")
library(ggrepel)
load("results/allsamples.withD180718.withD180913.withD181010.withD181017.withD181112.withD181206.summarizeOverlaps.exonbygene.withcolD.RData")
#loci = which(colData(se)[,3] == 3 | colData(se)[,3] == 4)
se = all
#all = se[, c(60:63)]
colnames(all) = colData(all)[,1]
colData(all) <- cbind(colData(all), condition = c("HPC","HSPC","HPC","HSPC"))
colData(all)[,3] <- as.factor(colData(all)[,3])
## construct dds object
dds <- DESeqDataSet(all, design = ~ id)

## pre filtering
nrow(dds)
keep <- rowSums(counts(dds)) >= 80
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

pdf("results/plots/D182106_vsdnormalize_cluster_distance.pdf", width = 16, height = 12)

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
pcaData <- plotPCA(vsd, intgroup = c( "condition.X"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition.X)) +
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
loci = seq(1,28)
mat.filter <- cutvar(assay(vsd),0.9)
res.pca <- prcomp(t(mat.filter))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)
g1 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = vsd$id)) +
  geom_point(size =6) + geom_label_repel(aes(label = colData(all)[,1]),
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

pdf("results/plots/D182106_vsdnormalize_PCA.pdf", width = 12, height = 12)
g1
dev.off()


dds <- DESeq(dds)
res <- results(dds, contrast = c("condition.X", "HSPC","HPC"), lfcThreshold = 0)
plotMA(res)
res.dat <- data.frame(res)
res.up <- res.dat[which(res.dat$log2FoldChange>=2 & res.dat$padj <= 0.05 ),]
res.up <- res.up[order( res.up$log2FoldChange, -log10(res.up$pvalue), decreasing = T),]
res.down <- res.dat[which(res.dat$log2FoldChange <= -2 & res.dat$padj <= 0.05),]
res.down <- res.down[order( abs(res.down$log2FoldChange), -log10(res.down$pvalue), decreasing = T),]

mat <- assay(vsd)
loci.up <- match(rownames(res.up), rownames(mat))
loci.down <- match(rownames(res.down), rownames(mat))
pdf("results/plots/HPC_VS_HSPC/T1T6_VS_T3T8_DEG_FC2_heatmap.pdf", width = 10, height = 8)
pheatmap(mat[c(loci.up, loci.down),], show_rownames = F, cluster_rows = F)
dev.off()
write.xlsx(res.up, file = "results/xlsx/DEG_T3T8_VS_T1T6_UPINHSPC_FC1.xlsx", colName = T, rowNames = T)                 

####### volcano plot
gene_list = res.dat[complete.cases(res.dat),]
#gene_list = gene_list[order(abs(gene_list$log2FoldChange), -log10(gene_list$pvalue), decreasing = T),]
gene_list$threshold = as.factor(abs(gene_list$log2FoldChange) > 1 & gene_list$padj <= 0.05)
gene_list$gene = rownames(gene_list)
loci.up = match(rownames(res.up), gene_list$gene)[1:20]
loci.down = match(rownames(res.down), gene_list$gene)[1:20]
g = ggplot(data=gene_list, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) +
  geom_point(size=1) +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
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
        plot.title = element_text(size=28, hjust = 0.5))+
  scale_colour_manual(values = c("FALSE"= "black", "TRUE"="red"))+
  geom_text_repel(data=gene_list[c(loci.up, loci.down),], aes(x=log2FoldChange, y=-log10(pvalue),label=gene))

pdf("results/DE-gene-volcano-CB-vs-BVSB-BVFSB-D18CW.pdf", width = 12, height = 10)
g
dev.off()


ws <- createWorkbook("matrix")

addWorksheet(ws, "gene-up-inCB")
writeData(ws, sheet = "gene-up-inCB",res.up , rowNames = T, colNames = T)

addWorksheet(ws, "gene-up-inBVSB-BVFSB-D18CW")
writeData(ws, sheet = "gene-up-inBVSB-BVFSB-D18CW",res.down, rowNames = T, colNames = T )
addWorksheet(ws, "allgenes-DE")
writeData(ws, sheet = "allgenes-DE",res.dat, rowNames = T, colNames = T )
addWorksheet(ws, "allgenes-rawcount")
writeData(ws, sheet = "allgenes-rawcount",counts(dds), rowNames = T, colNames = T )

saveWorkbook(ws, "results/DE-gene-table-CB-vs-BVSB-BVFSB-D18CW.xlsx", overwrite = T)


############ DE gene name
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene.up <- getBM(attributes=c('ensembl_gene_id'),
      filters=c('hgnc_symbol'),
      values = rownames(res.up),
      mart=ensembl)
gene.down <- getBM(attributes=c('ensembl_gene_id'),
                   filters=c('hgnc_symbol'),
                   values = rownames(res.down),
                   mart=ensembl)
write.table(gene.up, "results/DE-gene-name-upinCB-CB-vs-BVSB-BVFSB-D18CW.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(gene.down,  "results/DE-gene-name-downinCB-CB-vs-BVSB-BVFSB-D18CW.txt", row.names = F, col.names = F, sep = "\t", quote = F)


############### compare with public
public.up <- read.xlsx("public/nbt.3702-S3.xlsx", startRow = 2, sheet = 1)
comm.gene.up <- intersect(public.up$ILMN_Gene, rownames(res.up))
public.down <- read.xlsx("public/nbt.3702-S3.xlsx", startRow = 2, sheet = 2)
comm.gene.down <- intersect(public.down$ILMN_Gene, rownames(res.down))

loci.up <- match(comm.gene.up, gene_list$gene)
loci.down <- match(comm.gene.down, gene_list$gene)
g = ggplot(data=gene_list, aes(x=log2FoldChange, y=-log10(pvalue), colour=threshold)) +
  geom_point(size=1) +
  theme(legend.position = "none") +
  xlim(c(-15, 15)) + ylim(c(0, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+ggtitle("comm DE with NBT results")+
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
        plot.title = element_text(size=28, hjust = 0.5))+
  scale_colour_manual(values = c("FALSE"= "black", "TRUE"="red"))+
  geom_text_repel(data=gene_list[c(loci.up, loci.down),], aes(x=log2FoldChange, y=-log10(pvalue),label=gene))

pdf("results/DE-gene-volcano-CB-vs-BVSB-BVFSB-D18CW-comparewithNBT.pdf", width = 12, height = 10)
g
dev.off()
