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
loci = seq(1,47)
#rmat = read.xlsx("results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.xlsx", rowNames = T, colNames = T)
mat <- cutvar(assay(vsd)[,loci],0.95)
#mat <- cutvar(rmat, 0.95)
res.pca <- prcomp(t(mat))
eig <- (res.pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.decathlon2.active <- data.frame(eig = eig, variance = variance,
                                    cumvariance = cumvar)
head(eig.decathlon2.active)
gvsd.095.18 <- ggplot(data.frame(res.pca$x[,1:2]), aes(x = PC1, y = PC2, color = rld$type[loci])) +
  geom_point(size =6) + geom_label_repel(aes(label = colData(se)[loci,1]),
                                         box.padding   = 0.35, 
                                         point.padding = 0.5)+
  xlab(paste0("PC1 (", round(eig.decathlon2.active[1,2]), "% variance)")) +
  ylab(paste0("PC2 (", round(eig.decathlon2.active[2,2]), "% variance)")) +
  ggtitle(paste0(" PCA log2 (", nrow(mat), "genes) "))+
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

pdf("results/plots/PCA_nocombat_batch1234.pdf", width = 16, height = 12)
gvsd.0.47
gvsd.09.47
gvsd.095.47
#gvsd.0.18
#gvsd.09.18
#gvsd.095.18
dev.off()
