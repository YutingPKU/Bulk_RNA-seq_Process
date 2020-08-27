############## GO for genes in conserve loops and hgs loops
gob2a = read.table("results/tmp/HPC_VS_HSPC/chart_GO_DEG_T3T8_VS_T1T6_UPINHSPC_genename_FC1.txt", sep = "\t", header = T)
gob2a = gob2a[which(gob2a$Category != "Category"),]
gob2a$PValue = as.numeric(as.character(gob2a$PValue))
gob2a$Benjamini = as.numeric(as.character(gob2a$Benjamini))
gob2a = gob2a[order(gob2a$Benjamini, decreasing = F),]

gob2a.bp <- gob2a[which(gob2a$Category=="GOTERM_BP_DIRECT"),]
gob2a.bp = gob2a.bp[1:20,]
gob2a.bp$Term <- substr(gob2a.bp$Term, 12, 40)

pdf("results/plots/HPC_VS_HSPC/GO-term-BP_T3T8_VS_T1T6_UPINHSPC_genename_FC1.pdf", width = 18, height = 12)
par(mar = c(8,60,4,12))
barplot(-log10(gob2a.bp$Benjamini), main = "DEG up in HSPC BP ", horiz = T, names.arg = gob2a.bp$Term, las = 1, col = "royalblue4", 
        cex.names = 3.0, cex.axis = 3.0, xlab = "-log10(P.adjust value)", cex.main = 3, axes = F, cex.lab = 3)
axis(side = 1, at = seq(0,15,5), labels = seq(0,15,5), tick = T, lwd = 4, cex.axis = 2)

segments(x0 = -log10(0.05), y0 = 0, x1 = -log10(0.05), y1 = 25, col = 'red', lwd = 4)
dev.off()