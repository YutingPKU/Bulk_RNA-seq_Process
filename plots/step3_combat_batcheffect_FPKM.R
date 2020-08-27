## remove batch effect before clustering analysis
mat <- fread("results/allsamples-RNASeq-cqnnormalized-FPKM-log2-allgenes.txt", header = T)
rownames(mat) <- mat$id
mat$id <- NULL

se <- readRDS("results/allsamples.summarizeOverlaps.rds")
vec <- unlist(as.character(colData(se)[,1]))
sec <- lapply(vec[26:43], FUN = function(chr){
  ls <- unlist(strsplit(chr, "-"))[2]
})
sec <- unlist(sec)
vec <- c(vec[1:25], sec, vec[44:47])
sampleData = cbind(id = vec, 
                   type = c(substr(vec[1:25], 1, nchar(vec[1:25])-1), vec[26:43], "CD34","CD34CD38-","CD34","CD34CD38-"),
                   batch = c(rep(2,3),rep(4,4), 2, rep(4,3),rep(2,8), rep(4,3), rep(2,3), rep(3,18), rep(1,4)), 
                   replicate = c(substr(vec[1:25],nchar(vec[1:25]), nchar(vec[1:25])), rep(1,18),1,1,2,2))
                   

sampleData = data.frame(sampleData)
sampleData$type = as.factor(sampleData$type)
sampleData$replicate = as.factor(sampleData$replicate)
sampleData$batch = as.factor(sampleData$batch)

colData(se) = DataFrame(sampleData)

save(se, file = "results/allsamples.summarizeOverlaps.withcoldata.rds")
load("results/allsamples.summarizeOverlaps.withcoldat.RData")


pheno = data.frame(colData(se)[,1:3])
rownames(pheno) = pheno$id
batch = as.numeric(as.character(pheno$batch))
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=as.matrix(mat), batch=batch, mod=modcombat)
rownames(combat_edata) = rownames(mat)

write.table(combat_edata, "results/allsamples-RNASeq-cqnnormalized-combat-removebatch-FPKM-log2-rmd-allgenes.txt", 
            row.names = T, col.names = T, quote = F, sep = "\t")
