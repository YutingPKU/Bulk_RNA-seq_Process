
load("results/allsamples.backup.withD190128.withD190202.withD190303.withD190430.withD190724.summarizeOverlaps.exonbygene.withcolD.RData")
#se <- se[,-33]
#new <- readRDS("results/D180718.summarizeOverlaps.exonbygene.RData")
new <- readRDS("results/D190820.summarizeOverlaps.exonbygene.rds")

indir = "/lustre/user/liclab/liuyt/denglab/star/bamfiles/RNA-Seq_190820/"
filenames = list.files(indir, pattern = "bam",recursive = TRUE, full.names = T)
vec <- lapply(filenames, FUN=function(chr){
  list <- unlist(strsplit(basename(chr), "_"))[1]
})
vec <- unlist(vec)
sampleid = vec
#sampleid[1] <- "jiajia-1901"

sampleData <- cbind(id = sampleid, batch = rep(21,6))
sampleData <- DataFrame(sampleData)
colData(new) <- sampleData
colnames(new) <- sampleid

se <- all
se.count <- assay(se)
new.count <- assay(new)[-1,]
all.count <- cbind(se.count, new.count)

all.colD <- rbind(colData(se), colData(new))
colnames(all.count) <- rownames(all.colD)

all <- SummarizedExperiment(assays = all.count, colData = all.colD)

save(all,file =  "results/allsamples.backup.withD190128.withD190202.withD190303.withD190430.withD190724.withD190820.summarizeOverlaps.exonbygene.withcolD.RData")
