######## reads counting step for all samples

library(GenomicAlignments)
library(Rsamtools)
#library(GenomicFeatures)
#library(refGenome)
#library(BiocParallel)
#library(DESeq2)
library(data.table)
library(BiocParallel)

## step 1 definning gene models
#print("definning gene models")
#setwd("/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/")
#txdb <- makeTxDbFromGFF("genes.gtf", format="gtf")
#eByg <- exonsBy(txdb, by=c("gene"))
#ghs <- genes(txdb, columns = c("TXNAME","GENEID"))
#save(ghs, file = "/lustre/user/liclab/liuyt/denglab/backup/public/UCSC.hg19.genes.GRanges.RData")
load("/lustre/user/liclab/liuyt/denglab/backup/public/UCSC.hg19.genes.exonbygene.GRanges.RData")
#save(eByg, file = "/lustre/user/liclab/liuyt/denglab/backup/public/UCSC.hg19.genes.exonbygene.GRanges.RData")

## step 2 locating alinged files
#indir = "/lustre/user/liclab/liuyt/denglab/star/bamfiles/RNA-Seq_180718/"
#indir = "/lustre/user/liclab/liuyt/denglab/star/bamfiles/RNA-Seq_181017/"
indir = "/lustre/user/liclab/liuyt/denglab/star/bamfiles/RNA-Seq_200104/"
filenames = list.files(indir, pattern = "bam",recursive = TRUE, full.names = T)
file.exists(filenames)
bamfiles = BamFileList(filenames)
seqinfo(bamfiles[1])


## step 3 reads counting
registered()
register(MulticoreParam(workers = 6,RNGseed = 77394650, timeout = 144000000, log = F  ))
registered()
print("start counting reads:")
se <- summarizeOverlaps(features=eByg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
vec <- lapply(filenames, FUN=function(chr){
  list <- unlist(strsplit(basename(chr), "_"))[1]
})
vec <- unlist(vec)
sampleid = vec
#sampleid = c()
sampleData = cbind(id= sampleid)
                   #type = substr(vec,1, nchar(vec)-1),
                   #replicate = substr(vec,nchar(vec), nchar(vec)),
                   #batch = )
sampleData = data.frame(sampleData)
#sampleData$type = as.factor(sampleData$type)
#sampleData$replicate = as.factor(sampleData$replicate)

#colData(se) = DataFrame(sampleData)

#table = assay(se)
#table = data.frame(table)
#rownames(table) = mcols(my_gr)[,1]
#colnames(table) = colnames(se)
#write.table(table,"/lustre/user/liclab/liuyt/monkey-brain/human-brain/RNA-Seq/results/expression-mat.txt", sep = '\t', row.names = T, col.names = T, quote = F)
saveRDS(se, "/lustre/user/liclab/liuyt/denglab/star/results/D200104.summarizeOverlaps.exonbygene.rds")
print("Counting step is done!")
