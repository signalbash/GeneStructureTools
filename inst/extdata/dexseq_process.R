library(stringr)
library(DEXSeq)
options(stringsAsFactors = FALSE)

setwd("~/Downloads/GeneStructureTools_tuts/")

countFiles = list.files(full.names=TRUE, pattern=".dexseq.txt")
flattenedFile = "gencode.vM14.annotation.dexseq.gtf"

countFilesNames=basename(countFilesNames)

sampleTable=data.frame(row.names = countFilesNames,
                       condition=str_sub(countFilesNames, 2,3),
                       replicate=str_sub(countFilesNames, 1,1),
                       libType="paired-end", fileName=countFiles)

dxd = DEXSeqDataSetFromHTSeq(countFiles,
                             sampleTable,
                             design=~ sample + exon + condition:exon,
                             flattenedfile = flattenedFile)
colData(dxd)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )

plotDispEsts( dxd )

dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )

save(dxd, dxr1, sampleTable, file="dexseq_processed.Rdata")

signif_dex = as.data.frame(dxr1)
signif_dex = signif_dex[signif_dex$padj < 1e-12,]
signif_dex = signif_dex[which(abs(signif_dex$log2fold_21_01) > 2.34),]

write.table(signif_dex, file="dexseq_results_significant.txt", sep="\t", quote=FALSE)
