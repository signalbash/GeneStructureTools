#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
load(args[1])
m = match(introns$clusterID, clusters$clusterID)
introns$FDR = clusters$FDR[m]
write.table(introns, file=args[2], quote=FALSE, row.names=FALSE, sep="\t")

