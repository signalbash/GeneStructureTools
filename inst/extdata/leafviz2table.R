#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
load(args[1])
write.table(introns, file=args[2], quote=FALSE, row.names=FALSE, sep="\t")

