#!/usr/bin/env Rscript
# adapted from https://github.com/davidaknowles/leafcutter/blob/master/leafviz/prepare_results.R
library(optparse)

option_parser=OptionParser(
    usage="%prog [options] <name>_cluster_significance.txt <name>_effect_sizes.txt .",
    option_list=list(
        make_option( c("-o","--output"), default="per_intron_results.txt", help="The output file that will be created ready for isoform modelling by GeneStructureTools [%default]")
    )
 )

parsed_args <- parse_args(option_parser,  positional_arguments = 2)

cluster_significance_file <- parsed_args$args[1]
effect.sizes.file <- parsed_args$args[2]
results_file = parsed_args$options$output

library(data.table)
library(stringr)

effectSizes <- fread(effect.sizes.file, data.table=FALSE)
effectSizesSplit <-  as.data.frame(str_split_fixed(effectSizes$intron, ":", 4), stringsAsFactors = FALSE)
names(effectSizesSplit) <- c("chr","start","end","clusterID")

effectSizes <- cbind( effectSizes, effectSizesSplit)
effectSizes$cluster <- paste(effectSizesSplit$chr, effectSizesSplit$clusterID, sep = ":")

results <- fread(cluster_significance_file, data.table=FALSE)
results$FDR <- p.adjust( results$p, method = "fdr")

# Gather introns meeting the FDR threshold
all.introns <- merge(x = results, y = effectSizes, by = "cluster")

if( nrow(all.introns) == 0 ){
    stop("Merging the per-cluster results with the per-junction effect sizes produces an empty table. Please check your input files.")
}

all.introns <- all.introns[ order(all.introns$FDR),]

all.introns$start <- as.numeric(all.introns$start)
all.introns$end <- as.numeric(all.introns$end)
all.introns$cluster = NULL
all.introns$intron = NULL
all.introns = all.introns[,c(which(colnames(all.introns) == "clusterID"), which(!colnames(all.introns) == "clusterID"))]

write.table(all.introns, file = results_file, quote=FALSE, row.names = FALSE, sep="\t")

