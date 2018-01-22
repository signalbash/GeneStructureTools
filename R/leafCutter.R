#' Add set numbers to introns
#'
#' Converts a group of introns into non-overlapping sets
#' @param clusterGRanges Granges object with a cluster of intron locations
#' @return Granges object with a cluster of intron locations and corresponding set numbers
#' @keywords internal
#' @import GenomicRanges
#' @importFrom plyr desc
#' @author Beth Signal
addSets <- function(clusterGRanges){
    clusterGRanges$set <- 1

    ol <- as.data.frame(findOverlaps(clusterGRanges))
    ol <- ol[ol$queryHits != ol$subjectHits,]
    ol$setFrom <- clusterGRanges$set[ol$queryHits]
    ol$setTo <- clusterGRanges$set[ol$subjectHits]
    ol <- ol[ol$setFrom == ol$setTo,]

    while(nrow(ol) > 0){
        #find coord with most overlaps
        tab <- as.data.frame(table(ol$queryHits))
        tab <- tab[order(plyr::desc(tab$Freq)),]

        clusterGRanges$set[as.numeric(as.character(tab$Var1[1]))] <- max(clusterGRanges$set) + 1

        olSet <- findOverlaps(clusterGRanges[clusterGRanges$set == max(clusterGRanges$set)],
                              clusterGRanges[clusterGRanges$set != max(clusterGRanges$set)])

        line <- clusterGRanges[clusterGRanges$set!= max(clusterGRanges$set)][-olSet@to]
        if(length(line) > 0){
            line$set <- max(clusterGRanges$set)
            clusterGRanges <- c(clusterGRanges, line)
        }
        ol <- as.data.frame(findOverlaps(clusterGRanges))
        ol <- ol[ol$queryHits != ol$subjectHits,]
        ol$setFrom <- clusterGRanges$set[ol$queryHits]
        ol$setTo <- clusterGRanges$set[ol$subjectHits]
        ol <- ol[ol$setFrom == ol$setTo,]
    }
    return(clusterGRanges)
}

#' Remove exon duplicates
#'
#' Removes structural duplicates of exons in a GRanges object
#' @param exons GRanges object with exons
#' @return GRanges object with unique exons
#' @export
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' gtf.exons <- gtf[gtf$type=="exon"]
#' gtf.exons.duplicated <- c(gtf.exons[1:4], gtf.exons[1:4])
#' length(gtf.exons.duplicated)
#' gtf.exons.deduplicated <- removeSameExon(gtf.exons.duplicated)
#' length(gtf.exons.deduplicated)
removeSameExon <- function(exons){
    samesies <- findOverlaps(exons, type = "equal")
    samesies <- samesies[samesies@from > samesies@to]
    if(length(samesies) > 0){
        exons <- exons[-unique(samesies@from)]
    }
    return(exons)
}

#' Create transcripts with alternative intron usage
#'
#' Creates transcript isoforms from alternative intron usage tested by leafcutter
#' @param altIntronLocs data.frame containing information from the
#' per_intron_results.tab file output from leafcutter.
#' Note that only one cluster of alternative introns can be processed at a time.
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @return GRanges object with all potential alternative isoforms skipping the
#' introns specified in either the upregulated or downregulated locations
#' @export
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' leafcutterFiles <- list.files(system.file("extdata","leafcutter/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' leafcutterIntrons <- read.delim(leafcutterFiles[grep("intron_results",
#' leafcutterFiles)],stringsAsFactors=FALSE)
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' gtf.exons <- gtf[gtf$type=="exon"]
#' # single cluster processing
#' cluster <- leafcutterIntrons[leafcutterIntrons$cluster=="chr16:clu_1396",]
#' altIsoforms1396 <- alternativeIntronUsage(cluster, gtf.exons)
#' unique(altIsoforms1396$transcript_id)
#' cluster <- leafcutterIntrons[leafcutterIntrons$cluster=="chr16:clu_1395",]
#' altIsoforms1395 <- alternativeIntronUsage(cluster, gtf.exons)
#' unique(altIsoforms1395$transcript_id)
#' # multiple cluster processing
#' altIsoforms1396plus1395 <- alternativeIntronUsage(cluster, c(gtf.exons, altIsoforms1396))
#' unique(altIsoforms1396plus1395$transcript_id)

alternativeIntronUsage <- function(altIntronLocs, gtf.exons){
    clusterGRanges <- GRanges(seqnames=S4Vectors::Rle(altIntronLocs$chr),
                              ranges=IRanges::IRanges(start=as.numeric(altIntronLocs$start),
                                                      end=as.numeric(altIntronLocs$end)),
                              strand="*",
                              id=altIntronLocs$clusterID,
                              direction=ifelse(altIntronLocs$deltapsi >0, "+","-"))

    m <- match(altIntronLocs$ensemblID, gtf.exons$gene_id)
    strand(clusterGRanges)[which(!is.na(m))] <- strand(gtf.exons)[m][which(!is.na(m))]
    # maximum spanning region
    clusterGRanges.max <- clusterGRanges
    #start(clusterGRanges.max) <- min(start(clusterGRanges.max))
    #end(clusterGRanges.max) <- max(end(clusterGRanges.max))

    #find overlaps -- for when range overlaps multiple genes
    olExons <- as.data.frame(findOverlaps(clusterGRanges.max, gtf.exons))
    gtf.transcripts <- exonsToTranscripts(gtf.exons[gtf.exons$gene_id %in%
                                                        gtf.exons$gene_id[olExons$subjectHits]])
    # find transcripts which contain the cluster region
    olTrans <- as.data.frame(findOverlaps(clusterGRanges.max, gtf.transcripts, type = "within"))
    clusterTranscripts <- gtf.transcripts[unique(olTrans$subjectHits)]
    #transcript_exons <- gtf.exons[gtf.exons$transcript_id %in% clusterTranscripts$transcript_id,]
    clusterExons <- gtf.exons[gtf.exons$transcript_id %in% clusterTranscripts$transcript_id,]


    # add sets to cluster introns
    clusterGRanges.dnre <- addSets(clusterGRanges[clusterGRanges$direction=="-"])
    clusterGRanges.upre <- addSets(clusterGRanges[clusterGRanges$direction=="+"])
    clusterGRanges.upre$set <- clusterGRanges.upre$set + max(clusterGRanges.dnre$set)
    clusterGRanges <- c(clusterGRanges.upre, clusterGRanges.dnre)

    overlaps <- as.data.frame(findOverlaps(clusterGRanges, clusterExons))
    rmExons <- clusterExons[unique(overlaps$subjectHits)]

    clusterGRanges.intron <- clusterGRanges
    start(clusterGRanges.intron) <- start(clusterGRanges.intron) +1
    end(clusterGRanges.intron) <- end(clusterGRanges.intron) -1

    for(i in seq_along(unique(clusterGRanges$set))){

        clusterGRanges.max <- clusterGRanges.intron[clusterGRanges.intron$set==i]
        start(clusterGRanges.max) <- min(start(clusterGRanges.max))
        end(clusterGRanges.max) <- max(end(clusterGRanges.max))

        overlapsCluster <- findOverlaps(clusterGRanges.max, clusterExons)
        rmExons <- clusterExons[unique(overlapsCluster@to)]
        if(length(rmExons) > 0){
            clusterExonsBounding <- clusterExons[-unique(overlapsCluster@to)]
        }else{
            clusterExonsBounding <- clusterExons
        }

        #overlaps start of the intron
        clusterGRanges.start <- clusterGRanges[clusterGRanges$set==i]
        end(clusterGRanges.start) <- start(clusterGRanges.start)
        overlapsStart <- findOverlaps(clusterGRanges.start, rmExons, type = "end")
        #overlaps end of the intron
        clusterGRanges.end <- clusterGRanges[clusterGRanges$set==i]
        start(clusterGRanges.end) <- end(clusterGRanges.end)
        overlapsEnd <- findOverlaps(clusterGRanges.end, rmExons, type = "start")

        exonsStart <- rmExons[overlapsStart@to]
        exonsStart <- removeSameExon(exonsStart)

        exonsEnd <- rmExons[overlapsEnd@to]
        exonsEnd <- removeSameExon(exonsEnd)

        InternalExons <- removeSameExon(c(exonsStart,exonsEnd))

        # make sure exons are within the intronic region
        overlapsIntron <-
            as.data.frame(findOverlaps(InternalExons,
                                       clusterGRanges.max, type="within"))
        InternalExons <- InternalExons[unique(overlapsIntron$queryHits)]

        overlapsExonStart <- findOverlaps(clusterGRanges.start,
                                          clusterExonsBounding, type = "end")
        overlapsExonEnd <- findOverlaps(clusterGRanges.end,
                                        clusterExonsBounding, type = "start")

        keepTranscriptIds <- clusterExonsBounding$transcript_id[overlapsExonStart@to]
        keepTranscriptIds <- keepTranscriptIds[clusterExonsBounding$transcript_id[overlapsExonStart@to] %in%
                                                   clusterExonsBounding$transcript_id[overlapsExonEnd@to]]

        if(any(grepl("_dnre_", keepTranscriptIds) | grepl("_upre_", keepTranscriptIds))){
            # direction to to remove
            # if upre, remove dnre isforms
            removeDirection <- ifelse(clusterGRanges$direction[match(i, clusterGRanges$set)[1]] == "+", "dnre", "upre")
            keepTranscriptIds <- keepTranscriptIds[!grepl(removeDirection, keepTranscriptIds)]
        }

        clusterExonsBounding <- clusterExonsBounding[clusterExonsBounding$transcript_id %in% keepTranscriptIds]

        InternalExons.reps <- InternalExons[rep(seq_along(InternalExons), length(unique(keepTranscriptIds)))]
        InternalExons.reps$transcript_id <- rep(unique(keepTranscriptIds), each=length(InternalExons))

        clusterExonsBounding <- c(clusterExonsBounding, InternalExons.reps)
        clusterExonsBounding <- reorderExonNumbers(clusterExonsBounding)
        clusterExonsBounding$set <- as.numeric(i)

        if(i == 1){
        clusterExons.allSets <- clusterExonsBounding
        }else{
            clusterExons.allSets <- c(clusterExons.allSets, clusterExonsBounding)
        }
    }

    m <- match(clusterExons.allSets$set, clusterGRanges$set)
    #n <- which(!grepl("[+]", clusterExons.allSets$transcript_id))
    clusterExons.allSets$new_transcript_id <- NA
    clusterExons.allSets$new_transcript_id <-
        paste0(clusterExons.allSets$transcript_id, "+AS ",
               ifelse(clusterGRanges$direction[m] == "+", "upre","dnre"),
               "_",
               gsub("_", "", clusterGRanges$id[m]),
               "-", clusterExons.allSets$set)

    n <- which(grepl("[+]", clusterExons.allSets$transcript_id))
    clusterExons.allSets$new_transcript_id[n] <-
        paste0(clusterExons.allSets$transcript_id, ":",
               gsub("_", "", clusterGRanges$id[m]),
               "-", clusterExons.allSets$set)[n]

    clusterExons.allSets$transcript_id <-
        clusterExons.allSets$new_transcript_id
    clusterExons.allSets$new_transcript_id <- NULL
    clusterExons.allSets$set <- NULL

    clusterExons.allSets <-
        removeDuplicateTranscripts(clusterExons.allSets)
    return(clusterExons.allSets)
}

#' Convert an exon-level gtf annotation to a transcript-level gtf annotation
#'
#' @param gtf.exons GRanges object with exons
#' @return GRanges object with transcripts
#' @export
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' gtf.exons <- gtf[gtf$type=="exon" & gtf$transcript_id=="ENSMUST00000126412.1"]
#' gtf.exons
#' gtf.transcripts <- exonsToTranscripts(gtf.exons)
#' gtf.transcripts
exonsToTranscripts <- function(gtf.exons){
    gtf.transcripts <- gtf.exons[!duplicated(gtf.exons$transcript_id)]

    minStarts <- aggregate(start ~ transcript_id,
                           as.data.frame(gtf.exons), min)
    maxEnds <- aggregate(end ~ transcript_id,
                         as.data.frame(gtf.exons), max)


    ranges(gtf.transcripts) <-
        IRanges::IRanges(start=as.numeric(
            minStarts$start[match(gtf.transcripts$transcript_id,
                                  minStarts$transcript_id)]),
            end=as.numeric(
                maxEnds$end[match(gtf.transcripts$transcript_id,
                                  maxEnds$transcript_id)]))


    return(gtf.transcripts)
}
