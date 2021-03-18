#' Create transcripts with alternative intron usage
#'
#' Creates transcript isoforms from alternative intron usage tested by leafcutter
#' @param altIntronLocs data.frame containing information from the
#' per_intron_results.tab file output from leafcutter.
#' Note that only one cluster of alternative introns can be processed at a time.
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param replaceInternalExons replace any internal transcript exons with inferred exons from multi-intron leafcutter intron sets?
#' @param junctions GRanges object from readLeafcutterJunctions()
#' @return GRanges object with all potential alternative isoforms skipping the
#' introns specified in either the upregulated or downregulated locations
#' @export
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
#' @examples
#' leafcutterFiles <- list.files(system.file("extdata", "leaf_small/",
#'     package = "GeneStructureTools"
#' ), full.names = TRUE)
#' leafcutterIntrons <- read.delim(leafcutterFiles[grep(
#'     "intron_results",
#'     leafcutterFiles
#' )], stringsAsFactors = FALSE)
#' gtf <- rtracklayer::import(system.file("extdata", "gencode.vM25.small.gtf",
#'     package = "GeneStructureTools"
#' ))
#' exons <- gtf[gtf$type == "exon"]
#' # single cluster processing
#' cluster <- leafcutterIntrons[leafcutterIntrons$cluster == "chr17:clu_20975_+", ]
#' altIsoforms20975 <- alternativeIntronUsage(cluster, exons)
#' unique(altIsoforms20975$transcript_id)
alternativeIntronUsage <- function(altIntronLocs, exons, replaceInternalExons = TRUE, junctions = NULL) {
    clusterGRanges <-
        GRanges(
            seqnames = S4Vectors::Rle(altIntronLocs$chr),
            ranges = IRanges::IRanges(
                start = as.numeric(altIntronLocs$start),
                end = as.numeric(altIntronLocs$end)
            ),
            strand = altIntronLocs$strand,
            id = altIntronLocs$clusterID,
            direction = ifelse(altIntronLocs$deltapsi > 0, "+", "-"),
            verdict = altIntronLocs$verdict,
            deltapsi = altIntronLocs$deltapsi
        )

    # maximum spanning region
    clusterGRanges.max <- clusterGRanges
    # start(clusterGRanges.max) <- min(start(clusterGRanges.max))
    # end(clusterGRanges.max) <- max(end(clusterGRanges.max))

    # find overlaps -- for when range overlaps multiple genes
    olExons <- as.data.frame(findOverlaps(clusterGRanges.max, exons))
    if (nrow(olExons) > 0) {
        transcripts <-
            exonsToTranscripts(exons[exons$gene_id %in%
                exons$gene_id[olExons$subjectHits]])
        # find transcripts which contain the cluster region
        olTrans <- as.data.frame(findOverlaps(clusterGRanges.max,
            transcripts,
            type = "within"
        ))
        clusterTranscripts <- transcripts[unique(olTrans$subjectHits)]
        # transcript_exons <- exons[exons$transcript_id %in%
        # clusterTranscripts$transcript_id,]
        clusterExons <- exons[exons$transcript_id %in%
            clusterTranscripts$transcript_id, ]


        if (length(clusterExons) > 0) {

            # add sets to cluster introns
            if (any(clusterGRanges$direction == "-")) {
                clusterGRanges.dnre <-
                    addSets(clusterGRanges[clusterGRanges$direction == "-"])
                setlist <- unique(clusterGRanges.dnre$set)
                clusterGRanges.dnre$set <- match(clusterGRanges.dnre$set, setlist)
            }
            if (any(clusterGRanges$direction == "+")) {
                clusterGRanges.upre <-
                    addSets(clusterGRanges[clusterGRanges$direction == "+"])
                setlist <- unique(clusterGRanges.upre$set)
                clusterGRanges.upre$set <- match(clusterGRanges.upre$set, setlist)
            }
            if (any(clusterGRanges$direction == "+") & any(clusterGRanges$direction == "-")) {
                clusterGRanges.upre$set <-
                    clusterGRanges.upre$set + max(clusterGRanges.dnre$set)
                clusterGRanges <- c(clusterGRanges.upre, clusterGRanges.dnre)
            } else if (any(clusterGRanges$direction == "+")) {
                clusterGRanges <- clusterGRanges.upre
            } else if (any(clusterGRanges$direction == "-")) {
                clusterGRanges <- clusterGRanges.dnre
            }

            clusterGRanges$set <- clusterGRanges$set - min(clusterGRanges$set) + 1

            overlaps <- as.data.frame(findOverlaps(clusterGRanges, clusterExons))
            rmExons <- clusterExons[unique(overlaps$subjectHits)]

            clusterGRanges.intron <- clusterGRanges
            start(clusterGRanges.intron) <- start(clusterGRanges.intron) + 1
            end(clusterGRanges.intron) <- end(clusterGRanges.intron) - 1

            for (s in seq_along(unique(clusterGRanges$set))) {
                verdicts <- clusterGRanges$verdict[clusterGRanges$set == s]

                # make a new intermediate exon for funsies
                if (length(which(verdicts %in% c(c("cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored")))) > 1) {
                    clusterGRanges.newExon <- clusterGRanges[clusterGRanges$set == s]
                    clusterGRanges.newExon <- clusterGRanges.newExon[order(start(clusterGRanges.newExon))]

                    newExons <- clusterExons[rep(1, length(clusterGRanges.newExon) - 1)]
                    ranges(newExons) <- IRanges(
                        start = end(clusterGRanges.newExon)[-length(clusterGRanges.newExon)],
                        end = start(clusterGRanges.newExon)[-1]
                    )

                    newExons$exon_id <- paste0("newEXON_", c(seq_len(length(clusterGRanges.newExon) - 1)))

                    tids <- unique(clusterExons$transcript_id)
                    newExons <- rep(newExons, length(tids))
                    newExons$transcript_id <- rep(tids, each = length(clusterGRanges.newExon) - 1)

                    clusterExons <- c(clusterExons, newExons)
                    clusterExons <- reorderExonNumbers(clusterExons)
                }
            }


            move <- which(clusterGRanges$verdict %in% c("cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored"))
            clusterExons.novel <- NULL
            setTrack <- vector()
            for (m in seq_along(move)) {
                if (clusterGRanges$verdict[move[m]] == "cryptic_threeprime") {
                    clusterExons.alt <- makeNewLeafExons(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons, splice.type = 3)
                } else if (clusterGRanges$verdict[move[m]] == "cryptic_fiveprime") {
                    clusterExons.alt <- makeNewLeafExons(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons, splice.type = 5)
                } else if (clusterGRanges$verdict[move[m]] == "cryptic_unanchored") {
                    clusterExons.alt <- makeNewLeafExonsUnanchored(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons)
                }
                setTrack <- c(setTrack, rep(clusterGRanges.intron$set[move[m]], length(clusterExons.alt)))
                clusterExons.novel <- c(clusterExons.novel, clusterExons.alt)
            }
            if (length(move) > 1) {
                clusterExons.novel <- do.call("c", clusterExons.novel)
            } else if (length(move) == 1) {
                clusterExons.novel <- clusterExons.novel[[1]]
            }
            pairs <- which(clusterGRanges$verdict %in% c("novel annotated pair"))
            clusterExons.pairs <- NULL
            for (p in seq_along(pairs)) {
                clusterGRanges.start <- clusterGRanges[pairs[p]]
                end(clusterGRanges.start) <- start(clusterGRanges.start)
                ol <- findOverlaps(clusterGRanges.start, clusterExons, type = "end")
                clusterExons.start <- clusterExons[ol@to]
                rep.end <- length(clusterExons.start)
                start.ids <- clusterExons.start$transcript_id

                clusterGRanges.end <- clusterGRanges[pairs[p]]
                start(clusterGRanges.end) <- end(clusterGRanges.end)
                ol <- findOverlaps(clusterGRanges.end, clusterExons, type = "start")
                clusterExons.end <- clusterExons[ol@to]
                rep.start <- length(clusterExons.end)
                end.ids <- clusterExons.end$transcript_id

                clusterExons.start <- rep(clusterExons.start, rep.start)
                clusterExons.start$transcript_id <- rep(end.ids, each = rep.end)

                clusterExons.end <- rep(clusterExons.end, rep.end)
                clusterExons.end$transcript_id <- rep(start.ids, each = rep.start)

                clusterExons.pairs <- c(clusterExons.pairs, clusterExons.start, clusterExons.end)

                # remember to reorder exon numbers
                setTrack <- c(setTrack, rep(clusterGRanges.intron$set[pairs[p]], length(clusterExons.start) + length(clusterExons.end)))
            }
            if (is.list(clusterExons.pairs) & length(clusterExons.pairs) > 1) {
                clusterExons.pairs <- do.call("c", clusterExons.pairs)
            }
            if (!is.null(clusterExons.pairs) & !is.null(clusterExons.novel)) {
                clusterExons.novel <- c(clusterExons.novel, clusterExons.pairs)
            } else if (!is.null(clusterExons.pairs) & is.null(clusterExons.novel)) {
                clusterExons.novel <- clusterExons.pairs
            }

            # rm(clusterExons.allSets)
            clusterExons.allSets <- NULL

            if (!is.null(junctions)) {
                # find junctions which are within the cluster ranges,
                # but aren't the same junctions
                ol <- findOverlaps(clusterGRanges, junctions)
                jnc <- junctions[ol@to]

                newJunc <- which(!(paste0(start(jnc), "-", end(jnc) + 1) %in%
                    paste0(start(clusterGRanges), "-", end(clusterGRanges))))
                jnc <- jnc[newJunc]
            }

            for (i in seq_along(unique(clusterGRanges$set))) {
                if (i %in% setTrack) {
                    clusterExons.tid <- paste0(clusterExons$exon_id, clusterExons$transcript_id)
                    clusterExons.novel.tid <- paste0(clusterExons.novel$exon_id, clusterExons.novel$transcript_id)[setTrack == i]
                    keep <- which(!(clusterExons.tid %in% clusterExons.novel.tid))
                    clusterExons.replace <- c(clusterExons[keep], clusterExons.novel[setTrack == i])
                } else {
                    clusterExons.replace <- clusterExons
                }

                clusterExons.replace <- reorderExonNumbers(clusterExons.replace)

                clusterGRanges.max <-
                    clusterGRanges.intron[clusterGRanges.intron$set == i]
                start(clusterGRanges.max) <- min(start(clusterGRanges.max))
                end(clusterGRanges.max) <- max(end(clusterGRanges.max))

                overlapsCluster <- findOverlaps(clusterGRanges.max, clusterExons.replace)
                rmExons <- clusterExons.replace[unique(overlapsCluster@to)]
                if (length(rmExons) > 0) {
                    clusterExonsBounding <- clusterExons.replace[-unique(overlapsCluster@to)]
                } else {
                    clusterExonsBounding <- clusterExons.replace
                }

                # overlaps start of the intron
                clusterGRanges.start <- clusterGRanges[clusterGRanges$set == i]
                end(clusterGRanges.start) <- start(clusterGRanges.start)
                overlapsStart <-
                    findOverlaps(clusterGRanges.start, rmExons, type = "end")
                # overlaps end of the intron
                clusterGRanges.end <- clusterGRanges[clusterGRanges$set == i]
                start(clusterGRanges.end) <- end(clusterGRanges.end)
                overlapsEnd <-
                    findOverlaps(clusterGRanges.end, rmExons, type = "start")

                exonsStart <- rmExons[overlapsStart@to]
                exonsStart <- removeSameExon(exonsStart)

                exonsEnd <- rmExons[overlapsEnd@to]
                exonsEnd <- removeSameExon(exonsEnd)

                InternalExons <- removeSameExon(c(exonsStart, exonsEnd))

                # make sure exons are within the intronic region
                overlapsIntron <-
                    as.data.frame(findOverlaps(InternalExons,
                        clusterGRanges.max,
                        type = "within"
                    ))
                InternalExons <- InternalExons[unique(overlapsIntron$queryHits)]

                overlapsExonStart <- findOverlaps(clusterGRanges.start,
                    clusterExonsBounding,
                    type = "end"
                )
                overlapsExonEnd <- findOverlaps(clusterGRanges.end,
                    clusterExonsBounding,
                    type = "start"
                )

                keepTranscriptIds <-
                    clusterExonsBounding$transcript_id[overlapsExonStart@to]
                keepTranscriptIds <-
                    keepTranscriptIds[clusterExonsBounding$transcript_id[
                        overlapsExonStart@to
                    ] %in%
                        clusterExonsBounding$transcript_id[overlapsExonEnd@to]]

                if (length(keepTranscriptIds) > 0) {
                    if (any(grepl("_dnre_", keepTranscriptIds) |
                        grepl("_upre_", keepTranscriptIds))) {
                        # direction to to remove
                        # if upre, remove dnre isforms
                        removeDirection <- ifelse(clusterGRanges$direction[
                            match(i, clusterGRanges$set)[1]
                        ] == "+", "dnre", "upre")
                        keepTranscriptIds <- keepTranscriptIds[
                            !grepl(removeDirection, keepTranscriptIds)
                        ]
                    }

                    clusterExonsBounding <- clusterExonsBounding[
                        clusterExonsBounding$transcript_id %in% keepTranscriptIds
                    ]

                    if (replaceInternalExons == FALSE) {
                        InternalExons.reps <-
                            InternalExons[rep(
                                seq_along(InternalExons),
                                length(unique(keepTranscriptIds))
                            )]
                        InternalExons.reps$transcript_id <-
                            rep(unique(keepTranscriptIds), each = length(InternalExons))

                        clusterExonsBounding <- c(clusterExonsBounding, InternalExons.reps)
                    } else {
                        if (!is.null(junctions)) {
                            replaceExons <- leafcutterIntronsToExons(clusterGRanges[clusterGRanges$set == i])
                            # keep junctions which overlap any exons created from the cluster junctions
                            if (!is.null(replaceExons)) {
                                keep <- findOverlaps(jnc, replaceExons, type = "within")
                                jnc.set <- jnc[unique(keep@from)]
                                if (length(jnc.set) > 0) {
                                    # annotate the same as clusterGRanges
                                    jnc.set$id <- clusterGRanges$id[clusterGRanges$set == i][1]
                                    jnc.set$direction <- clusterGRanges$direction[clusterGRanges$set == i][1]
                                    jnc.set$verdict <- "junction"
                                    jnc.set$deltapsi <- 0
                                    jnc.set$set <- i
                                    jnc.set$count <- NULL

                                    # remove any overlapping junctions
                                    jnc.set <- jnc.set[rev(order(width(jnc.set)))]
                                    ol <- findOverlaps(jnc.set)
                                    rm <- as.data.frame(table(ol@from))
                                    while (any(rm$Freq > 1)) {
                                        rm <- rm[rev(order(rm$Freq)), ]
                                        jnc.set <- jnc.set[-as.numeric(rm$Var1[1])]
                                        ol <- findOverlaps(jnc.set)
                                        rm <- as.data.frame(table(ol@from))
                                    }
                                    replaceExons <- leafcutterIntronsToExons(c(clusterGRanges[clusterGRanges$set == i], jnc.set))
                                }
                            }
                        } else {
                            replaceExons <- leafcutterIntronsToExons(clusterGRanges[clusterGRanges$set == i])
                        }
                        if (!is.null(replaceExons)) {
                            replaceExonsFormat <- clusterExonsBounding[rep(1, length(replaceExons))]
                            ranges(replaceExonsFormat) <- ranges(replaceExons)

                            replaceExonsFormat.reps <-
                                replaceExonsFormat[rep(
                                    seq_along(replaceExonsFormat),
                                    length(unique(keepTranscriptIds))
                                )]
                            replaceExonsFormat.reps$transcript_id <-
                                rep(unique(keepTranscriptIds), each = length(replaceExonsFormat))

                            clusterExonsBounding <- c(clusterExonsBounding, replaceExonsFormat.reps)
                        }
                    }
                    clusterExonsBounding <- reorderExonNumbers(clusterExonsBounding)
                    clusterExonsBounding$set <- as.numeric(i)

                    if (!exists("clusterExons.allSets")) {
                        clusterExons.allSets <- clusterExonsBounding
                    } else {
                        clusterExons.allSets <- c(
                            clusterExons.allSets,
                            clusterExonsBounding
                        )
                    }
                }
            }

            if (!is.null(clusterExons.allSets)) {
                if (length(clusterExons.allSets) == 1) {
                    clusterExons.allSets <- clusterExons.allSets[[1]]
                } else if (length(clusterExons.allSets) > 1) {
                    clusterExons.allSets <- do.call("c", clusterExons.allSets)
                }

                m <- match(clusterExons.allSets$set, clusterGRanges$set)
                # n <- which(!grepl("[+]", clusterExons.allSets$transcript_id))
                clusterExons.allSets$new_transcript_id <- NA
                clusterExons.allSets$new_transcript_id <-
                    paste0(
                        clusterExons.allSets$transcript_id, "+AS ",
                        ifelse(clusterGRanges$direction[m] == "+", "upre", "dnre"),
                        "_",
                        gsub("_", "", clusterGRanges$id[m]),
                        "-", clusterExons.allSets$set
                    )

                n <- which(grepl("[+]", clusterExons.allSets$transcript_id))
                clusterExons.allSets$new_transcript_id[n] <-
                    paste0(
                        clusterExons.allSets$transcript_id, ":",
                        gsub("_", "", clusterGRanges$id[m]),
                        "-", clusterExons.allSets$set
                    )[n]

                clusterExons.allSets$transcript_id <-
                    clusterExons.allSets$new_transcript_id
                clusterExons.allSets$new_transcript_id <- NULL
                clusterExons.allSets$set <- NULL

                dnre <- clusterExons.allSets[grep("dnre", clusterExons.allSets$transcript_id)]
                upre <- clusterExons.allSets[grep("upre", clusterExons.allSets$transcript_id)]
                if (length(dnre) > 0 & length(upre) > 0) {
                    clusterExons.allSets <- c(
                        removeDuplicateTranscripts(dnre),
                        removeDuplicateTranscripts(upre)
                    )
                    return(clusterExons.allSets)
                } else {
                    return(NULL)
                }
            } else {
                return(NULL)
            }
        } else {
            return(NULL)
        }
    } else {
        return(NULL)
    }
}

#' Create junction ranges from STAR junction files
#'
#' Create junction ranges from STAR junction files
#' @param junctionFiles list of files with STAR junctions
#' @param minReads minimum number of reads a junction requires to be kept
#' @return GRanges object with junctions
#' @export
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
#' @examples
#' junction_files <- list.files(system.file("extdata","leaf_small",
#' package="GeneStructureTools"), full.names=TRUE, pattern=".junc")
#' leaf_junc <- readLeafcutterJunctions(junction_files)

readLeafcutterJunctions <- function(junctionFiles, minReads = 10) {
    junctions.all <- NULL
    for (f in seq_along(junctionFiles)) {
        junctions <- utils::read.delim(junctionFiles[f], header = FALSE)

        index <- match(
            paste(junctions$V1, junctions$V2, junctions$V3, junctions$V6, sep = "_"),
            paste(junctions.all$V1, junctions.all$V2, junctions.all$V3, junctions.all$V6, sep = "_")
        )

        junctions.all$V5[index[which(!is.na(index))]] <- junctions.all$V5[index[which(!is.na(index))]] + junctions$V5[which(!is.na(index))]

        junctions.all <- rbind(junctions.all, junctions[which(is.na(index)), ])
    }
    junctions.all <- junctions.all[junctions.all$V5 >= minReads, ]

    junctionGrange <- GRanges(
        seqnames = junctions.all$V1,
        ranges = IRanges(
            start = junctions.all$V2,
            end = junctions.all$V3
        ),
        strand = junctions.all$V6,
        count = junctions.all$V5
    )

    return(junctionGrange)
}
