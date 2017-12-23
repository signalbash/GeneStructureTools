#' Given the location of a whole retained intron, find transcripts which splice out this intron
#' @param eventCoords GRanges object with ranges for introns
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param match what type of matching to perform? exact = only exons which bound the intron exactly,
#' introns = any exon pairs which overlap the intron,
#' all = any exon pairs AND single exons which overlap the intron
#' @return data.frame with all flanking exon pairs
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/", package = "GeneStructureTools"), full.names = TRUE)
#' diffFiles <- whippetFiles[grep(".diff", whippetFiles)]
#' whippetDiffSplice <- readWhippetDIFFfiles(diffFiles)
#' whippetCoords <- formatWhippetEvents(whippetDiffSplice)
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' gtf.exons <- gtf[gtf$type=="exon"]
#' gtf.transcripts <- gtf[gtf$type=="transcript"]
#' event.intronRetention <- whippetDiffSplice[which(whippetDiffSplice$type=="RI")[1],]
#' coords.intronRetention <- whippetCoords[whippetCoords$id %in% event.intronRetention$coord]
#' exons.intronRetention <- findIntronContainingTranscripts(coords.intronRetention, gtf.exons)
findIntronContainingTranscripts <- function(eventCoords, gtf.exons, match="exact"){

    moved <- FALSE
    eventCoords = ranges.ri
    # remove any duplicates
    overlaps <- GenomicRanges::findOverlaps(eventCoords, type="equal")
    overlaps <- overlaps[which(overlaps@from != overlaps@to)]
    if(length(overlaps) > 0){
        overlaps <- overlaps[which(overlaps@from < overlaps@to)]
        if(length(overlaps) > 0){
            duplicates <- unique(overlaps@to)
            eventCoords <- eventCoords[-duplicates]
        }
    }

    # start of intron // end of exon a
    rangeRI.start <- eventCoords
    end(rangeRI.start) <- start(rangeRI.start)

    overlaps <- GenomicRanges::findOverlaps(rangeRI.start, gtf.exons, type="end")

    # catch if intron coords dont overlap the 1nt exon start/end
    if(length(overlaps) == 0){
        start(rangeRI.start) <- start(rangeRI.start) -1
        end(rangeRI.start) <- start(rangeRI.start)
        overlaps <- GenomicRanges::findOverlaps(rangeRI.start, gtf.exons, type="end")
        # fix original
        start(eventCoords) <- start(eventCoords) -1
        end(eventCoords) <- end(eventCoords) +1
        moved <- TRUE
    }

    gtf.fromA <- gtf.exons[overlaps@to]
    gtf.fromA$from <- overlaps@from
    gtf.fromA$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf.fromA)),
                             paste0(transcript_id, "_",from))

    # end of intron // start of exon b
    rangeRI.end <- eventCoords
    start(rangeRI.end) <- end(rangeRI.end)

    overlaps <- GenomicRanges::findOverlaps(rangeRI.end, gtf.exons, type="start")
    if(length(overlaps) == 0){
        end(rangeRI.end) <- end(rangeRI.end) +1
        start(rangeRI.end) <- end(rangeRI.end)
        overlaps <- GenomicRanges::findOverlaps(rangeRI.end, gtf.exons, type="start")
        moved <- TRUE
    }

    gtf.toA <- gtf.exons[overlaps@to]
    gtf.toA$from <- overlaps@from
    gtf.toA$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf.toA)),
                           paste0(transcript_id, "_",from))

    keep.from <- which(gtf.fromA@elementMetadata$new_id %in% gtf.toA@elementMetadata$new_id)
    gtf.fromA <- gtf.fromA[keep.from]

    m <- match(gtf.fromA$new_id, gtf.toA$new_id)
    gtf.toA <- gtf.toA[m]

    mcols(gtf.toA)$from_exon_number <- as.numeric(gtf.fromA$exon_number)
    mcols(gtf.toA)$to_exon_number <- as.numeric(gtf.toA$exon_number)

    mcols(gtf.toA)$intron_exon_number <- apply(
        GenomicRanges::mcols(gtf.toA)[,c('to_exon_number','from_exon_number')],
        1, mean)
    if(length(gtf.toA) > 0){
        mcols(gtf.toA)$overlaps <- "intron"
    }else{
        mcols(gtf.toA) <- cbind(mcols(gtf.toA), DataFrame(overlaps=character()))
    }
    if(match == "introns" | match=="all"){
        # non-exact overlaps
        ol <- findOverlaps(eventCoords, gtf.exons)
        gtf.overlaps <- gtf.exons[ol@to]
        gtf.overlaps$from <- ol@from
        gtf.overlaps$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf.overlaps)),
                                    paste0(transcript_id, "_",from))
        gtf.overlaps$exon_number <- as.numeric(gtf.overlaps$exon_number)
        # remove perfect match pairs
        gtf.overlaps <- gtf.overlaps[which(!(gtf.overlaps$new_id %in% gtf.toA$new_id))]

        transcriptTable <- as.data.frame(table(gtf.overlaps$new_id))
        gtf.overlapsPairs <- gtf.overlaps[gtf.overlaps$new_id %in%
                                              transcriptTable$Var1[transcriptTable$Freq == 2]]
        startPair <- gtf.overlapsPairs$new_id[start(gtf.overlapsPairs) <
                                                  start(eventCoords[gtf.overlapsPairs$from])]
        endPair <- gtf.overlapsPairs$new_id[end(gtf.overlapsPairs) >
                                                end(eventCoords[gtf.overlapsPairs$from])]
        keep <- startPair[startPair %in% endPair]
        gtf.overlapsPairs <- gtf.overlapsPairs[gtf.overlapsPairs$new_id %in% keep]

        if(length(gtf.overlapsPairs) > 0){
            exon_numbers <- aggregate(exon_number ~ new_id, mcols(gtf.overlapsPairs), mean)
            gtf.overlapsPairs$intron_exon_number <-
                exon_numbers$exon_number[match(gtf.overlapsPairs$new_id,
                                               exon_numbers$new_id)]
            gtf.overlapsPairs <- gtf.overlapsPairs[!duplicated(gtf.overlapsPairs$new_id)]
            gtf.overlapsPairs$from_exon_number <- gtf.overlapsPairs$intron_exon_number - 0.5
            gtf.overlapsPairs$to_exon_number <- gtf.overlapsPairs$intron_exon_number + 0.5
            gtf.overlapsPairs$overlaps <- "nonexact"
            mcols(gtf.overlapsPairs) <-
                mcols(gtf.overlapsPairs[,match(colnames(mcols(gtf.toA)),
                                               colnames(mcols(gtf.overlapsPairs)))])

            gtf.toA <- c(gtf.toA, gtf.overlapsPairs)
        }
    }
    if(match=="all"){
        # overlaps an exon
        gtf.overlaps <- gtf.overlaps[which(!(gtf.overlaps$new_id %in%
                                                 gtf.overlapsPairs$new_id))]
        keep <- which(start(gtf.overlaps) <= start(eventCoords[gtf.overlaps$from]) &
                          end(gtf.overlaps) >= end(eventCoords[gtf.overlaps$from]))
        gtf.overlaps <- gtf.overlaps[keep]
        if(length(gtf.overlaps) > 0){
            gtf.overlaps$intron_exon_number <- gtf.overlaps$exon_number
            gtf.overlaps$from_exon_number <- gtf.overlaps$exon_number
            gtf.overlaps$to_exon_number <- gtf.overlaps$exon_number
            gtf.overlaps$overlaps <- "exon"
            mcols(gtf.overlaps) <-
                mcols(gtf.overlaps[,match(colnames(mcols(gtf.toA)),
                                          colnames(mcols(gtf.overlaps)))])
            gtf.toA <- c(gtf.toA, gtf.overlaps)
        }
    }

    if(length(gtf.toA) > 0){
        flankingExons <-
            as.data.frame(GenomicRanges::mcols(gtf.toA)[,c('gene_id','transcript_id',
                                                           'transcript_type',
                                                           'from','from_exon_number',
                                                           'intron_exon_number',
                                                           'to_exon_number','overlaps')])

        flankingExons$from <- eventCoords$id[flankingExons$from]
        flankingExons$moved <- moved
        return(flankingExons)
    }else{
        return(NULL)
    }
}

#' Add a retained intron to the transcripts it is skipped by
#' @param eventCoords GRanges object with ranges for introns
#' @param flankingExons data.frame generataed by findIntronContainingTranscripts()
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param glueExons Join together exons that are not seperated by introns?
#' @param match what type of match replacement should be done?
#' exact: exact matches to the intron only
#' retain: keep non-exact intron match coordinates in spliced sets, and retain them in retained sets
#' replace: replace non-exact intron match coordinates with event coordinates in spliced sets, and retain in retained sets
#' @return GRanges with transcripts containing retained introns
#' @export
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/", package = "GeneStructureTools"), full.names = TRUE)
#' diffFiles <- whippetFiles[grep(".diff", whippetFiles)]
#' whippetDiffSplice <- readWhippetDIFFfiles(diffFiles)
#' whippetCoords <- formatWhippetEvents(whippetDiffSplice)
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' gtf.exons <- gtf[gtf$type=="exon"]
#' gtf.transcripts <- gtf[gtf$type=="transcript"]
#' event.intronRetention <- whippetDiffSplice[which(whippetDiffSplice$type=="RI")[1],]
#' coords.intronRetention <- whippetCoords[whippetCoords$id %in% event.intronRetention$coord]
#' exons.intronRetention <- findIntronContainingTranscripts(coords.intronRetention, gtf.exons)
#' IntronRetentionTranscripts <- addIntronInTranscript(coords.intronRetention, exons.intronRetention, gtf.exons)
addIntronInTranscript <- function(eventCoords,
                                  flankingExons,
                                  gtf.exons,
                                  glueExons=TRUE,
                                  match="exact"){

    move <- which(flankingExons$moved == TRUE)
    move.RIindex <- which(eventCoords$id %in% flankingExons$from[move])
    start(eventCoords)[move.RIindex] <- start(eventCoords)[move.RIindex] -1
    end(eventCoords)[move.RIindex] <- end(eventCoords)[move.RIindex] +1

    if(!(match %in% c("exact","retain","replace"))){
        message("match must be 'exact', 'retain', or 'replace'")
        message("using default match = 'exact'")
        match <- "exact"
    }

    if(match == "exact"){
        keep <- which(flankingExons$overlaps == "intron")
        flankingExons <- flankingExons[keep,]
        eventCoords <- eventCoords[eventCoords$id %in% flankingExons$from]
    }

    eventCoords <- eventCoords[match(flankingExons$from, eventCoords$id)]
    eventCoords$exon_number <- flankingExons$intron_exon_number
    eventCoords$transcript_id <- flankingExons$transcript_id
    eventCoords$transcript_type <- flankingExons$transcript_type
    eventCoords$gene_id <- flankingExons$gene_id
    eventCoords$exon_id <- eventCoords$id

    transcripts <- as.data.frame(table(flankingExons$transcript_id))
    gtfTranscripts <- gtf.exons[gtf.exons$transcript_id %in% transcripts$Var1]
    m <- match(gtfTranscripts$transcript_id, eventCoords$transcript_id)
    mcols(gtfTranscripts) <-
        cbind(mcols(gtfTranscripts),
              DataFrame(new_transcript_id=paste0(gtfTranscripts$transcript_id,
                                                 "+AS ", eventCoords$exon_id[m])))
    mcols(eventCoords) <-
        cbind(mcols(eventCoords),
              DataFrame(new_transcript_id = paste0(eventCoords$transcript_id,
                                                   "+AS ", eventCoords$exon_id)))
    flankingExons$new_transcript_id <- paste0(flankingExons$transcript_id,
                                              "+AS ", flankingExons$from)

    #mcols(gtfTranscripts)$new_transcript_id <- paste0(gtfTranscripts$transcript_id,"+AS ",eventCoords$exon_id[m])
    #mcols(eventCoords)$new_transcript_id <- paste0(eventCoords$transcript_id,"+AS ",eventCoords$exon_id)

    mcols(gtfTranscripts) <- mcols(gtfTranscripts)[,c('gene_id','transcript_id',
                                                      'transcript_type','exon_id',
                                                      'exon_number','new_transcript_id')]
    mcols(eventCoords) <- mcols(eventCoords)[,c('gene_id','transcript_id',
                                                'transcript_type','exon_id',
                                                'exon_number','new_transcript_id')]

    needsDuplicated <- which(!(eventCoords$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts.add <- gtfTranscripts[gtfTranscripts$transcript_id %in%
                                                 eventCoords$transcript_id[needsDuplicated]]
    }

    while(length(needsDuplicated) > 0){
        gtfTranscripts.add <-
            gtfTranscripts.add[gtfTranscripts.add$transcript_id %in%
                                   eventCoords$transcript_id[needsDuplicated]]
        m <- match(gtfTranscripts.add$transcript_id,
                   eventCoords$transcript_id[needsDuplicated])
        gtfTranscripts.add$new_transcript_id <-
            paste0(gtfTranscripts.add$transcript_id,
                   "+AS ",eventCoords$exon_id[needsDuplicated][m])
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts.add)
        needsDuplicated <- which(!(eventCoords$new_transcript_id
                                   %in% gtfTranscripts$new_transcript_id))
    }

    # fix starts/ends of introns
    gtfTranscripts$new_id_ex <-
        with(mcols(gtfTranscripts), paste0(new_transcript_id, "_", exon_number))
    flankingExons$from_ex <-
        with(flankingExons, paste0(new_transcript_id, "_", from_exon_number))
    flankingExons$to_ex <-
        with(flankingExons, paste0(new_transcript_id, "_", to_exon_number))

    if(match == "replace"){
        wi <- which(flankingExons$overlaps == "nonexact")
        wi.pos <-
            wi[as.character(strand(gtfTranscripts[match(flankingExons$from_ex[wi],
                                                        gtfTranscripts$new_id_ex)])) == "+"]
        wi.neg <-
            wi[as.character(strand(gtfTranscripts[match(flankingExons$from_ex[wi],
                                                        gtfTranscripts$new_id_ex)])) == "-"]

        if(length(wi.pos) > 0){
            end(gtfTranscripts)[match(flankingExons$from_ex[wi.pos],
                                      gtfTranscripts$new_id_ex)] <-
                start(eventCoords)[match(flankingExons$new_transcript_id[wi.pos],
                                         eventCoords$new_transcript_id)]
            start(gtfTranscripts)[match(flankingExons$to_ex[wi.pos],
                                        gtfTranscripts$new_id_ex)] <-
                end(eventCoords)[match(flankingExons$new_transcript_id[wi.pos],
                                       eventCoords$new_transcript_id)]
        }
        if(length(wi.neg) > 0){
            end(gtfTranscripts)[match(flankingExons$to_ex[wi.neg],
                                      gtfTranscripts$new_id_ex)] <-
                start(eventCoords)[match(flankingExons$new_transcript_id[wi.neg],
                                         eventCoords$new_transcript_id)]
            start(gtfTranscripts)[match(flankingExons$from_ex[wi.neg],
                                        gtfTranscripts$new_id_ex)] <-
                end(eventCoords)[match(flankingExons$new_transcript_id[wi.neg],
                                       eventCoords$new_transcript_id)]
        }
    }

    # create an intron in 'exons'
    we <- which(flankingExons$overlaps == "exon")
    if(length(we) > 0){
        replace <- match(flankingExons$from_ex[we],gtfTranscripts$new_id_ex)
        start.replacement <- gtfTranscripts[replace]
        end(start.replacement) <- start(eventCoords)[match(flankingExons$new_transcript_id[we],
                                                           eventCoords$new_transcript_id)]
        end.replacement <- gtfTranscripts[replace]
        start(end.replacement) <- end(eventCoords)[match(flankingExons$new_transcript_id[we],
                                                         eventCoords$new_transcript_id)]
        gtfTranscripts <- c(gtfTranscripts[-replace], start.replacement, end.replacement)
    }
    gtfTranscripts <- reorderExonNumbers(gtfTranscripts, by="new_transcript_id")
    gtfTranscripts$new_id_ex <- NULL

    gtfTranscripts.withIntron <- c(gtfTranscripts, eventCoords)
    gtfTranscripts.withIntron <- reorderExonNumbers(gtfTranscripts.withIntron,
                                                    by="new_transcript_id")

    mcols(gtfTranscripts.withIntron) <-
        mcols(gtfTranscripts.withIntron)[,c('gene_id','new_transcript_id',
                                            'transcript_type','exon_id','exon_number')]
    colnames(mcols(gtfTranscripts.withIntron))[2] <- "transcript_id"

    mcols(gtfTranscripts) <-
        mcols(gtfTranscripts)[,c('gene_id','new_transcript_id',
                                 'transcript_type','exon_id','exon_number')]
    colnames(mcols(gtfTranscripts))[2] <- "transcript_id"

    gtfTranscripts.withIntron$exon_number <- as.numeric(gtfTranscripts.withIntron$exon_number)
    order <- order(gtfTranscripts.withIntron$transcript_id,
                   gtfTranscripts.withIntron$exon_number)
    gtfTranscripts.withIntron <- gtfTranscripts.withIntron[order]

    gtfTranscripts$set <- "spliced_intron"
    gtfTranscripts.withIntron$set <- "retained_intron"

    gtfTranscripts.withIntron <- c(gtfTranscripts.withIntron, gtfTranscripts)

    # rename retained/spliced isoforms
    gtfTranscripts.withIntron$transcript_id[which(
        gtfTranscripts.withIntron$set=="spliced_intron")] <-
        gsub("AS", "ASSI", gtfTranscripts.withIntron$transcript_id[which(
            gtfTranscripts.withIntron$set=="spliced_intron")])
    gtfTranscripts.withIntron$transcript_id[which(
        gtfTranscripts.withIntron$set=="retained_intron")] <-
        gsub("AS", "ASRI", gtfTranscripts.withIntron$transcript_id[which(
            gtfTranscripts.withIntron$set=="retained_intron")])

    #join together exons that are not seperated by an intron
    if(glueExons==TRUE){

        # split pos/neg ordering
        gtfTranscripts.withIntron.neg <-
            gtfTranscripts.withIntron[which(
                as.logical(strand(gtfTranscripts.withIntron) == "-"))]
        order <- order(gtfTranscripts.withIntron.neg$transcript_id,
                       plyr::desc(gtfTranscripts.withIntron.neg$exon_number))
        gtfTranscripts.withIntron.neg <- gtfTranscripts.withIntron.neg[order]

        gtfTranscripts.withIntron.pos <-
            gtfTranscripts.withIntron[which(
                as.logical(strand(gtfTranscripts.withIntron) == "+"))]
        gtfTranscripts.withIntron <- c(gtfTranscripts.withIntron.pos,
                                       gtfTranscripts.withIntron.neg)

        #extend starts <---<---<---
        w <- which(end(ranges(gtfTranscripts.withIntron))[
            -length(gtfTranscripts.withIntron)] ==
                start(ranges(gtfTranscripts.withIntron[-1])))
        gtfTranscripts.withIntron <- gtfTranscripts.withIntron
        GenomicRanges::start(GenomicRanges::ranges(gtfTranscripts.withIntron))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(gtfTranscripts.withIntron))[w]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtfTranscripts.withIntron, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to - 1]
        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtfTranscripts.withIntron <- gtfTranscripts.withIntron[-rm]
        }

        #extend ends --->--->--->
        overlaps <- findOverlaps(gtfTranscripts.withIntron)
        overlaps <- overlaps[overlaps@from==overlaps@to +1]

        GenomicRanges::end(GenomicRanges::ranges(gtfTranscripts.withIntron))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(gtfTranscripts.withIntron))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtfTranscripts.withIntron, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtfTranscripts.withIntron <- gtfTranscripts.withIntron[-rm]
        }

        #order <- order(gtfTranscripts.withIntron$transcript_id, gtfTranscripts.withIntron$exon_number)
        #gtfTranscripts.withIntron <- gtfTranscripts.withIntron[order]

    }

    gtfTranscripts.withIntron$whippet_id <- unlist(lapply(stringr::str_split(
        gtfTranscripts.withIntron$transcript_id, " "),"[[",2))
    gtfTranscripts.withIntron$overlaps <- NULL

    return(gtfTranscripts.withIntron)
}

