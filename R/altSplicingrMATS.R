#' Generate isoforms with and without a skipped exon (or mututally exclusive exons)
#' @param rmatsEvents data.frame containing RMATS SE or MXE events
#' @param eventType type of event to skip exons for. "SE" - skipped exons, or "MXE" - mutally exclusive exons
#' @param exons reference exons GRanges
#' @return data.frame with overlapping event/exons
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
skipExonByJunction <- function(rmatsEvents,
                               eventType="SE",
                               exons){

    # find reference exons that overlap the up/downstream exons (junctions)
    # check for MXE
    granges.upstream <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$upstreamES+1, end=rmatsEvents$upstreamEE), strand=rmatsEvents$strand,
                               id=rmatsEvents$ID, event_id=make.unique(paste0(rmatsEvents$exonStart_0base+1, "-", rmatsEvents$exonEnd), sep="_"))
    granges.downstream <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$downstreamES+1, end=rmatsEvents$downstreamEE), strand=rmatsEvents$strand,
                                 id=rmatsEvents$ID, event_id=make.unique(paste0(rmatsEvents$exonStart_0base+1, "-", rmatsEvents$exonEnd), sep="_"))
    if(eventType == "MXE"){
        granges.upstream$event_id <- make.unique(paste0(rmatsEvents$`1stExonStart_0base`+1, "-", rmatsEvents$`2ndExonEnd`), sep="_")
        granges.downstream$event_id <- make.unique(paste0(rmatsEvents$`1stExonStart_0base`+1, "-", rmatsEvents$`2ndExonEnd`), sep="_")
    }

    seqlevelsStyle(granges.upstream) <- seqlevelsStyle(exons)[1]
    seqlevelsStyle(granges.downstream) <- seqlevelsStyle(exons)[1]
    #
    ol.upstream <- annotateOverlapRmats(granges.upstream, exons, exon_number=1)
    ol.upstream$new_end <- end(granges.upstream)[ol.upstream$queryHits]
    ol.downstream <- annotateOverlapRmats(granges.downstream, exons, exon_number=2)
    ol.downstream$new_start <- start(granges.downstream)[ol.downstream$queryHits]

    betweenExons <- ol.upstream[which(paste0(ol.upstream$from_id, "_", ol.upstream$transcript_id) %in% paste0(ol.downstream$from_id, "_", ol.downstream$transcript_id)),-c(1,2)]
    m.between <- match(paste0(betweenExons$from_id, "_", betweenExons$transcript_id),paste0(ol.downstream$from_id, "_", ol.downstream$transcript_id))
    betweenExons$exon_id2 <- ol.downstream$exon_id[m.between]
    betweenExons$exon_number2 <- ol.downstream$exon_number[m.between]
    betweenExons$new_start <- ol.downstream$new_start[m.between]
    betweenExons$new_transcript_id <- paste0(betweenExons$transcript_id, "+AS ", betweenExons$from_id, "-", betweenExons$event_id)

    betweenExons <- removeDuplicatePairs(betweenExons)
    gtfTranscripts <- duplicateReference(betweenExons, exons)
    gtfTranscripts.rm <- removeExonsBetween(betweenExons, gtfTranscripts)

    # make sure long exons are split
    x <- splitLongExons(betweenExons, gtfTranscripts.rm)
    betweenExons <- x$between
    gtfTranscripts.rm <- x$ranges

    # replace boundries of up/downstream exons
    m1 <- match(paste0(betweenExons$exon_id1, ":", betweenExons$new_transcript_id), paste0(gtfTranscripts.rm$exon_id, ":", gtfTranscripts.rm$new_transcript_id))
    m2 <- match(paste0(betweenExons$exon_id2, ":", betweenExons$new_transcript_id), paste0(gtfTranscripts.rm$exon_id, ":", gtfTranscripts.rm$new_transcript_id))

    end(gtfTranscripts.rm)[m1] <- betweenExons$new_end
    start(gtfTranscripts.rm)[m2] <- betweenExons$new_start

    #### Add in alternatively skipped exons

    # make eventCoords
    # check that rmatsEvents has upstream/downstream EE/ES
    if(eventType == "CE" | eventType =="SE"){
        eventCoords <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$exonStart_0base+1, end=rmatsEvents$exonEnd), strand=rmatsEvents$strand, id=rmatsEvents$ID)
        seqlevelsStyle(eventCoords) <- seqlevelsStyle(exons)[1]
    }else if(eventType == "MXE"){
        eventCoords <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$`1stExonStart_0base`+1, end=rmatsEvents$`1stExonEnd`), strand=rmatsEvents$strand, id=rmatsEvents$ID)
        eventCoords.mxe2 <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$`2ndExonStart_0base`+1, end=rmatsEvents$`2ndExonEnd`), strand=rmatsEvents$strand, id=rmatsEvents$ID)
        seqlevelsStyle(eventCoords) <- seqlevelsStyle(exons)[1]
        seqlevelsStyle(eventCoords.mxe2) <- seqlevelsStyle(exons)[1]
    }

    eventCoords <- annotateEventCoords(eventCoords, betweenExons, exons)

    # create ranges for skipped exon 2 (in MXE pairs)
    if(eventType == "MXE"){
        eventCoords.mxe2 <- annotateEventCoords(eventCoords.mxe2, betweenExons, exons)
    }

    mcols(gtfTranscripts.rm) <- mcols(
        gtfTranscripts.rm)[,c('gene_id','new_transcript_id',
                              'transcript_type','exon_id',
                              'exon_number')]
    colnames(mcols(gtfTranscripts.rm))[2] <- "transcript_id"

    mcols(eventCoords) <- mcols(
        eventCoords)[,c('gene_id','new_transcript_id',
                        'transcript_type','exon_id',
                        'exon_number')]
    colnames(mcols(eventCoords))[2] <- "transcript_id"

    if(eventType == "MXE"){
        mcols(eventCoords.mxe2) <- mcols(
            eventCoords.mxe2)[,c('gene_id','new_transcript_id',
                                 'transcript_type','exon_id',
                                 'exon_number')]
        colnames(mcols(eventCoords.mxe2))[2] <- "transcript_id"

        # add skipped exon back in
        gtfTranscripts.withExon1 <- gtfTranscripts.rm
        gtfTranscripts.withExon1 <- c(gtfTranscripts.withExon1, eventCoords)
        gtfTranscripts.withExon1 <- reorderExonNumbers(gtfTranscripts.withExon1)

        gtfTranscripts.withExon2 <- gtfTranscripts.rm
        gtfTranscripts.withExon2 <- c(gtfTranscripts.withExon2, eventCoords.mxe2)
        gtfTranscripts.withExon2 <- reorderExonNumbers(gtfTranscripts.withExon2)

        gtfTranscripts.withExon1$set <- "included_exon1"
        gtfTranscripts.withExon2$set <- "included_exon2"

        altSplicedTranscripts <- c(gtfTranscripts.withExon1, gtfTranscripts.withExon2)

        # rename included isoforms
        altSplicedTranscripts$transcript_id[which(
            altSplicedTranscripts$set=="included_exon1")] <-
            gsub("AS", "ASMXE1", altSplicedTranscripts$transcript_id[
                which(altSplicedTranscripts$set=="included_exon1")])
        # rename included isoforms
        altSplicedTranscripts$transcript_id[which(
            altSplicedTranscripts$set=="included_exon2")] <-
            gsub("AS", "ASMXE2", altSplicedTranscripts$transcript_id[
                which(altSplicedTranscripts$set=="included_exon2")])

    }else{

        # add skipped exon back in
        gtfTranscripts.withExon <- gtfTranscripts.rm
        gtfTranscripts.withExon <- c(gtfTranscripts.withExon, eventCoords)

        gtfTranscripts.withExon <- reorderExonNumbers(gtfTranscripts.withExon)
        gtfTranscripts.rm <- reorderExonNumbers(gtfTranscripts.rm)

        gtfTranscripts.rm$set <- "skipped_exon"
        gtfTranscripts.withExon$set <- "included_exon"

        altSplicedTranscripts <- c(gtfTranscripts.rm, gtfTranscripts.withExon)

        # rename skipped/included isoforms
        altSplicedTranscripts$transcript_id[which(
            altSplicedTranscripts$set=="skipped_exon")] <-
            gsub("AS", "ASSE", altSplicedTranscripts$transcript_id[
                which(altSplicedTranscripts$set=="skipped_exon")])
        altSplicedTranscripts$transcript_id[which(
            altSplicedTranscripts$set=="included_exon")] <-
            gsub("AS", "ASIE", altSplicedTranscripts$transcript_id[
                which(altSplicedTranscripts$set=="included_exon")])
    }


    altSplicedTranscripts <- reorderExonNumbers(altSplicedTranscripts)

    altSplicedTranscripts$event_id <-
        unlist(lapply(str_split(lapply(stringr::str_split(altSplicedTranscripts$transcript_id, " "),"[[",2), "-"), "[[" , 1))

    mcols(altSplicedTranscripts) <-
        mcols(altSplicedTranscripts)[,c('gene_id','transcript_id',
                                        'transcript_type','exon_id',
                                        'exon_number',
                                        'set', 'event_id')]

    return(altSplicedTranscripts)
}

#' Generate isoforms with and without a retain intron from RMATS data
#' shortcut function: uses the same base function for intron retention
#' @param rmatsEvents data.frame containing RMATS RI events
#' @param exons reference exons GRanges
#' @return GRanges retained and skipped intron isoforms
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#'
altIntronRmats <- function(rmatsEvents, exons){
    events.RI <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$upstreamEE+1,
                                                               end=rmatsEvents$downstreamES),
                        strand=rmatsEvents$strand)
    seqlevelsStyle(events.RI) <- seqlevelsStyle(exons)[1]
    events.RI$id <- paste0(rmatsEvents$ID, "-",
                          as.character(seqnames(events.RI)),":",
                          rmatsEvents$riExonStart_0base+1, "-", rmatsEvents$riExonEnd)

    exons.intronRetention <- findIntronContainingTranscripts(rmatsEvents=events.RI, exons, match="exact")
    isoforms.RI <- addIntronInTranscript(flankingExons=exons.intronRetention, exons, match="retain")

    return(isoforms.RI)

}

#' Generate isoforms with different 5' or 3' splice site usage from RMATS data
#' @param rmatsEvents data.frame containing RMATS RI events
#' @param exons reference exons GRanges
#' @param eventType type of event. "A5E" - alternative 5', or "A3E" - alternative 3'
#' @return GRanges of isoforms
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
altSpliceSiteRmats <- function(rmatsEvents,
                              exons,
                              eventType){

    # new event id.. .to match whippet
    if(eventType %in% c("A5E", "A5SS")){
        rmatsEvents$event_range <- ifelse(rmatsEvents$strand == "+", paste0(rmatsEvents$shortEE, "-", rmatsEvents$longExonEnd), paste0(rmatsEvents$longExonStart_0base+1, "-", rmatsEvents$shortES+1))
    }else{
        rmatsEvents$event_range <- ifelse(rmatsEvents$strand == "-", paste0(rmatsEvents$shortEE, "-", rmatsEvents$longExonEnd), paste0(rmatsEvents$longExonStart_0base+1, "-", rmatsEvents$shortES+1))
    }

    granges.longExon <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$longExonStart_0base+1, end=rmatsEvents$longExonEnd), strand=rmatsEvents$strand, id=rmatsEvents$ID, event_id=rmatsEvents$event_range)
    granges.shortExon <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$shortES+1, end=rmatsEvents$shortEE), strand=rmatsEvents$strand, id=rmatsEvents$ID, event_id=rmatsEvents$event_range)
    granges.flankingExon <- GRanges(seqnames=rmatsEvents$chr, ranges=IRanges(start=rmatsEvents$flankingES+1, end=rmatsEvents$flankingEE), strand=rmatsEvents$strand, id=rmatsEvents$ID, event_id=rmatsEvents$event_range)


    seqlevelsStyle(granges.longExon) <- seqlevelsStyle(exons)[1]
    seqlevelsStyle(granges.shortExon) <- seqlevelsStyle(exons)[1]
    seqlevelsStyle(granges.flankingExon) <- seqlevelsStyle(exons)[1]


    ol.longExon <- annotateOverlapRmats(granges.longExon, exons, exon_number=1)
    ol.longExon$long_start <- start(granges.longExon)[ol.longExon$queryHits]
    ol.longExon$long_end <- end(granges.longExon)[ol.longExon$queryHits]
    ol.longExon$short_start <- start(granges.shortExon)[ol.longExon$queryHits]
    ol.longExon$short_end <- end(granges.shortExon)[ol.longExon$queryHits]

    ol.flankingExon <- annotateOverlapRmats(granges.flankingExon, exons, exon_number=2)
    ol.flankingExon$flanking_start <- start(granges.flankingExon)[ol.flankingExon$queryHits]
    ol.flankingExon$flanking_end <- end(granges.flankingExon)[ol.flankingExon$queryHits]


    betweenExons <- ol.longExon[which(paste0(ol.longExon$from_id, "_", ol.longExon$transcript_id) %in% paste0(ol.flankingExon$from_id, "_", ol.flankingExon$transcript_id)),-c(1,2)]
    m.between <- match(paste0(betweenExons$from_id, "_", betweenExons$transcript_id),paste0(ol.flankingExon$from_id, "_", ol.flankingExon$transcript_id))
    betweenExons$exon_id2 <- ol.flankingExon$exon_id[m.between]
    betweenExons$exon_number2 <- ol.flankingExon$exon_number[m.between]
    betweenExons$flanking_start <- ol.flankingExon$flanking_start[m.between]
    betweenExons$flanking_end <- ol.flankingExon$flanking_end[m.between]
    betweenExons$new_transcript_id <- paste0(betweenExons$transcript_id, "+AS ", betweenExons$from_id, "-", betweenExons$event_id)
    betweenExons <- removeDuplicatePairs(betweenExons)

    gtfTranscripts <- duplicateReference(betweenExons, exons)
    gtfTranscripts.rm <- removeExonsBetween(betweenExons, gtfTranscripts)

    # make sure long exons are split
    x <- splitLongExons(betweenExons, gtfTranscripts.rm)
    betweenExons <- x$between
    gtfTranscripts.rm <- x$ranges

    # replace boundaries of up/downstream exons
    # m1 <- match to the LONG/SHORT exon
    m1 <- match(paste0(betweenExons$exon_id1, ":", betweenExons$new_transcript_id), paste0(gtfTranscripts.rm$exon_id, ":", gtfTranscripts.rm$new_transcript_id))
    # m2 <- match to the FLANKING exon
    m2 <- match(paste0(betweenExons$exon_id2, ":", betweenExons$new_transcript_id), paste0(gtfTranscripts.rm$exon_id, ":", gtfTranscripts.rm$new_transcript_id))
    strand <- as.character(strand(gtfTranscripts.rm)[m1])

    if(eventType %in% c("A5E", "A5SS")){
        long_5 <- which(strand == "+")
        long_3 <- which(strand == "-")
    }else{
        long_5 <- which(strand == "-")
        long_3 <- which(strand == "+")
    }

    # replace the non-alt junction first
    start(gtfTranscripts.rm[m2][long_5]) <- betweenExons$flanking_start[long_5]
    end(gtfTranscripts.rm[m2][long_3]) <- betweenExons$flanking_end[long_3]

    # extend exons (((IF))) they aren't covered by both >>>short<<< & long junctions, otherwise keep them the same

    extend3 <- which(end(gtfTranscripts.rm[m1][long_3]) < betweenExons$short_start[long_3])
    if(length(extend3) > 0){
        end(gtfTranscripts.rm[m1][long_3][extend3]) <- betweenExons$short_end[long_3][extend3]
    }
    extend5 <- which(start(gtfTranscripts.rm[m1][long_5]) > betweenExons$short_end[long_5])
    if(length(extend5) > 0){
        start(gtfTranscripts.rm[m1][long_5][extend5]) <- betweenExons$short_start[long_5][extend5]
    }
    # make a long/short version

    altJunctionLong <- gtfTranscripts.rm
    altJunctionShort <- gtfTranscripts.rm

    end(altJunctionLong[m1][long_5]) <- betweenExons$long_end[long_5]
    end(altJunctionShort[m1][long_5]) <- betweenExons$short_end[long_5]

    start(altJunctionLong[m1][long_3]) <- betweenExons$long_start[long_3]
    start(altJunctionShort[m1][long_3]) <- betweenExons$short_start[long_3]

    altJunctionLong$set <- ifelse(eventType %in% c("A5E", "A5SS"), "alt5_splicesite_long", "alt3_splicesite_long")
    altJunctionLong$new_transcript_id[altJunctionLong$set == "alt5_splicesite_long"] <- gsub("AS", "ASA5L", altJunctionLong$new_transcript_id[altJunctionLong$set == "alt5_splicesite_long"])
    altJunctionLong$new_transcript_id[altJunctionLong$set == "alt3_splicesite_long"] <- gsub("AS", "ASA3L", altJunctionLong$new_transcript_id[altJunctionLong$set == "alt3_splicesite_long"])
    altJunctionShort$set <- ifelse(eventType %in% c("A5E", "A5SS"), "alt5_splicesite_short","alt3_splicesite_short")
    altJunctionShort$new_transcript_id[altJunctionShort$set == "alt5_splicesite_short"] <- gsub("AS", "ASA5S", altJunctionShort$new_transcript_id[altJunctionShort$set == "alt5_splicesite_short"])
    altJunctionShort$new_transcript_id[altJunctionShort$set == "alt3_splicesite_short"] <- gsub("AS", "ASA3S", altJunctionShort$new_transcript_id[altJunctionShort$set == "alt3_splicesite_short"])

    altSplicedTranscripts <- c(altJunctionLong, altJunctionShort)

    altSplicedTranscripts$transcript_id <- altSplicedTranscripts$new_transcript_id
    altSplicedTranscripts <- reorderExonNumbers(altSplicedTranscripts)

    altSplicedTranscripts$event_id <-
        unlist(lapply(str_split(lapply(stringr::str_split(altSplicedTranscripts$transcript_id, " "),"[[",2), "-"), "[[" , 1))

    mcols(altSplicedTranscripts) <-
        mcols(altSplicedTranscripts)[,c('gene_id','transcript_id',
                                        'transcript_type','exon_id',
                                        'exon_number',
                                        'set', 'event_id')]

    return(altSplicedTranscripts)
}

