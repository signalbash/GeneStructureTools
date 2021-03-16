#' Find overlaps from a rmats GRanges to reference exons and annotate with ids
#' @param rmatsGRanges rmats event GRanges object
#' @param exons reference exons GRanges
#' @param exon_number which exon in the rmats event this is for (only used for annotating output)
#' @return data.frame with overlapping event/exons
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
annotateOverlapRmats <- function(rmatsGRanges, exons, exon_number=1){
    ol <- as.data.frame(findOverlaps(rmatsGRanges, exons))
    ol$from_id <- rmatsGRanges$id[ol$queryHits]
    ol$transcript_id <- exons$transcript_id[ol$subjectHits]
    ol$exon_id <- exons$exon_id[ol$subjectHits]
    ol$exon_number <- exons$exon_number[ol$subjectHits]
    ol$event_id <- paste0(seqnames(rmatsGRanges)[ol$queryHits], ":", rmatsGRanges$event_id[ol$queryHits])
    colnames(ol)[match(c('exon_id', 'exon_number') ,colnames(ol))] <- paste0(c('exon_id', 'exon_number'), exon_number)
    return(ol)
}

#' remove any duplicate pairs of events/reference transcripts (i.e. long event range which overlaps 2+ exons)
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @return data.frame with related differential splicing event ids and reference transcript_ids
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
removeDuplicatePairs <- function(betweenExons){

    hasDups <- which(duplicated(paste0(betweenExons$new_transcript_id)))
    if(length(hasDups) > 0){
        betweenExons.duplicates <- betweenExons[betweenExons$new_transcript_id %in% betweenExons$new_transcript_id[hasDups],]
        betweenExons.duplicates$exon_num_range <- abs(as.numeric(betweenExons.duplicates$exon_number1) - as.numeric(betweenExons.duplicates$exon_number2))
        betweenExons.duplicates <- arrange(betweenExons.duplicates, new_transcript_id, plyr::desc(exon_num_range))
        betweenExons.duplicates <- betweenExons.duplicates[!duplicated(betweenExons.duplicates$new_transcript_id),]
        betweenExons.duplicates$exon_num_range <- NULL

        betweenExons <- rbind(betweenExons[which(!(betweenExons$new_transcript_id %in% betweenExons$new_transcript_id[hasDups])),],
                              betweenExons.duplicates)
    }

    return(betweenExons)

}

#' Duplicate a reference Granges (exon-level) for each diff-splicing event/transcript combination required.
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param exons reference exons GRanges
#' @return
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
duplicateReference <- function(betweenExons, exons){
    # transcripts containing the exon pairs
    transcripts <- as.data.frame(table(betweenExons$transcript_id))
    gtfTranscripts <- exons[exons$transcript_id %in% transcripts$Var1]

    m <- match(gtfTranscripts$transcript_id, betweenExons$transcript_id)

    mcols(gtfTranscripts) <- cbind(mcols(gtfTranscripts),
                                   DataFrame(new_transcript_id=paste0(
                                       gtfTranscripts$transcript_id, "+AS ",
                                       betweenExons$from_id[m], "-",
                                       betweenExons$event_id[m])))

    needsDuplicated <- which(!(betweenExons$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts[
            gtfTranscripts$transcript_id %in%
                betweenExons$transcript_id[needsDuplicated]]
    }

    while(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts_add[
            gtfTranscripts_add$transcript_id %in%
                betweenExons$transcript_id[needsDuplicated]]
        m <- match(gtfTranscripts_add$transcript_id,
                   betweenExons$transcript_id[needsDuplicated])
        gtfTranscripts_add$new_transcript_id <- paste0(
            gtfTranscripts_add$transcript_id, "+AS ",
            betweenExons$from_id[needsDuplicated][m], "-",
            betweenExons$event_id[needsDuplicated][m])
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts_add)
        needsDuplicated <- which(!(betweenExons$new_transcript_id %in%
                                       gtfTranscripts$new_transcript_id))
    }

    return(gtfTranscripts)
}

#' Generate vector of integers between two numbers (non-inclusive)
#' @param a first integer
#' @param b second integer
#' @return vector of integers
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
betweenNumbers <- function(a, b){
    ab.range <- seq(as.numeric(a),as.numeric(b))
    ab.range <- ab.range[!(ab.range %in% c(a,b))]
    return(ab.range)
}


#' Remove any exons in a transcript within an event range
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param gtfTranscripts Granges with altered transcript structures
#' @return Granges with altered transcript structures
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
removeExonsBetween <- function(betweenExons, gtfTranscripts){
    remove <- apply(betweenExons[,c('exon_number1', 'exon_number2', 'new_transcript_id')], 1, function(x) paste0(x[3], " ", betweenNumbers(x[1], x[2])))
    gtfTranscripts.rm <- gtfTranscripts[which(!(paste0(gtfTranscripts$new_transcript_id, " ", gtfTranscripts$exon_number) %in% unlist(remove)))]
    return(gtfTranscripts.rm)
}

#' Split long exons in two if they overlap an event
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param gtfTranscripts Granges with altered transcript structures
#' @return Granges with altered transcript structures
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
splitLongExons <- function(betweenExons, gtfTranscripts){
    # make sure there is > 1 exon in the pair (i.e. split a single long exon into two)
    duplicate.index <- which(betweenExons$exon_number1 == betweenExons$exon_number2)
    if(length(duplicate.index) > 0){
        duplicate <- paste0(betweenExons$exon_id1[duplicate.index], "_", betweenExons$new_transcript_id[duplicate.index])
        betweenExons$exon_id1[duplicate.index] <- paste0(betweenExons$exon_id1[duplicate.index], "_1")
        betweenExons$exon_id2[duplicate.index] <- paste0(betweenExons$exon_id2[duplicate.index], "_2")
        betweenExons$exon_number1[duplicate.index] <- as.numeric(betweenExons$exon_number1[duplicate.index]) - 0.5
        betweenExons$exon_number2[duplicate.index] <- as.numeric(betweenExons$exon_number2[duplicate.index]) + 0.5

        duplicate.granges.1 <- gtfTranscripts[paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate]
        duplicate.granges.1$exon_number <- as.character(as.numeric(duplicate.granges.1$exon_number) - 0.5)
        duplicate.granges.1$exon_id <- paste0(duplicate.granges.1$exon_id, "_1")
        duplicate.granges.2 <- gtfTranscripts[paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate]
        duplicate.granges.2$exon_number <- as.character(as.numeric(duplicate.granges.2$exon_number) + 0.5)
        duplicate.granges.2$exon_id <- paste0(duplicate.granges.2$exon_id, "_2")

        gtfTranscripts <- c(gtfTranscripts[which(!(paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate))],
                           duplicate.granges.1, duplicate.granges.2)
    }

    splitReturn <- list(ranges=gtfTranscripts, between=betweenExons)
    return(splitReturn)
}


#' annotate event coordinates with exon/transcript/gene ids.
#' @param eventCoords Granges with event coordinates
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param exons reference exons GRanges
#' @return eventCoords Granges with annotations
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
annotateEventCoords <- function(eventCoords, betweenExons, exons){

    #eventCoords <- GRanges(seqnames=input$chr, ranges=IRanges(start=input$exonStart_0base+1, end=input$exonEnd), strand=input$strand, id=input$ID)
    #seqlevelsStyle(eventCoords)=seqlevelsStyle(exons)[1]

    # create ranges for skipped exons
    eventCoords <- eventCoords[match(betweenExons$from_id, eventCoords$id)]
    eventCoords$exon_number <- NA
    eventCoords$transcript_id <- betweenExons$transcript_id
    eventCoords$transcript_type <- exons$transcript_type[match(eventCoords$transcript_id, exons$transcript_id)]
    eventCoords$gene_id <- exons$gene_id[match(eventCoords$transcript_id, exons$transcript_id)]
    eventCoords$exon_id <- unlist(lapply(str_split(betweenExons$new_transcript_id, "[+]AS[ ]"), "[[", 2))
    eventCoords$exon_number <- NA
    mcols(eventCoords) <- cbind(mcols(eventCoords),
                                DataFrame(new_transcript_id=paste0(
                                    eventCoords$transcript_id,"+AS ",
                                    eventCoords$exon_id)))
    return(eventCoords)

}

#' reformat a exons GRanges with first/last annotations and force 'biotype' to 'type'
#' @param exons reference exons GRanges
#' @return reference exons GRanges
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
reformatExons = function(exons){

    # add first/last annotation (speeds up later steps)
    if(!("first_last" %in% colnames(mcols(exons)))){
        t <- as.data.frame(table(exons$transcript_id))
        exons$first_last <- NA
        exons$first_last[exons$exon_number == 1] <- "first"
        exons$first_last[exons$exon_number == t$Freq[match(exons$transcript_id, t$Var1)]] <- "last"
    }

    colnames(mcols(exons)) <- gsub("biotype", "type", colnames(mcols(exons)))
    return(exons)
}

#' Generate a transcripts Granges from an exons Granges
#' @param exons reference exons GRanges
#' @return reference transcripts GRanges
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
transcriptsFromExons = function(exons){

    exons.df <- data.frame(t_id=exons$transcript_id, start=start(exons), end=end(exons))
    txRanges <- aggregate(start ~ t_id, exons.df, min)
    txRanges$end <- aggregate(end ~ t_id, exons.df, max)[,2]

    transcripts <- exons[!duplicated(exons$transcript_id),]
    transcripts$exon_id <- NULL
    transcripts$exon_number <- NULL
    transcripts$first_last <- NULL

    start(transcripts) <- txRanges$start[match(transcripts$transcript_id, txRanges$t_id)]
    end(transcripts) <- txRanges$end[match(transcripts$transcript_id, txRanges$t_id)]

    return(transcripts)

}
#' Generate a introns Granges from an exons Granges
#' @param exons reference exons GRanges
#' @return reference introns GRanges
#' @export
#' @import methods
#' @importFrom dplyr lead
#' @importFrom dplyr left_join
#' @family rmats data processing
#' @author Beth Signal
#' @examples
exonsToIntrons = function(exons){

    exons_df <- as.data.frame(exons)
    exons_df <- exons_df[,c("seqnames", "start", "end", "strand", "transcript_id", "exon_number")]
    exons_df <- arrange(exons_df, transcript_id, start, end)

    exons_df$intron_start <- exons_df$end
    exons_df$intron_end <- dplyr::lead(exons_df$start)

    rm <- which(dplyr::lead(exons_df$transcript_id) != exons_df$transcript_id)

    exons_df <- exons_df[-rm,]
    exons_df <- exons_df[-nrow(exons_df),]
    min_exon_n <- aggregate(exon_number ~ transcript_id, exons_df, min)
    if(!all(min_exon_n$exon_number == 1)){
        exons_df$exon_number[exons_df$strand=="-"] <- as.numeric(exons_df$exon_number[exons_df$strand=="-"]) - 1
    }

    introns <- GRanges(seqnames=exons_df$seqnames, ranges=IRanges(start=exons_df$intron_start, end=exons_df$intron_end),
                      strand=exons_df$strand, transcript_id=exons_df$transcript_id, exon_number=exons_df$exon_number)
    m <- match(introns$transcript_id, exons$transcript_id)
    introns$gene_id <- exons$gene_id[m]
    introns$gene_name <- exons$gene_name[m]

    mcols(introns) <- DataFrame(dplyr::left_join(as.data.frame(mcols(introns)),  as.data.frame(mcols(exons))))
    return(introns)

}

#' Find overlaps where the start/end coordinates are the same
#' @param query a GRanges object
#' @param subject a GRanges object
#' @return Hits object
#' @export
#' @import methods
#' @family data processing
#' @author Beth Signal
#' @examples
findOverlaps.junc = function(query, subject, type=c("start", "end")){

    if(any(type == "start")){

        query.start <- query
        end(query.start) <- start(query.start)
        subject.start <- subject
        end(subject.start) <- start(subject.start)

        ol.start <- findOverlaps(query.start, subject.start, type="equal")
    }

    if(any(type == "end")){

        query.end <- query
        start(query.end) <- end(query.end)
        subject.end <- subject
        start(subject.end) <- end(subject.end)

        ol.end <- findOverlaps(query.end, subject.end, type="equal")
    }

    if("start" %in% type & "end" %in% type){
        ol.df <- rbind(as.data.frame(ol.start), as.data.frame(ol.end))
        ol.df <- arrange(ol.df, queryHits, subjectHits)
        ol <- Hits(from=ol.df$queryHits, to=ol.df$subjectHits, nLnode=nLnode(ol.start), nRnode=nRnode(ol.start), sort.by.query=TRUE)
    }else if("start" %in% type){
        ol <- ol.start
    }else if("end" %in% type){
        ol <- ol.end
    }

    return(ol)

}
