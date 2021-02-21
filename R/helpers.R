## commonly used RMATs functions

annotateOverlapRmats <- function(granges, exons, exon_number=1){
    ol <- as.data.frame(findOverlaps(granges, exons))
    ol$from_id <- granges$id[ol$queryHits]
    ol$transcript_id <- exons$transcript_id[ol$subjectHits]
    ol$exon_id <- exons$exon_id[ol$subjectHits]
    ol$exon_number <- exons$exon_number[ol$subjectHits]
    ol$event_id <- paste0(seqnames(granges)[ol$queryHits], ":", granges$event_id[ol$queryHits])
    #ol$event_id <- make.unique(ol$event_id, sep="_")
    colnames(ol)[match(c('exon_id', 'exon_number') ,colnames(ol))] <- paste0(c('exon_id', 'exon_number'), exon_number)
    return(ol)
}

duplicateReference <- function(betweenExons, exons){
    # transcripts containing the exon pairs
    transcripts <- as.data.frame(table(betweenExons$transcript_id))
    gtfTranscripts <- exons[exons$transcript_id %in% transcripts$Var1]

    m <- match(gtfTranscripts$transcript_id, betweenExons$transcript_id)

    mcols(gtfTranscripts) <- cbind(mcols(gtfTranscripts),
                                   DataFrame(new_transcript_id=paste0(
                                       gtfTranscripts$transcript_id,"+AS ",
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
            gtfTranscripts_add$transcript_id,"+AS ",
            betweenExons$from_id[needsDuplicated][m], "-",
            betweenExons$event_id[needsDuplicated][m])
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts_add)
        needsDuplicated <- which(!(betweenExons$new_transcript_id %in%
                                       gtfTranscripts$new_transcript_id))
    }

    return(gtfTranscripts)
}

# remove any exons inbetween up/downstream exons
betweenNumbers <- function(a, b){
    ab.range <- seq(as.numeric(a),as.numeric(b))
    ab.range <- ab.range[!(ab.range %in% c(a,b))]
    return(ab.range)
}

removeExonsBetween <- function(betweenExons, gtfTranscripts){
    remove <- apply(betweenExons[,c('exon_number1', 'exon_number2', 'new_transcript_id')], 1, function(x) paste0(x[3], " ", betweenNumbers(x[1], x[2])))
    gtfTranscripts.rm <- gtfTranscripts[which(!(paste0(gtfTranscripts$new_transcript_id, " ", gtfTranscripts$exon_number) %in% unlist(remove)))]
    return(gtfTranscripts.rm)
}


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

    splitReturn <- list(ranges=gtfTranscripts, between =betweenExons)
    return(splitReturn)
}


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

reformatExons = function(exons){

    # add first/last annotation (speeds up later steps)
    if(!("first_last" %in% colnames(mcols(exons)))){
        t <- as.data.frame(table(exons$transcript_id))
        exons$first_last <- NA
        exons$first_last[exons$exon_number == 1] <- "first"
        exons$first_last[exons$exon_number ==
                             t$Freq[match(exons$transcript_id, t$Var1)]] <- "last"
    }

    colnames(mcols(exons)) = gsub("biotype", "type", colnames(mcols(exons)))
    return(exons)
}

transcriptsFromExons = function(exons){

    exons.df = data.frame(t_id = exons$transcript_id, start = start(exons), end = end(exons))
    txRanges = aggregate(start ~ t_id, exons.df, min)
    txRanges$end = aggregate(end ~ t_id, exons.df, max)[,2]

    transcripts = exons[!duplicated(exons$transcript_id),]
    transcripts$exon_id = NULL
    transcripts$exon_number = NULL
    transcripts$first_last = NULL

    start(transcripts) = txRanges$start[match(transcripts$transcript_id, txRanges$t_id)]
    end(transcripts) = txRanges$end[match(transcripts$transcript_id, txRanges$t_id)]

    return(transcripts)

}
