#' Given the location of a whole retained exon, find transcripts which can splice out this exon
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param exons GRanges object made from a GTF containing exon coordinates
#' @param variableWidth How many nts overhang is allowed for finding matching exons
#' (default = 0, i.e. complete match)
#' @param findIntrons Find transcripts where the event occurs within the intron?
#' @param transcripts GRanges object made from a GTF containing transcript coordinates
#' (only required if findIntrons=TRUE)
#' @return data.frame with all overlapping exons
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @author Beth Signal
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#' wds <- filterWhippetEvents(wds)
#'
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' transcripts <- gtf[gtf$type=="transcript"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' wds.exonSkip <- filterWhippetEvents(wds, eventTypes="CE",psiDelta = 0.2)
#' exons.exonSkip <- findExonContainingTranscripts(wds.exonSkip, exons,
#' variableWidth=0, findIntrons=FALSE, transcripts)
findExonContainingTranscripts <- function(whippetDataSet,
                                          exons,
                                          variableWidth=0,
                                          findIntrons=FALSE,
                                          transcripts){
    # check all are CE
    whippetDataSet <- filterWhippetEvents(whippetDataSet,
                                                     probability = 0,
                                                     psiDelta = 0,
                                                     eventTypes="CE")

    eventCoords <- coordinates(whippetDataSet)

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

    # whole match
    if(variableWidth == 0){
        overlaps <- GenomicRanges::findOverlaps(eventCoords, exons, type="equal")
        gtf.equal <- exons[overlaps@to]
        gtf.equal$from <- overlaps@from
    }else{
        # find all overlaps
        overlaps <- GenomicRanges::findOverlaps(eventCoords, exons)
        gtf.equal <- exons[overlaps@to]
        gtf.equal$from <- overlaps@from
        startDiff <- abs(start(gtf.equal) - start(eventCoords[gtf.equal$from]))
        endDiff <- abs(end(gtf.equal) - end(eventCoords[gtf.equal$from]))
        totalDiff <- startDiff + endDiff
        # remove overlaps with change in start+end greater than specified
        keep <- which(totalDiff <= variableWidth)
        gtf.equal <- gtf.equal[keep]
    }

    gtf.equal$new_id <- paste(gtf.equal$transcript_id,gtf.equal$from, sep="_")
    #gtf.equal$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf.equal)),
    # paste0(transcript_id, "_",from))
    gtf.equal$exon_number <- as.numeric(gtf.equal$exon_number)
    skippedExons <- as.data.frame(GenomicRanges::mcols(gtf.equal)[,c('gene_id',
                                                                    'transcript_id',
                                                                    'transcript_type',
                                                                    'from',
                                                                    'exon_number')])
    skippedExons$alt_id <- eventCoords$id[skippedExons$from]
    skippedExons$from <- NULL
    skippedExons$start <- start(gtf.equal)
    skippedExons$end <- end(gtf.equal)
    skippedExons$overlaps <- "exon"


    if(findIntrons == TRUE){
        # overlaps a transcript (i.e. can overlap an intron)
        overlaps <- GenomicRanges::findOverlaps(eventCoords, transcripts)
        overlapsDF <- as.data.frame(overlaps)
        overlapsDF$from <- eventCoords$id[overlapsDF$queryHits]
        overlapsDF$to <- transcripts$transcript_id[overlapsDF$subjectHits]

        # annotate first/last exons
        # (takes ~ 2.5 sec, please do before running this function multiple times)
         if(!("first_last" %in% colnames(mcols(exons)))){
             t <- as.data.frame(table(exons$transcript_id))
             exons$first_last <- NA
             exons$first_last[exons$exon_number == 1] <- "first"
             exons$first_last[exons$exon_number == t$Freq[match(exons$transcript_id, t$Var1)]] <- "last"
         }

        # check that skipped exon doesn't overlap the first/last exon
        overlapsExons <- GenomicRanges::findOverlaps(
            eventCoords, exons[which(exons$first_last %in% c("first","last"))])
        removeTranscripts <- exons$transcript_id[
            which(exons$first_last %in% c("first","last"))][overlapsExons@to]
        overlapsDF <- overlapsDF[which(!(overlapsDF$to %in% removeTranscripts)),]

        # gtf with ALL exons where there is an intron overlap
        gtf.within <- exons[
            which(exons$transcript_id %in% overlapsDF$to)]
        gtf.within$from <- overlapsDF$from[
            match(gtf.within$transcript_id, overlapsDF$to)]
        gtf.within <- gtf.within[
            which(!(gtf.within$transcript_id %in% c(gtf.equal$transcript_id)))]

        if(length(gtf.within) > 0){

            # check for non-eact exon matches
            overlappingExons <- as.data.frame(findOverlaps(eventCoords, gtf.within))
            overlappingExons$from_id <- eventCoords$id[overlappingExons$queryHits]
            overlappingExons <- cbind(overlappingExons,
                                      as.data.frame(gtf.within[overlappingExons$subjectHits]))
            overlappingExons$to_id <- gtf.within$transcript_id[overlappingExons$subjectHits]


            gtf.within$new_id <- paste0(gtf.within$transcript_id, "_",gtf.within$from)
            gtf.within$exon_number <- 0

            overlappingIntrons <- as.data.frame(GenomicRanges::mcols(gtf.within)[,c('gene_id',
                                                                            'transcript_id',
                                                                            'transcript_type',
                                                                            'from',
                                                                            'exon_number')])


            overlappingIntrons <- overlappingIntrons[
                which(!(duplicated(paste0(overlappingIntrons$transcript_id, "-",
                                          overlappingIntrons$from)))),]
            overlappingIntrons$alt_id <- overlappingIntrons$from
            overlappingIntrons$from <- NULL
            overlappingIntrons$start <- start(eventCoords[match(overlappingIntrons$alt_id,
                                                                eventCoords$id)])
            overlappingIntrons$end <- end(eventCoords[match(overlappingIntrons$alt_id,
                                                            eventCoords$id)])
            overlappingIntrons$overlaps <- "intron"



            # non exact matches (i.e. that have a partial intron overlap)
            nonexact <- which(overlappingIntrons$transcript_id %in% overlappingExons$to_id)
            if(length(nonexact) > 0){
                overlappingIntrons$overlaps[nonexact] <- "nonexact"
                m <- match(overlappingIntrons$transcript_id, overlappingExons$transcript_id)
                # replace coordinates with known coordinates
                overlappingIntrons$start <- overlappingExons$start[m]
                overlappingIntrons$end <- overlappingExons$end[m]
                overlappingIntrons$exon_number <- overlappingExons$exon_number[m]
            }

            skippedExons <- rbind(skippedExons, overlappingIntrons)

        }
    }


    return(skippedExons)
}

#' Remove and include a skipped exon from the transcripts it overlaps
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param skippedExons data.frame generataed by findExonContainingTranscripts()
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param glueExons Join together exons that are not seperated by exons?
#' @param match what type of match replacement should be done?
#' exact: exact matches to the skipped event only, also removes any intron overlaps
#' skip: keep non-exact exon match coordinates in included sets, and skip them in skipped sets
#' replace: replace non-exact exon match coordinates with event coordinates in included sets,
#' and skip them in skipped sets
#' @return GRanges with transcripts skipping exons
#' @export
#' @import GenomicRanges
#' @importFrom plyr desc
#' @importFrom rtracklayer import
#' @author Beth Signal
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#' wds <- filterWhippetEvents(wds)
#'
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' transcripts <- gtf[gtf$type=="transcript"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' wds.exonSkip <- filterWhippetEvents(wds, eventTypes="CE",psiDelta = 0.2)
#' exons.exonSkip <- findExonContainingTranscripts(wds.exonSkip, exons,
#' variableWidth=0, findIntrons=FALSE, transcripts)
#' ExonSkippingTranscripts <- skipExonInTranscript(wds.exonSkip, exons.exonSkip, exons)
skipExonInTranscript <- function(whippetDataSet,
                                 skippedExons,
                                 exons,
                                 glueExons=TRUE,
                                 match="exact"){

    if(!(match %in% c("exact","skip","replace"))){
        message("match must be 'exact', 'skip', or 'replace'")
        message("using default match = 'exact'")
        match <- "exact"
    }
    # check all are CE
    whippetDataSet <- filterWhippetEvents(whippetDataSet,
                                                     probability = 0,
                                                     psiDelta = 0,
                                                     eventTypes="CE")

    eventCoords <- coordinates(whippetDataSet)

    # remove non-exact matches
    if(match == "exact"){
        m <- match(skippedExons$alt_id, eventCoords$id)
        keep <- which((skippedExons$start) == start(eventCoords)[m] &
                          (skippedExons$end == end(eventCoords)[m]) &
                          skippedExons$overlaps=="exon")
        skippedExons <- skippedExons[keep,]
    }

    eventCoords <- eventCoords[match(skippedExons$alt_id, eventCoords$id)]
    eventCoords$exon_number <- skippedExons$exon_number
    eventCoords$transcript_id <- skippedExons$transcript_id
    eventCoords$transcript_type <- skippedExons$transcript_type
    eventCoords$gene_id <- skippedExons$gene_id
    eventCoords$exon_id <- eventCoords$id
    eventCoords$overlaps <- skippedExons$overlaps

    # replace exon coordinates with the event coordinates
    if(match == "replace"){
        oldStarts <- start(eventCoords)
        oldEnds <- end(eventCoords)
    }
    start(eventCoords) <- skippedExons$start
    end(eventCoords) <- skippedExons$end


    # transcripts containing the exon
    transcripts <- as.data.frame(table(skippedExons$transcript_id))
    gtfTranscripts <- exons[exons$transcript_id %in% transcripts$Var1]
    m <- match(gtfTranscripts$transcript_id, eventCoords$transcript_id)
    mcols(gtfTranscripts) <- cbind(mcols(gtfTranscripts),
                                   DataFrame(new_transcript_id=paste0(
                                       gtfTranscripts$transcript_id,"+AS ",
                                       eventCoords$exon_id[m])))
    mcols(eventCoords) <- cbind(mcols(eventCoords),
                               DataFrame(new_transcript_id = paste0(
                                   eventCoords$transcript_id,"+AS ",
                                   eventCoords$exon_id)))

    mcols(gtfTranscripts) <- mcols(gtfTranscripts)[,c('gene_id','transcript_id',
                                                      'transcript_type','exon_id',
                                                      'exon_number', 'new_transcript_id')]
    mcols(eventCoords) <- mcols(eventCoords)[,c('gene_id','transcript_id',
                                                'transcript_type','exon_id',
                                                'exon_number','new_transcript_id',
                                                'overlaps')]

    needsDuplicated <- which(!(eventCoords$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts[gtfTranscripts$transcript_id %in%
                                                 eventCoords$transcript_id[needsDuplicated]]
    }

    while(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts_add[
            gtfTranscripts_add$transcript_id %in% eventCoords$transcript_id[needsDuplicated]]
        m <- match(gtfTranscripts_add$transcript_id, eventCoords$transcript_id[needsDuplicated])
        gtfTranscripts_add$new_transcript_id <- paste0(
            gtfTranscripts_add$transcript_id,"+AS ",
            eventCoords$exon_id[needsDuplicated][m])
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts_add)
        needsDuplicated <- which(!(eventCoords$new_transcript_id %in%
                                       gtfTranscripts$new_transcript_id))
    }

    exon_names <- with(as.data.frame(gtfTranscripts), paste0(seqnames, ":", start,"-",end))

    ol <- findOverlaps(gtfTranscripts, eventCoords, type="equal")
    ol <- as.data.frame(ol)
    ol$gtf_trans_id <- gtfTranscripts$new_transcript_id[ol$queryHits]
    ol$eventCoords_id <- eventCoords$new_transcript_id[ol$subjectHits]
    ol <- ol[ol$gtf_trans_id == ol$eventCoords_id,]

    # add/remove exons from skipped isoforms if introns are used
    if(match!="exact" & any(eventCoords$overlaps!="exon")){
        ol.var <- findOverlaps(gtfTranscripts, eventCoords[eventCoords$overlaps != "exon"])
        ol.var <- as.data.frame(ol.var)
        ol.var$gtf_trans_id <- gtfTranscripts$new_transcript_id[ol.var$queryHits]
        ol.var$eventCoords_id <- eventCoords$new_transcript_id[
            which(eventCoords$overlaps != "exon")][ol.var$subjectHits]
        ol.var <- ol.var[ol.var$gtf_trans_id == ol.var$eventCoords_id,]
        if(nrow(ol.var) > 0){
            ol <- rbind(ol, ol.var)
        }
    }

    # remove the skipped exon
    rm <- unique(ol$queryHits)
    gtfTranscripts.rm <- gtfTranscripts[-rm]

    mcols(gtfTranscripts.rm) <- mcols(gtfTranscripts.rm)[,c('gene_id','new_transcript_id',
                                                            'transcript_type','exon_id',
                                                            'exon_number')]
    colnames(mcols(gtfTranscripts.rm))[2] <- "transcript_id"

    mcols(eventCoords) <- mcols(eventCoords)[,c('gene_id','new_transcript_id',
                                                'transcript_type','exon_id',
                                                'exon_number','overlaps')]
    colnames(mcols(eventCoords))[2] <- "transcript_id"

    gtfTranscripts.rm$exon_number <- as.numeric(gtfTranscripts.rm$exon_number)
    order <- order(gtfTranscripts.rm$transcript_id, gtfTranscripts.rm$exon_number)
    gtfTranscripts.rm <- gtfTranscripts.rm[order]
    gtfTranscripts.rm$overlaps <- eventCoords$overlaps[
        match(gtfTranscripts.rm$transcript_id, eventCoords$new_transcript_id)]

    # add skipped exon back in
    gtfTranscripts.withExon <- gtfTranscripts.rm
    # replace with event coordinates
    if(match == "replace"){
        start(eventCoords) <- oldStarts
        end(eventCoords) <- oldEnds
    }
    gtfTranscripts.withExon <- c(gtfTranscripts.withExon, eventCoords)
    gtfTranscripts.withExon <- reorderExonNumbers(gtfTranscripts.withExon)

    gtfTranscripts.rm$set <- "skipped_exon"
    gtfTranscripts.withExon$set <- "included_exon"

    gtfTranscripts.rm <- c(gtfTranscripts.rm, gtfTranscripts.withExon)

    # rename skipped/included isoforms
    gtfTranscripts.rm$transcript_id[which(gtfTranscripts.rm$set=="skipped_exon")] <-
        gsub("AS", "ASSE", gtfTranscripts.rm$transcript_id[
            which(gtfTranscripts.rm$set=="skipped_exon")])
    gtfTranscripts.rm$transcript_id[which(gtfTranscripts.rm$set=="included_exon")] <-
        gsub("AS", "ASIE", gtfTranscripts.rm$transcript_id[
            which(gtfTranscripts.rm$set=="included_exon")])

    #join together exons that are not seperated by an exon
    if(glueExons==TRUE){

        # split pos/neg ordering
        gtfTranscripts.rm.neg <-
            gtfTranscripts.rm[which(as.logical(strand(gtfTranscripts.rm) == "-"))]
        order <- order(gtfTranscripts.rm.neg$transcript_id,
                       plyr::desc(gtfTranscripts.rm.neg$exon_number))
        gtfTranscripts.rm.neg <- gtfTranscripts.rm.neg[order]

        gtfTranscripts.rm.pos <-
            gtfTranscripts.rm[which(as.logical(strand(gtfTranscripts.rm) == "+"))]
        gtfTranscripts.rm <- c(gtfTranscripts.rm.pos,gtfTranscripts.rm.neg)

        #extend starts <---<---<---
        w <- which(end(ranges(gtfTranscripts.rm))[-length(gtfTranscripts.rm)] ==
                       start(ranges(gtfTranscripts.rm[-1])))
        gtfTranscripts.rm <- gtfTranscripts.rm
        GenomicRanges::start(GenomicRanges::ranges(gtfTranscripts.rm))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(gtfTranscripts.rm))[w]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtfTranscripts.rm, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to - 1]
        rm <- unique(overlaps@from)
        if(length(rm) > 0){
            gtfTranscripts.rm <- gtfTranscripts.rm[-rm]
        }

        #extend ends --->--->--->
        overlaps <- findOverlaps(gtfTranscripts.rm)
        overlaps <- overlaps[overlaps@from==overlaps@to +1]

        GenomicRanges::end(GenomicRanges::ranges(gtfTranscripts.rm))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(gtfTranscripts.rm))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtfTranscripts.rm, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtfTranscripts.rm <- gtfTranscripts.rm[-rm]
        }

        #order <- order(gtfTranscripts.rm$transcript_id, gtfTranscripts.rm$exon_number)
        #gtfTranscripts.rm <- gtfTranscripts.rm[order]

    }
    gtfTranscripts.rm$whippet_id <- unlist(lapply(stringr::str_split(
        gtfTranscripts.rm$transcript_id, " "),"[[",2))
    gtfTranscripts.rm$overlaps <- NULL
    return(gtfTranscripts.rm)
}

#' Reorder the exon numbers in a gtf annotation
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param by what column are the transcripts grouped by?
#' @return The same input GRanges, but with exon numbers reordered.
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' exons <- reorderExonNumbers(exons)
reorderExonNumbers <- function(exons, by="transcript_id"){
    n <- which(colnames(mcols(exons)) == by)

    order <- order(mcols(exons)[,n], start(exons))

    exons <- exons[order]

    transcriptTable <- as.data.frame(table(mcols(exons)[,n]))
    transcriptTable$strand <- as.character(strand(exons[match(transcriptTable$Var1,
                                                                  mcols(exons)[,n])]))

    exons$exon_number <- as.numeric(unlist(apply(transcriptTable, 1,
                                          function(x) if(x[3] == "+"){
                                              c(1:(x[2]))}else{c((x[2]:1))})))
    return(exons)
}

