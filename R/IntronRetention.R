#' Given the location of a whole retained intron, find transcripts which splice out this intron
#' @param input whippetDataSet generated from \code{readWhippetDataSet()} or a Granges of intron coordinates
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param match what type of matching to perform? exact = only exons which bound the intron exactly,
#' introns = any exon pairs which overlap the intron,
#' all = any exon pairs AND single exons which overlap the intron
#' @return data.frame with all flanking exon pairs
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family whippet splicing isoform creation
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
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' wds.intronRetention <- filterWhippetEvents(wds, eventTypes="RI")
#' exons.intronRetention <- findIntronContainingTranscripts(input=wds.intronRetention, exons)
#'
#' exonsFromGRanges <- exons[exons$transcript_id=="ENSMUST00000139129.8" &
#' exons$exon_number %in% c(3,4)]
#' intronFromGRanges <- exonsFromGRanges[1]
#' GenomicRanges::start(intronFromGRanges) <-
#' GenomicRanges::end(exonsFromGRanges[exonsFromGRanges$exon_number==3])
#' GenomicRanges::end(intronFromGRanges) <-
#' GenomicRanges::start(exonsFromGRanges[exonsFromGRanges$exon_number==4])
#' exons.intronRetention <- findIntronContainingTranscripts(intronFromGRanges, exons)
findIntronContainingTranscripts <- function(input,
                                            exons,
                                            match="exact"){
    moved <- FALSE

    if(class(input)=="whippetDataSet"){
        whippetDataSet <- filterWhippetEvents(input,
                                              probability = 0,
                                              psiDelta = 0,
                                              eventTypes="RI")

        eventCoords <- coordinates(whippetDataSet)

    }else if(class(input) == "GRanges"){
        eventCoords <- input
        if(!("id" %in% names(mcols(eventCoords))) &
           "exon_id" %in% names(mcols(eventCoords))){
            eventCoords$id <- eventCoords$exon_id
        }else{
            stop("please specify \"id\" or \"exon_id\" in the input")
        }
    }
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

    overlaps <- GenomicRanges::findOverlaps(rangeRI.start, exons, type="end")

    # catch if intron coords dont overlap the 1nt exon start/end
    # less than 1/3 to catch stuff
    if(length(overlaps) < length(rangeRI.start)/3){
        GenomicRanges::start(rangeRI.start) <-
            GenomicRanges::start(rangeRI.start) -1
        GenomicRanges::end(rangeRI.start) <-
            GenomicRanges::start(rangeRI.start)
        overlaps <- GenomicRanges::findOverlaps(rangeRI.start,
                                                exons, type="end")
        # fix original
        GenomicRanges::start(eventCoords) <-
            GenomicRanges::start(eventCoords) -1
        GenomicRanges::end(eventCoords) <-
            GenomicRanges::end(eventCoords) +1
        moved <- TRUE
    }

    gtf.fromA <- exons[overlaps@to]
    gtf.fromA$from <- overlaps@from
    gtf.fromA$new_id <- paste(gtf.fromA$transcript_id, gtf.fromA$from, sep="_")

    # end of intron // start of exon b
    rangeRI.end <- eventCoords
    GenomicRanges::start(rangeRI.end) <- GenomicRanges::end(rangeRI.end)

    overlaps <- GenomicRanges::findOverlaps(rangeRI.end, exons, type="start")
    if(length(overlaps) < length(rangeRI.end)/3){
        GenomicRanges::end(rangeRI.end) <-
            GenomicRanges::end(rangeRI.end) +1
        GenomicRanges::start(rangeRI.end) <-
            GenomicRanges::end(rangeRI.end)
        overlaps <- GenomicRanges::findOverlaps(rangeRI.end,
                                                exons, type="start")
        moved <- TRUE
    }

    gtf.toA <- exons[overlaps@to]
    gtf.toA$from <- overlaps@from
    gtf.toA$new_id <- paste0(gtf.toA$transcript_id, "_",gtf.toA$from)

    keep.from <- which(gtf.fromA@elementMetadata$new_id %in%
                           gtf.toA@elementMetadata$new_id)
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
        mcols(gtf.toA) <- cbind(mcols(gtf.toA),
                                S4Vectors::DataFrame(overlaps=character()))
    }
    if(match == "introns" | match=="all"){
        # non-exact overlaps
        ol <- findOverlaps(eventCoords, exons)
        gtf.overlaps <- exons[ol@to]
        gtf.overlaps$from <- ol@from
        gtf.overlaps$new_id <- paste0(gtf.overlaps$transcript_id,
                                      "_",gtf.overlaps$from)
        gtf.overlaps$exon_number <- as.numeric(gtf.overlaps$exon_number)
        # remove perfect match pairs
        gtf.overlaps <- gtf.overlaps[which(!(gtf.overlaps$new_id %in%
                                                 gtf.toA$new_id))]

        transcriptTable <- as.data.frame(table(gtf.overlaps$new_id))
        gtf.overlapsPairs <-
            gtf.overlaps[gtf.overlaps$new_id %in%
                             transcriptTable$Var1[transcriptTable$Freq == 2]]
        startPair <-
            gtf.overlapsPairs$new_id[
                start(gtf.overlapsPairs) <
                    start(eventCoords[gtf.overlapsPairs$from])]
        endPair <-
            gtf.overlapsPairs$new_id[
                end(gtf.overlapsPairs) >
                    end(eventCoords[gtf.overlapsPairs$from])]
        keep <- startPair[startPair %in% endPair]
        gtf.overlapsPairs <-
            gtf.overlapsPairs[gtf.overlapsPairs$new_id %in% keep]

        if(length(gtf.overlapsPairs) > 0){
            exon_numbers <- aggregate(exon_number ~ new_id,
                                      mcols(gtf.overlapsPairs), mean)
            gtf.overlapsPairs$intron_exon_number <-
                exon_numbers$exon_number[match(gtf.overlapsPairs$new_id,
                                               exon_numbers$new_id)]
            gtf.overlapsPairs <-
                gtf.overlapsPairs[!duplicated(gtf.overlapsPairs$new_id)]
            gtf.overlapsPairs$from_exon_number <-
                gtf.overlapsPairs$intron_exon_number - 0.5
            gtf.overlapsPairs$to_exon_number <-
                gtf.overlapsPairs$intron_exon_number + 0.5
            gtf.overlapsPairs$overlaps <- "nonexact"
            mcols(gtf.overlapsPairs) <-
                mcols(gtf.overlapsPairs[,match(colnames(mcols(gtf.toA)),
                                               colnames(mcols(
                                                   gtf.overlapsPairs)))])

            gtf.toA <- c(gtf.toA, gtf.overlapsPairs)
        }
    }
    if(match=="all"){
        # overlaps an exon
        gtf.overlaps <- gtf.overlaps[which(!(gtf.overlaps$new_id %in%
                                                 gtf.overlapsPairs$new_id))]
        keep <- which(start(gtf.overlaps) <=
                          start(eventCoords[gtf.overlaps$from]) &
                          end(gtf.overlaps) >=
                          end(eventCoords[gtf.overlaps$from]))
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
            as.data.frame(GenomicRanges::mcols(
                gtf.toA)[,c('gene_id','transcript_id',
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
#' @param flankingExons data.frame generataed by findIntronContainingTranscripts()
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param match what type of match replacement should be done?
#' exact: exact matches to the intron only
#' retain: keep non-exact intron match coordinates in spliced sets, and retain them in retained sets
#' replace: replace non-exact intron match coordinates with event coordinates in spliced sets,
#' and retain in retained sets
#' @param glueExons Join together exons that are not seperated by introns?
#' @return GRanges with transcripts containing retained introns
#' @export
#' @import GenomicRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom plyr desc
#' @importFrom rtracklayer import
#' @family whippet splicing isoform creation
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
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' wds.intronRetention <- filterWhippetEvents(wds, eventTypes="RI")
#' exons.intronRetention <- findIntronContainingTranscripts(wds.intronRetention, exons)
#' IntronRetentionTranscripts <- addIntronInTranscript(exons.intronRetention, exons,
#' whippetDataSet=wds.intronRetention)
#'
#' exonsFromGRanges <- exons[exons$transcript_id=="ENSMUST00000139129.8" &
#' exons$exon_number %in% c(3,4)]
#' intronFromGRanges <- exonsFromGRanges[1]
#' GenomicRanges::start(intronFromGRanges) <-
#' GenomicRanges::end(exonsFromGRanges[exonsFromGRanges$exon_number==3])
#' GenomicRanges::end(intronFromGRanges) <-
#' GenomicRanges::start(exonsFromGRanges[exonsFromGRanges$exon_number==4])
#' exons.intronRetention <- findIntronContainingTranscripts(intronFromGRanges, exons)
#'
#' IntronRetentionTranscripts <-
#' addIntronInTranscript(exons.intronRetention, exons, match="retain")
addIntronInTranscript <- function(flankingExons,
                                  exons,
                                  whippetDataSet = NULL,
                                  match="exact",
                                  glueExons=TRUE){


    if(!(match %in% c("exact","retain","replace"))){
        message("match must be 'exact', 'retain', or 'replace'")

        if(!is.null(whippetDataSet)){
            message("using match = 'exact'")
            match <- "exact"
        }else{
            message("using match = 'retain'")
            match <- "retain"
        }
    }

    if(is.null(whippetDataSet) & match=="exact"){
        message("cannot use match = 'exact' without a whippetDataSet")
        message("using match = 'retain'")
        match <- "retain"
    }

    if(!is.null(whippetDataSet)){
        whippetDataSet <- filterWhippetEvents(whippetDataSet,
                                              probability = 0,
                                              psiDelta = 0,
                                              eventTypes="RI")

        eventCoords <- coordinates(whippetDataSet)

        if(match == "exact"){
            keep <- which(flankingExons$overlaps == "intron")
            flankingExons <- flankingExons[keep,]
            eventCoords <- eventCoords[eventCoords$id %in% flankingExons$from]
        }

        move <- which(flankingExons$moved == TRUE)
        move.RIindex <- which(eventCoords$id %in% flankingExons$from[move])
        GenomicRanges::start(eventCoords)[move.RIindex] <-
            GenomicRanges::start(eventCoords)[move.RIindex] -1
        GenomicRanges::end(eventCoords)[move.RIindex] <-
            GenomicRanges::end(eventCoords)[move.RIindex] +1

    }else{
        m1 <- match(paste0(flankingExons$transcript_id, flankingExons$from_exon_number),
                   paste0(exons$transcript_id, exons$exon_number))
        m2 <- match(paste0(flankingExons$transcript_id, flankingExons$to_exon_number),
                    paste0(exons$transcript_id, exons$exon_number))

        eventCoords <- GRanges(seqnames=seqnames(exons[c(m1)]),
                               ranges = IRanges(start=end(exons[c(m1)]),
                                                end=start(exons[c(m2)])),
                               strand=strand(exons[c(m1)]),
                               id=flankingExons$from)
    }

    eventCoords <- eventCoords[match(flankingExons$from, eventCoords$id)]
    eventCoords$exon_number <- flankingExons$intron_exon_number
    eventCoords$transcript_id <- flankingExons$transcript_id
    eventCoords$transcript_type <- flankingExons$transcript_type
    eventCoords$gene_id <- flankingExons$gene_id
    eventCoords$exon_id <- eventCoords$id

    transcripts <- as.data.frame(table(flankingExons$transcript_id))
    gtfTranscripts <- exons[exons$transcript_id %in% transcripts$Var1]
    m <- match(gtfTranscripts$transcript_id, eventCoords$transcript_id)
    mcols(gtfTranscripts) <-
        cbind(mcols(gtfTranscripts),
              S4Vectors::DataFrame(
                  new_transcript_id=paste0(gtfTranscripts$transcript_id,
                                           "+AS ", eventCoords$exon_id[m])))
    mcols(eventCoords) <-
        cbind(mcols(eventCoords),
              S4Vectors::DataFrame(
                  new_transcript_id = paste0(eventCoords$transcript_id,
                                             "+AS ", eventCoords$exon_id)))
    flankingExons$new_transcript_id <- paste0(flankingExons$transcript_id,
                                              "+AS ", flankingExons$from)

    mcols(gtfTranscripts) <- mcols(
        gtfTranscripts)[,c('gene_id','transcript_id',
                           'transcript_type','exon_id',
                           'exon_number','new_transcript_id')]
    mcols(eventCoords) <- mcols(
        eventCoords)[,c('gene_id','transcript_id',
                        'transcript_type','exon_id',
                        'exon_number','new_transcript_id')]

    needsDuplicated <- which(!(eventCoords$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts.add <- gtfTranscripts[
            gtfTranscripts$transcript_id %in%
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

    ##########################
    # fix starts/ends of introns
    gtfTranscripts$new_id_ex <-
        paste0(gtfTranscripts$new_transcript_id, "_",
               gtfTranscripts$exon_number)
    flankingExons$from_ex <-
        paste0(flankingExons$new_transcript_id, "_",
               flankingExons$from_exon_number)
    flankingExons$to_ex <-
        paste0(flankingExons$new_transcript_id, "_",
               flankingExons$to_exon_number)

    if(match == "replace"){
        wi <- which(flankingExons$overlaps == "nonexact")
        wi.pos <-
            wi[as.character(strand(gtfTranscripts[
                match(flankingExons$from_ex[wi],
                      gtfTranscripts$new_id_ex)])) == "+"]
        wi.neg <-
            wi[as.character(strand(gtfTranscripts[
                match(flankingExons$from_ex[wi],
                      gtfTranscripts$new_id_ex)])) == "-"]

        if(length(wi.pos) > 0){
            end(gtfTranscripts)[match(flankingExons$from_ex[wi.pos],
                                      gtfTranscripts$new_id_ex)] <-
                start(eventCoords)[
                    match(flankingExons$new_transcript_id[wi.pos],
                          eventCoords$new_transcript_id)]
            start(gtfTranscripts)[match(flankingExons$to_ex[wi.pos],
                                        gtfTranscripts$new_id_ex)] <-
                end(eventCoords)[match(flankingExons$new_transcript_id[wi.pos],
                                       eventCoords$new_transcript_id)]
        }
        if(length(wi.neg) > 0){
            end(gtfTranscripts)[match(flankingExons$to_ex[wi.neg],
                                      gtfTranscripts$new_id_ex)] <-
                start(eventCoords)[
                    match(flankingExons$new_transcript_id[wi.neg],
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
        end(start.replacement) <- start(eventCoords)[
            match(flankingExons$new_transcript_id[we],
                  eventCoords$new_transcript_id)]
        end.replacement <- gtfTranscripts[replace]
        start(end.replacement) <- end(eventCoords)[
            match(flankingExons$new_transcript_id[we],
                  eventCoords$new_transcript_id)]
        gtfTranscripts <- c(gtfTranscripts[-replace],
                            start.replacement, end.replacement)
    }
    gtfTranscripts <- reorderExonNumbers(gtfTranscripts,
                                         by="new_transcript_id")
    gtfTranscripts$new_id_ex <- NULL

    gtfTranscripts.withIntron <- c(gtfTranscripts, eventCoords)
    gtfTranscripts.withIntron <- reorderExonNumbers(gtfTranscripts.withIntron,
                                                    by="new_transcript_id")

    mcols(gtfTranscripts.withIntron) <-
        mcols(gtfTranscripts.withIntron)[,c(
            'gene_id','new_transcript_id',
            'transcript_type','exon_id','exon_number')]
    colnames(mcols(gtfTranscripts.withIntron))[2] <- "transcript_id"

    mcols(gtfTranscripts) <-
        mcols(gtfTranscripts)[,c('gene_id','new_transcript_id',
                                 'transcript_type','exon_id','exon_number')]
    colnames(mcols(gtfTranscripts))[2] <- "transcript_id"

    gtfTranscripts.withIntron$exon_number <-
        as.numeric(gtfTranscripts.withIntron$exon_number)
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
        GenomicRanges::start(GenomicRanges::ranges(
            gtfTranscripts.withIntron))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(
                gtfTranscripts.withIntron))[w]

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

        GenomicRanges::end(GenomicRanges::ranges(
            gtfTranscripts.withIntron))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(
                gtfTranscripts.withIntron))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtfTranscripts.withIntron, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtfTranscripts.withIntron <- gtfTranscripts.withIntron[-rm]
        }

        #order <- order(gtfTranscripts.withIntron$transcript_id,
        #gtfTranscripts.withIntron$exon_number)
        #gtfTranscripts.withIntron <- gtfTranscripts.withIntron[order]

    }

    gtfTranscripts.withIntron$whippet_id <- unlist(lapply(stringr::str_split(
        gtfTranscripts.withIntron$transcript_id, " "),"[[",2))
    gtfTranscripts.withIntron$overlaps <- NULL

    return(gtfTranscripts.withIntron)
}

