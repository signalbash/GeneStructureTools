#' Find transcripts overlapping first/last exons and alter their transcription start site or transcription end site
#'
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param exons GRanges object made from a GTF containing exon coordinates
#' @param type type of Whippet event (TS/TE).
#' Note only one event type should be processed at a time.
#' @param unfilteredWDS unfiltered whippetDataSet generated from \code{readWhippetDataSet()}
#' Note that this should contain ALL TS/TE events, regardless of significance, as these are used to contruct the entire exon range.
#' @return GRanges object with transcripts containing alternative junctions.
#' @export
#' @importFrom rtracklayer import
#' @import GenomicRanges
#' @importFrom dplyr inner_join
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#' @family whippet splicing isoform creation
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","gencode.vM25.small.gtf", package="GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' whippetFiles <- system.file("extdata","whippet_small/",
#' package="GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#'
#' wds.TS <- filterWhippetEvents(wds, eventTypes="TS", psiDelta=0.1, probability=0.95)
#' alterTranscriptStartEnds(whippetDataSet=wds.TS, exons=exons, unfilteredWDS=wds, type="TS")
alterTranscriptStartEnds <- function(whippetDataSet,
                                     exons,
                                     unfilteredWDS,
                                     type=NA){

    if(is.na(type)[1]){
        stop("Please specify event type as TS or TE")
    }else{
        type <- toupper(type)
        if(type != "TS" & type != "TE"){
            stop("Please specify event type as TS or TE")
        }
    }

    colnames(mcols(exons))[which(colnames(mcols(exons)) == "transcript_biotype")] <- "transcript_type"
    # check all are TS
    whippetDataSet <- filterWhippetEvents(whippetDataSet,
                                          probability=0,
                                          psiDelta=0,
                                          eventTypes=type)

    #wds.t_events=filterWhippetEvents(wds, eventTypes="TS", probability=0, psiDelta=0, idList=diffSplicingResults(whippetDataSet)$gene)
    # filter out ALL ranges involved in each diff spliced cluster - these form a contigous range, which we need to characterise the true start/ends of the exon of interest
    allEventRanges <- filterWhippetEvents(unfilteredWDS, eventTypes=type, probability=0, psiDelta=0, idList=diffSplicingResults(whippetDataSet)$gene)

    ## Group TS/E events into whole ranges ##
    allEventCoords <- coordinates(allEventRanges)

    ## find ranges which join together
    # make 2nd range with only the range starts so equal ends == exon starts at end +1
    coordEnds <- allEventCoords
    ranges(coordEnds) <- IRanges(start=start(coordEnds)-1, end=start(coordEnds)-1)
    ol.ends <- findOverlaps(allEventCoords, coordEnds, type="end")
    # make 2nd range with only the range ends so equal starts == exon ends at start-1
    coordStarts <- allEventCoords
    ranges(coordStarts) <- IRanges(start=end(coordStarts)+1, end=end(coordStarts)+1)
    ol.starts <- findOverlaps(allEventCoords, coordStarts, type="start")

    ol.coords <- rbind(as.data.frame(ol.ends), as.data.frame(ol.starts))

    # find all overlapping combinations (as in above will only be pairs)
    combinations <- unique(apply(ol.coords, 1, function(x) paste(sort(x), collapse=',')))
    totalInCombinations <- as.data.frame(table(unlist(stringr::str_split(combinations, ","))))

    # loop through each with >1 combination
    while(any(totalInCombinations$Freq > 1)){
        # index number of the first instance in >1 combination
        index.n <- as.numeric(as.character(totalInCombinations[which(totalInCombinations$Freq > 1)[1],]$Var1))
        # index for those needed to recombine
        index.combine <- which(unlist(lapply(stringr::str_split(combinations, ","), function(x) any(x==index.n))))
        # paste together
        combinations.new <- paste(sort(as.numeric(unique(unlist(stringr::str_split(combinations[index.combine], ","))))), collapse=",")
        # replace old (non-joined) with new (joined)
        combinations <- c(combinations[-index.combine], combinations.new)
        totalInCombinations <- as.data.frame(table(unlist(stringr::str_split(combinations, ","))))
    }

    combinations <- combinations[order(as.numeric(unlist(lapply(stringr::str_split(combinations, ","), "[[", 1))))]
    # convert strings to vectors of indices
    combinations.index <- lapply(combinations, function(x) as.numeric(as.character(unlist(stringr::str_split(x, ",")))))

    combo.start <- unlist(lapply(combinations.index, function(x) min(start(allEventCoords[x]))))
    combo.end <- unlist(lapply(combinations.index, function(x) max(end(allEventCoords[x]))))
    combo.chrom <- unlist(lapply(combinations.index, function(x) as.character(seqnames(allEventCoords[x])[1])))
    combo.strand <- unlist(lapply(combinations.index, function(x) as.character(strand(allEventCoords[x])[1])))

    combo.replacementRanges <- paste0(combo.chrom, ":", combo.start, "-", combo.end)

    # give each diff splicing event a 'group name' which has the entire range that the event should cover
    # will need to remove all coord groups for which there is no contigous group (no CLUE WTF these ranges are supposed to represent, but they are in the minority)
    allEventCoords$group_name[unlist(combinations.index)] <- unlist(mapply(function(x, y) rep(x, length(y)), x=combo.replacementRanges, y=combinations.index))


    combinations.granges <- GRanges(seqnames=combo.chrom,
                                   ranges=IRanges(start=combo.start, end= combo.end),
                                   strand=combo.strand,
                                   group_name=combo.replacementRanges)

    # add in alternative up/down sites if not signif. for both directions
    ### TODO: add option to just compare to reference if no alt. up/down?
    signifEvents <- diffSplicingResults(whippetDataSet)
    signifEvents$group_name <- allEventCoords$group_name[match(signifEvents$coord, allEventCoords$id)]
    # will need to remove all coord groups for which there is no contigous group (no CLUE WTF these ranges are supposed to represent, but they are in the minority)
    rm <- which(is.na(signifEvents$group_name))
    if(length(rm) > 0){signifEvents <- signifEvents[-rm,]}

    signifEvents$direction <- ifelse(signifEvents$psi_delta > 0, "up", "down")
    signif.groups <- unique(signifEvents$group_name)
    # add in an alternative up/down regulated if one doesn't exist in the original filtering
    add.up <- signif.groups[which(!(signif.groups %in% signifEvents$group_name[signifEvents$direction == "up"]))]
    add.down <- signif.groups[which(!(signif.groups %in% signifEvents$group_name[signifEvents$direction == "down"]))]

    allSignifEvents <- diffSplicingResults(allEventRanges)
    allSignifEvents$group_name <- allEventCoords$group_name[match(allSignifEvents$coord, allEventCoords$id)]
    allSignifEvents <- allSignifEvents[which(!(allSignifEvents$coord %in% signifEvents$coord)),]
    allSignifEvents$direction <- ifelse(allSignifEvents$psi_delta > 0, "up", "down")
    # add event with next largest psi_delta (in the correct direction)
    allSignifEvents <- plyr::arrange(allSignifEvents, psi_delta)
    signifEvents <- plyr::rbind.fill(signifEvents, allSignifEvents[match(add.down, allSignifEvents$group_name),])
    allSignifEvents <- plyr::arrange(allSignifEvents, plyr::desc(psi_delta))
    signifEvents <- plyr::rbind.fill(signifEvents, allSignifEvents[match(add.up, allSignifEvents$group_name),])

    ## find/replace transcription start sites
    # rep overlaps for each up/down TSS

    ol.combo <- as.data.frame(GenomicRanges::findOverlaps(combinations.granges, exons))
    ol.combo <- cbind(ol.combo,
                      transcript_id=exons$transcript_id[ol.combo$subjectHits],
                      group_name=combinations.granges$group_name[ol.combo$queryHits])

    ol.combo <- dplyr::inner_join(ol.combo, signifEvents[,c('group_name', 'coord', "direction", 'strand')], by="group_name")
    ol.combo$new_transcript_id <- paste0(ol.combo$transcript_id, "+AS", type, ifelse(ol.combo$direction=="up", "U", "D"), " ", ol.combo$group_name)

    ## table of transcripts overlapping the site
    # tid: transcript id
    tids <- unique(ol.combo$transcript_id)

    ## all transcripts for structural altercations
    gtfTranscripts <- exons[exons$transcript_id %in% tids]
    mcols(gtfTranscripts) <-
        mcols(gtfTranscripts)[,c('gene_id','transcript_id',
                                 'transcript_type','exon_id',
                                 'exon_number')]
    m <- match(gtfTranscripts$transcript_id, ol.combo$transcript_id)
    # add new transcript id
    gtfTranscripts$new_transcript_id <- ol.combo$new_transcript_id[m]

    # duplicate core transcripts if needed
    needsDuplicated <- which(!(ol.combo$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts.add <- gtfTranscripts[
            gtfTranscripts$transcript_id %in%
                ol.combo$transcript_id[needsDuplicated]]
    }

    while(length(needsDuplicated) > 0){
        gtfTranscripts.add <- gtfTranscripts.add[
            gtfTranscripts.add$transcript_id %in%
                ol.combo$transcript_id[needsDuplicated]]
        m <- match(gtfTranscripts.add$transcript_id,
                   ol.combo$transcript_id[needsDuplicated])
        gtfTranscripts.add$new_transcript_id <- ol.combo$new_transcript_id[needsDuplicated][m]
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts.add)
        needsDuplicated <- which(!(ol.combo$new_transcript_id %in%
                                       gtfTranscripts$new_transcript_id))
    }

    gtfTranscripts$from <- unlist(lapply(stringr::str_split(
        gtfTranscripts$new_transcript_id, "AS[A-Z]* "),"[[",2))
    gtfTranscripts <- gtfTranscripts[order(gtfTranscripts$transcript_id,
                                           start(gtfTranscripts))]


    # find which exons to alter in the new GTF annotation
    replaceWith <- as.data.frame(GenomicRanges::findOverlaps(combinations.granges, gtfTranscripts))
    replaceWith$from_id <- combinations.granges$group_name[replaceWith$queryHits]
    replaceWith$to_id <- gtfTranscripts$from[replaceWith$subjectHits]
    replaceWith$new_transcript_id <- gtfTranscripts$new_transcript_id[replaceWith$subjectHits]
    replaceWith <- replaceWith[replaceWith$from_id == replaceWith$to_id,]
    replaceWith$strand <- as.character(strand(combinations.granges[replaceWith$queryHits]))


    # if(type == "TS"){
    #     rm=which((end(combinations.granges[replaceWith$queryHits]) != end(gtfTranscripts[replaceWith$subjectHits])) & replaceWith$strand == "+")
    #     if(length(rm)>0){replaceWith=replaceWith[-rm,]}
    #     rm=which((start(combinations.granges[replaceWith$queryHits]) != start(gtfTranscripts[replaceWith$subjectHits])) & replaceWith$strand == "-")
    #     if(length(rm)>0){replaceWith=replaceWith[-rm,]}
    # }else if(type == "TE"){
    #     rm=which((end(combinations.granges[replaceWith$queryHits]) != end(gtfTranscripts[replaceWith$subjectHits])) & replaceWith$strand == "-")
    #     if(length(rm)>0){replaceWith=replaceWith[-rm,]}
    #     rm=which((start(combinations.granges[replaceWith$queryHits]) != start(gtfTranscripts[replaceWith$subjectHits])) & replaceWith$strand == "+")
    #     if(length(rm)>0){replaceWith=replaceWith[-rm,]}
    # }

    signifEvents$new_id <- paste0("AS", type, ifelse(signifEvents$direction=="up", "U", "D"), " ", signifEvents$group_name)
    replaceWith$coord <- signifEvents$coord[match(unlist(lapply(stringr::str_split(replaceWith$new_transcript_id, "[+]"), "[[", 2)), signifEvents$new_id)]

    replaceWith$newstart <- unlist(lapply(stringr::str_split(lapply(stringr::str_split(replaceWith$coord, "[-]"), "[[", 1 ), "[:]"), "[[" , 2))
    replaceWith$newend <- unlist(lapply(stringr::str_split(replaceWith$coord, "[-]"), "[[", 2 ))


    if(type == "TS"){
        # check that the exon to be replaced is the first exon
        overlaps_n_exon <- data.frame(new_transcript_id=gtfTranscripts$new_transcript_id[replaceWith$subjectHits], exon_number=as.numeric(as.character(gtfTranscripts$exon_number[replaceWith$subjectHits])))
        cant_replace <- overlaps_n_exon$new_transcript_id[which(overlaps_n_exon$exon_number != 1)]

        # check that the new start doesn't occur after the end of the exon
        compare <- as.numeric(as.character(replaceWith$newstart[which(replaceWith$strand == "+")])) < end(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "+")])
        cant_replace <- c(cant_replace, replaceWith$new_transcript_id[which(replaceWith$strand == "+")][compare==FALSE])
        start(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "+")])[which(compare)] <- as.numeric(as.character(replaceWith$newstart[which(replaceWith$strand == "+")]))[which(compare)]

        # check that the new end doesn't occur before the start of the exon
        compare <- as.numeric(as.character(replaceWith$newend[which(replaceWith$strand == "-")])) > start(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "-")])
        cant_replace <- c(cant_replace, replaceWith$new_transcript_id[which(replaceWith$strand == "-")][compare==FALSE])
        end(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "-")])[which(compare)] <- as.numeric(as.character(replaceWith$newend[which(replaceWith$strand == "-")]))[which(compare)]

        # remove any transcripts which were unable to be altered
        gtfTranscripts <- gtfTranscripts[which(!(gsub("TS[UD]", "TS", gtfTranscripts$new_transcript_id) %in% gsub("TS[UD]", "TS", cant_replace)))]
        gtfTranscripts$set <- ifelse(unlist(lapply(stringr::str_split(lapply(stringr::str_split(gtfTranscripts$new_transcript_id, "[+]AS"), "[[", 2), "[ ]"), "[[", 1))=="TSU", "tss_upregulated", "tss_downregulated")

    }else if(type == "TE"){

        # check that the exon to be replaced is the terminal exon
        n_exons <- as.data.frame(table(gtfTranscripts$new_transcript_id))
        overlaps_n_exon <- data.frame(new_transcript_id=gtfTranscripts$new_transcript_id[replaceWith$subjectHits], exon_number=as.numeric(as.character(gtfTranscripts$exon_number[replaceWith$subjectHits])))
        n_exons <- dplyr::left_join(n_exons, overlaps_n_exon, by=c('Var1'='new_transcript_id'))
        cant_replace <- n_exons$Var1[which(n_exons$exon_number != n_exons$Freq)]

        # check that the new end doesn't occur before the start of the exon
        compare <- as.numeric(as.character(replaceWith$newend[which(replaceWith$strand == "+")])) > start(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "+")])
        cant_replace <- c(cant_replace, replaceWith$new_transcript_id[which(replaceWith$strand == "+")][compare==FALSE])
        end(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "+")])[which(compare)] <- as.numeric(as.character(replaceWith$newend[which(replaceWith$strand == "+")]))[which(compare)]

        # check that the new start doesn't occur after the end of the exon
        compare <- as.numeric(as.character(replaceWith$newstart[which(replaceWith$strand == "-")])) < end(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "-")])
        cant_replace <- c(cant_replace, replaceWith$new_transcript_id[which(replaceWith$strand == "-")][compare==FALSE])
        start(gtfTranscripts[replaceWith$subjectHits][which(replaceWith$strand == "-")])[which(compare)] <- as.numeric(as.character(replaceWith$newstart[which(replaceWith$strand == "-")]))[which(compare)]

        # remove any transcripts which were unable to be altered
        gtfTranscripts <- gtfTranscripts[which(!(gsub("TE[UD]", "TE", gtfTranscripts$new_transcript_id) %in% gsub("TE[UD]", "TE", cant_replace)))]
        gtfTranscripts$set <- ifelse(unlist(lapply(stringr::str_split(lapply(stringr::str_split(gtfTranscripts$new_transcript_id, "[+]AS"), "[[", 2), "[ ]"), "[[", 1))=="TEU", "tes_upregulated", "tes_downregulated")

    }

    gtfTranscripts$event_id <- ol.combo$coord[match(gtfTranscripts$new_transcript_id, ol.combo$new_transcript_id)]

    gtfTranscripts$from <- NULL
    gtfTranscripts$transcript_id <- gtfTranscripts$new_transcript_id
    gtfTranscripts$new_transcript_id <- NULL


    return(gtfTranscripts)
}
