#' Filter out significant events from a whippet diff comparison
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param probability minimum probability required to call event as significant
#' @param psiDelta minimum change in psi required to call an event as significant
#' @param eventTypes which event type to filter for? default = \code{"all"}
#' @param idList (optional) list of gene ids to filter for
#' @param minCounts minumum number of counts for all replicates
#' in at least one condition to call an event as significant
#' @param medianCounts median count for all replicates
#' in at least one condition to call an event as significant
#' @param sampleTable data.frame with sample names and conditions.
#' Only needed if filtering with counts.
#' @return filtered whippet differential comparison data.frame
#' @export
#' @importFrom stats median
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
#' wds <- filterWhippetEvents(wds)
filterWhippetEvents <- function(whippetDataSet,
                                probability = 0.95,
                                psiDelta = 0.1,
                                eventTypes = "all",
                                idList = NA,
                                minCounts = NA,
                                medianCounts = NA,
                                sampleTable){

    if(eventTypes[1] == "all"){
        eventTypes <- unique(diffSplicingResults(whippetDataSet)$type)
    }

    if(is.na(idList[1])){
        significantEventsIndex <-
            which(diffSplicingResults(whippetDataSet)$probability >
                      probability &
                      abs(diffSplicingResults(whippetDataSet)$psi_delta) >
                      psiDelta &
                      diffSplicingResults(whippetDataSet)$type %in%
                      eventTypes)
    }else{
        significantEventsIndex <-
            which(diffSplicingResults(whippetDataSet)$probability >
                      probability &
                      abs(diffSplicingResults(whippetDataSet)$psi_delta) >
                      psiDelta &
                      diffSplicingResults(whippetDataSet)$type %in%
                      eventTypes &
                      diffSplicingResults(whippetDataSet)$gene %in% idList)
    }

    slot(whippetDataSet, "diffSplicingResults") <-
        diffSplicingResults(whippetDataSet)[significantEventsIndex,]

    m <- match(diffSplicingResults(whippetDataSet)$coord,
               coordinates(whippetDataSet)$id)

    slot(whippetDataSet, "coordinates") <-
        coordinates(whippetDataSet)[unique(m),]

    keep <- which(readCounts(whippetDataSet)$Gene %in%
                      diffSplicingResults(whippetDataSet)$gene)

    slot(whippetDataSet, "readCounts") <-
        readCounts(whippetDataSet)[keep,]

    keep <- which(junctions(whippetDataSet)$gene %in%
                      diffSplicingResults(whippetDataSet)$gene)

    slot(whippetDataSet, "junctions") <-
        junctions(whippetDataSet)[keep,]


    #filter by read counts?
    if(nrow(slot(whippetDataSet, "readCounts")) > 0 & (!is.na(minCounts) |
                                                       !is.na(medianCounts))){
        m <- match(diffSplicingResults(whippetDataSet)$unique_name,
                   readCounts(whippetDataSet)$unique_name)

        if(!all(sampleTable$sample %in% colnames(readCounts(whippetDataSet)))){
            notFound <- sampleTable$sample[which(!(sampleTable$sample %in% colnames(readCounts(whippetDataSet))))]
            message(paste0("Can't find ", paste0(notFound, collapse = ", "), " in the read counts file column names"))
            maybe <- colnames(readCounts(whippetDataSet))
            maybe <- maybe[which(!(tolower(maybe) %in% c('gene', "node", "coord", "strand", "type", "na_count", "unique_name")))]
            message(paste0("Did you mean to specifiy 'sample' in sampleTable as ", paste0(maybe, collapse = ", "), " ?"))
            message("skipping filtering based on read counts...")
        }else{

            n1 <- match(sampleTable$sample[
                sampleTable$condition %in%
                    unique(diffSplicingResults(whippetDataSet)$condition_1)],
                colnames(readCounts(whippetDataSet)))
            n2 <- match(sampleTable$sample[
                sampleTable$condition %in%
                    unique(diffSplicingResults(whippetDataSet)$condition_2)],
                colnames(readCounts(whippetDataSet)))

            diffSplicingResultsTemp <- diffSplicingResults(whippetDataSet)

            diffSplicingResultsTemp$condition_1_counts <-
                apply(readCounts(whippetDataSet)[m,n1], 1, stats::median)
            diffSplicingResultsTemp$condition_2_counts <-
                apply(readCounts(whippetDataSet)[m,n2], 1, stats::median)


            if(!is.na(minCounts)){
                keep <- which(apply(readCounts(whippetDataSet)[m,n1], 1,
                                    function(x) all(x > minCounts)) |
                                  apply(readCounts(whippetDataSet)[m,n2], 1,
                                        function(x) all(x > minCounts)))
                slot(whippetDataSet, "diffSplicingResults") <-
                    diffSplicingResultsTemp[keep,]
            }

            if(!is.na(medianCounts)){
                keep <- which(diffSplicingResultsTemp$condition_1_counts >=
                                  medianCounts |
                                  diffSplicingResultsTemp$condition_2_counts >=
                                  medianCounts)
                slot(whippetDataSet, "diffSplicingResults") <-
                    diffSplicingResultsTemp[keep,]
            }
        }

    }

    return(whippetDataSet)

}
#' Compare open reading frames for whippet differentially spliced events
#' @param whippetDataSet whippetDataSet generated from \code{readWhippetDataSet()}
#' @param unfilteredWDS unfiltered whippetDataSet generated from \code{readWhippetDataSet()}
#' Note that this should contain ALL TS/TE events, regardless of significance, as these are used to construct the entire exon range.
#' @param eventTypes which event type to filter for? default = "all"
#' @param exons GRanges gtf annotation of exons
#' @param transcripts GRanges gtf annotation of transcripts
#' @param gtf.all GRanges gtf annotation (can be used instead of specifying exons and transcripts)
#' @param BSgenome BSGenome object containing the genome for the species analysed
#' @param NMD Use NMD predictions? (Note: notNMD must be installed to use this feature)
#' @param exportGTF file name to export alternative isoform GTFs (default=NULL)
#' @param rearrangeXY should PSI directionality be taken into account?
#' @param uniprotData data.frame of uniprot sequence information
#' @param uniprotSeqFeatures data.frame of uniprot sequence features
#' @param selectLongest passed to getORFs()
#' @return data.frame containing signficant whippet diff data and ORF change summaries
#' @export
#' @importFrom rtracklayer import
#' @import GenomicRanges
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds.all <- readWhippetDataSet(whippetFiles)
#' wds.signif <- filterWhippetEvents(wds.all)
#'
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' whippetTranscriptChangeSummary(whippetDataSet = wds.signif, unfilteredWDS=wds.all ,BSgenome = g)
#'
#'
whippetTranscriptChangeSummary <- function(whippetDataSet,
                                           unfilteredWDS,
                                           BSgenome,
                                           eventTypes = "all",
                                           exons,
                                           NMD = FALSE,
                                           exportGTF = NULL,
                                           selectLongest=1){


    if(eventTypes[1] == "all"){
        eventTypes <- unique(diffSplicingResults(whippetDataSet)$type)
    }
    if(!is.null(exportGTF)){
        exportedTranscripts <- list()
    }
    if(exists("SignificantEvents.withORF")){
        # somehow everything screws up if this exists globally...
        suppressWarnings(rm(SignificantEvents.withORF))
    }

    if(any(eventTypes == "CE") &
       any(diffSplicingResults(whippetDataSet)$type == "CE")){


        whippetDataSet.ce <- filterWhippetEvents(whippetDataSet,
                                                 probability = 0,
                                                 psiDelta = 0,
                                                 eventTypes="CE")

        significantEvents.ce <- diffSplicingResults(whippetDataSet.ce)

        exons.ce <- findExonContainingTranscripts(input=whippetDataSet.ce,
                                                  exons=exons,
                                                  variableWidth=0,
                                                  findIntrons=FALSE)
        # make skipped exon transcripts
        skippedExonTranscripts <- skipExonInTranscript(whippetDataSet=whippetDataSet.ce,
                                                       skippedExons = exons.ce,
                                                       exons = exons,
                                                       glueExons = TRUE)

        orfChanges.ce <- transcriptChangeSummary(
            transcriptsX = skippedExonTranscripts[skippedExonTranscripts$set=="included_exon"],
            transcriptsY = skippedExonTranscripts[skippedExonTranscripts$set=="skipped_exon"],
            BSgenome = BSgenome,NMD = NMD, dataSet=whippetDataSet.ce,
            selectLongest = selectLongest)

        add = which(!(significantEvents.ce$coord %in% orfChanges.ce$id))
        if(length(add) > 0){
            significantEvents.ce <- plyr::rbind.fill(orfChanges.ce, significantEvents.ce[add,])
        }else{
            significantEvents.ce <- orfChanges.ce
        }

        if(exists("SignificantEvents.withORF")){
            SignificantEvents.withORF <- rbind(SignificantEvents.withORF,
                                               significantEvents.ce)
        }else{
            SignificantEvents.withORF <- significantEvents.ce
        }
        if(!is.null(exportGTF)){
            exportedTranscripts <- c(exportedTranscripts,
                                     skippedExonTranscripts)
        }
    }
    if(any(eventTypes == "RI") & any(diffSplicingResults(whippetDataSet)$type ==
                                     "RI")){

        whippetDataSet.ri <- filterWhippetEvents(whippetDataSet,
                                                 probability = 0,
                                                 psiDelta = 0,
                                                 eventTypes="RI")

        significantEvents.ri <- diffSplicingResults(whippetDataSet.ri)

        exons.ri <- findIntronContainingTranscripts(input=whippetDataSet.ri,
                                                    exons=exons)

        # find introns in the gtf that overlap whippet introns
        significantEvents.ri <-
            diffSplicingResults(whippetDataSet)[which(
                diffSplicingResults(whippetDataSet)$type=="RI"),]
        # add the intron into transcripts
        retainedIntronTranscripts <- addIntronInTranscript(whippetDataSet=whippetDataSet.ri,
                                                           flankingExons = exons.ri,
                                                           exons = exons,
                                                           glueExons = TRUE)


        orfChanges.ri <- transcriptChangeSummary(
            transcriptsX = retainedIntronTranscripts[retainedIntronTranscripts$set==
                                                         "spliced_intron"],
            transcriptsY = retainedIntronTranscripts[retainedIntronTranscripts$set==
                                                         "retained_intron"],
            BSgenome = BSgenome,NMD = NMD, dataSet=whippetDataSet.ri,
            selectLongest = selectLongest)

        add = which(!(significantEvents.ri$coord %in% orfChanges.ri$id))
        if(length(add) > 0){
            significantEvents.ri <- plyr::rbind.fill(orfChanges.ri, significantEvents.ri[add,])
        }else{
            significantEvents.ri <- orfChanges.ri
        }
        if(exists("SignificantEvents.withORF")){
            SignificantEvents.withORF <- rbind(SignificantEvents.withORF,
                                               significantEvents.ri)
        }else{
            SignificantEvents.withORF <- significantEvents.ri
        }
        if(!is.null(exportGTF)){
            exportedTranscripts <- c(exportedTranscripts,
                                     retainedIntronTranscripts)
        }
    }
    if(any(eventTypes %in% c("TS", "TE")) & any(diffSplicingResults(whippetDataSet)$type %in% c("TE","TS"))){

        events.t <- eventTypes[eventTypes %in% c("TE","TS")]
        events.significant <- unique(diffSplicingResults(whippetDataSet)$type)
        events.significant <- events.significant[events.significant %in% events.t]

        for(e in seq_along(events.significant)){
            event <- events.significant[e]
            whippetDataSet.t <- filterWhippetEvents(whippetDataSet,
                                                    eventTypes=event)
            significantEvents.t <- diffSplicingResults(whippetDataSet.t)

            transcripts.altT <- alterTranscriptStartEnds(whippetDataSet.t, exons, unfilteredWDS, type=event)

            orfChanges.t <- transcriptChangeSummary(
                transcriptsX = transcripts.altT[transcripts.altT$set== "tss_downregulated" | transcripts.altT$set== "tes_downregulated"],
                transcriptsY = transcripts.altT[transcripts.altT$set== "tss_upregulated" | transcripts.altT$set== "tes_upregulated"],
                BSgenome = BSgenome,NMD = NMD, dataSet=whippetDataSet.t,
                selectLongest = selectLongest)

            # add to significantEvents.t
            # this is different to other events due to how TE/TS coords are given.
            # only ONE orf chance per 'group' of events, instead of duplicating the orf change in the opposite direction
            m <- match(significantEvents.t$coord, orfChanges.t$coord)
            significantEvents.t <-  orfChanges.t[m[which(!is.na(m))],]


            if(exists("SignificantEvents.withORF")){
                SignificantEvents.withORF <- rbind(SignificantEvents.withORF,
                                                   significantEvents.t)
            }else{
                SignificantEvents.withORF <- significantEvents.t
            }
            if(!is.null(exportGTF)){
                exportedTranscripts <- c(exportedTranscripts,
                                         transcripts.altT)
            }
        }
    }
    if(any(eventTypes %in% c("AA","AD","AF","AL")) &
       any(diffSplicingResults(whippetDataSet)$type %in%
           c("AA","AD","AF","AL"))){

        events.junctions <- eventTypes[eventTypes %in% c("AA","AD","AF","AL")]
        events.significant <- unique(diffSplicingResults(whippetDataSet)$type)
        events.significant <- events.significant[events.significant %in%
                                                     events.junctions]

        for(e in seq_along(events.significant)){

            event <- events.significant[e]

            whippetDataSet.jnc <- filterWhippetEvents(whippetDataSet,
                                                      probability = 0,
                                                      psiDelta = 0,
                                                      eventTypes=event)

            significantEvents.jnc <- diffSplicingResults(whippetDataSet.jnc)

            junctionPairs <- findJunctionPairs(whippetDataSet.jnc, type=event)

            # check for pairs
            ids.x <- unique(junctionPairs$event_id[junctionPairs$set=="X"])
            ids.x <- ids.x[ids.x %in% unique(
                junctionPairs$event_id[junctionPairs$set=="Y"])]

            significantEvents.jnc <-significantEvents.jnc[
                which(significantEvents.jnc$coord %in% ids.x),]
            junctionPairs <- junctionPairs[
                which(junctionPairs$event_id %in% ids.x),]

            if(nrow(significantEvents.jnc) > 0){
                # make transcripts with alternative junction usage
                altTranscripts <- replaceJunction(whippetDataSet.jnc,
                                                  junctionPairs,
                                                  exons, type=event)
                orfChanges.jnc <- transcriptChangeSummary(
                    transcriptsX = altTranscripts[
                        altTranscripts$set==paste0(event, "_X")],
                    transcriptsY = altTranscripts[
                        altTranscripts$set==paste0(event, "_Y")],
                    BSgenome = BSgenome,NMD = NMD,
                    dataSet=whippetDataSet.jnc,
                    selectLongest = selectLongest)

                # add to significantEvents
                add = which(!(significantEvents.jnc$coord %in% orfChanges.jnc$id))
                if(length(add) > 0){
                    significantEvents.jnc <- plyr::rbind.fill(orfChanges.jnc, significantEvents.jnc[add,])
                }else{
                    significantEvents.jnc <- orfChanges.jnc
                }

                if(exists("SignificantEvents.withORF")){
                    SignificantEvents.withORF <-
                        rbind(SignificantEvents.withORF,
                              significantEvents.jnc)
                }else{
                    SignificantEvents.withORF <- significantEvents.jnc
                }
                if(!is.null(exportGTF)){
                    exportedTranscripts <- c(exportedTranscripts,
                                             altTranscripts)
                }
            }
        }

    }
    if(!is.null(exportGTF)){
        exportedTranscripts <- do.call("c", exportedTranscripts)
        rtracklayer::export.gff(exportedTranscripts, con=exportGTF,
                                format="gtf")
    }
    if(exists("SignificantEvents.withORF")){
        return(SignificantEvents.withORF)
    }else{
        return(NULL)
    }
}
