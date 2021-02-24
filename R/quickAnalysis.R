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
#' Compare open reading frames for two sets of paired transcripts
#' @param transcriptsX GRanges object with exon annotations for
#' all transcripts to be compared for the 'normal' condition
#' @param transcriptsY GRanges object with exon annotations for
#' all transcripts to be compared for the 'alternative' condition
#' @param BSgenome BSGenome object containing the genome for the species analysed
#' @param exons GRanges object made from a GTF containing exon coordinates
#' @param NMD Use NMD predictions? (Note: notNMD must be installed to use this feature)
#' @param NMDModel Use the "base" or "lncRNA" NMD model?
#' @param orfPrediction What type of orf predictions to return. default= \code{"allFrames"}
#' @param compareBy compare isoforms by 'transcript' id, or aggregate all changes occuring by 'gene'
#' @param compareToGene compare alternative isoforms to all normal gene isoforms (in exons)
#' @param dataSet whippetDataSet/rMATSDataSet generated from \code{readWhippetDataSet() or \code{readrMATSDataSet()}}
#' Use if PSI directionality should be taken into account when comparing isoforms.
#' @param exportGTF file name to export alternative isoform GTFs (default=\code{NULL})
#' @param uniprotData data.frame of uniprot sequence information
#' @param uniprotSeqFeatures data.frame of uniprot sequence features
#' @param selectLongest passed to getORFs()
#' @return Summarised ORF changes data.frame
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @importFrom utils installed.packages
#' @family transcript isoform comparisons
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
#' wds.exonSkip <- filterWhippetEvents(wds, eventTypes="CE",psiDelta = 0.2)
#'
#' exons.exonSkip <- findExonContainingTranscripts(wds.exonSkip, exons,
#' variableWidth=0, findIntrons=FALSE, transcripts)
#' ExonSkippingTranscripts <- skipExonInTranscript(exons.exonSkip, exons, whippetDataSet=wds.exonSkip)
#' transcriptChangeSummary(ExonSkippingTranscripts[ExonSkippingTranscripts$set=="included_exon"],
#' ExonSkippingTranscripts[ExonSkippingTranscripts$set=="skipped_exon"],
#' BSgenome=g,exons)

transcriptChangeSummary <- function(transcriptsX,
                                    transcriptsY,
                                    BSgenome,
                                    exons,
                                    dataSet = NULL,
                                    exportGTF = NULL,

                                    NMD = FALSE,

                                    compareBy="gene",
                                    orfPrediction = "allFrames",
                                    compareToGene = FALSE,

                                    selectLongest=1){


    if(!is.null(dataSet)){
        if(class(dataSet)[1] == "whippetDataSet"){
            whippetEvents <- diffSplicingResults(dataSet)

            allTranscripts <- c(transcriptsX, transcriptsY)
            txEvents <- unlist(lapply(str_split(lapply(str_split(allTranscripts$transcript_id, "[+]AS"), "[[", 2), "[ ]"), "[[", 1))

            if(any(txEvents %in% c("TEU", "TED", "TSU", "TSD"))){
                transcriptsX <- allTranscripts[txEvents %in% c("TED", "TSD")]
                transcriptsY <- allTranscripts[txEvents %in% c("TEU", "TSU")]
                allTranscripts <- allTranscripts[-which(txEvents %in% c("TED", "TSD","TEU", "TSU"))]
            }else{
                # clear transcripts X/Y
                transcriptsX <- transcriptsX[-(1:length(transcriptsX))]
                transcriptsY <- transcriptsY[-(1:length(transcriptsY))]
            }

            # for non-TS/TE events
            if(length(allTranscripts) > 0){

                type <- gsub("_X","",gsub("_Y","", allTranscripts$set))
                type <- gsub("skipped_exon", "CE", gsub("included_exon","CE", type))
                type <- gsub("retained_intron", "RI", gsub("spliced_intron","RI", type))

                m <- match(paste0(allTranscripts$event_id,"_",type),
                           paste0(whippetEvents$coord,"_",
                                  whippetEvents$type))
                # A -- psi in condition 1 (A) is higher (i.e. included -- > skipped)
                normA <- which(whippetEvents$psi_a > whippetEvents$psi_b)
                # B -- psi in condition 2 (B) is higher (i.e. skipped -- > included)
                normB <- which(whippetEvents$psi_a < whippetEvents$psi_b)

                #sets for X (+A)
                setsX <- c(paste0(unique(type), "_Y"), "included_exon","retained_intron")
                #sets for Y (+B)
                setsY <- c(paste0(unique(type), "_X"), "skipped_exon","spliced_intron")

                transcriptsX <- c(transcriptsX, allTranscripts[
                    which((m %in% normA & allTranscripts$set %in% setsX) |
                              (m %in% normB & allTranscripts$set %in% setsY))])

                transcriptsY <- c(transcriptsY, allTranscripts[
                    which((m %in% normA & allTranscripts$set %in% setsY) |
                                 (m %in% normB & allTranscripts$set %in% setsX))])
            }


        }

        if(class(dataSet)[1] == "rMATSDataSet"){

            # combine all events (only PSI direction/eventID)
            allEvents = NULL
            for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
                psiDiffs = slot(rds, eventType)[,c("ID","GeneID","geneSymbol", "PValue","FDR","IncLevelDifference")]
                psiDiffs$type = eventType
                allEvents = rbind(allEvents, psiDiffs)
            }

            allTranscripts <- c(transcriptsX, transcriptsY)
            allTranscripts$type = NA
            allTranscripts$type[allTranscripts$set %in% c("included_exon","skipped_exon")] = "SE"
            allTranscripts$type[allTranscripts$set %in% c("included_exon1","included_exon2")] = "MXE"
            allTranscripts$type[allTranscripts$set %in% c("spliced_intron","retained_intron")] = "RI"
            allTranscripts$type[allTranscripts$set %in% c("alt3_splicesite_long","alt3_splicesite_short")] = "A3SS"
            allTranscripts$type[allTranscripts$set %in% c("alt5_splicesite_long","alt5_splicesite_short")] = "A5SS"

            eventId = unlist(lapply(str_split(lapply(str_split(allTranscripts$transcript_id, "[ ]"),"[[", 2), "[-]"),"[[", 1))

            allEvents$event_id = unlist(lapply(str_split(allTranscripts$transcript_id[match(paste0(allEvents$ID, "_", allEvents$type),
                                                                                            paste0(eventId, "_", allTranscripts$type))],
                                                         "[ ]"), "[[" ,2))

            m = match(paste0(eventId, "_", allTranscripts$type),
                      paste0(allEvents$ID, "_", allEvents$type))

            normA = which(allEvents$IncLevelDifference > 0)
            normB = which(allEvents$IncLevelDifference < 0)

            setsA = c("included_exon", "included_exon2", "retained_intron","alt3_splicesite_long","alt5_splicesite_long")
            setsB = c("Skipped_exon", "included_exon1", "spliced_intron","alt3_splicesite_short","alt5_splicesite_short")

            transcriptsX <- allTranscripts[
                which((m %in% normA & allTranscripts$set %in% setsA) |
                          (m %in% normB & allTranscripts$set %in% setsB))]

            transcriptsY <- allTranscripts[
                which((m %in% normA & allTranscripts$set %in% setsB) |
                          (m %in% normB & allTranscripts$set %in% setsA))]

        }
    }

    if(!is.null(exportGTF)){
        transcriptsX$comp_set <- "X"
        transcriptsY$comp_set <- "Y"
        transcriptsXY <- c(transcriptsX, transcriptsY)

        rtracklayer::export.gff(transcriptsXY, con=exportGTF, format="gtf")

        transcriptsX$comp_set <- NULL
        transcriptsY$comp_set <- NULL
    }

    if(orfPrediction == "allFrames"){
        orfsX <- getOrfs(transcriptsX, BSgenome, returnLongestOnly = FALSE,
                         allFrames = TRUE, uORFs = TRUE,
                         selectLongest = selectLongest)
        orfsY <- getOrfs(transcriptsY, BSgenome, returnLongestOnly = FALSE,
                         allFrames = TRUE, uORFs = TRUE,
                         selectLongest = selectLongest)
    }else{
        orfsX <- getOrfs(transcriptsX, BSgenome,returnLongestOnly = TRUE,
                         uORFs = TRUE, selectLongest = selectLongest)
        orfsY <- getOrfs(transcriptsY, BSgenome,returnLongestOnly = TRUE,
                         uORFs = TRUE, selectLongest = selectLongest)
    }

    if(all(!grepl("[+]", orfsX$id))){

        if(orfPrediction == "allFrames"){
            Yid.withFrame <- paste0(unlist(lapply(
                str_split(orfsY$id, "[+]"),"[[",1)),"_", orfsY$frame)
            Xid.withFrame <- paste0(orfsX$id,"_", orfsX$frame)
            m <- match(Yid.withFrame, Xid.withFrame)
        }else{
            m <- match(unlist(lapply(str_split(orfsY$id, "[+]"),"[[",1)),
                       orfsX$id)
        }

        orfsX<- orfsX[m,]
        orfsX$id <- orfsY$id
        #orfsX <- orfsX[which(!duplicated(orfsX$id)),]
    }

    orfsX <- orfsX[which(!is.na(orfsX$orf_length)),]
    orfsY <- orfsY[which(!is.na(orfsY$orf_length)),]

    if(NMD == TRUE){
        orfsX <- manualNMD(orfsX)
        orfsY <- manualNMD(orfsY)
    }

    if(compareToGene == TRUE){
        orfAllGenes <- getOrfs(exons[exons$gene_id %in% unique(c(transcriptsX$gene_id, transcriptsY$gene_id)) & exons$transcript_type=="protein_coding"],
                               BSgenome=BSgenome, returnLongestOnly = TRUE)

        if(NMD == TRUE){
            orfAllGenes <- manualNMD(orfAllGenes)
        }

        orfChange <- orfDiff(orfsX, orfsY, filterNMD = NMD, compareBy = "gene",
                             geneSimilarity = TRUE,
                             allORFs = orfAllGenes, compareUTR = TRUE)
    }else{

        orfChange <- orfDiff(orfsX, orfsY, filterNMD = NMD, compareBy = "gene",
                             compareUTR = TRUE)
    }


    if(NMD == TRUE){
        nmdChangeMan <- attrChangeAltSpliced(orfsX,
                                             orfsY,
                                             attribute="nmd_prob_manual",
                                             compareBy="gene",
                                             useMax=FALSE)
        m <- match(orfChange$id, nmdChangeMan$id)
        orfChange <- cbind(orfChange, nmdChangeMan[m,-1])
    }

    if(!is.null(dataSet)){
        if(class(dataSet)[1] == "whippetDataSet"){
            m <- match(orfChange$id, whippetEvents$coord)

            # for TS/TE events
            # need to collapse/re-match to a 'group' id, otherwise you get lots of doubleups
            if(all(is.na(m))){

                tidToEvent = data.frame(tid=c(transcriptsX$transcript_id, transcriptsY$transcript_id),
                                        event_id=c(transcriptsX$event_id, transcriptsY$event_id))
                tidToEvent = dplyr::distinct(tidToEvent)

                #m = match(orfChange$id, unlist(lapply(str_split(tidToEvent$tid, "[ ]"), "[[", 2)))
                m = match(whippetEvents$coord, tidToEvent$event_id)
                whippetEvents = whippetEvents[which(!is.na(m)),]
                m = match(whippetEvents$coord, tidToEvent$event_id)
                whippetEvents$group_id = unlist(lapply(str_split(tidToEvent$tid[m], "[ ]"), "[[", 2))

                whippetEvents = arrange(whippetEvents, group_id, plyr::desc(probability), plyr::desc(abs(psi_delta)))
                whippetEvents = whippetEvents[!duplicated(whippetEvents$group_id),]

                m = match(orfChange$id, whippetEvents$group_id)

                orfChange <- cbind(whippetEvents[m,], orfChange)
                orfChange$group_id = NULL

            }else{
                orfChange <- cbind(whippetEvents[m,], orfChange)
            }
        }else if(class(dataSet)[1] == "rMATSDataSet"){
            m <- match(orfChange$id, allEvents$event_id)
            orfChange = cbind(allEvents[m,], orfChange[,-1])
        }
    }
    return(orfChange)

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
#' wds <- readWhippetDataSet(whippetFiles)
#' wds <- filterWhippetEvents(wds)
#'
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' whippetTranscriptChangeSummary(wds, gtf.all=gtf,BSgenome = g)
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
                                                     probability = 0.99,
                                                     psiDelta = 0.4,
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
#' Compare open reading frames for whippet differentially spliced events
#' @param significantEvents  data.frame containing information from the
#' per_intron_results.tab file output from leafcutter.
#' @param combineGeneEvents combine clusters occuring in the same gene?
#' Currently not reccomended.
#' @param exons GRanges gtf annotation of exons
#' @param BSgenome BSGenome object containing the genome for the species analysed
#' @param NMD Use NMD predictions? (Note: notNMD must be installed to use this feature)
#' @param showProgressBar show a progress bar of alternative isoform generation?
#' @param junctions junctions GRanges object from readLeafcutterJunctions()
#' @param exportGTF file name to export alternative isoform GTFs (default=NULL)
#' @param uniprotData data.frame of uniprot sequence information
#' @param uniprotSeqFeatures data.frame of uniprot sequence features
#' @return data.frame containing signficant whippet diff data and ORF change summaries
#' @export
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom stringr str_split
#' @importFrom rtracklayer import
#' @import GenomicRanges
#' @family leafcutter data processing
#' @author Beth Signal
#' @examples
#' leafcutterFiles <- list.files(system.file("extdata","leafcutter/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' leafcutterIntrons <- read.delim(leafcutterFiles[
#' grep("intron_results", leafcutterFiles)],stringsAsFactors=FALSE)
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' leafcutterTranscriptChangeSummary(significantEvents = leafcutterIntrons,
#' exons=exons,BSgenome = g,NMD=FALSE)

leafcutterTranscriptChangeSummary <- function(significantEvents,
                                              combineGeneEvents=FALSE,
                                              exons,
                                              BSgenome,
                                              NMD = FALSE,
                                              showProgressBar=TRUE,
                                              junctions=NULL,
                                              exportGTF=NULL){


    geneEvents <- as.data.frame(table(significantEvents$ensemblID,
                                      significantEvents$clusterID))
    geneEvents <- geneEvents[geneEvents$Freq!=0,]

    if(combineGeneEvents == FALSE){
        if(showProgressBar){
            message(paste0("Generating alternative isoforms for ",
                           nrow(geneEvents), " clusters:"))
            pb <- utils::txtProgressBar(min = 0, max = nrow(geneEvents),
                                        style = 3)
        }

        altIso <- alternativeIntronUsage(altIntronLocs = significantEvents[significantEvents$clusterID == geneEvents$Var2[1],],
            exons, junctions=junctions)

        if(showProgressBar){utils::setTxtProgressBar(pb, 1)}

        if(nrow(geneEvents) > 1){
            for(i in 2:nrow(geneEvents)){
                altIntronLocs = significantEvents[
                    significantEvents$clusterID == geneEvents$Var2[i],]
                if(!("deltapsi" %in% colnames(altIntronLocs))){
                    altIntronLocs$deltapsi <- altIntronLocs$PSI_a - altIntronLocs$PSI_b
                }
                if(all(altIntronLocs$deltapsi < 0) | all(altIntronLocs$deltapsi > 0)){
                    altIntronLocs <- altIntronLocs[-(1:nrow(altIntronLocs)),]
                }

                if(nrow(altIntronLocs) > 1){
                    altIso1 <- alternativeIntronUsage(altIntronLocs, exons, junctions=junctions)
                    if(!is.null(altIso1)){
                        altIso <- c(altIso, altIso1)
                    }
                }
                if(showProgressBar){utils::setTxtProgressBar(pb, i)}

            }
        }
    }else{
        genes <- unique(geneEvents$Var1)
        if(showProgressBar){
            message(paste0("Generating alternative isoforms for ",
                           nrow(geneEvents), " genes:"))
            pb <- utils::txtProgressBar(min = 0, max = length(genes), style = 3)
        }
        for(j in seq_along(genes)){
            clusters <- geneEvents[geneEvents$Var1==genes[j],]

            altIso <- alternativeIntronUsage(significantEvents[
                significantEvents$clusterID == clusters$Var2[1],],
                exons)

            if(nrow(clusters) > 1){
                for(i in 2:nrow(clusters)){
                    altIntronLocs = significantEvents[
                        significantEvents$clusterID == clusters$Var2[i],]
                    altIntronLocs <- altIntronLocs[altIntronLocs$verdict==
                                                       "annotated",]
                    if(nrow(altIntronLocs) > 1){
                        altIso1 <- alternativeIntronUsage(altIntronLocs,
                                                          c(exons, altIso))
                        altIso <- c(altIso, altIso1)
                    }

                }
            }
            if(showProgressBar){utils::setTxtProgressBar(pb, j)}

        }


    }
    if(is.list(altIso)){
        if(length(altIso) > 1){
            altIso <- do.call("c", altIso)
        }
    }

    altIso$spliced_id <- unlist(lapply(
        stringr::str_split(altIso$transcript_id, " "),"[[",2))

    transcriptsX <- altIso[grep("dnre", altIso$transcript_id)]
    transcriptsY <- altIso[grep("upre", altIso$transcript_id)]

    if(!is.null(exportGTF)){
        rtracklayer::export.gff(altIso, con=exportGTF, format="gtf")
    }

    orfDiff <- transcriptChangeSummary(transcriptsX,
                                       transcriptsY,
                                       BSgenome = BSgenome,
                                       NMD = NMD)
    m <- match(gsub("_","",significantEvents$clusterID), orfDiff$id)
    significantEvents.withORF <- cbind(significantEvents, orfDiff[m,-1])
    #significantEvents.withORF <- significantEvents.withORF[!duplicated(m),]

    return(significantEvents.withORF)
}
################ FUCK

readRMATS = function(directory,
                     type = "JC"){

    RMATSEventTypes = c("SE", "MXE", "RI", "A3SS", "A5SS")
    RMATSFileList = paste0(RMATSEventTypes, ".MATS.", type, ".txt")

    allFiles = list.files(directory, full.names = TRUE)
    diffSpliceFiles = allFiles[basename(allFiles) %in% RMATSFileList]

    if(length(diffSpliceFiles)==0){
        stop("no RMATs files in the specified directory")
    }else if(!all(RMATSFileList %in% basename(diffSpliceFiles))){
        for(f in RMATSFileList[which(!(RMATSFileList %in% basename(diffSpliceFiles)))]){
            message(paste0("Can't find file: ", f, " please check if this file should exist in the directory"))
        }
    }

    diffSplice.SE = fread(diffSpliceFiles[basename(diffSpliceFiles) == RMATSFileList[RMATSEventTypes=="SE"]], header=TRUE, data.table=FALSE)
    diffSplice.MXE = fread(diffSpliceFiles[basename(diffSpliceFiles) == RMATSFileList[RMATSEventTypes=="MXE"]], header=TRUE, data.table=FALSE)
    diffSplice.RI = fread(diffSpliceFiles[basename(diffSpliceFiles) == RMATSFileList[RMATSEventTypes=="RI"]], header=TRUE, data.table=FALSE)
    diffSplice.A3SS = fread(diffSpliceFiles[basename(diffSpliceFiles) == RMATSFileList[RMATSEventTypes=="A3SS"]], header=TRUE, data.table=FALSE)
    diffSplice.A5SS = fread(diffSpliceFiles[basename(diffSpliceFiles) == RMATSFileList[RMATSEventTypes=="A5SS"]], header=TRUE, data.table=FALSE)

    rMATSEvents = plyr::rbind.fill(cbind(diffSplice.SE, event = "SE"),
                                   cbind(diffSplice.MXE, event = "MXE"),
                                   cbind(diffSplice.RI, event = "RI"),
                                   cbind(diffSplice.A3SS, event = "A3SS"),
                                   cbind(diffSplice.A5SS, event = "A5SS"))

}
rMATSEvents = readRMATS(directory = "~/Projects/UTAS/IsoformModeller/rmats_00_48/", type="JC")

rMATSEvents.signif = rMATSEvents[rMATSEvents$FDR < 0.01 & abs(rMATSEvents$IncLevelDifference) > 0.2,]
table(rMATSEvents.signif$event)


whippetTranscriptChangeSummary <- function(whippetDataSet,
                                           unfilteredWDS,
                                           gtf.all=NULL,
                                           BSgenome,
                                           eventTypes = "all",
                                           exons=NULL,
                                           transcripts=NULL,
                                           NMD = FALSE,
                                           exportGTF = NULL,
                                           uniprotData=NULL,
                                           uniprotSeqFeatures=NULL,
                                           selectLongest=1){


    # diffSplicingResults(whippetDataSet)

    if(any(eventTypes == "CE") &
       any(diffSplicingResults(whippetDataSet)$type == "CE")){

        diffSplice.SE.signif = rMATSEvents.signif[which(rMATSEvents.signif$event == "SE"),]
        isoforms.SE = skipExonByJunction(input = diffSplice.SE.signif, eventType = "SE", exons=exons)

        cond1_isoforms = isoforms.SE[(isoforms.SE$rmats_id %in%
                                         diffSplice.SE.signif$ID[diffSplice.SE.signif$IncLevelDifference >0] &
                                         isoforms.SE$set=="included_exon") |
                                         (isoforms.SE$rmats_id %in%
                                              diffSplice.SE.signif$ID[diffSplice.SE.signif$IncLevelDifference <0] &
                                              isoforms.SE$set=="skipped_exon")]
        cond2_isoforms = isoforms.SE[(isoforms.SE$rmats_id %in%
                                          diffSplice.SE.signif$ID[diffSplice.SE.signif$IncLevelDifference >0] &
                                          isoforms.SE$set=="skipped_exon") |
                                         (isoforms.SE$rmats_id %in%
                                              diffSplice.SE.signif$ID[diffSplice.SE.signif$IncLevelDifference <0] &
                                              isoforms.SE$set=="included_exon")]

        orfChanges.se <- transcriptChangeSummary(
            transcriptsX=cond1_isoforms,
            transcriptsY=cond2_isoforms,
            BSgenome = BSgenome,
            NMD = NMD,
            rearrangeXY = FALSE,
            uniprotData = uniprotData,
            uniprotSeqFeatures = uniprotSeqFeatures,
            selectLongest = selectLongest)


        whippetDataSet.ce <- filterWhippetEvents(whippetDataSet,
                                                 probability = 0,
                                                 psiDelta = 0,
                                                 eventTypes="CE")

        significantEvents.ce <- diffSplicingResults(whippetDataSet.ce)

        # make skipped exon transcripts
        skippedExonTranscripts <- skipExonInTranscript(whippetDataSet=whippetDataSet.ce,
                                                       skippedExons = exons.ce,
                                                       exons = exons,
                                                       glueExons = TRUE)

        orfChanges.ce <- transcriptChangeSummary(
            skippedExonTranscripts[skippedExonTranscripts$set=="included_exon"],
            skippedExonTranscripts[skippedExonTranscripts$set=="skipped_exon"],
            BSgenome = BSgenome,NMD = NMD, whippetDataSet=whippetDataSet.ce,
            rearrangeXY = rearrangeXY,
            uniprotData = uniprotData,
            uniprotSeqFeatures = uniprotSeqFeatures,
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
            retainedIntronTranscripts[retainedIntronTranscripts$set==
                                          "spliced_intron"],
            retainedIntronTranscripts[retainedIntronTranscripts$set==
                                          "retained_intron"],
            BSgenome = BSgenome,NMD = NMD, whippetDataSet=whippetDataSet.ri,
            rearrangeXY = rearrangeXY,
            uniprotData = uniprotData,
            uniprotSeqFeatures = uniprotSeqFeatures,
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
                                                    probability = 0.99,
                                                    psiDelta = 0.4,
                                                    eventTypes=event)
            significantEvents.t <- diffSplicingResults(whippetDataSet.t)

            transcripts.altT <- alterTranscriptStartEnds(whippetDataSet.t, exons, unfilteredWDS, type=event)

            orfChanges.t <- transcriptChangeSummary(
                transcriptsX = transcripts.altT[transcripts.altT$set== "tss_upregulated" | transcripts.altT$set== "tes_upregulated"],
                transcriptsY = transcripts.altT[transcripts.altT$set== "tss_downregulated" | transcripts.altT$set== "tes_downregulated"],
                BSgenome = BSgenome,NMD = NMD, whippetDataSet=whippetDataSet.t,
                rearrangeXY = FALSE,
                uniprotData = uniprotData,
                uniprotSeqFeatures = uniprotSeqFeatures,
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
            ids.x <- unique(junctionPairs$whippet_id[junctionPairs$set=="X"])
            ids.x <- ids.x[ids.x %in% unique(
                junctionPairs$whippet_id[junctionPairs$set=="Y"])]

            significantEvents.jnc <-significantEvents.jnc[
                which(significantEvents.jnc$coord %in% ids.x),]
            junctionPairs <- junctionPairs[
                which(junctionPairs$whippet_id %in% ids.x),]

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
                    whippetDataSet=whippetDataSet.jnc,
                    rearrangeXY = rearrangeXY,
                    uniprotData = uniprotData,
                    uniprotSeqFeatures = uniprotSeqFeatures,
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
