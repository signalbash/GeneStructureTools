#' Import RMATS Turbo results files as a rmatsDataSet
#' @param filePath path to RMATS differential splicing output files
#' @param type type of counts to use. "JC" or "JCEC"
#' @return rmatsDataSet
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' rmats_filePath <- system.file("extdata","rmats_small/", package="GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_filePath)
readRmatsDataSet <- function(filePath, type="JC"){

    rds <- new("rmatsDataSet", filePath=filePath)
    rmatsEventTypes <- c("SE", "MXE", "RI", "A3SS", "A5SS")
    rmatsFileList <- paste0(rmatsEventTypes, ".MATS.", type, ".txt")
    allFiles <- list.files(filePath, full.names=TRUE)
    diffSpliceFiles <- allFiles[basename(allFiles) %in% rmatsFileList]

    if(length(diffSpliceFiles) == 0){
        stop("no rmats files in the specified filePath")
    }else if(!all(rmatsFileList %in% basename(diffSpliceFiles))){
        for(f in rmatsFileList[which(!(rmatsFileList %in% basename(diffSpliceFiles)))]){
            message(paste0("Can't find file: ", f, " please check if this file should exist in the filePath"))
        }
    }

    for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
        if(any(basename(diffSpliceFiles) == paste0(eventType, ".MATS.", type, ".txt"))){
            slot(rds, eventType) <- fread(diffSpliceFiles[basename(diffSpliceFiles) == rmatsFileList[rmatsEventTypes == eventType]], header=TRUE, data.table=FALSE)
        }
    }
    return(rds)
}



#' Filter out significant events from a RMATS dataset
#' @param rmatsDataSet rmatsDataSet generated from \code{readRmatsDataSet()}
#' @param FDR maximum FDR required to call event as significant
#' @param psiDelta minimum change in psi required to call an event as significant
#' @param idList (optional) list of gene ids to filter for
#' @param minCounts minumum number of counts for all replicates
#' in at least one condition to call an event as significant
#' @param medianCounts median count for all replicates
#' in at least one condition to call an event as significant
#' @param sampleTable data.frame with sample names and conditions.
#' Only needed if filtering with counts.
#' @return filtered rmatsDataSet
#' @export
#' @importFrom stats median
#' @import stringr
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' rmats_directory <- system.file("extdata","rmats_small/", package="GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_directory)
#' rds.filtered <- filterRmatsEvents(rds, FDR=0.01, psiDelta=0.1)
#' # filter by gene name/id
#' rds.Tmem208 <- filterRmatsEvents(rds, idList="Tmem208", FDR=1, psiDelta=0)
filterRmatsEvents <- function(rmatsDataSet,
                                FDR=0.05,
                                psiDelta=0.1,
                                idList=NA,
                                minCounts=NA,
                                medianCounts=NA,
                                sampleTable){

    #set FDR/psiDelta if NA
    if(is.na(FDR[1])){
        FDR <- 1
    }
    if(is.na(psiDelta[1])){
        psiDelta <- 0
    }

    for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
        tmp <- slot(rmatsDataSet, eventType)
        if(nrow(tmp) > 0){
            significantEventsIndex <- which(tmp$FDR <= FDR & abs(tmp$IncLevelDifference) > psiDelta)
            slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex,]
        }
    }

    if(!is.na(idList[1])){
        for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
            tmp <- slot(rmatsDataSet, eventType)
            if(nrow(tmp) > 0){
                significantEventsIndex <- which(tmp$GeneID %in% idList | tmp$geneSymbol %in% idList)
                slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex,]
            }
        }
    }

    if(!is.na(minCounts[1])){
        for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
            tmp <- slot(rmatsDataSet, eventType)
            if(nrow(tmp) > 0){
                cond1 <- mapply(function(x,y) as.numeric(x) + as.numeric(y), x=stringr::str_split(tmp$IJC_SAMPLE_1, ","), y=stringr::str_split(tmp$SJC_SAMPLE_1, ","))
                cond2 <- mapply(function(x,y) as.numeric(x) + as.numeric(y), x=stringr::str_split(tmp$IJC_SAMPLE_2, ","), y=stringr::str_split(tmp$SJC_SAMPLE_2, ","))
                significantEventsIndex <- which(apply(rbind(cond1, cond2), 2, function(x) all(x >= minCounts)))
                slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex,]
            }
        }
    }
    if(!is.na(medianCounts[1])){
        for(eventType in c("SE", "MXE", "RI", "A3SS", "A5SS")){
            tmp <- slot(rmatsDataSet, eventType)
            if(nrow(tmp) > 0){
                cond1 <- mapply(function(x,y) as.numeric(x) + as.numeric(y), x=stringr::str_split(tmp$IJC_SAMPLE_1, ","), y=stringr::str_split(tmp$SJC_SAMPLE_1, ","))
                cond2 <- mapply(function(x,y) as.numeric(x) + as.numeric(y), x=stringr::str_split(tmp$IJC_SAMPLE_2, ","), y=stringr::str_split(tmp$SJC_SAMPLE_2, ","))
                significantEventsIndex <- which(apply(cond1, 2, median) >= medianCounts | apply(cond2, 2, median) >= medianCounts)
                slot(rmatsDataSet, eventType) <- tmp[significantEventsIndex,]
            }
        }
    }

    return(rmatsDataSet)

}


#' Compare open reading frames for RMATS differentially spliced events
#' @param rmatsDataSet rmatsDataSet generated from \code{readRmatsDataSet()}
#' @param eventTypes which event type to filter for? default="all"
#' @param exons GRanges gtf annotation of exons
#' @param BSgenome BSGenome object containing the genome for the species analysed
#' @param NMD Use NMD predictions? (Note: notNMD must be installed to use this feature)
#' @param exportGTF file name to export alternative isoform GTFs (default=NULL)
#' @return data.frame containing significant RMATS differential splicing data and ORF change summaries
#' @export
#' @importFrom rtracklayer export.gff
#' @import GenomicRanges
#' @family rmats data processing
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","gencode.vM25.small.gtf", package="GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'
#' rmats_directory <- system.file("extdata","rmats_small/", package="GeneStructureTools")
#' rds <- readRmatsDataSet(rmats_directory)
#' rds.filtered <- filterRmatsEvents(rds, FDR=0.01, psiDelta=0.1)
#' rmats_summary <- rmatsTranscriptChangeSummary(rmatsDataSet=rds.filtered, exons, BSgenome=g)

rmatsTranscriptChangeSummary <- function(rmatsDataSet,
                                         exons=NULL,
                                         eventTypes="all",
                                         BSgenome,
                                         NMD=TRUE,
                                         exportGTF=NULL){

    if(eventTypes[1] == "all"){
        eventTypes <- c("SE", "MXE", "RI", "A3SS", "A5SS")
    }

    allTranscripts <- GRanges()
    orfChanges <- NULL

    diffSplice.SE.signif <- extractEvent(rmatsDataSet, "SE")
    if(nrow(diffSplice.SE.signif) > 0 & "SE" %in% eventTypes){
        message(paste0("Creating isoforms for ", nrow(diffSplice.SE.signif), " SE events"))
        isoforms.SE <- skipExonByJunction(diffSplice.SE.signif, eventType="SE", exons=exons)
        orfChanges.SE <- transcriptChangeSummary(isoforms.SE[isoforms.SE$set == "included_exon"],
                                                 isoforms.SE[isoforms.SE$set == "skipped_exon"],
                                                 BSgenome=BSgenome, NMD=NMD, exportGTF=NULL, dataSet=rmatsDataSet)
        m <- match(unlist(lapply(stringr::str_split(orfChanges.SE$id, "[-]"), "[[", 1)), diffSplice.SE.signif$ID)
        orfChanges <- rbind(orfChanges, cbind(diffSplice.SE.signif[m,c('ID', 'GeneID', 'geneSymbol', 'PValue','FDR', 'IncLevelDifference')], type="SE", orfChanges.SE))
        allTranscripts <- c(allTranscripts, isoforms.SE)
    }

    diffSplice.MXE.signif <- extractEvent(rmatsDataSet, "MXE")
    if(nrow(diffSplice.MXE.signif) > 0 & "MXE" %in% eventTypes){
        message(paste0("Creating isoforms for ", nrow(diffSplice.MXE.signif), " MXE events"))
        isoforms.MXE <- skipExonByJunction(diffSplice.MXE.signif, eventType="MXE", exons=exons)
        orfChanges.MXE <- transcriptChangeSummary(isoforms.MXE[isoforms.MXE$set == "included_exon1"],
                                                 isoforms.MXE[isoforms.MXE$set == "included_exon2"],
                                                 BSgenome=BSgenome, NMD=NMD, exportGTF=NULL, dataSet=rmatsDataSet)
        m <- match(unlist(lapply(stringr::str_split(orfChanges.MXE$id, "[-]"), "[[", 1)), diffSplice.MXE.signif$ID)
        orfChanges <- rbind(orfChanges, cbind(diffSplice.MXE.signif[m,c('ID', 'GeneID', 'geneSymbol', 'PValue','FDR', 'IncLevelDifference')], type="MXE", orfChanges.MXE))
        allTranscripts <- c(allTranscripts, isoforms.MXE)

    }

    diffSplice.RI.signif <- extractEvent(rmatsDataSet, "RI")
    if(nrow(diffSplice.RI.signif) > 0 & "RI" %in% eventTypes){
        message(paste0("Creating isoforms for ", nrow(diffSplice.RI.signif), " RI events"))
        isoforms.RI <- altIntronRmats(diffSplice.RI.signif, exons=exons)
        orfChanges.RI <- transcriptChangeSummary(isoforms.RI[isoforms.RI$set == "spliced_intron"],
                                                 isoforms.RI[isoforms.RI$set == "retained_intron"],
                                                 BSgenome=BSgenome, NMD=NMD, exportGTF=NULL, dataSet=rmatsDataSet)
        m <- match(unlist(lapply(stringr::str_split(orfChanges.RI$id, "[-]"), "[[", 1)), diffSplice.RI.signif$ID)
        orfChanges <- rbind(orfChanges, cbind(diffSplice.RI.signif[m,c('ID', 'GeneID', 'geneSymbol', 'PValue','FDR', 'IncLevelDifference')], type="RI", orfChanges.RI))
        allTranscripts <- c(allTranscripts, isoforms.RI)
    }

    diffSplice.A3SS.signif <- extractEvent(rmatsDataSet, "A3SS")
    if(nrow(diffSplice.A3SS.signif) > 0 & "A3SS" %in% eventTypes){
        message(paste0("Creating isoforms for ", nrow(diffSplice.A3SS.signif), " A3SS events"))
        isoforms.A3SS <- altSpliceSiteRmats(diffSplice.A3SS.signif, eventType="A3SS", exons=exons)
        orfChanges.A3SS <- transcriptChangeSummary(isoforms.A3SS[isoforms.A3SS$set == "alt3_splicesite_long"],
                                                 isoforms.A3SS[isoforms.A3SS$set == "alt3_splicesite_short"],
                                                 BSgenome=BSgenome, NMD=NMD, exportGTF=NULL, dataSet=rmatsDataSet)
        m <- match(unlist(lapply(stringr::str_split(orfChanges.A3SS$id, "[-]"), "[[", 1)), diffSplice.A3SS.signif$ID)
        orfChanges <- rbind(orfChanges, cbind(diffSplice.A3SS.signif[m,c('ID', 'GeneID', 'geneSymbol', 'PValue','FDR', 'IncLevelDifference')], type="A3SS", orfChanges.A3SS))
        allTranscripts <- c(allTranscripts, isoforms.A3SS)
    }

    diffSplice.A5SS.signif <- extractEvent(rmatsDataSet, "A5SS")
    if(nrow(diffSplice.A5SS.signif) > 0 & "A5SS" %in% eventTypes){
        message(paste0("Creating isoforms for ", nrow(diffSplice.A5SS.signif), " A5SS events"))
        isoforms.A5SS <- altSpliceSiteRmats(diffSplice.A5SS.signif, eventType="A5SS", exons=exons)
        orfChanges.A5SS <- transcriptChangeSummary(isoforms.A5SS[isoforms.A5SS$set == "alt5_splicesite_long"],
                                                 isoforms.A5SS[isoforms.A5SS$set == "alt5_splicesite_short"],
                                                 BSgenome=BSgenome, NMD=NMD, exportGTF=NULL, dataSet=rmatsDataSet)
        m <- match(unlist(lapply(stringr::str_split(orfChanges.A5SS$id, "[-]"), "[[", 1)), diffSplice.A5SS.signif$ID)
        orfChanges <- rbind(orfChanges, cbind(diffSplice.A5SS.signif[m,c('ID', 'GeneID', 'geneSymbol', 'PValue','FDR', 'IncLevelDifference')], type="A5SS", orfChanges.A5SS))
        allTranscripts <- c(allTranscripts, isoforms.A5SS)
    }


    if(!is.null(exportGTF)){
        rtracklayer::export.gff(allTranscripts, con=exportGTF,
                                format="gtf")
    }

    return(orfChanges)
}

