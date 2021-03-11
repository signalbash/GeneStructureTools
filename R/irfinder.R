readIrfDataSet <- function(path){

    irf <- new("irfDataSet", filePath=path)

    irfResultsFile = fread(path, data.table = F)
    irfResultsFile$FDR = p.adjust(irfResultsFile$`p-diff`, "fdr")
    irfResultsFile$psi_diff = irfResultsFile$`A-IRratio` - irfResultsFile$`B-IRratio`
    irfResultsFile$gene_name = unlist(lapply(str_split(irfResultsFile$`Intron-GeneName/GeneID`, "/"), "[[", 1))
    irfResultsFile$gene_id = unlist(lapply(str_split(irfResultsFile$`Intron-GeneName/GeneID`, "/"), "[[", 2))
    irfResultsFile$status = unlist(lapply(str_split(irfResultsFile$`Intron-GeneName/GeneID`, "/"), "[[", 3))
    irfResultsFile$intron_id = paste0(irfResultsFile$Chr, ":", irfResultsFile$Start, "-", irfResultsFile$End)


    events.RI = GRanges(seqnames = irfResultsFile$Chr,
                        ranges=IRanges(start=irfResultsFile$Start+1, end = irfResultsFile$End),
                        strand = irfResultsFile$Direction,
                        id = irfResultsFile$intron_id)

    slot(irf, "IRFresults") = irfResultsFile
    slot(irf, "coordinates") = events.RI

    return(irf)
}
#' Filter out significant events from a whippet diff comparison
#' @param irfDataSet whippetDataSet generated from \code{readWhippetDataSet()}
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
filterIrfEvents <- function(irfDataSet,
                              FDR = 0.05,
                              psiDelta = 0.1,
                              idList = NA){

    #set FDR/psiDelta if NA
    if(is.na(FDR[1])){
        FDR <- 1
    }
    if(is.na(psiDelta[1])){
        psiDelta <- 0
    }

    tmp <- slot(irfDataSet, "IRFresults")
    significantEventsIndex <- which(tmp$FDR <= FDR & abs(tmp$psi_diff) > psiDelta)
    slot(irfDataSet, "IRFresults") <- tmp[significantEventsIndex,]


    if(!is.na(idList[1])){
        tmp <- slot(irfDataSet, "IRFresults")
        significantEventsIndex <- which(tmp$gene_name %in% idList | tmp$gene_id %in% idList)
        slot(irfDataSet, "IRFresults") <- tmp[significantEventsIndex,]
    }

    tmpCoords = slot(irfDataSet, "coordinates")
    tmpCoords = tmpCoords[which(tmpCoords$id %in% slot(irfDataSet, "IRFresults")$intron_id),]
    slot(irfDataSet, "coordinates") = tmpCoords

    return(irfDataSet)

}


irfTranscriptChangeSummary <- function(irfDataSet = irf,
                                         BSgenome,
                                        intronMatchType = "exact",
                                         exons=NULL,
                                         NMD=TRUE,
                                         exportGTF=NULL){

    irfCoords = slot(irfDataSet, "coordinates")

    seqlevelsStyle(irfCoords) <- seqlevelsStyle(exons)[1]

    exons.intronRetention <- findIntronContainingTranscripts(input=irfCoords, exons, match=intronMatchType )
    isoforms.RI <- addIntronInTranscript(flankingExons=exons.intronRetention, exons, match="retain")

    orfChanges.RI <- transcriptChangeSummary(transcriptsX = isoforms.RI[isoforms.RI$set == "retained_intron"],
                                             transcriptsY = isoforms.RI[isoforms.RI$set == "spliced_intron"],
                                             BSgenome=g, NMD=NMD, exportGTF = NULL, dataSet=irfDataSet)

    irf.signif = irfResults(irfDataSet)
    m = match(orfChanges.RI$id, irf.signif$intron_id)

    orfChanges = cbind(irf.signif[m,c('intron_id', 'gene_id', 'gene_name', 'p-diff','FDR', 'psi_diff')], type = "RI", orfChanges.RI)


    if(!is.null(exportGTF)){
        rtracklayer::export.gff(isoforms.RI, con=exportGTF,
                                format="gtf")
    }

    return(orfChanges)
}




#' Class irfDataSet
#'
#' Class \code{irfDataSet} contains information read from irf output files
#'
#' @name irfDataSet-class
#' @rdname irfDataSet-class
#' @exportClass irfDataSet
#' @imports methods
setClass("irfDataSet", slots=list(coordinates="GRanges",
                                    IRFresults="data.frame",
                                    filePath="character"))

#' Method irfResults
#' @name irfResults
#' @rdname irfResults-methods
#' @exportMethod irfResults
#' @imports methods
#' @param irfDataSet irfDataSet generated from \code{readIRFDataSet()}
setGeneric("irfResults",
           def=function(irfDataSet)
           {
               standardGeneric("irfResults")
           }
)

#' @rdname irfResults-methods
#' @return differential splicing results data.frame
#' (originally from a whippet .diff file)
#' @family irf data processing
setMethod("irfResults", signature="irfDataSet",
          definition=function(irfDataSet)
          {
              return(slot(irfDataSet, "IRFresults"))
          }
)
