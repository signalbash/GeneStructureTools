#' Annotate UTRs from Gencode GTF as 5' or 3'
#'
#' Annotate UTRs from Gencode GTF as 5' or 3'
#' @param gtf GRanges object of the GTF
#' @return gtf annotation GRanges object
#' @export
#' @import GenomicRanges
#' @examples
#' gtfFile <- system.file("extdata","gencode.vM14.annotation.small.gtf",
#' package = "GeneStructureTools")
#' gtf <- rtracklayer::import(gtfFile)
#' gtf <- UTR2UTR53(gtf)
#' table(gtf$type)
#' @author Beth Signal
UTR2UTR53 <- function(gtf){

    #transcripts with CDS/UTR annotated
    gtf.cdsutr <- gtf[gtf@elementMetadata$type %in% c("UTR", "CDS")]
    gtf.cdsutr <- gtf.cdsutr[order(gtf.cdsutr@elementMetadata$transcript_id,
                                   GenomicRanges::start(GenomicRanges::ranges(gtf.cdsutr)))]

    #UTRs
    UTRTranscripts <- which(gtf.cdsutr@elementMetadata$type == "UTR")
    id <- gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts]

    #5'UTRs (+)
    pos5 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts] ==
                      gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts +1] &
                      as.character(gtf.cdsutr@strand[UTRTranscripts]) == "+" &
                      gtf.cdsutr@elementMetadata$type[UTRTranscripts +1 ] %in%
                      c("CDS", "exon"))
    #3'UTRs (+)
    pos3 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                      gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                      as.character(gtf.cdsutr@strand[UTRTranscripts[-1]]) == "+" &
                      gtf.cdsutr@elementMetadata$type[UTRTranscripts[-1] -1] %in%
                      c("CDS","exon")) + 1

    #3'UTRs (-)
    neg3 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts] ==
                      gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts +1] &
                      as.character(gtf.cdsutr@strand[UTRTranscripts]) == "-" &
                      gtf.cdsutr@elementMetadata$type[UTRTranscripts +1 ] %in%
                      c("CDS","exon"))
    #5'UTRs (-)
    neg5 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                      gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                      as.character(gtf.cdsutr@strand[UTRTranscripts[-1]]) == "-" &
                      gtf.cdsutr@elementMetadata$type[UTRTranscripts[-1] -1] %in%
                      c("CDS","exon")) + 1

    #new type var
    gtf.cdsutr@elementMetadata$type2 <- as.character(gtf.cdsutr@elementMetadata$type)
    gtf.cdsutr@elementMetadata$type2[UTRTranscripts][pos5] <- "UTR5"
    gtf.cdsutr@elementMetadata$type2[UTRTranscripts][neg5] <- "UTR5"
    gtf.cdsutr@elementMetadata$type2[UTRTranscripts][pos3] <- "UTR3"
    gtf.cdsutr@elementMetadata$type2[UTRTranscripts][neg3] <- "UTR3"

    while(any(gtf.cdsutr@elementMetadata$type2[UTRTranscripts] == "UTR")){

        UTRTranscripts <- which(gtf.cdsutr@elementMetadata$type2 == "UTR")

        #5'UTRs (+)
        pos5 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts] ==
                          gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts +1] &
                          as.character(gtf.cdsutr@strand[UTRTranscripts]) == "+" &
                          gtf.cdsutr@elementMetadata$type2[UTRTranscripts +1 ] %in%
                          c("CDS", "exon","UTR5"))
        #3'UTRs (+)
        pos3 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                          gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                          as.character(gtf.cdsutr@strand[UTRTranscripts[-1]]) == "+" &
                          gtf.cdsutr@elementMetadata$type2[UTRTranscripts[-1] -1] %in%
                          c("CDS","exon","UTR3")) + 1

        #3'UTRs (-)
        neg3 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts] ==
                          gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts +1] &
                          as.character(gtf.cdsutr@strand[UTRTranscripts]) == "-" &
                          gtf.cdsutr@elementMetadata$type2[UTRTranscripts +1 ] %in%
                          c("CDS","exon","UTR3"))
        #5'UTRs (-)
        neg5 <- which(gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                          gtf.cdsutr@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                          as.character(gtf.cdsutr@strand[UTRTranscripts[-1]]) == "-" &
                          gtf.cdsutr@elementMetadata$type2[UTRTranscripts[-1] -1] %in%
                          c("CDS","exon","UTR5")) + 1

        gtf.cdsutr@elementMetadata$type2[UTRTranscripts][pos5] <- "UTR5"
        gtf.cdsutr@elementMetadata$type2[UTRTranscripts][neg5] <- "UTR5"
        gtf.cdsutr@elementMetadata$type2[UTRTranscripts][pos3] <- "UTR3"
        gtf.cdsutr@elementMetadata$type2[UTRTranscripts][neg3] <- "UTR3"

        if(all(gtf.cdsutr@elementMetadata$type2[UTRTranscripts] == "UTR")){
            gtf.cdsutr@elementMetadata$type2[UTRTranscripts] <- "UTR_NA"
        }else{
            UTRTranscripts <- which(gtf.cdsutr@elementMetadata$type2 == "UTR")
        }
    }

    gtf.cdsutrDF <- as.data.frame(gtf.cdsutr)

    gtf.cdsutrNames <- paste(gtf.cdsutrDF$seqnames,
                             gtf.cdsutrDF$start,
                             gtf.cdsutrDF$end,
                             gtf.cdsutrDF$transcript_id, sep="_")
    #gtf.cdsutrNames <- with(gtf.cdsutrDF, paste(seqnames,start,end,transcript_id, sep="_"))
    gtfDF <- as.data.frame(gtf)
    gtfNames <- paste(gtfDF$seqnames,gtfDF$start,gtfDF$end,gtfDF$transcript_id, sep="_")
    #gtfNames <- with(gtfDF, paste(seqnames,start,end,transcript_id, sep="_"))

    m5 <- which(gtfNames %in% gtf.cdsutrNames[gtf.cdsutr@elementMetadata$type2 == "UTR5"])
    m3 <- which(gtfNames %in% gtf.cdsutrNames[gtf.cdsutr@elementMetadata$type2 == "UTR3"])

    gtf@elementMetadata$type <- as.character(gtf@elementMetadata$type)

    gtf@elementMetadata$type[m5] <- "UTR5"
    gtf@elementMetadata$type[m3] <- "UTR3"

    return(gtf)
}
