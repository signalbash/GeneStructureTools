#' Annotate UTRs from Gencode GTF as 5' or 3'
#'
#' Annotate UTRs from Gencode GTF as 5' or 3'
#' @param gtf GRanges object of the GTF
#' @return gtf annotation GRanges object
#' @export
#' @import GenomicRanges
#' @examples
#' gtf_file <- system.file("extdata","gencode.vM14.annotation.small.gtf",
#' package = "GeneStructureTools")
#' gtf <- rtracklayer::import(gtf_file)
#' gtf <- UTR2UTR53(gtf)
#' table(gtf$type)
#' @author Beth Signal
UTR2UTR53 <- function(gtf){

    #transcripts with CDS/UTR annotated
    gtf_sub <- gtf[gtf@elementMetadata$type %in% c("UTR", "CDS")]
    gtf_sub <- gtf_sub[order(gtf_sub@elementMetadata$transcript_id, GenomicRanges::start(GenomicRanges::ranges(gtf_sub)))]

    #UTRs
    UTRTranscripts <- which(gtf_sub@elementMetadata$type == "UTR")
    id <- gtf_sub@elementMetadata$transcript_id[UTRTranscripts]

    #5'UTRs (+)
    pos5 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts] ==
                      gtf_sub@elementMetadata$transcript_id[UTRTranscripts +1] &
                      as.character(gtf_sub@strand[UTRTranscripts]) == "+" &
                      gtf_sub@elementMetadata$type[UTRTranscripts +1 ] %in% c("CDS", "exon"))
    #3'UTRs (+)
    pos3 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                      gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                      as.character(gtf_sub@strand[UTRTranscripts[-1]]) == "+" &
                      gtf_sub@elementMetadata$type[UTRTranscripts[-1] -1] %in% c("CDS","exon")) + 1

    #3'UTRs (-)
    neg3 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts] ==
                      gtf_sub@elementMetadata$transcript_id[UTRTranscripts +1] &
                      as.character(gtf_sub@strand[UTRTranscripts]) == "-" &
                      gtf_sub@elementMetadata$type[UTRTranscripts +1 ] %in% c("CDS","exon"))
    #5'UTRs (-)
    neg5 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1]] == gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                      as.character(gtf_sub@strand[UTRTranscripts[-1]]) == "-" &
                      gtf_sub@elementMetadata$type[UTRTranscripts[-1] -1] %in% c("CDS","exon")) + 1

    #new type var
    gtf_sub@elementMetadata$type2 <- as.character(gtf_sub@elementMetadata$type)
    gtf_sub@elementMetadata$type2[UTRTranscripts][pos5] <- "UTR5"
    gtf_sub@elementMetadata$type2[UTRTranscripts][neg5] <- "UTR5"
    gtf_sub@elementMetadata$type2[UTRTranscripts][pos3] <- "UTR3"
    gtf_sub@elementMetadata$type2[UTRTranscripts][neg3] <- "UTR3"

    while(any(gtf_sub@elementMetadata$type2[UTRTranscripts] == "UTR")){

        UTRTranscripts <- which(gtf_sub@elementMetadata$type2 == "UTR")

        #5'UTRs (+)
        pos5 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts] ==
                          gtf_sub@elementMetadata$transcript_id[UTRTranscripts +1] &
                          as.character(gtf_sub@strand[UTRTranscripts]) == "+" &
                          gtf_sub@elementMetadata$type2[UTRTranscripts +1 ] %in% c("CDS", "exon","UTR5"))
        #3'UTRs (+)
        pos3 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                          gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                          as.character(gtf_sub@strand[UTRTranscripts[-1]]) == "+" &
                          gtf_sub@elementMetadata$type2[UTRTranscripts[-1] -1] %in% c("CDS","exon","UTR3")) + 1

        #3'UTRs (-)
        neg3 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts] ==
                          gtf_sub@elementMetadata$transcript_id[UTRTranscripts +1] &
                          as.character(gtf_sub@strand[UTRTranscripts]) == "-" &
                          gtf_sub@elementMetadata$type2[UTRTranscripts +1 ] %in% c("CDS","exon","UTR3"))
        #5'UTRs (-)
        neg5 <- which(gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1]] ==
                          gtf_sub@elementMetadata$transcript_id[UTRTranscripts[-1] -1] &
                          as.character(gtf_sub@strand[UTRTranscripts[-1]]) == "-" &
                          gtf_sub@elementMetadata$type2[UTRTranscripts[-1] -1] %in% c("CDS","exon","UTR5")) + 1

        gtf_sub@elementMetadata$type2[UTRTranscripts][pos5] <- "UTR5"
        gtf_sub@elementMetadata$type2[UTRTranscripts][neg5] <- "UTR5"
        gtf_sub@elementMetadata$type2[UTRTranscripts][pos3] <- "UTR3"
        gtf_sub@elementMetadata$type2[UTRTranscripts][neg3] <- "UTR3"

        if(all(gtf_sub@elementMetadata$type2[UTRTranscripts] == "UTR")){
            gtf_sub@elementMetadata$type2[UTRTranscripts] <- "UTR_NA"
        }else{
            UTRTranscripts <- which(gtf_sub@elementMetadata$type2 == "UTR")
        }
    }

    gtf_sub_df <- data.frame(chr=gtf_sub@seqnames,
                             gtf_sub@ranges,
                             transcript_id=gtf_sub@elementMetadata$transcript_id)
    gtf_sub_names <- with(gtf_sub_df, paste(chr,start,end,transcript_id, sep="_"))
    gtf_df <- data.frame(chr=gtf@seqnames,
                         gtf@ranges,
                         transcript_id=gtf@elementMetadata$transcript_id)
    gtf_names <- with(gtf_df, paste(chr,start,end,transcript_id, sep="_"))

    m5 <- which(gtf_names %in% gtf_sub_names[gtf_sub@elementMetadata$type2 == "UTR5"])
    m3 <- which(gtf_names %in% gtf_sub_names[gtf_sub@elementMetadata$type2 == "UTR3"])

    gtf@elementMetadata$type <- as.character(gtf@elementMetadata$type)

    gtf@elementMetadata$type[m5] <- "UTR5"
    gtf@elementMetadata$type[m3] <- "UTR3"

    return(gtf)
}
