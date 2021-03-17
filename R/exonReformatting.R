#' reformat a exons GRanges with first/last annotations and force 'biotype' to 'type'
#' @param exons reference exons GRanges
#' @return reference exons GRanges
#' @export
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
#' @examples
reformatExons = function(exons){

    # add first/last annotation (speeds up later steps)
    if(!("first_last" %in% colnames(mcols(exons)))){
        t <- as.data.frame(table(exons$transcript_id))
        exons$first_last <- NA
        exons$first_last[exons$exon_number == 1] <- "first"
        exons$first_last[exons$exon_number == t$Freq[match(exons$transcript_id, t$Var1)]] <- "last"
    }

    colnames(mcols(exons)) <- gsub("biotype", "type", colnames(mcols(exons)))
    return(exons)
}


#' Generate a introns Granges from an exons Granges
#' @param exons reference exons GRanges
#' @return reference introns GRanges
#' @export
#' @import methods
#' @importFrom dplyr lead
#' @importFrom dplyr left_join
#' @importFrom rlang .data
#' @family rmats data processing
#' @author Beth Signal
#' @examples
exonsToIntrons = function(exons){

    exons_df <- as.data.frame(exons)
    exons_df <- exons_df[,c("seqnames", "start", "end", "strand", "transcript_id", "exon_number")]
    exons_df <- arrange(exons_df, transcript_id, start, end)

    exons_df$intron_start <- exons_df$end
    exons_df$intron_end <- dplyr::lead(exons_df$start)

    rm <- which(dplyr::lead(exons_df$transcript_id) != exons_df$transcript_id)

    exons_df <- exons_df[-rm,]
    exons_df <- exons_df[-nrow(exons_df),]
    min_exon_n <- aggregate(exon_number ~ transcript_id, exons_df, min)
    if(!all(min_exon_n$exon_number == 1)){
        exons_df$exon_number[exons_df$strand=="-"] <- as.numeric(exons_df$exon_number[exons_df$strand=="-"]) - 1
    }

    introns <- GRanges(seqnames=exons_df$seqnames, ranges=IRanges(start=exons_df$intron_start, end=exons_df$intron_end),
                       strand=exons_df$strand, transcript_id=exons_df$transcript_id, exon_number=exons_df$exon_number)
    m <- match(introns$transcript_id, exons$transcript_id)
    introns$gene_id <- exons$gene_id[m]
    introns$gene_name <- exons$gene_name[m]

    mcols(introns) <- DataFrame(dplyr::left_join(as.data.frame(mcols(introns)),  as.data.frame(mcols(exons))))
    return(introns)

}
#' Convert an exon-level gtf annotation to a transcript-level gtf annotation
#'
#' @param exons GRanges object with exons
#' @return GRanges object with transcripts
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import
#' @family gtf manipulation
#' @author Beth Signal
#' @examples
exonsToTranscripts <- function(exons){

    exons.df <- data.frame(t_id=exons$transcript_id, start=start(exons), end=end(exons))
    txRanges <- aggregate(start ~ t_id, exons.df, min)
    txRanges$end <- aggregate(end ~ t_id, exons.df, max)[,2]

    transcripts <- exons[!duplicated(exons$transcript_id),]
    transcripts$exon_id <- NULL
    transcripts$exon_number <- NULL
    transcripts$first_last <- NULL

    start(transcripts) <- txRanges$start[match(transcripts$transcript_id, txRanges$t_id)]
    end(transcripts) <- txRanges$end[match(transcripts$transcript_id, txRanges$t_id)]

    return(transcripts)


    return(transcripts)
}

#' Remove version number from ensembl gene/transcript ids
#'
#' @param ids vector of ensembl ids
#' @return vector of ensembl ids without the version number
#' @import stringr
#' @export
#' @author Beth Signal
#' @examples
#' removeVersion("ENSMUSG00000001017.15")

removeVersion <- function(ids){
    return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[", 1)))
}
