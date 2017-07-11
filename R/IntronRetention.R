#' Given the location of a whole retained intron, find transcripts which splice out this intron
#' @param intronRanges GRanges object with ranges for introns
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @return data.frame with all flanking exon pairs
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
findIntronContainingTranscripts <- function(intronRanges, gtf.exons){

    # start of intron // end of exon a
    IR_range_start <- intronRanges
    end(IR_range_start) <- start(IR_range_start)

    overlaps <- findOverlaps(IR_range_start, gtf.exons, type="end")
    gtf_from_a <- gtf.exons[overlaps@to]
    gtf_from_a$from <- overlaps@from
    gtf_from_a$new_id <- with(mcols(gtf_from_a), paste0(transcript_id, "_",from))

    # end of intron // start of exon b
    IR_range_end <- intronRanges
    start(IR_range_end) <- end(IR_range_end)

    overlaps <- findOverlaps(IR_range_end, gtf.exons, type="start")
    gtf_to_a <- gtf.exons[overlaps@to]
    gtf_to_a$from <- overlaps@from
    gtf_to_a$new_id <- with(mcols(gtf_to_a), paste0(transcript_id, "_",from))

    keep_from <- which(gtf_from_a@elementMetadata$new_id %in% gtf_to_a@elementMetadata$new_id)
    gtf_from_a <- gtf_from_a[keep_from]

    m <- match(gtf_from_a$new_id, gtf_to_a$new_id)
    gtf_to_a <- gtf_to_a[m]

    gtf_to_a$from_exon_number <- as.numeric(gtf_from_a$exon_number)
    gtf_to_a$to_exon_number <- as.numeric(gtf_to_a$exon_number)

    gtf_to_a$intron_exon_number <- apply(mcols(gtf_to_a)[,c('to_exon_number','from_exon_number')], 1, mean)

    flanking_exons <- as.data.frame(mcols(gtf_to_a)[,c('gene_id','transcript_id','transcript_type',
                                                       'from','from_exon_number',
                                                       'intron_exon_number','to_exon_number')])

    flanking_exons$from <- intronRanges$id[flanking_exons$from]
    return(flanking_exons)
}

#' Add a retained intron to the transcripts it is skipped by
#' @param intronRanges GRanges object with ranges for introns
#' @param flanking_exons data.frame generataed by findIntronContainingTranscripts()
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @return GRanges with transcripts containing retained introns
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
addIntronInTranscript <- function(intronRanges = IR_range, flanking_exons, gtf.exons){

    intronRanges <- intronRanges[match(flanking_exons$from, intronRanges$id)]
    intronRanges$exon_number <- flanking_exons$intron_exon_number
    intronRanges$transcript_id <- flanking_exons$transcript_id
    intronRanges$transcript_type <- flanking_exons$transcript_type
    intronRanges$gene_id <- flanking_exons$gene_id
    intronRanges$exon_id <- intronRanges$id

    transcripts <- as.data.frame(table(flanking_exons$transcript_id))
    gtf_transcripts <- gtf.exons[gtf.exons$transcript_id %in% transcripts$Var1]
    m <- match(gtf_transcripts$transcript_id, intronRanges$transcript_id)
    mcols(gtf_transcripts) <- cbind(mcols(gtf_transcripts), DataFrame(new_transcript_id=paste0(gtf_transcripts$transcript_id,"+INTRON ",intronRanges$exon_id[m])))
    mcols(intronRanges) <- cbind(mcols(intronRanges), DataFrame(new_transcript_id = paste0(intronRanges$transcript_id,"+INTRON ",intronRanges$exon_id)))

    #mcols(gtf_transcripts)$new_transcript_id <- paste0(gtf_transcripts$transcript_id,"+INTRON ",intronRanges$exon_id[m])
    #mcols(intronRanges)$new_transcript_id <- paste0(intronRanges$transcript_id,"+INTRON ",intronRanges$exon_id)

    mcols(gtf_transcripts) <- mcols(gtf_transcripts)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number', 'new_transcript_id')]
    mcols(intronRanges) <- mcols(intronRanges)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number','new_transcript_id')]

    needs_dup <- which(!(intronRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))

    while(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts[gtf_transcripts$transcript_id %in% intronRanges$transcript_id[needs_dup]]
        m <- match(gtf_transcripts_add$transcript_id, intronRanges$transcript_id[needs_dup])
        gtf_transcripts_add$new_transcript_id <- paste0(gtf_transcripts_add$transcript_id,"+INTRON ",intronRanges$exon_id[needs_dup][m])
        gtf_transcripts <- c(gtf_transcripts, gtf_transcripts_add)
        needs_dup <- which(!(intronRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))
    }

    gtf_transcripts_all <- c(gtf_transcripts, intronRanges)
    mcols(gtf_transcripts_all) <- mcols(gtf_transcripts_all)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts_all))[2] <- "transcript_id"
    return(gtf_transcripts_all)
}


#' Find the largest distance between two vectors of numbers
#' Helper function for get_orfs
#' @param start_site vector of start sites - i.e Met amino acid positions
#' @param stop_site vector of stop sites - i.e Stop (*) amino acid positions
#' @param longest 1
#' @return sequential start site and end site with the greatest difference
#' @export
#' @import plyr
#' @examples
#' @author Beth Signal
maxLocation <- function(start_site, stop_site, longest = 1){
    if(length(start_site) > 0){
        diffs <- unlist(lapply(start_site, function(x) stop_site[which(stop_site > x)[1]] - x))
        order <- order(diffs, decreasing=TRUE)
        max_loc <- order[longest]
        start <- start_site[max_loc]
        stop <- stop_site[which(stop_site > start)[1]]
        return(c(start,stop))
    }else{
        return(c(NA,NA))
    }
}

#' get distance to second exon-exon junction
#' @param diffs distances to junctions
#' @return distance to closest second exon-exon junction
#' @export
#' @examples
#' @author Beth Signal
distance_to_junction_b <- function(diffs){
    keep <- which(diffs <= 0)
    if(length(keep) > 1){
        return(min(abs(diffs[diffs <= 0])))
    }else{
        return(NA)
    }
}

#' Get open reading frames for transcripts
#' @param transcripts GRanges object with ONLY exon annotations (no gene, transcript, CDS etc.) with all transcripts for orf retrevial
#' @param BSgenome BSgenome object
#' @param returnLongestOnly TRUE - currently only returns longest ORF
#' @return data.frame with longest orf details
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @examples
#' @author Beth Signal
get_orfs <- function(transcripts, BSgenome = g, returnLongestOnly=TRUE){

    transcripts$exon_number <- as.numeric(transcripts$exon_number)
    order <- order(transcripts$transcript_id, transcripts$exon_number)
    transcripts <- transcripts[order]
    transcripts$seq <- as.character(Biostrings::getSeq(g, transcripts))

    seq_cat <- aggregate(seq ~ transcript_id, mcols(transcripts), function(x) (paste(x, collapse="")))
    ids <- as.character(seq_cat$transcript_id)
    seq_cat <- seq_cat$seq

    # 3 frames
    seq_cat <- c(seq_cat, str_sub(seq_cat, 2), str_sub(seq_cat, 3))
    frames <- rep(c(1,2,3), each = length(ids))
    ids <- c(ids,ids,ids)

    orf <- suppressWarnings(unlist(lapply(seq_cat, function(x) as.character(Biostrings::translate(Biostrings::DNAString(x))))))

    orf_df <- data.frame(id = ids,
                         aa_sequence = orf,
                         frame = frames,
                         stringsAsFactors = FALSE
    )

    orf_df$seq_length <- nchar(orf_df$aa_sequence)

    start_sites <- str_locate_all(orf_df$aa_sequence, "M")
    start_sites <-
        lapply(start_sites, function(x) as.numeric(x[,2]))
    stop_sites <- str_locate_all(orf_df$aa_sequence, "[*]")
    stop_sites <-
        mapply(function(x,y) c(as.numeric(x[,2]), nchar(y)), stop_sites, orf_df$aa_sequence)

    max_loc <- mapply(function(x,y) maxLocation(x,y), start_sites, stop_sites)

    orf_df$start_site <- max_loc[1,]
    orf_df$stop_site <- max_loc[2,]
    orf_df$orf_sequence <- str_sub(orf_df$aa_sequence, orf_df$start_site, orf_df$stop_site - 1)
    orf_df$orf_length <- nchar(orf_df$orf_sequence)

    widths <- data.frame(w = width(transcripts),
                         id = transcripts$transcript_id)
    widths <- aggregate(w ~ id, widths, function(x) cumsum(c(1,x))[-1])

    m <- match(orf_df$id, widths$id)
    diffs <- mapply(function(x , y) x*3 - y, orf_df$stop_site, widths$w[m])

    orf_df$min_dist_to_junction_a <- suppressWarnings(unlist(lapply(diffs, function(x) min(x[x >= 0]))))
    orf_df$min_dist_to_junction_a[is.infinite(orf_df$min_dist_to_junction_a)] <- NA
    orf_df$min_dist_to_junction_b <- unlist(lapply(diffs, distance_to_junction_b)) - orf_df$frame
    orf_df$aa_sequence <- NULL

    if(returnLongestOnly == TRUE){
        orf_df <- plyr::arrange(orf_df, plyr::desc(orf_length))
        orf_df <- orf_df[!duplicated(orf_df$id),]
    }

    orf_df <- plyr::arrange(orf_df, id)
    return(orf_df)

}
