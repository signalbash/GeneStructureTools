#' Given the location of a whole retained intron, find transcripts which splice out this intron
#' @param intronRanges GRanges object with ranges for introns
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @return data.frame with all flanking exon pairs
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
findIntronContainingTranscripts <- function(intronRanges, gtf.exons){

    # remove any duplicates
    overlaps <- GenomicRanges::findOverlaps(intronRanges, type="equal")
    overlaps <- overlaps[which(overlaps@from != overlaps@to)]
    if(length(overlaps) > 0){
        overlaps <- overlaps[which(overlaps@from < overlaps@to)]
        if(length(overlaps) > 0){
            duplicates <- unique(overlaps@to)
            intronRanges <- intronRanges[-duplicates]
        }
    }

    # start of intron // end of exon a
    IR_range_start <- intronRanges
    end(IR_range_start) <- start(IR_range_start)

    overlaps <- GenomicRanges::findOverlaps(IR_range_start, gtf.exons, type="end")
    gtf_from_a <- gtf.exons[overlaps@to]
    gtf_from_a$from <- overlaps@from
    gtf_from_a$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_from_a)), paste0(transcript_id, "_",from))

    # end of intron // start of exon b
    IR_range_end <- intronRanges
    start(IR_range_end) <- end(IR_range_end)

    overlaps <- GenomicRanges::findOverlaps(IR_range_end, gtf.exons, type="start")
    gtf_to_a <- gtf.exons[overlaps@to]
    gtf_to_a$from <- overlaps@from
    gtf_to_a$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_to_a)), paste0(transcript_id, "_",from))

    keep_from <- which(gtf_from_a@elementMetadata$new_id %in% gtf_to_a@elementMetadata$new_id)
    gtf_from_a <- gtf_from_a[keep_from]

    m <- match(gtf_from_a$new_id, gtf_to_a$new_id)
    gtf_to_a <- gtf_to_a[m]

    gtf_to_a$from_exon_number <- as.numeric(gtf_from_a$exon_number)
    gtf_to_a$to_exon_number <- as.numeric(gtf_to_a$exon_number)

    gtf_to_a$intron_exon_number <- apply(GenomicRanges::mcols(gtf_to_a)[,c('to_exon_number','from_exon_number')], 1, mean)

    flanking_exons <- as.data.frame(GenomicRanges::mcols(gtf_to_a)[,c('gene_id','transcript_id','transcript_type',
                                                       'from','from_exon_number',
                                                       'intron_exon_number','to_exon_number')])

    flanking_exons$from <- intronRanges$id[flanking_exons$from]

    return(flanking_exons)
}

#' Add a retained intron to the transcripts it is skipped by
#' @param intronRanges GRanges object with ranges for introns
#' @param flanking_exons data.frame generataed by findIntronContainingTranscripts()
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param glueExons Join together exons that are not seperated by introns?
#' @return GRanges with transcripts containing retained introns
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
addIntronInTranscript <- function(intronRanges = IR_range, flanking_exons, gtf.exons, glueExons=TRUE){

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

    if(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts[gtf_transcripts$transcript_id %in% intronRanges$transcript_id[needs_dup]]
    }

    while(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts_add[gtf_transcripts_add$transcript_id %in% intronRanges$transcript_id[needs_dup]]
        m <- match(gtf_transcripts_add$transcript_id, intronRanges$transcript_id[needs_dup])
        gtf_transcripts_add$new_transcript_id <- paste0(gtf_transcripts_add$transcript_id,"+INTRON ",intronRanges$exon_id[needs_dup][m])
        gtf_transcripts <- c(gtf_transcripts, gtf_transcripts_add)
        needs_dup <- which(!(intronRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))
    }

    gtf_transcripts_all <- c(gtf_transcripts, intronRanges)
    mcols(gtf_transcripts_all) <- mcols(gtf_transcripts_all)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts_all))[2] <- "transcript_id"

    gtf_transcripts_all$exon_number <- as.numeric(gtf_transcripts_all$exon_number)
    order <- order(gtf_transcripts_all$transcript_id, gtf_transcripts_all$exon_number)
    gtf_transcripts_all <- gtf_transcripts_all[order]

    #join together exons that are not seperated by an intron
    if(glueExons==TRUE){

        # split pos/neg ordering
        gtf_transcripts_all.neg <-
            gtf_transcripts_all[which(as.logical(strand(gtf_transcripts_all) == "-"))]
        order <- order(gtf_transcripts_all.neg$transcript_id,
                       plyr::desc(gtf_transcripts_all.neg$exon_number))
        gtf_transcripts_all.neg <- gtf_transcripts_all.neg[order]

        gtf_transcripts_all.pos <-
            gtf_transcripts_all[which(as.logical(strand(gtf_transcripts_all) == "+"))]
        gtf_transcripts_all <- c(gtf_transcripts_all.pos,gtf_transcripts_all.neg)

        #extend starts <---<---<---
        w <- which(end(ranges(gtf_transcripts_all))[-length(gtf_transcripts_all)] == start(ranges(gtf_transcripts_all[-1])))
        gtf_transcripts_all <- gtf_transcripts_all
        GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_all))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_all))[w]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_all, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to - 1]
        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_all <- gtf_transcripts_all[-rm]
        }

        #extend ends --->--->--->
        overlaps <- findOverlaps(gtf_transcripts_all)
        overlaps <- overlaps[overlaps@from==overlaps@to +1]

        GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_all))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_all))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_all, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_all <- gtf_transcripts_all[-rm]
        }

        #order <- order(gtf_transcripts_all$transcript_id, gtf_transcripts_all$exon_number)
        #gtf_transcripts_all <- gtf_transcripts_all[order]

    }
    return(gtf_transcripts_all)
}

