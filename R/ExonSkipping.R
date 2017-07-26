#' Given the location of a whole retained exon, find transcripts which can splice out this exon
#' @param exonRanges GRanges object with ranges for exons
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @return data.frame with all overlapping exons
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
findExonContainingTranscripts <- function(exonRanges, gtf.exons){

    # remove any duplicates
    overlaps <- GenomicRanges::findOverlaps(exonRanges, type="equal")
    overlaps <- overlaps[which(overlaps@from != overlaps@to)]
    if(length(overlaps) > 0){
        overlaps <- overlaps[which(overlaps@from < overlaps@to)]
        if(length(overlaps) > 0){
            duplicates <- unique(overlaps@to)
            exonRanges <- exonRanges[-duplicates]
        }
    }


    # whole match

    overlaps <- GenomicRanges::findOverlaps(exonRanges, gtf.exons, type="equal")
    gtf_equal <- gtf.exons[overlaps@to]
    gtf_equal$from <- overlaps@from
    gtf_equal$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_equal)), paste0(transcript_id, "_",from))

    gtf_equal$exon_number <- as.numeric(gtf_equal$exon_number)

    equal_exons <- as.data.frame(GenomicRanges::mcols(gtf_equal)[,c('gene_id',
                                                                       'transcript_id',
                                                                       'transcript_type',
                                                                       'from',
                                                                       'exon_number')])

    equal_exons$from <- exonRanges$id[equal_exons$from]

    return(equal_exons)
}

#' Remove a skipped exon from the transcripts it is contained in
#' @param exonRanges GRanges object with ranges for exons
#' @param equal_exons data.frame generataed by findExonContainingTranscripts()
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param glueExons Join together exons that are not seperated by exons?
#' @return GRanges with transcripts skipping exons
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
removeExonInTranscript <- function(exonRanges, equal_exons, gtf.exons, glueExons=TRUE){

    exonRanges <- exonRanges[match(equal_exons$from, exonRanges$id)]
    exonRanges$exon_number <- equal_exons$exon_number
    exonRanges$transcript_id <- equal_exons$transcript_id
    exonRanges$transcript_type <- equal_exons$transcript_type
    exonRanges$gene_id <- equal_exons$gene_id
    exonRanges$exon_id <- exonRanges$id

    transcripts <- as.data.frame(table(equal_exons$transcript_id))
    gtf_transcripts <- gtf.exons[gtf.exons$transcript_id %in% transcripts$Var1]
    m <- match(gtf_transcripts$transcript_id, exonRanges$transcript_id)
    mcols(gtf_transcripts) <- cbind(mcols(gtf_transcripts), DataFrame(new_transcript_id=paste0(gtf_transcripts$transcript_id,"-EXON ",exonRanges$exon_id[m])))
    mcols(exonRanges) <- cbind(mcols(exonRanges), DataFrame(new_transcript_id = paste0(exonRanges$transcript_id,"-EXON ",exonRanges$exon_id)))

    mcols(gtf_transcripts) <- mcols(gtf_transcripts)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number', 'new_transcript_id')]
    mcols(exonRanges) <- mcols(exonRanges)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number','new_transcript_id')]

    needs_dup <- which(!(exonRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))

    if(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts[gtf_transcripts$transcript_id %in% exonRanges$transcript_id[needs_dup]]
    }

    while(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts_add[gtf_transcripts_add$transcript_id %in% exonRanges$transcript_id[needs_dup]]
        m <- match(gtf_transcripts_add$transcript_id, exonRanges$transcript_id[needs_dup])
        gtf_transcripts_add$new_transcript_id <- paste0(gtf_transcripts_add$transcript_id,"-EXON ",exonRanges$exon_id[needs_dup][m])
        gtf_transcripts <- c(gtf_transcripts, gtf_transcripts_add)
        needs_dup <- which(!(exonRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))
    }

    exon_names <- with(as.data.frame(gtf_transcripts), paste0(seqnames, ":", start,"-",end))
    rm <- which(unlist(lapply(stringr::str_split(gtf_transcripts$new_transcript_id, "-EXON "), "[[", 2)) == exon_names)
    gtf_transcripts_rm <- gtf_transcripts[-rm]

    mcols(gtf_transcripts_rm) <- mcols(gtf_transcripts_rm)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts_rm))[2] <- "transcript_id"
    gtf_transcripts_rm$exon_number <- as.numeric(gtf_transcripts_rm$exon_number)
    order <- order(gtf_transcripts_rm$transcript_id, gtf_transcripts_rm$exon_number)
    gtf_transcripts_rm <- gtf_transcripts_rm[order]

    #join together exons that are not seperated by an exon
    if(glueExons==TRUE){

        # split pos/neg ordering
        gtf_transcripts_rm.neg <-
            gtf_transcripts_rm[which(as.logical(strand(gtf_transcripts_rm) == "-"))]
        order <- order(gtf_transcripts_rm.neg$transcript_id,
                       plyr::desc(gtf_transcripts_rm.neg$exon_number))
        gtf_transcripts_rm.neg <- gtf_transcripts_rm.neg[order]

        gtf_transcripts_rm.pos <-
            gtf_transcripts_rm[which(as.logical(strand(gtf_transcripts_rm) == "+"))]
        gtf_transcripts_rm <- c(gtf_transcripts_rm.pos,gtf_transcripts_rm.neg)

        #extend starts <---<---<---
        w <- which(end(ranges(gtf_transcripts_rm))[-length(gtf_transcripts_rm)] == start(ranges(gtf_transcripts_rm[-1])))
        gtf_transcripts_rm <- gtf_transcripts_rm
        GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_rm))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_rm))[w]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_rm, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to - 1]
        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_rm <- gtf_transcripts_rm[-rm]
        }

        #extend ends --->--->--->
        overlaps <- findOverlaps(gtf_transcripts_rm)
        overlaps <- overlaps[overlaps@from==overlaps@to +1]

        GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_rm))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_rm))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_rm, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_rm <- gtf_transcripts_rm[-rm]
        }

        #order <- order(gtf_transcripts_rm$transcript_id, gtf_transcripts_rm$exon_number)
        #gtf_transcripts_rm <- gtf_transcripts_rm[order]

    }
    return(gtf_transcripts_rm)
}

