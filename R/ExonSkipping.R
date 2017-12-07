#' Given the location of a whole retained exon, find transcripts which can splice out this exon
#' @param exonRanges GRanges object with ranges for exons
#' @param gtf GRanges object made from a GTF
#' @param variableWidth How many nts overhang is allowed for finding matching exons (default = 0, i.e. complete match)
#' @param findIntrons Find transcripts where the event occurs within the intron?
#' @return data.frame with all overlapping exons
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
findExonContainingTranscripts <- function(exonRanges, gtf, variableWidth=0, findIntrons=FALSE){

    gtf.exons <- gtf[gtf$type=="exon"]

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
    if(variableWidth == 0){
        overlaps <- GenomicRanges::findOverlaps(exonRanges, gtf.exons, type="equal")
        gtf_equal <- gtf.exons[overlaps@to]
        gtf_equal$from <- overlaps@from
    }else{
        # find all overlaps
        overlaps <- GenomicRanges::findOverlaps(exonRanges, gtf.exons)
        gtf_equal <- gtf.exons[overlaps@to]
        gtf_equal$from <- overlaps@from
        start_diff <- abs(start(gtf_equal) - start(exonRanges[gtf_equal$from]))
        end_diff <- abs(end(gtf_equal) - end(exonRanges[gtf_equal$from]))
        total_diff <- start_diff + end_diff
        # remove overlaps with change in start+end greater than specified
        keep <- which(total_diff <= variableWidth)
        gtf_equal <- gtf_equal[keep]
    }


    gtf_equal$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_equal)), paste0(transcript_id, "_",from))
    gtf_equal$exon_number <- as.numeric(gtf_equal$exon_number)
    equal_exons <- as.data.frame(GenomicRanges::mcols(gtf_equal)[,c('gene_id',
                                                                    'transcript_id',
                                                                    'transcript_type',
                                                                    'from',
                                                                    'exon_number')])
    equal_exons$alt_id <- exonRanges$id[equal_exons$from]
    equal_exons$from <- NULL
    equal_exons$start <- start(gtf_equal)
    equal_exons$end <- end(gtf_equal)
    equal_exons$overlaps <- "exon"

    if(findIntrons == TRUE){
        # overlaps a transcript (i.e. can overlap an intron)
        overlaps <- GenomicRanges::findOverlaps(exonRanges, gtf[gtf$type=="transcript"])
        overlaps_df <- as.data.frame(overlaps)
        overlaps_df$from <- exonRanges$id[overlaps_df$queryHits]
        overlaps_df$to <- gtf[gtf$type=="transcript"]$transcript_id[overlaps_df$subjectHits]

        # annotate first/last exons (plz move elsewhere)
        if(!("first_last" %in% colnames(mcols(gtf.exons)))){
            t <- as.data.frame(table(gtf.exons$transcript_id))
            gtf.exons$first_last <- NA
            gtf.exons$first_last[gtf.exons$exon_number == 1] <- "first"
            gtf.exons$first_last[gtf.exons$exon_number == t$Freq[match(gtf.exons$transcript_id, t$Var1)]] <- "last"
        }
        # check that skipped exon doesn't overlap the first/last exon
        overlaps_exons <- GenomicRanges::findOverlaps(exonRanges, gtf.exons[which(gtf.exons$first_last %in% c("first","last"))])
        rm_transcripts <- gtf.exons$transcript_id[which(gtf.exons$first_last %in% c("first","last"))][overlaps_exons@to]
        overlaps_df <- overlaps_df[which(!(overlaps_df$to %in% rm_transcripts)),]

        # gtf with ALL exons where there is an intron overlap
        gtf_within <- gtf.exons[which(gtf.exons$transcript_id %in% overlaps_df$to)]
        gtf_within$from <- overlaps_df$from[match(gtf_within$transcript_id, overlaps_df$to)]
        gtf_within <- gtf_within[which(!(gtf_within$transcript_id %in% c(gtf_equal$transcript_id)))]

        if(length(gtf_within) > 0){
            gtf_within$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_within)), paste0(transcript_id, "_",from))
            gtf_within$exon_number <- 0

            overlaping_introns <- as.data.frame(GenomicRanges::mcols(gtf_within)[,c('gene_id',
                                                                            'transcript_id',
                                                                            'transcript_type',
                                                                            'from',
                                                                            'exon_number')])
            overlaping_introns <- overlaping_introns[which(!(duplicated(paste0(overlaping_introns$transcript_id, "-", overlaping_introns$from)))),]
            overlaping_introns$alt_id <- overlaping_introns$from
            overlaping_introns$from <- NULL
            overlaping_introns$start <- start(exonRanges[match(overlaping_introns$alt_id, exonRanges$id)])
            overlaping_introns$end <- end(exonRanges[match(overlaping_introns$alt_id, exonRanges$id)])
            overlaping_introns$overlaps <- "intron"
            equal_exons <- rbind(equal_exons, overlaping_introns)
        }
    }

    return(equal_exons)
}

#' Remove and iclude a skipped exon from the transcripts it overlaps
#' @param exonRanges GRanges object with ranges for exons
#' @param equal_exons data.frame generataed by findExonContainingTranscripts()
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param glueExons Join together exons that are not seperated by exons?
#' @param replaceVariableExons if overlapping exons are not exact matches for the event, should they be replaced by the event?
#' @return GRanges with transcripts skipping exons
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
skipExonInTranscript <- function(exonRanges, equal_exons, gtf.exons, glueExons=TRUE, replaceVariableExons=TRUE){

    exonRanges <- exonRanges[match(equal_exons$alt_id, exonRanges$id)]
    exonRanges$exon_number <- equal_exons$exon_number
    exonRanges$transcript_id <- equal_exons$transcript_id
    exonRanges$transcript_type <- equal_exons$transcript_type
    exonRanges$gene_id <- equal_exons$gene_id
    exonRanges$exon_id <- exonRanges$id
    exonRanges$overlaps <- equal_exons$overlaps

    if(replaceVariableExons == FALSE){
        old_starts <- start(exonRanges)
        old_ends <- end(exonRanges)
    }
    start(exonRanges) <- equal_exons$start
    end(exonRanges) <- equal_exons$end

    # transcripts containing the exon
    transcripts <- as.data.frame(table(equal_exons$transcript_id))
    gtf_transcripts <- gtf.exons[gtf.exons$transcript_id %in% transcripts$Var1]
    m <- match(gtf_transcripts$transcript_id, exonRanges$transcript_id)
    mcols(gtf_transcripts) <- cbind(mcols(gtf_transcripts), DataFrame(new_transcript_id=paste0(gtf_transcripts$transcript_id,"+AS ",exonRanges$exon_id[m])))
    mcols(exonRanges) <- cbind(mcols(exonRanges), DataFrame(new_transcript_id = paste0(exonRanges$transcript_id,"+AS ",exonRanges$exon_id)))

    mcols(gtf_transcripts) <- mcols(gtf_transcripts)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number', 'new_transcript_id')]
    mcols(exonRanges) <- mcols(exonRanges)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number','new_transcript_id', 'overlaps')]

    needs_dup <- which(!(exonRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))

    if(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts[gtf_transcripts$transcript_id %in% exonRanges$transcript_id[needs_dup]]
    }

    while(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts_add[gtf_transcripts_add$transcript_id %in% exonRanges$transcript_id[needs_dup]]
        m <- match(gtf_transcripts_add$transcript_id, exonRanges$transcript_id[needs_dup])
        gtf_transcripts_add$new_transcript_id <- paste0(gtf_transcripts_add$transcript_id,"+AS ",exonRanges$exon_id[needs_dup][m])
        gtf_transcripts <- c(gtf_transcripts, gtf_transcripts_add)
        needs_dup <- which(!(exonRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))
    }

    exon_names <- with(as.data.frame(gtf_transcripts), paste0(seqnames, ":", start,"-",end))

    ol <- findOverlaps(gtf_transcripts, exonRanges, type="equal")
    ol <- as.data.frame(ol)
    ol$gtf_trans_id <- gtf_transcripts$new_transcript_id[ol$queryHits]
    ol$exonRanges_id <- exonRanges$new_transcript_id[ol$subjectHits]
    ol <- ol[ol$gtf_trans_id == ol$exonRanges_id,]

    rm <- unique(ol$queryHits)
    gtf_transcripts_rm <- gtf_transcripts[-rm]

    mcols(gtf_transcripts_rm) <- mcols(gtf_transcripts_rm)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts_rm))[2] <- "transcript_id"

    mcols(exonRanges) <- mcols(exonRanges)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number','overlaps')]
    colnames(mcols(exonRanges))[2] <- "transcript_id"

    gtf_transcripts_rm$exon_number <- as.numeric(gtf_transcripts_rm$exon_number)
    order <- order(gtf_transcripts_rm$transcript_id, gtf_transcripts_rm$exon_number)
    gtf_transcripts_rm <- gtf_transcripts_rm[order]
    gtf_transcripts_rm$overlaps <- exonRanges$overlaps[match(gtf_transcripts_rm$transcript_id, exonRanges$new_transcript_id)]

    gtf_transcripts_withExon <- gtf_transcripts_rm
    if(replaceVariableExons == FALSE){
        start(exonRanges) <- old_starts
        end(exonRanges) <- old_ends
    }
    gtf_transcripts_withExon <- c(gtf_transcripts_withExon, exonRanges)
    gtf_transcripts_withExon <- reorderExonNumbers(gtf_transcripts_withExon)

    gtf_transcripts_rm$set <- "skipped_exon"
    gtf_transcripts_withExon$set <- "included_exon"

    gtf_transcripts_rm <- c(gtf_transcripts_rm, gtf_transcripts_withExon)

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

#' Reorder the exon numbers in a gtf annotation
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param by what column are the transcripts grouped by?
#' @return The same input GRanges, but with exon numbers reordered.
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
reorderExonNumbers <- function(gtf.exons, by="transcript_id"){
    n <- which(colnames(mcols(gtf.exons)) == by)

    order <- order(mcols(gtf.exons)[,n], start(gtf.exons))

    gtf.exons <- gtf.exons[order]

    transcript_table <- as.data.frame(table(mcols(gtf.exons)[,n]))
    transcript_table$strand <- as.character(strand(gtf.exons[match(transcript_table$Var1, mcols(gtf.exons)[,n])]))

    gtf.exons$exon_number <- unlist(apply(transcript_table, 1, function(x) if(x[3] == "+"){c(1:(x[2]))}else{c((x[2]:1))}))
    return(gtf.exons)
}

