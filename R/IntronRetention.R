#' Given the location of a whole retained intron, find transcripts which splice out this intron
#' @param intronRanges GRanges object with ranges for introns
#' @param gtf.exons GRanges object made from a GTF with ONLY exon annotations (no gene, transcript, CDS etc.)
#' @param match what type of matching to perform? perfect = only exons which bound the intron exactly,
#' introns = any exon pairs which overlap the intron,
#' all = any exon pairs AND single exons which overlap the intron
#' @return data.frame with all flanking exon pairs
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
findIntronContainingTranscripts <- function(intronRanges, gtf.exons, match="perfect"){

    #intronRanges = ranges.ri
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

    mcols(gtf_to_a)$from_exon_number <- as.numeric(gtf_from_a$exon_number)
    mcols(gtf_to_a)$to_exon_number <- as.numeric(gtf_to_a$exon_number)

    mcols(gtf_to_a)$intron_exon_number <- apply(GenomicRanges::mcols(gtf_to_a)[,c('to_exon_number','from_exon_number')], 1, mean)
    mcols(gtf_to_a)$match <- "perfect"

    if(match == "introns" | match=="all"){
        # non-exact overlaps
        ol <- findOverlaps(intronRanges, gtf.exons)
        gtf_over <- gtf.exons[ol@to]
        gtf_over$from <- ol@from
        gtf_over$new_id <- with(as.data.frame(GenomicRanges::mcols(gtf_over)), paste0(transcript_id, "_",from))
        gtf_over$exon_number <- as.numeric(gtf_over$exon_number)
        # remove perfect match pairs
        gtf_over <- gtf_over[which(!(gtf_over$new_id %in% gtf_to_a$new_id))]

        transcript_table <- as.data.frame(table(gtf_over$new_id))
        gtf_over_pairs <- gtf_over[gtf_over$new_id %in% transcript_table$Var1[transcript_table$Freq == 2]]
        start_pair <- gtf_over_pairs$new_id[start(gtf_over_pairs) < start(intronRanges[gtf_over_pairs$from])]
        end_pair <- gtf_over_pairs$new_id[end(gtf_over_pairs) > end(intronRanges[gtf_over_pairs$from])]
        keep <- start_pair[start_pair %in% end_pair]
        gtf_over_pairs <- gtf_over_pairs[gtf_over_pairs$new_id %in% keep]

        if(length(gtf_over_pairs) > 0){
            exon_numbers <- aggregate(exon_number ~ new_id, mcols(gtf_over_pairs), mean)
            gtf_over_pairs$intron_exon_number <- exon_numbers$exon_number[match(gtf_over_pairs$new_id, exon_numbers$new_id)]
            gtf_over_pairs <- gtf_over_pairs[!duplicated(gtf_over_pairs$new_id)]
            gtf_over_pairs$from_exon_number <- gtf_over_pairs$intron_exon_number - 0.5
            gtf_over_pairs$to_exon_number <- gtf_over_pairs$intron_exon_number + 0.5
            gtf_over_pairs$match <- "intron"
            mcols(gtf_over_pairs) <-  mcols(gtf_over_pairs[,match(colnames(mcols(gtf_to_a)),colnames(mcols(gtf_over_pairs)))])

            gtf_to_a <- c(gtf_to_a, gtf_over_pairs)
        }
    }
    if(match=="all"){
        # overlaps an exon
        gtf_over <- gtf_over[which(!(gtf_over$new_id %in% gtf_over_pairs$new_id))]
        keep <- which(start(gtf_over) <= start(intronRanges[gtf_over$from]) &
                          end(gtf_over) >= end(intronRanges[gtf_over$from]))
        gtf_over <- gtf_over[keep]
        if(length(gtf_over) > 0){
            gtf_over$intron_exon_number <- gtf_over$exon_number
            gtf_over$from_exon_number <- gtf_over$exon_number
            gtf_over$to_exon_number <- gtf_over$exon_number
            gtf_over$match <- "exon"
            mcols(gtf_over) <-  mcols(gtf_over[,match(colnames(mcols(gtf_to_a)),colnames(mcols(gtf_over)))])
            gtf_to_a <- c(gtf_to_a, gtf_over)
        }
    }

    flanking_exons <- as.data.frame(GenomicRanges::mcols(gtf_to_a)[,c('gene_id','transcript_id','transcript_type',
                                                       'from','from_exon_number',
                                                       'intron_exon_number','to_exon_number','match')])

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
addIntronInTranscript <- function(intronRanges, flanking_exons, gtf.exons, glueExons=TRUE){

    #intronRanges <- ranges.ri
    #flanking_exons <- exons.ri

    intronRanges <- intronRanges[match(flanking_exons$from, intronRanges$id)]
    intronRanges$exon_number <- flanking_exons$intron_exon_number
    intronRanges$transcript_id <- flanking_exons$transcript_id
    intronRanges$transcript_type <- flanking_exons$transcript_type
    intronRanges$gene_id <- flanking_exons$gene_id
    intronRanges$exon_id <- intronRanges$id

    transcripts <- as.data.frame(table(flanking_exons$transcript_id))
    gtf_transcripts <- gtf.exons[gtf.exons$transcript_id %in% transcripts$Var1]
    m <- match(gtf_transcripts$transcript_id, intronRanges$transcript_id)
    mcols(gtf_transcripts) <- cbind(mcols(gtf_transcripts), DataFrame(new_transcript_id=paste0(gtf_transcripts$transcript_id,"+AS ",intronRanges$exon_id[m])))
    mcols(intronRanges) <- cbind(mcols(intronRanges), DataFrame(new_transcript_id = paste0(intronRanges$transcript_id,"+AS ",intronRanges$exon_id)))
    flanking_exons$new_transcript_id <-  paste0(flanking_exons$transcript_id,"+AS ",flanking_exons$from)

    #mcols(gtf_transcripts)$new_transcript_id <- paste0(gtf_transcripts$transcript_id,"+AS ",intronRanges$exon_id[m])
    #mcols(intronRanges)$new_transcript_id <- paste0(intronRanges$transcript_id,"+AS ",intronRanges$exon_id)

    mcols(gtf_transcripts) <- mcols(gtf_transcripts)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number', 'new_transcript_id')]
    mcols(intronRanges) <- mcols(intronRanges)[,c('gene_id','transcript_id','transcript_type','exon_id','exon_number','new_transcript_id')]

    needs_dup <- which(!(intronRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))

    if(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts[gtf_transcripts$transcript_id %in% intronRanges$transcript_id[needs_dup]]
    }

    while(length(needs_dup) > 0){
        gtf_transcripts_add <- gtf_transcripts_add[gtf_transcripts_add$transcript_id %in% intronRanges$transcript_id[needs_dup]]
        m <- match(gtf_transcripts_add$transcript_id, intronRanges$transcript_id[needs_dup])
        gtf_transcripts_add$new_transcript_id <- paste0(gtf_transcripts_add$transcript_id,"+AS ",intronRanges$exon_id[needs_dup][m])
        gtf_transcripts <- c(gtf_transcripts, gtf_transcripts_add)
        needs_dup <- which(!(intronRanges$new_transcript_id %in% gtf_transcripts$new_transcript_id))
    }

    # fix starts/ends of introns
    gtf_transcripts$new_id_ex <- with(mcols(gtf_transcripts), paste0(new_transcript_id, "_", exon_number))
    flanking_exons$from_ex <- with(flanking_exons, paste0(new_transcript_id, "_", from_exon_number))
    flanking_exons$to_ex <- with(flanking_exons, paste0(new_transcript_id, "_", to_exon_number))

    wi <- which(flanking_exons$match == "intron")
    wi_p <- wi[as.character(strand(gtf_transcripts[match(flanking_exons$from_ex[wi], gtf_transcripts$new_id_ex)])) == "+"]
    wi_n <- wi[as.character(strand(gtf_transcripts[match(flanking_exons$from_ex[wi], gtf_transcripts$new_id_ex)])) == "-"]

    if(length(wi_p) > 0){
        end(gtf_transcripts)[match(flanking_exons$from_ex[wi_p], gtf_transcripts$new_id_ex)] <-
            start(intronRanges)[match(flanking_exons$new_transcript_id[wi_p], intronRanges$new_transcript_id)]
        start(gtf_transcripts)[match(flanking_exons$to_ex[wi_p], gtf_transcripts$new_id_ex)] <-
            end(intronRanges)[match(flanking_exons$new_transcript_id[wi_p], intronRanges$new_transcript_id)]
    }
    if(length(wi_n) > 0){
        end(gtf_transcripts)[match(flanking_exons$to_ex[wi_n], gtf_transcripts$new_id_ex)] <-
            start(intronRanges)[match(flanking_exons$new_transcript_id[wi_n], intronRanges$new_transcript_id)]
        start(gtf_transcripts)[match(flanking_exons$from_ex[wi_n], gtf_transcripts$new_id_ex)] <-
            end(intronRanges)[match(flanking_exons$new_transcript_id[wi_n], intronRanges$new_transcript_id)]
    }
    # create an intron in 'exons'

    we <- which(flanking_exons$match == "exon")
    replace <- match(flanking_exons$from_ex[we],gtf_transcripts$new_id_ex)
    replacement_start <- gtf_transcripts[replace]
    end(replacement_start) <- start(intronRanges)[match(flanking_exons$new_transcript_id[we], intronRanges$new_transcript_id)]
    replacement_end <- gtf_transcripts[replace]
    start(replacement_end) <- end(intronRanges)[match(flanking_exons$new_transcript_id[we], intronRanges$new_transcript_id)]
    gtf_transcripts <- c(gtf_transcripts[-replace], replacement_start, replacement_end)

    gtf_transcripts <- reorderExonNumbers(gtf_transcripts, by="new_transcript_id")
    gtf_transcripts$new_id_ex <- NULL

    gtf_transcripts_withIntron <- c(gtf_transcripts, intronRanges)
    gtf_transcripts_withIntron <- reorderExonNumbers(gtf_transcripts_withIntron, by="new_transcript_id")

    mcols(gtf_transcripts_withIntron) <- mcols(gtf_transcripts_withIntron)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts_withIntron))[2] <- "transcript_id"

    mcols(gtf_transcripts) <- mcols(gtf_transcripts)[,c('gene_id','new_transcript_id','transcript_type','exon_id','exon_number')]
    colnames(mcols(gtf_transcripts))[2] <- "transcript_id"

    gtf_transcripts_withIntron$exon_number <- as.numeric(gtf_transcripts_withIntron$exon_number)
    order <- order(gtf_transcripts_withIntron$transcript_id, gtf_transcripts_withIntron$exon_number)
    gtf_transcripts_withIntron <- gtf_transcripts_withIntron[order]

    #join together exons that are not seperated by an intron
    if(glueExons==TRUE){

        # split pos/neg ordering
        gtf_transcripts_withIntron.neg <-
            gtf_transcripts_withIntron[which(as.logical(strand(gtf_transcripts_withIntron) == "-"))]
        order <- order(gtf_transcripts_withIntron.neg$transcript_id,
                       plyr::desc(gtf_transcripts_withIntron.neg$exon_number))
        gtf_transcripts_withIntron.neg <- gtf_transcripts_withIntron.neg[order]

        gtf_transcripts_withIntron.pos <-
            gtf_transcripts_withIntron[which(as.logical(strand(gtf_transcripts_withIntron) == "+"))]
        gtf_transcripts_withIntron <- c(gtf_transcripts_withIntron.pos,gtf_transcripts_withIntron.neg)

        #extend starts <---<---<---
        w <- which(end(ranges(gtf_transcripts_withIntron))[-length(gtf_transcripts_withIntron)] == start(ranges(gtf_transcripts_withIntron[-1])))
        gtf_transcripts_withIntron <- gtf_transcripts_withIntron
        GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_withIntron))[w+1] <-
            GenomicRanges::start(GenomicRanges::ranges(gtf_transcripts_withIntron))[w]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_withIntron, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to - 1]
        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_withIntron <- gtf_transcripts_withIntron[-rm]
        }

        #extend ends --->--->--->
        overlaps <- findOverlaps(gtf_transcripts_withIntron)
        overlaps <- overlaps[overlaps@from==overlaps@to +1]

        GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_withIntron))[overlaps@to] <-
            GenomicRanges::end(GenomicRanges::ranges(gtf_transcripts_withIntron))[overlaps@from]

        # remove exons that cover now redundant regions
        overlaps <- findOverlaps(gtf_transcripts_withIntron, type="within")
        overlaps <- overlaps[overlaps@from == overlaps@to + 1]

        rm <- unique(overlaps@from)
        if(length(rm) >0){
            gtf_transcripts_withIntron <- gtf_transcripts_withIntron[-rm]
        }

        #order <- order(gtf_transcripts_withIntron$transcript_id, gtf_transcripts_withIntron$exon_number)
        #gtf_transcripts_withIntron <- gtf_transcripts_withIntron[order]

    }
    gtf_transcripts$set <- "spliced_intron"
    gtf_transcripts_withIntron$set <- "retained_intron"

    return(c(gtf_transcripts_withIntron, gtf_transcripts))
}

