#' Find the largest distance between two vectors of numbers
#' Helper function for get_orfs
#' @param start_site vector of start sites - i.e Met amino acid positions
#' @param stop_site vector of stop sites - i.e Stop (*) amino acid positions
#' @param longest how many pairs to return
#' @return sequential start site and end site with the greatest difference
#' @export
#' @import plyr
#' @examples
#' @author Beth Signal
maxLocation <- function(start_site, stop_site, longest = 1){
    if(length(start_site) > 0){
        #diffs <- unlist(lapply(start_site, function(x) stop_site[which(stop_site > x)[1]] - x))
        #order <- order(diffs, decreasing=TRUE)
        #max_loc <- order[longest]
        #start <- start_site[max_loc]
        #stop <- stop_site[which(stop_site > start)[1]]
        #return(c(start,stop))

        # make start / stop pairs
        stop_pair_index <- unlist(lapply(start_site, function(x) which.max(1/(stop_site - x))))
        pairs <- data.frame(start=start_site, stop=stop_site[stop_pair_index])
        pairs$len <- pairs$stop - pairs$start
        pairs <- plyr::arrange(pairs, plyr::desc(len))
        pairs <- pairs[!duplicated(pairs$stop),]
        return(as.numeric(pairs[longest,1:2]))    }else{
            return(c(NA,NA))
        }
}

#' Get open reading frames for transcripts
#' @param transcripts GRanges object with ONLY exon annotations (no gene, transcript, CDS etc.) with all transcripts for orf retrevial
#' @param BSgenome BSgenome object
#' @param returnLongestOnly only return longest ORF?
#' @param all_frames return longest ORF for all 3 frames?
#' @param longest return x longest ORFs (regardless of frames)
#' @return data.frame with longest orf details
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @import stringr
#' @examples
#' @author Beth Signal
getOrfs <- function(transcripts, BSgenome = g, returnLongestOnly=TRUE, all_frames=FALSE, longest=1){

    if (all_frames == TRUE) {
        returnLongestOnly = FALSE
        longest = 1
    }


    # check -ve dist to junction b calls
    # check exon number ORF start/ends

    transcripts$exon_number <-
        as.numeric(transcripts$exon_number)
    order <-
        order(transcripts$transcript_id, transcripts$exon_number)
    transcripts <- transcripts[order]
    transcripts$seq <-
        as.character(Biostrings::getSeq(g, transcripts))

    seq_cat <-
        aggregate(seq ~ transcript_id, mcols(transcripts), function(x)
            (paste(x, collapse = "")))
    ids <- as.character(seq_cat$transcript_id)
    seq_cat <- seq_cat$seq
    rm <- which(grepl("N", seq_cat))

    if (length(rm) > 0) {
        seq_cat <- seq_cat[-rm]
        remove_id <- ids[rm]
        ids <- ids[-rm]
        transcripts <-
            transcripts[-which(transcripts$transcript_id %in% remove_id)]
    }

    # 3 frames
    seq_cat <-
        c(seq_cat, str_sub(seq_cat, 2), str_sub(seq_cat, 3))
    frames <- rep(c(1, 2, 3), each = length(ids))
    ids <- c(ids, ids, ids)

    orf <-
        suppressWarnings(unlist(lapply(seq_cat, function(x)
            as.character(Biostrings::translate(Biostrings::DNAString(x))))))

    orf_df <- data.frame(
        id = ids,
        aa_sequence = orf,
        frame = frames,
        stringsAsFactors = FALSE
    )

    orf_df$seq_length <- nchar(orf_df$aa_sequence)
    orf_df$seq_length_nt <- nchar(seq_cat) + orf_df$frame -1

    start_sites <-
        stringr::str_locate_all(orf_df$aa_sequence, "M")
    start_sites <-
        lapply(start_sites, function(x)
            as.numeric(x[, 2]))
    stop_sites <- str_locate_all(orf_df$aa_sequence, "[*]")
    stop_sites <-
        mapply(function(x, y)
            c(as.numeric(x[, 2]), nchar(y)),
            stop_sites,
            orf_df$aa_sequence)

    max_loc <-
        mapply(function(x, y)
            maxLocation(x, y), start_sites, stop_sites)

    if (longest >= 2 & returnLongestOnly == FALSE) {
        orf_df_longest <- orf_df

        for (i in 2:longest) {
            max_loc <- cbind(max_loc,
                             mapply(
                                 function(x, y)
                                     maxLocation(x, y, longest = i),
                                 start_sites,
                                 stop_sites
                             ))
            orf_df_longest <- rbind(orf_df_longest, orf_df)
        }

        o <- order(max_loc[2, ] - max_loc[1, ], decreasing = TRUE)

        orf_df_longest$start_site <- max_loc[1, ]
        orf_df_longest$stop_site <- max_loc[2, ]

        orf_df_longest <- orf_df_longest[o, ]
        keep <- which(!duplicated(orf_df_longest$id))
        orf_df <- orf_df_longest[keep, ]
        orf_df_longest <- orf_df_longest[-keep, ]

        for (i in 2:longest) {
            keep <- which(!duplicated(orf_df_longest$id))
            orf_df <- rbind(orf_df, orf_df_longest[keep, ])
            orf_df_longest <- orf_df_longest[-keep, ]
        }

    } else{
        orf_df$start_site <- max_loc[1, ]
        orf_df$stop_site <- max_loc[2, ]
    }

    orf_df$orf_sequence <-
        stringr::str_sub(orf_df$aa_sequence, orf_df$start_site, orf_df$stop_site - 1)
    orf_df$orf_length <- nchar(orf_df$orf_sequence)

    orf_df$start_site_nt <-
        (orf_df$start_site * 3)- 3 + orf_df$frame
    orf_df$stop_site_nt <- (orf_df$orf_length * 3) + orf_df$start_site_nt + 3
    widths <- data.frame(w = width(transcripts),
                         id = transcripts$transcript_id)
    pad <- max(table(widths$id))
    if (pad > 1) {
        if (length(unique(transcripts$transcript_id)) == 1) {
            w <- cumsumANDpad(widths$w, pad)
            diffs <-
                lapply(orf_df$stop_site_nt, function(x)
                    x - w)
            diffs <-
                matrix(unlist(diffs), ncol = length(orf_df$id))
        } else{
            #widths_w <- aggregate(w ~ id, widths, function(x) cumsum(c(1,x))[-1])
            widths_w2 <-
                aggregate(w ~ id, widths, function(x)
                    cumsumANDpad(x, pad_length = pad))
            m <- match(orf_df$id, widths_w2$id)
            widths_w2 <- widths_w2[m, -1]
            widths_w2  <- split(widths_w2, seq(nrow(widths_w2)))
            diffs <-
                mapply(function(x , y)
                    x - y, orf_df$stop_site_nt, widths_w2)
        }

        orf_df$min_dist_to_junction_a <-
            suppressWarnings(apply(diffs, 2, function(x)
                min(x[x > 0 & !is.na(x)])))
        orf_df$min_dist_to_junction_a[which(is.infinite(orf_df$min_dist_to_junction_a))] <-
            NA
        orf_df$exon_a_from_start <-
            (apply(diffs, 2, function(x)
                length(x[x > 0 & !is.na(x)]))) - 1

        orf_df$min_dist_to_junction_b <-
            suppressWarnings((apply(diffs, 2, function(x)
                max(x[x <= 0 & !is.na(x)])) * -1) + 1)
        orf_df$min_dist_to_junction_b[which(is.infinite(orf_df$min_dist_to_junction_b))] <-
            NA
        orf_df$exon_b_from_final <-
            (apply(diffs, 2, function(x)
                length(x[x <= 0 & !is.na(x)]))) - 1

        exon_num <-
            apply(diffs, 2, function(x)
                length(which(!is.na(x))))
        orf_df$exon_a_from_start[exon_num == 1] <- NA
        orf_df$exon_b_from_final[exon_num == 1] <- NA
    } else{
        # all single exon transcripts -- therefore no junctions
        orf_df$min_dist_to_junction_a <- NA
        orf_df$exon_a_from_start <- NA
        orf_df$min_dist_to_junction_b <- NA
        orf_df$exon_b_from_final <- NA
    }
    orf_df$utr3_length <-
        (orf_df$seq_length_nt - orf_df$stop_site_nt) + 1

    orf_df$aa_sequence <- NULL

    if (returnLongestOnly == TRUE) {
        orf_df <- plyr::arrange(orf_df, plyr::desc(orf_length))
        orf_df <- orf_df[!duplicated(orf_df$id), ]
    }

    orf_df <- plyr::arrange(orf_df, id)
    m <- match(orf_df$id, transcripts$transcript_id)
    orf_df$gene_id <- transcripts$gene_id[m]
    orf_df <- orf_df[,c(1, ncol(orf_df), 2:(ncol(orf_df)-1))]

    return(orf_df)
}

#' Cumulative sum of a sequence of numbers, padded with NA
#' @param x input numeric vector
#' @param pad_length length to pad output to
#' @return vector with cumulative sum, padded with NA
#' @export
#' @examples
#' @author Beth Signal
cumsumANDpad <- function(x, pad_length){
    y <- cumsum(c(1,x))[-1]
    if(length(y) < pad_length){
        y <- c(y, rep(NA, pad_length - length(y)))
    }
    return(y)
}

