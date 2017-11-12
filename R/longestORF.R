#' Find the largest distance between two vectors of numbers
#' Helper function for get_orfs
#' @param startSite vector of start sites - i.e Met amino acid positions
#' @param stopSite vector of stop sites - i.e Stop (*) amino acid positions
#' @param longest which pair to return (1 = longest pair, 2= 2nd longest pair etc.)
#' @return sequential start site and end site with the greatest difference
#' @export
#' @import plyr
#' @examples
#' @author Beth Signal
#' @examples
#' starts <- c(1,10,15,25)
#' stops <- c(4,16,50,55)
#' # longest start site = 25, longest stop site = 50
#' maxLocation(starts, stops, longest = 1)
#' starts <- c(1,10,15,25)
#' stops <- c(4,14,50,55)
#' # longest start site = 15, longest stop site = 50
#' maxLocation(starts, stops, longest = 1)
#' # 2nd longest start site = 10, 2nd longest stop site = 14
#' maxLocation(starts, stops, longest = 2)
maxLocation <- function(startSite, stopSite, longest = 1){
    if(length(startSite) > 0){
        #diffs <- unlist(lapply(startSite, function(x) stopSite[which(stopSite > x)[1]] - x))
        #order <- order(diffs, decreasing=TRUE)
        #max_loc <- order[longest]
        #start <- startSite[max_loc]
        #stop <- stopSite[which(stopSite > start)[1]]
        #return(c(start,stop))

        # make start / stop pairs
        stop_pair_index <- unlist(lapply(startSite, function(x) which.max(1/(stopSite - x))))
        pairs <- data.frame(start=startSite, stop=stopSite[stop_pair_index])
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
#' @param allFrames return longest ORF for all 3 frames?
#' @param longest return x longest ORFs (regardless of frames)
#' @return data.frame with longest orf details
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @import stringr
#' @examples
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "gencode.vM14.neurl1a.gtf", package="GeneStructureTools"))
#' transcript <- gtf[gtf$type=="exon"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' # longest ORF for each transcripts
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = TRUE)
#' # longest ORF in all 3 frames for each transcript
#' orfs <- getOrfs(transcript, BSgenome = g, allFrames = TRUE)
#' # longest 3 ORFS in eacht transcript
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = FALSE, longest=3)
getOrfs <- function(transcripts, BSgenome = g, returnLongestOnly=TRUE, allFrames=FALSE, longest=1){

    if (allFrames == TRUE) {
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
    # add first site as potential start (if no M)
    start_sites <-
        lapply(start_sites, function(x)
            if(length(x) == 0){1}else{
                as.numeric(x[, 2])})

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
    orf_df$utr3_length <-
        (orf_df$seq_length_nt - orf_df$stop_site_nt) + 1

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
                    cumsumANDpad(x, pad))
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
            orf_df$start_site_nt[which(is.infinite(orf_df$min_dist_to_junction_a))]
        orf_df$exon_a_from_start <-
            (apply(diffs, 2, function(x)
                length(x[x > 0 & !is.na(x)]))) - 1

        orf_df$min_dist_to_junction_b <-
            suppressWarnings((apply(diffs, 2, function(x)
                max(x[x <= 0 & !is.na(x)])) * -1) + 1)
        orf_df$min_dist_to_junction_b[which(is.infinite(orf_df$min_dist_to_junction_b))] <-
            orf_df$utr3_length[which(is.infinite(orf_df$min_dist_to_junction_b))]
        orf_df$exon_b_from_final <-
            (apply(diffs, 2, function(x)
                length(x[x <= 0 & !is.na(x)]))) - 1

        exon_num <-
            apply(diffs, 2, function(x)
                length(which(!is.na(x))))
        orf_df$exon_a_from_start[exon_num == 1] <- 0
        orf_df$exon_b_from_final[exon_num == 1] <- 0
    } else{
        # all single exon transcripts -- therefore no junctions
        orf_df$min_dist_to_junction_a <- orf_df$start_site_nt
        orf_df$exon_a_from_start <- 0
        orf_df$min_dist_to_junction_b <- orf_df$utr3_length
        orf_df$exon_b_from_final <- 0
    }

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
#' @param padLength length to pad output to
#' @return vector with cumulative sum, padded with NA
#' @export
#' @examples
#' @author Beth Signal
#' @examples
#' x <- c(1,4,7,2,5)
#' cumsumANDpad(x, 10)
cumsumANDpad <- function(x, padLength){
    y <- cumsum(c(1,x))[-1]
    if(length(y) < padLength){
        y <- c(y, rep(NA, padLength - length(y)))
    }
    return(y)
}

