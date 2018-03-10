#' Find the largest distance between two vectors of numbers
#' Helper function for get_orfs
#' @param startSite vector of start sites - i.e Met amino acid positions
#' @param stopSite vector of stop sites - i.e Stop (*) amino acid positions
#' @param longest which pair to return (1 = longest pair, 2= 2nd longest pair etc.)
#' @return sequential start site and end site with the greatest difference
#' @export
#' @import plyr
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
        # make start / stop pairs
        stopPairIndex <-
            unlist(lapply(startSite, function(x) which.max(1/(stopSite - x))))
        pairs <- data.frame(start=startSite, stop=stopSite[stopPairIndex])
        pairs$len <- pairs$stop - pairs$start
        pairs <- pairs[order(plyr::desc(pairs$len)),]
        #pairs <- plyr::arrange(pairs, plyr::desc(len))
        pairs <- pairs[!duplicated(pairs$stop),]
        return(as.numeric(pairs[longest,1:2]))    }else{
            return(c(NA,NA))
        }
}

#' Get open reading frames for transcripts
#' @param transcripts GRanges object with ONLY exon annotations
#' (no gene, transcript, CDS etc.) with all transcripts for orf retrevial
#' @param BSgenome BSgenome object
#' @param returnLongestOnly only return longest ORF?
#' @param allFrames return longest ORF for all 3 frames?
#' @param longest return x longest ORFs (regardless of frames)
#' @param exportFasta export a .fa.gz file with nucleotide sequences for each transcript?
#' @param fastaFile file name for .fa.gz export
#' @param uORFs get uORF summaries?
#' @return data.frame with longest orf details
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @import stringr
#' @importFrom plyr arrange
#' @importFrom plyr desc
#' @importFrom stats aggregate
#' @importFrom rtracklayer import
#' @importFrom utils write.table
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "example_gtf.gtf",
#' package="GeneStructureTools"))
#' transcript <- gtf[gtf$type=="exon" & gtf$gene_name=="Neurl1a"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' # longest ORF for each transcripts
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = TRUE)
#' # longest ORF in all 3 frames for each transcript
#' orfs <- getOrfs(transcript, BSgenome = g, allFrames = TRUE)
#' # longest 3 ORFS in eacht transcript
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = FALSE, longest=3)
getOrfs <- function(transcripts,
                    BSgenome = NULL,
                    returnLongestOnly=TRUE,
                    allFrames=FALSE,
                    longest=1,
                    exportFasta=FALSE,
                    fastaFile=NULL,
                    uORFs=FALSE){

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
        as.character(Biostrings::getSeq(BSgenome, transcripts))

    seqCat <-
        aggregate(seq ~ transcript_id, mcols(transcripts), function(x)
            (paste(x, collapse = "")))
    ids <- as.character(seqCat$transcript_id)

    if(exportFasta & !is.null(fastaFile)){
        fastaFile <- ifelse(stringr::str_sub(fastaFile, -3,-1) == ".gz",
                            fastaFile, paste0(fastaFile,".gz"))
        fa <- seqCat
        fa$transcript_id <- paste0("> ", fa$transcript_id)

        # gzip files to save some space
        gz <- gzfile(fastaFile, "w")
        utils::write.table(fa, gz, col.names = FALSE,
                           row.names = FALSE, quote=FALSE, sep="\n")
        close(gz)
    }else if(exportFasta & is.null(fastaFile)){
        message("skipping writing .fa file")
        message("please specify a file name for export")
    }

    seqCat <- seqCat$seq
    rm <- which(grepl("N", seqCat))

    if (length(rm) > 0) {
        seqCat <- seqCat[-rm]
        removeId <- ids[rm]
        ids <- ids[-rm]
        transcripts <-
            transcripts[-which(transcripts$transcript_id %in% removeId)]
    }

    # 3 frames
    seqCat <-
        c(seqCat, stringr::str_sub(seqCat, 2), stringr::str_sub(seqCat, 3))
    frames <- rep(c(1, 2, 3), each = length(ids))
    ids <- c(ids, ids, ids)

    orf <-
        suppressWarnings(unlist(lapply(seqCat, function(x)
            as.character(Biostrings::translate(Biostrings::DNAString(x))))))

    orfDF <- data.frame(
        id = ids,
        aa_sequence = orf,
        frame = frames,
        stringsAsFactors = FALSE
    )

    orfDF$seq_length <- nchar(orfDF$aa_sequence)
    orfDF$seq_length_nt <- nchar(seqCat) + orfDF$frame -1

    startSites <-
        stringr::str_locate_all(orfDF$aa_sequence, "M")
    # add first site as potential start (if no M)
    startSites <-
        lapply(startSites, function(x)
            if(length(x) == 0){1}else{
                as.numeric(x[, 2])})

    stopSites <- str_locate_all(orfDF$aa_sequence, "[*]")
    stopSites <-
        mapply(function(x, y)
            c(as.numeric(x[, 2]), nchar(y)),
            stopSites,
            orfDF$aa_sequence)

    maxLoc <-
        mapply(function(x, y)
            maxLocation(x, y), startSites, stopSites)

    if (longest >= 2 & returnLongestOnly == FALSE) {
        orfDF.longest <- orfDF

        for (i in 2:longest) {
            maxLoc <- cbind(maxLoc,
                            mapply(
                                function(x, y)
                                    maxLocation(x, y, longest = i),
                                startSites,
                                stopSites
                            ))
            orfDF.longest <- rbind(orfDF.longest, orfDF)
        }

        o <- order(maxLoc[2, ] - maxLoc[1, ], decreasing = TRUE)

        orfDF.longest$start_site <- maxLoc[1, ]
        orfDF.longest$stop_site <- maxLoc[2, ]

        orfDF.longest <- orfDF.longest[o, ]
        keep <- which(!duplicated(orfDF.longest$id))
        orfDF <- orfDF.longest[keep, ]
        orfDF.longest <- orfDF.longest[-keep, ]

        for (i in 2:longest) {
            keep <- which(!duplicated(orfDF.longest$id))
            orfDF <- rbind(orfDF, orfDF.longest[keep, ])
            orfDF.longest <- orfDF.longest[-keep, ]
        }

    } else{
        orfDF$start_site <- maxLoc[1, ]
        orfDF$stop_site <- maxLoc[2, ]
    }

    orfDF$orf_sequence <-
        stringr::str_sub(orfDF$aa_sequence, orfDF$start_site,
                         orfDF$stop_site - 1)
    orfDF$orf_length <- nchar(orfDF$orf_sequence)

    if(returnLongestOnly==TRUE){
        orfDF <- plyr::arrange(orfDF, plyr::desc(orf_length))
        orfDF <- orfDF[!duplicated(orfDF$id),]
    }

    orfDF$start_site_nt <-
        (orfDF$start_site * 3)- 3 + orfDF$frame
    orfDF$stop_site_nt <- (orfDF$orf_length * 3) + orfDF$start_site_nt + 3
    orfDF$utr3_length <-
        (orfDF$seq_length_nt - orfDF$stop_site_nt) + 1

    widths <- data.frame(w = width(transcripts),
                         id = transcripts$transcript_id)
    pad <- max(table(widths$id))
    if (pad > 1) {
        if (length(unique(transcripts$transcript_id)) == 1) {
            w <- cumsumANDpad(widths$w, pad)
            diffs <-
                lapply(orfDF$stop_site_nt, function(x)
                    x - w)
            diffs <-
                matrix(unlist(diffs), ncol = length(orfDF$id))
        } else{
            w2 <-
                aggregate(w ~ id, widths, function(x)
                    cumsumANDpad(x, pad))
            m <- match(orfDF$id, w2$id)
            w2 <- w2[m, -1]
            w2  <- split(w2, seq(nrow(w2)))
            diffs <-
                mapply(function(x , y)
                    x - y, orfDF$stop_site_nt, w2)
        }

        orfDF$min_dist_to_junction_a <-
            suppressWarnings(apply(diffs, 2, function(x)
                min(x[x > 0 & !is.na(x)])))
        orfDF$min_dist_to_junction_a[
            which(is.infinite(orfDF$min_dist_to_junction_a))] <-
            orfDF$start_site_nt[
                which(is.infinite(orfDF$min_dist_to_junction_a))]
        orfDF$exon_a_from_start <-
            (apply(diffs, 2, function(x)
                length(x[x > 0 & !is.na(x)])))

        orfDF$min_dist_to_junction_b <-
            suppressWarnings((apply(diffs, 2, function(x)
                max(x[x <= 0 & !is.na(x)])) * -1) + 1)
        orfDF$min_dist_to_junction_b[
            which(is.infinite(orfDF$min_dist_to_junction_b))] <-
            orfDF$utr3_length[which(is.infinite(orfDF$min_dist_to_junction_b))]
        orfDF$exon_b_from_final <-
            (apply(diffs, 2, function(x)
                length(x[x <= 0 & !is.na(x)]))) - 1

        exonNumber <-
            apply(diffs, 2, function(x)
                length(which(!is.na(x))))
        orfDF$exon_a_from_start[exonNumber == 1] <- 0
        orfDF$exon_b_from_final[exonNumber == 1] <- 0
    } else{
        # all single exon transcripts -- therefore no junctions
        orfDF$min_dist_to_junction_a <- orfDF$start_site_nt
        orfDF$exon_a_from_start <- 0
        orfDF$min_dist_to_junction_b <- orfDF$utr3_length
        orfDF$exon_b_from_final <- 0
    }

    orfDF$aa_sequence <- NULL

    orfDF <- plyr::arrange(orfDF, id)
    m <- match(orfDF$id, transcripts$transcript_id)
    orfDF$gene_id <- transcripts$gene_id[m]
    orfDF <- orfDF[,c(1, ncol(orfDF), 2:(ncol(orfDF)-1))]

    if(uORFs == TRUE){
        upstreamORFs <- getUOrfs(transcripts=transcripts,
                                 BSgenome = BSgenome,
                                 orfs = orfDF,
                                 findExonB = TRUE)

        if(nrow(upstreamORFs) > 0){
            uORFS.bytranscript <- aggregate(
                overlaps_main_ORF ~ id+frame, upstreamORFs,
                function(x) length(x))
            colnames(uORFS.bytranscript)[3] <- "total_uorfs"

            uORFS.bytranscript.newVal <- aggregate(
                overlaps_main_ORF ~ id+frame,
                upstreamORFs, function(x)
                    length(x[which(x=="upstream")]))
            uORFS.bytranscript$upstream_count <-
                uORFS.bytranscript.newVal[match(paste0(uORFS.bytranscript$id,
                                                       uORFS.bytranscript$frame),
                                                paste0(uORFS.bytranscript.newVal$id,
                                                       uORFS.bytranscript.newVal$frame)),3]

            uORFS.bytranscript.newVal <-
                aggregate(overlaps_main_ORF ~ id+frame,
                          upstreamORFs, function(x)
                              length(x[which(x=="downstream")]))
            uORFS.bytranscript$downstream_count <-
                uORFS.bytranscript.newVal[match(paste0(uORFS.bytranscript$id,
                                                       uORFS.bytranscript$frame),
                                                paste0(uORFS.bytranscript.newVal$id,
                                                       uORFS.bytranscript.newVal$frame)),3]

            uORFS.bytranscript.newVal <-
                aggregate(uorf_length ~ id+frame,upstreamORFs, function(x) max(x))
            uORFS.bytranscript$max_uorf <-
                uORFS.bytranscript.newVal[match(paste0(uORFS.bytranscript$id,
                                                       uORFS.bytranscript$frame),
                                                paste0(uORFS.bytranscript.newVal$id,
                                                       uORFS.bytranscript.newVal$frame)),3]

            if(any(upstreamORFs$exon_b_from_final != 0)){
                uORFS.bytranscript.newVal <-
                    aggregate(min_dist_to_junction_b ~ id+frame,
                              upstreamORFs[upstreamORFs$exon_b_from_final !=0,],
                              function(x) max(x))
                uORFS.bytranscript$uorf_maxb <-
                    uORFS.bytranscript.newVal[match(paste0(uORFS.bytranscript$id,
                                                           uORFS.bytranscript$frame),
                                                    paste0(uORFS.bytranscript.newVal$id,
                                                           uORFS.bytranscript.newVal$frame)),3]
            }else{
                uORFS.bytranscript$uorf_maxb <- NA
            }
            m <- match(paste0(orfDF$id, orfDF$frame), paste0(uORFS.bytranscript$id,
                                                             uORFS.bytranscript$frame))
            orfDF <- cbind(orfDF, uORFS.bytranscript[m,-c(1:2)])
        }else{
            orfDF$total_uorfs <- NA
            orfDF$upstream_count <- NA
            orfDF$downstream_count <- NA
            orfDF$max_uorf <- NA
            orfDF$uorf_maxb <- NA
        }

        # replace any non-matching uorf summaries with 0
        for(i in 17:ncol(orfDF)){
            orfDF[which(is.na(orfDF[,i])),i] <- 0
        }


    }

    rownames(orfDF) <- NULL
    return(orfDF)
}

#' Cumulative sum of a sequence of numbers, padded with NA
#' @param x input numeric vector
#' @param padLength length to pad output to
#' @return vector with cumulative sum, padded with NA
#' @export
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

#' Get upstream open reading frames for transcripts with annotated main ORFs
#' @param transcripts GRanges object with ONLY exon annotations
#' (no gene, transcript, CDS etc.) with all transcripts for orf retrevial
#' @param BSgenome BSgenome object
#' @param orfs orf annotation for the transcripts object. Generated by getOrfs(transcripts, ...)
#' @param findExonB find the distance to and exon number of the downstream (B) junction?
#' @return data.frame with all upstream ORF details.
#' @export
#' @import GenomicRanges
#' @import Biostrings
#' @import stringr
#' @importFrom rtracklayer import
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "example_gtf.gtf",
#' package="GeneStructureTools"))
#' transcript <- gtf[gtf$type=="exon" & gtf$gene_name=="Neurl1a"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' # longest ORF for each transcripts
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = FALSE)
#' uORFS <- getUOrfs(transcript, BSgenome = g, orfs = orfs, findExonB = TRUE)
getUOrfs <- function(transcripts,
                     BSgenome = NULL,
                     orfs, findExonB=FALSE){

    transcripts$exon_number <-
        as.numeric(transcripts$exon_number)
    order <-
        order(transcripts$transcript_id, transcripts$exon_number)
    transcripts <- transcripts[order]
    transcripts$seq <-
        as.character(Biostrings::getSeq(BSgenome, transcripts))

    seqCat <-
        aggregate(seq ~ transcript_id, mcols(transcripts), function(x)
            (paste(x, collapse = "")))
    ids <- as.character(seqCat$transcript_id)
    seqCat <- seqCat$seq
    rm <- which(grepl("N", seqCat))

    if (length(rm) > 0) {
        seqCat <- seqCat[-rm]
        removeId <- ids[rm]
        ids <- ids[-rm]
        transcripts <-
            transcripts[-which(transcripts$transcript_id %in% removeId)]
    }

    # 3 frames
    seqCat <-
        c(seqCat, stringr::str_sub(seqCat, 2), stringr::str_sub(seqCat, 3))
    frames <- rep(c(1, 2, 3), each = length(ids))
    ids <- c(ids, ids, ids)

    orf <-
        suppressWarnings(unlist(lapply(seqCat, function(x)
            as.character(Biostrings::translate(Biostrings::DNAString(x))))))

    orfDF <- data.frame(
        id = ids,
        aa_sequence = orf,
        frame = frames,
        stringsAsFactors = FALSE
    )

    orfDF$seq_length <- nchar(orfDF$aa_sequence)
    orfDF$seq_length_nt <- nchar(seqCat) + orfDF$frame -1

    startSites <-
        stringr::str_locate_all(orfDF$aa_sequence, "M")
    # NOTE THAT M MUST BE FIRST CODON
    startSites <-
        lapply(startSites, function(x)
            as.numeric(x[, 2]))

    # startSites <- mapply(function(x, y)
    #     x[which(x < y)],
    #     startSites,
    #     orfs$start_site[match(orfDF$id, orfs$id)])

    stopSites <- str_locate_all(orfDF$aa_sequence, "[*]")
    stopSites <-
        mapply(function(x, y)
            c(as.numeric(x[, 2]), nchar(y)),
            stopSites,
            orfDF$aa_sequence)

    maxLoc1 <-
        mapply(function(x, y)
            maxLocation(x, y), startSites, stopSites)
    maxLoc <- maxLoc1
    maxLoc.id <- orfDF$id
    maxLoc.frame <- orfDF$frame

    redo <- which(!(apply(maxLoc1,2, function(x) all(is.na(x)))))
    n=2
    while(length(redo) > 0){
        maxLoc1 <- mapply(function(x, y)
            maxLocation(x, y, n), startSites[redo], stopSites[redo])
        maxLoc.id <- c(maxLoc.id, orfDF$id[redo])
        maxLoc.frame <- c(maxLoc.frame, orfDF$frame[redo])
        maxLoc <- cbind(maxLoc, maxLoc1)
        n=n+1
        redo <- which(!(apply(maxLoc1,2, function(x) all(is.na(x)))))
    }

    upstreamORFs <- data.frame(id=maxLoc.id, frame=maxLoc.frame, t(maxLoc))
    colnames(upstreamORFs)[3:4] <- c("start","stop")

    if(any(duplicated(orfs$id))){
        orfs$id <- paste0(orfs$id, "_frame",orfs$frame)

        upstreamORFs <- rbind(upstreamORFs, upstreamORFs, upstreamORFs)
        upstreamORFs$id <-
            paste0(upstreamORFs$id, "_frame",
                   rep(1:3, each=nrow(upstreamORFs)/3))
    }

    m <- match(upstreamORFs$id, orfs$id)
    upstreamORFs <- upstreamORFs[which((upstreamORFs$start - orfs$start_site[m]) < 0),]
    m <- match(upstreamORFs$id, orfs$id)
    upstreamORFs$dist_to_start <- orfs$start_site[m] - upstreamORFs$stop
    upstreamORFs$overlaps_main_ORF <- ifelse(upstreamORFs$dist_to_start > 0,
                                             "upstream", "downstream")
    upstreamORFs <- upstreamORFs[!is.na(upstreamORFs$start),]
    upstreamORFs$uorf_length <- upstreamORFs$stop - upstreamORFs$start

    upstreamORFs$start_site_nt <-
        (upstreamORFs$start * 3)- 3 + upstreamORFs$frame
    upstreamORFs$stop_site_nt <- (upstreamORFs$uorf_length * 3) +
        upstreamORFs$start_site_nt + 3
    m <- match(upstreamORFs$id,orfs$id, orfs$frame)
    upstreamORFs$dist_to_start_nt <- orfs$start_site_nt[m] -
        upstreamORFs$stop_site_nt
    upstreamORFs <- upstreamORFs[which((upstreamORFs$start_site_nt -
                                            orfs$start_site_nt[m]) < 0),]

    if(findExonB == TRUE & nrow(upstreamORFs) > 0){

        upstreamORFs$utr3_length <-
            (orfs$seq_length_nt[match(upstreamORFs$id,orfs$id,)] -
                 upstreamORFs$stop_site_nt) + 1

        widths <- data.frame(w = width(transcripts),
                             id = transcripts$transcript_id)

        pad <- max(table(widths$id))
        if (pad > 1) {
            if (length(unique(transcripts$transcript_id)) == 1) {
                w <- cumsumANDpad(widths$w, pad)
                diffs <-
                    lapply(upstreamORFs$stop_site_nt, function(x)
                        x - w)
                diffs <-
                    matrix(unlist(diffs), ncol = length(upstreamORFs$id))
            } else{
                w2 <-
                    aggregate(w ~ id, widths, function(x)
                        cumsumANDpad(x, pad))
                id <- gsub("_frame1","",
                           gsub("_frame2","",
                                gsub("_frame3","", upstreamORFs$id)))
                m <- match(id, w2$id)
                w2 <- w2[m, -1]
                w2  <- split(w2, seq(nrow(w2)))
                diffs <-
                    mapply(function(x , y)
                        x - y, upstreamORFs$stop_site_nt, w2)
            }

            upstreamORFs$min_dist_to_junction_b <-
                suppressWarnings((apply(diffs, 2, function(x)
                    max(x[x <= 0 & !is.na(x)])) * -1) + 1)
            upstreamORFs$min_dist_to_junction_b[
                which(is.infinite(upstreamORFs$min_dist_to_junction_b))] <-
                upstreamORFs$utr3_length[
                    which(is.infinite(upstreamORFs$min_dist_to_junction_b))]
            upstreamORFs$exon_b_from_final <-
                (apply(diffs, 2, function(x)
                    length(x[x <= 0 & !is.na(x)]))) - 1

            exonNumber <-
                apply(diffs, 2, function(x)
                    length(which(!is.na(x))))
            upstreamORFs$exon_b_from_final[exonNumber == 1] <- 0
        } else{

            upstreamORFs$min_dist_to_junction_b <- upstreamORFs$utr3_length
            upstreamORFs$exon_b_from_final <- 0
        }
    }

    upstreamORFs$start <- NULL
    upstreamORFs$stop <- NULL
    upstreamORFs$dist_to_start <- NULL
    upstreamORFs$utr3_length <- NULL

    m <- match(upstreamORFs$id, orfs$id)
    upstreamORFs$frame <- orfs$frame[m]
    replaceName <- which(stringr::str_sub(upstreamORFs$id, -6,-1) %in% c("frame1", "frame2","frame3"))

    upstreamORFs$id[replaceName] <- gsub("_frame1", "",
                                         gsub("_frame2", "",
                                              gsub("_frame3", "",
                                                   upstreamORFs$id[replaceName])))

    upstreamORFs <- plyr::arrange(upstreamORFs, id, frame, start_site_nt)
    return(upstreamORFs)
}
