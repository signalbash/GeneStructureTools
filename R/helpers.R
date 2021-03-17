#' Find overlaps from a rmats GRanges to reference exons and annotate with ids
#' @param rmatsGRanges rmats event GRanges object
#' @param exons reference exons GRanges
#' @param exon_number which exon in the rmats event this is for (only used for annotating output)
#' @return data.frame with overlapping event/exons
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
annotateOverlapRmats <- function(rmatsGRanges, exons, exon_number=1){
    ol <- as.data.frame(findOverlaps(rmatsGRanges, exons))
    ol$from_id <- rmatsGRanges$id[ol$queryHits]
    ol$transcript_id <- exons$transcript_id[ol$subjectHits]
    ol$exon_id <- exons$exon_id[ol$subjectHits]
    ol$exon_number <- exons$exon_number[ol$subjectHits]
    ol$event_id <- paste0(seqnames(rmatsGRanges)[ol$queryHits], ":", rmatsGRanges$event_id[ol$queryHits])
    colnames(ol)[match(c('exon_id', 'exon_number') ,colnames(ol))] <- paste0(c('exon_id', 'exon_number'), exon_number)
    return(ol)
}

#' remove any duplicate pairs of events/reference transcripts (i.e. long event range which overlaps 2+ exons)
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @return data.frame with related differential splicing event ids and reference transcript_ids
#' @keywords internal
#' @import methods
#' @importFrom rlang .data
#' @family rmats data processing
#' @author Beth Signal
removeDuplicatePairs <- function(betweenExons){

    hasDups <- which(duplicated(paste0(betweenExons$new_transcript_id)))
    if(length(hasDups) > 0){
        betweenExons.duplicates <- betweenExons[betweenExons$new_transcript_id %in% betweenExons$new_transcript_id[hasDups],]
        betweenExons.duplicates$exon_num_range <- abs(as.numeric(betweenExons.duplicates$exon_number1) - as.numeric(betweenExons.duplicates$exon_number2))
        betweenExons.duplicates <- arrange(betweenExons.duplicates, new_transcript_id, plyr::desc(exon_num_range))
        betweenExons.duplicates <- betweenExons.duplicates[!duplicated(betweenExons.duplicates$new_transcript_id),]
        betweenExons.duplicates$exon_num_range <- NULL

        betweenExons <- rbind(betweenExons[which(!(betweenExons$new_transcript_id %in% betweenExons$new_transcript_id[hasDups])),],
                              betweenExons.duplicates)
    }

    return(betweenExons)

}

#' Duplicate a reference Granges (exon-level) for each diff-splicing event/transcript combination required.
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param exons reference exons GRanges
#' @return
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
duplicateReference <- function(betweenExons, exons){
    # transcripts containing the exon pairs
    transcripts <- as.data.frame(table(betweenExons$transcript_id))
    gtfTranscripts <- exons[exons$transcript_id %in% transcripts$Var1]

    m <- match(gtfTranscripts$transcript_id, betweenExons$transcript_id)

    mcols(gtfTranscripts) <- cbind(mcols(gtfTranscripts),
                                   DataFrame(new_transcript_id=paste0(
                                       gtfTranscripts$transcript_id, "+AS ",
                                       betweenExons$from_id[m], "-",
                                       betweenExons$event_id[m])))

    needsDuplicated <- which(!(betweenExons$new_transcript_id %in%
                                   gtfTranscripts$new_transcript_id))

    if(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts[
            gtfTranscripts$transcript_id %in%
                betweenExons$transcript_id[needsDuplicated]]
    }

    while(length(needsDuplicated) > 0){
        gtfTranscripts_add <- gtfTranscripts_add[
            gtfTranscripts_add$transcript_id %in%
                betweenExons$transcript_id[needsDuplicated]]
        m <- match(gtfTranscripts_add$transcript_id,
                   betweenExons$transcript_id[needsDuplicated])
        gtfTranscripts_add$new_transcript_id <- paste0(
            gtfTranscripts_add$transcript_id, "+AS ",
            betweenExons$from_id[needsDuplicated][m], "-",
            betweenExons$event_id[needsDuplicated][m])
        gtfTranscripts <- c(gtfTranscripts, gtfTranscripts_add)
        needsDuplicated <- which(!(betweenExons$new_transcript_id %in%
                                       gtfTranscripts$new_transcript_id))
    }

    return(gtfTranscripts)
}

#' Generate vector of integers between two numbers (non-inclusive)
#' @param a first integer
#' @param b second integer
#' @return vector of integers
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
betweenNumbers <- function(a, b){
    ab.range <- seq(as.numeric(a),as.numeric(b))
    ab.range <- ab.range[!(ab.range %in% c(a,b))]
    return(ab.range)
}


#' Remove any exons in a transcript within an event range
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param gtfTranscripts Granges with altered transcript structures
#' @return Granges with altered transcript structures
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
removeExonsBetween <- function(betweenExons, gtfTranscripts){
    remove <- apply(betweenExons[,c('exon_number1', 'exon_number2', 'new_transcript_id')], 1, function(x) paste0(x[3], " ", betweenNumbers(x[1], x[2])))
    gtfTranscripts.rm <- gtfTranscripts[which(!(paste0(gtfTranscripts$new_transcript_id, " ", gtfTranscripts$exon_number) %in% unlist(remove)))]
    return(gtfTranscripts.rm)
}

#' Split long exons in two if they overlap an event
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param gtfTranscripts Granges with altered transcript structures
#' @return Granges with altered transcript structures
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
splitLongExons <- function(betweenExons, gtfTranscripts){
    # make sure there is > 1 exon in the pair (i.e. split a single long exon into two)
    duplicate.index <- which(betweenExons$exon_number1 == betweenExons$exon_number2)
    if(length(duplicate.index) > 0){
        duplicate <- paste0(betweenExons$exon_id1[duplicate.index], "_", betweenExons$new_transcript_id[duplicate.index])
        betweenExons$exon_id1[duplicate.index] <- paste0(betweenExons$exon_id1[duplicate.index], "_1")
        betweenExons$exon_id2[duplicate.index] <- paste0(betweenExons$exon_id2[duplicate.index], "_2")
        betweenExons$exon_number1[duplicate.index] <- as.numeric(betweenExons$exon_number1[duplicate.index]) - 0.5
        betweenExons$exon_number2[duplicate.index] <- as.numeric(betweenExons$exon_number2[duplicate.index]) + 0.5

        duplicate.granges.1 <- gtfTranscripts[paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate]
        duplicate.granges.1$exon_number <- as.character(as.numeric(duplicate.granges.1$exon_number) - 0.5)
        duplicate.granges.1$exon_id <- paste0(duplicate.granges.1$exon_id, "_1")
        duplicate.granges.2 <- gtfTranscripts[paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate]
        duplicate.granges.2$exon_number <- as.character(as.numeric(duplicate.granges.2$exon_number) + 0.5)
        duplicate.granges.2$exon_id <- paste0(duplicate.granges.2$exon_id, "_2")

        gtfTranscripts <- c(gtfTranscripts[which(!(paste0(gtfTranscripts$exon_id, "_", gtfTranscripts$new_transcript_id) %in% duplicate))],
                           duplicate.granges.1, duplicate.granges.2)
    }

    splitReturn <- list(ranges=gtfTranscripts, between=betweenExons)
    return(splitReturn)
}


#' annotate event coordinates with exon/transcript/gene ids.
#' @param eventCoords Granges with event coordinates
#' @param betweenExons data.frame with related differential splicing event ids and reference transcript_ids
#' @param exons reference exons GRanges
#' @return eventCoords Granges with annotations
#' @keywords internal
#' @import methods
#' @family rmats data processing
#' @author Beth Signal
annotateEventCoords <- function(eventCoords, betweenExons, exons){

    #eventCoords <- GRanges(seqnames=input$chr, ranges=IRanges(start=input$exonStart_0base+1, end=input$exonEnd), strand=input$strand, id=input$ID)
    #seqlevelsStyle(eventCoords)=seqlevelsStyle(exons)[1]

    # create ranges for skipped exons
    eventCoords <- eventCoords[match(betweenExons$from_id, eventCoords$id)]
    eventCoords$exon_number <- NA
    eventCoords$transcript_id <- betweenExons$transcript_id
    eventCoords$transcript_type <- exons$transcript_type[match(eventCoords$transcript_id, exons$transcript_id)]
    eventCoords$gene_id <- exons$gene_id[match(eventCoords$transcript_id, exons$transcript_id)]
    eventCoords$exon_id <- unlist(lapply(str_split(betweenExons$new_transcript_id, "[+]AS[ ]"), "[[", 2))
    eventCoords$exon_number <- NA
    mcols(eventCoords) <- cbind(mcols(eventCoords),
                                DataFrame(new_transcript_id=paste0(
                                    eventCoords$transcript_id,"+AS ",
                                    eventCoords$exon_id)))
    return(eventCoords)

}

#' Find overlaps where the start/end coordinates are the same
#' @param query a GRanges object
#' @param subject a GRanges object
#' @return Hits object
#' @keywords internal
#' @import methods
#' @importFrom rlang .data
#' @family data processing
#' @author Beth Signal
findOverlaps.junc = function(query, subject, type=c("start", "end")){

    if(any(type == "start")){

        query.start <- query
        end(query.start) <- start(query.start)
        subject.start <- subject
        end(subject.start) <- start(subject.start)

        ol.start <- findOverlaps(query.start, subject.start, type="equal")
    }

    if(any(type == "end")){

        query.end <- query
        start(query.end) <- end(query.end)
        subject.end <- subject
        start(subject.end) <- end(subject.end)

        ol.end <- findOverlaps(query.end, subject.end, type="equal")
    }

    if("start" %in% type & "end" %in% type){
        ol.df <- rbind(as.data.frame(ol.start), as.data.frame(ol.end))
        ol.df <- arrange(ol.df, queryHits, subjectHits)
        ol <- S4Vectors::Hits(from=ol.df$queryHits, to=ol.df$subjectHits, nLnode=S4Vectors::nLnode(ol.start), nRnode=S4Vectors::nRnode(ol.start), sort.by.query=TRUE)
    }else if("start" %in% type){
        ol <- ol.start
    }else if("end" %in% type){
        ol <- ol.end
    }

    return(ol)

}
#' Create a 0/1 matrix from a Ranges overlap
#' @param ol ranges overlaps from findOverlaps()
#' @return matrix of all ranges by all ranges, with 0/1 coded overlaps
#' @keywords internal
#' @import GenomicRanges
#' @author Beth Signal
overlap2matrix <- function(ol, maxn=7){
    mat <- matrix(nrow = maxn, ncol = maxn, data = 1)

    for(i in 1:nrow(mat)){
        mat[i,i] <- 0
        mat[i, ol$subjectHits[ol$queryHits == i]] <- 0
    }

    return(mat)
}
#' Get a list of all potential range combinations, covering the full length of the region.
#' @param mat matrix of all ranges by all ranges, with 0/1 coded overlaps
#' @return list of range combinations possible
#' @keywords internal
#' @import GenomicRanges
#' @author Beth Signal
matrix2combinations <-function(mat){
    maxn <- nrow(mat)
    ignorecols <- vector()
    newMat <- matrix(nrow=1, ncol=maxn, data=NA)
    newMatLine <- newMat

    for(j in 1:maxn){

        if(all(is.na(newMatLine[,j]))){
            newMatLine[,j] <- 1

            if(any(mat[-c(j, ignorecols),j] == 0) & !all(mat[-c(j, ignorecols),j] == 0)){
                n <- which(mat[,j] == 0)
                n <- n[which(!(n %in% c(j, ignorecols)))]

                if(length(n) <= 1){
                    newMatLineadd <- newMatLine
                    if(n > j){
                        newMatLineadd[,j] <- NA
                        newMatLineadd[,n] <- 1

                        newMatLine <- rbind(newMatLine, newMatLineadd)
                    }else{
                        newMatLineadd[,j][which(!is.na(newMatLineadd[,n]))] <- NA
                        newMatLine <- newMatLineadd
                    }
                }else{
                    newMatLine[,j] <- NA
                    newMatLineadd <- newMatLine

                    n.low <- n[n < j]
                    for(nx in n.low){
                        newMatLineadd <- newMatLineadd[is.na(newMatLineadd[,nx]),]
                        if(!is.matrix(newMatLineadd)){
                            newMatLineadd <- t(as.matrix(newMatLineadd))
                        }
                    }
                    n.high <- n[n > j]

                    newMatLineadd[,n.high] <- NA
                    newMatLineadd[,j] <- 1
                    newMatLine <- rbind(newMatLine, newMatLineadd)

                }
            }else if(all(mat[-j,j] == 0)){
                ignorecols <- append(ignorecols, j)
                newMatLine[,-j] <- NA
                newMatLine[,j] <- 1

                if(exists("newMat.singles")){
                    newMat.singles <- rbind(newMat.singles,newMatLine[1,])
                }else{
                    newMat.singles <- newMatLine[1,]
                }
                newMatLine <- newMat
            }
            newMat <- newMatLine
            newMatLine <- newMat
        }

    }
    if(exists("newMat.singles")){
        newMat <- rbind(newMat, newMat.singles)
    }

    rownames(newMat) <- NULL
    names(newMat) <- NULL

    # lens <- vector()
    # for(i in 1:nrow(newMat)){
    #     colz <- which(!is.na(newMat[i,]))
    #     if(length(colz) > 1){
    #         lens[i] <- length(which(apply(newMat[,colz],1, function(x) all(x == 1))))
    #     }else{
    #         lens[i] <- length(which(newMat[,colz] == 1))
    #     }
    # }
    # rm <- which(lens > 1)
    # if(length(rm) > 0){
    #     newMat <- newMat[-rm,]
    # }

    for(x in 1:nrow(mat)){

        dontCombine <- which(mat[x,] == 0)
        dontCombine <- dontCombine[dontCombine != x]

        for(y in seq_along(dontCombine)){
            rm <- which(newMat[,x] == 1 & newMat[,dontCombine[y]] == 1)
            if(length(rm) > 0){
                newMat <- newMat[-rm,]
            }
        }
    }


    combos <- (apply(newMat, 1, function(x) which(!is.na(x))))
    if(is.matrix(combos)){
        combinationList <- list()
        for(i in 1:ncol(combos)){
            combinationList[[i]] <- combos[,i]
        }
        return(combinationList)
    }else{
        return(combos)
    }

}
#' Add set numbers to introns
#'
#' Converts a group of introns into all non-overlapping sets
#' @param clusterGRanges.noset Granges object with a cluster of intron locations
#' @return Granges object with a cluster of intron locations and corresponding set numbers
#' @keywords internal
#' @import GenomicRanges
#' @importFrom plyr desc
#' @author Beth Signal
addSets <- function(clusterGRanges.noset){

    clusterGRanges.noset$set <- 1

    ol <- as.data.frame(findOverlaps(clusterGRanges.noset))
    ol <- ol[ol$queryHits != ol$subjectHits,]
    ol$setFrom <- clusterGRanges.noset$set[ol$queryHits]
    ol$setTo <- clusterGRanges.noset$set[ol$subjectHits]
    ol <- ol[ol$setFrom == ol$setTo,]

    if(nrow(ol) > 0){
        combinationList <- matrix2combinations(overlap2matrix(ol, maxn=length(clusterGRanges.noset)))

        sets <- unlist(mapply(function(x,y) rep(x,y), 1:length(combinationList), lapply(combinationList, length)))

        clusterGRanges.sets <- clusterGRanges.noset[unlist(combinationList)]
        clusterGRanges.sets$set <- sets
    }else{
        clusterGRanges.sets <- clusterGRanges.noset
        clusterGRanges.sets$set <- 1
    }

    setlist <- unique(clusterGRanges.sets$set)
    newSetList <- 1:length(setlist)

    clusterGRanges.sets$set <- newSetList[match(clusterGRanges.sets$set, setlist)]

    return(clusterGRanges.sets)
}

#' Remove exon duplicates
#'
#' Removes structural duplicates of exons in a GRanges object
#' @param exons GRanges object with exons
#' @return GRanges object with unique exons
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family gtf manipulation
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' exons.duplicated <- c(exons[1:4], exons[1:4])
#' length(exons.duplicated)
#' exons.deduplicated <- removeSameExon(exons.duplicated)
#' length(exons.deduplicated)
removeSameExon <- function(exons){
    samesies <- findOverlaps(exons, type = "equal")
    samesies <- samesies[samesies@from > samesies@to]
    if(length(samesies) > 0){
        exons <- exons[-unique(samesies@from)]
    }
    return(exons)
}
#' Create exon ranges from leafcutter intron ranges
#'
#' Create exon ranges from leafcutter intron ranges
#' @param clusterGRanges GRanges of intron clusters with a 'set' column
#' @return GRanges object with intron clusters converted to exons
#' @keywords internal
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
leafcutterIntronsToExons <- function(clusterGRanges, junctions){
    clusterGRanges <- clusterGRanges[order(start(clusterGRanges))]
    sets <- as.data.frame(table(clusterGRanges$set))
    sets <- sets[which(sets$Freq > 1),]

    for(s in seq_along(nrow(sets))){
        clusterGRanges.set <- clusterGRanges[clusterGRanges$set==sets$Var1]
        clusterGRanges.set <- clusterGRanges.set[order(start(clusterGRanges.set))]
        newStarts <- end(clusterGRanges.set)[-length(clusterGRanges.set)]
        newEnds <- start(clusterGRanges.set)[-1]
        clusterExons <- clusterGRanges.set[-1]
        ranges(clusterExons) <- IRanges(start=newStarts, end=newEnds)

        if(s == 1){
            clusterExons.all <- clusterExons
        }else{
            clusterExons.all <- c(clusterExons.all, clusterExons)
        }
    }

    if(nrow(sets) > 0){
        return(clusterExons.all)
    }else{
        return(NULL)
    }

}
#' Create new exon ranges for a cryptic splice junction
#'
#' Create new exon ranges for a cryptic splice junction
#' @param junctionGRanges GRanges of splice junction (must be intron version)
#' @param clusterExons GRanges of annotated exons from the target gene
#' @param splice.type which end of the junciton is cryptic? 3 or 5.
#' @return GRanges object with new exons
#' @keywords internal
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
makeNewLeafExons <- function(altIntronLocs, junctionRanges, clusterExons, splice.type){

    strand.var <- as.character(strand(junctionRanges))

    # find anything that overlaps the 'cryptic' end of the junction
    junctionRanges.end <- junctionRanges
    if((strand.var == "+" & splice.type==3) | (strand.var == "-" & splice.type==5)){
        start(junctionRanges.end) <- end(junctionRanges.end)
    }else{
        end(junctionRanges.end) <- start(junctionRanges.end)
    }
    ol <- findOverlaps(junctionRanges.end, clusterExons)

    # if junction is not a smaller version of known exon
    if(length(ol@to) == 0){

        #next exon following/preceding on the cryptic side
        if(splice.type==3){
            near <- precede(junctionRanges, clusterExons)
        }else{
            near <- follow(junctionRanges, clusterExons)
        }
        # if there's nothing (i.e. novel first/last junction)
        if(is.na(near)){
            near <- nearest(junctionRanges, clusterExons)
            clusterExons.alt <- clusterExons[near]
            #make an exon with length==100
            if((strand.var == "+" & splice.type==3) | (strand.var == "-" & splice.type==5)){
                ranges(clusterExons.alt) <- IRanges(start=end(junctionRanges)+1, width=100)
            }else{
                ranges(clusterExons.alt) <- IRanges(end=start(junctionRanges)-1, width=100)
            }
            if(splice.type == 3){
                clusterExons.alt$exon_number <- as.numeric(clusterExons.alt$exon_number)+1
            }else{
                clusterExons.alt$exon_number <- as.numeric(clusterExons.alt$exon_number)-1
            }
        }else{
            clusterExons.alt <- clusterExons[near]
            # all exons w/ same exon start/end
            if((strand.var == "+" & splice.type==3) | (strand.var == "-" & splice.type==5)){
                same <- findOverlaps(clusterExons.alt, clusterExons, type="start")
            }else{
                same <- findOverlaps(clusterExons.alt, clusterExons, type="end")
            }
            clusterExons.alt <- clusterExons[same@to]
        }
        #paired exon -- for transcript id matching
        junctionRanges.start <- junctionRanges
        junctionRanges.cluster <- junctionRanges
        ranges(junctionRanges.cluster) <- IRanges(start=min(altIntronLocs$start),
                                                  end=max(altIntronLocs$end))

        if(splice.type == 3){
            if(strand.var == "+"){
                start(junctionRanges.start) <- start(junctionRanges.start) -1
                end(junctionRanges.start) <- start(junctionRanges.start)
                ol <- findOverlaps(junctionRanges.start, clusterExons, type="end")@to
                end(junctionRanges.cluster) <- start(junctionRanges.cluster)
                ol <- c(ol,findOverlaps(junctionRanges.cluster, clusterExons, type="end")@to)
            }else{
                end(junctionRanges.start) <- end(junctionRanges.start) +1
                start(junctionRanges.start) <- end(junctionRanges.start)
                ol <- findOverlaps(junctionRanges.start, clusterExons, type="start")@to
                start(junctionRanges.cluster) <- end(junctionRanges.cluster)
                ol <- c(ol,findOverlaps(junctionRanges.cluster, clusterExons, type="start")@to)
            }
            ol <- c(ol, precede(junctionRanges, clusterExons))
        }else{
            if(strand.var == "+"){
                end(junctionRanges.start) <- end(junctionRanges.start) +1
                start(junctionRanges.start) <- end(junctionRanges.start)
                ol <- findOverlaps(junctionRanges.start, clusterExons, type="start")@to
                start(junctionRanges.cluster) <- end(junctionRanges.cluster)
                ol <- c(ol,findOverlaps(junctionRanges.cluster, clusterExons, type="start")@to)
            }else{
                start(junctionRanges.start) <- start(junctionRanges.start) -1
                end(junctionRanges.start) <- start(junctionRanges.start)
                ol <- findOverlaps(junctionRanges.start, clusterExons, type="end")@to
                end(junctionRanges.cluster) <- start(junctionRanges.cluster)
                ol <- c(ol,findOverlaps(junctionRanges.cluster, clusterExons, type="end")@to)
            }
            ol <- c(ol, precede(junctionRanges, clusterExons))
        }
        ol <- unique(ol)
        ol <- ol[which(!is.na(ol))]

        if((strand.var == "+" & splice.type==3) | (strand.var == "-" & splice.type==5)){
            ol <- findOverlaps(clusterExons[ol], clusterExons, type="end")
        }else{
            ol <- findOverlaps(clusterExons[ol], clusterExons, type="start")
        }
        ol <- ol[!duplicated(ol@to)]

        ids <- clusterExons$transcript_id[ol@to]
        ce.alt.length <- length(clusterExons.alt)

        #rep ids so transcript id will match across junctions
        clusterExons.alt <- rep(clusterExons.alt, length(ids)+1)
        clusterExons.alt$transcript_id[-(1:ce.alt.length)] <-
            rep(ids, each=ce.alt.length)

        exon.ids <- paste(start(clusterExons.alt), end(clusterExons.alt), clusterExons.alt$transcript_id)
        clusterExons.alt <- clusterExons.alt[!duplicated(exon.ids)]

    }else{
        clusterExons.alt <- clusterExons[ol@to]

        #check for a F*&$#@^& pair
        junctionRanges.cluster <- junctionRanges
        ranges(junctionRanges.cluster) <- IRanges(start=min(altIntronLocs$start),
                                                  end=max(altIntronLocs$end))


        if((strand.var == "+" & splice.type == 3)|(strand.var == "-" & splice.type == 5)){
            end(junctionRanges.cluster) <- start(junctionRanges.cluster)
            ol <- findOverlaps(junctionRanges.cluster, clusterExons, type="end")@to
        }else{
            start(junctionRanges.cluster) <- end(junctionRanges.cluster)
            ol <- findOverlaps(junctionRanges.cluster, clusterExons, type="start")@to
        }

        ids <- clusterExons$transcript_id[ol]

        ce.alt.length <- length(clusterExons.alt)

        #rep ids so transcript id will match across junctions
        clusterExons.alt <- rep(clusterExons.alt, length(ids)+1)
        clusterExons.alt$transcript_id[-(1:ce.alt.length)] <-
            rep(ids, each=ce.alt.length)

        exon.ids <- paste(start(clusterExons.alt), end(clusterExons.alt), clusterExons.alt$transcript_id)
        clusterExons.alt <- clusterExons.alt[!duplicated(exon.ids)]

    }

    #make sure all the ends are the same
    if((strand.var == "+" & splice.type==3) | (strand.var == "-" & splice.type==5)){
        start(clusterExons.alt) <- end(junctionRanges)+1
    }else{
        end(clusterExons.alt) <- start(junctionRanges)-1
    }

    return(clusterExons.alt)
}
#' Create new exon ranges for a cryptic unachored splice junction
#'
#' Create new exon ranges for a cryptic unachored splice junction
#' @param junctionGRanges GRanges of splice junction (must be intron version)
#' @param clusterExons GRanges of annotated exons from the target gene
#' @return GRanges object with new exons
#' @keywords internal
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
makeNewLeafExonsUnanchored <- function(altIntronLocs, junctionRanges, clusterExons){
    clusterExons.alt3 <- makeNewLeafExons(altIntronLocs, junctionRanges, clusterExons, splice.type = 3)
    clusterExons.alt5 <- makeNewLeafExons(altIntronLocs, junctionRanges, clusterExons, splice.type = 5)

    ids <- clusterExons.alt5$transcript_id
    ce.alt.length <- length(clusterExons.alt3)
    #rep ids so transcript id will match across junctions
    clusterExons.alt3 <- rep(clusterExons.alt3, length(ids)+1)
    clusterExons.alt3$transcript_id[-(1:ce.alt.length)] <-
        rep(ids, each=ce.alt.length)
    exon.ids <- paste(start(clusterExons.alt3), end(clusterExons.alt3), clusterExons.alt3$transcript_id)
    clusterExons.alt3 <- clusterExons.alt3[!duplicated(exon.ids)]

    ids <- clusterExons.alt3$transcript_id
    ce.alt.length <- length(clusterExons.alt5)
    #rep ids so transcript id will match across junctions
    clusterExons.alt5 <- rep(clusterExons.alt5, length(ids)+1)
    clusterExons.alt5$transcript_id[-(1:ce.alt.length)] <-
        rep(ids, each=ce.alt.length)
    exon.ids <- paste(start(clusterExons.alt5), end(clusterExons.alt5), clusterExons.alt5$transcript_id)
    clusterExons.alt5 <- clusterExons.alt5[!duplicated(exon.ids)]

    clusterExons <- c(clusterExons.alt3, clusterExons.alt5)
    return(clusterExons)
}
