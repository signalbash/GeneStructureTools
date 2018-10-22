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

#' Create transcripts with alternative intron usage
#'
#' Creates transcript isoforms from alternative intron usage tested by leafcutter
#' @param altIntronLocs data.frame containing information from the
#' per_intron_results.tab file output from leafcutter.
#' Note that only one cluster of alternative introns can be processed at a time.
#' @param exons GRanges object made from a GTF with ONLY exon annotations
#' (no gene, transcript, CDS etc.)
#' @param replaceInternalExons replace any internal transcript exons with inferred exons from multi-intron leafcutter intron sets?
#' @return GRanges object with all potential alternative isoforms skipping the
#' introns specified in either the upregulated or downregulated locations
#' @export
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
#' @examples
#' leafcutterFiles <- list.files(system.file("extdata","leafcutter/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' leafcutterIntrons <- read.delim(leafcutterFiles[grep("intron_results",
#' leafcutterFiles)],stringsAsFactors=FALSE)
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon"]
#' # single cluster processing
#' cluster <- leafcutterIntrons[leafcutterIntrons$cluster=="chr16:clu_1396",]
#' altIsoforms1396 <- alternativeIntronUsage(cluster, exons)
#' unique(altIsoforms1396$transcript_id)
#' cluster <- leafcutterIntrons[leafcutterIntrons$cluster=="chr16:clu_1395",]
#' altIsoforms1395 <- alternativeIntronUsage(cluster, exons)
#' unique(altIsoforms1395$transcript_id)
#' # multiple cluster processing
#' altIsoforms1396plus1395 <- alternativeIntronUsage(cluster, c(exons, altIsoforms1396))
#' unique(altIsoforms1396plus1395$transcript_id)
alternativeIntronUsage <- function(altIntronLocs, exons,replaceInternalExons=TRUE, junctions=NULL){

    altIntronLocs$deltapsi <- altIntronLocs$PSI_a - altIntronLocs$PSI_b

    clusterGRanges <-
        GRanges(seqnames=S4Vectors::Rle(altIntronLocs$chr),
                ranges=IRanges::IRanges(start=as.numeric(altIntronLocs$start),
                                        end=as.numeric(altIntronLocs$end)),
                strand="*",
                id=altIntronLocs$clusterID,
                direction=ifelse(altIntronLocs$deltapsi >0, "+","-"),
                verdict=altIntronLocs$verdict,
                deltapsi=altIntronLocs$deltapsi)


    m <- match(altIntronLocs$gene, exons$gene_name)
    if(all(!is.na(m))){
    strand(clusterGRanges)[which(!is.na(m))] <-
        strand(exons)[m][which(!is.na(m))]
    }else{
        near <- nearest(clusterGRanges, exons)
        strand(clusterGRanges) <-
            strand(exons)[near]
    }

    # maximum spanning region
    clusterGRanges.max <- clusterGRanges
    #start(clusterGRanges.max) <- min(start(clusterGRanges.max))
    #end(clusterGRanges.max) <- max(end(clusterGRanges.max))

    #find overlaps -- for when range overlaps multiple genes
    olExons <- as.data.frame(findOverlaps(clusterGRanges.max, exons))
    if(nrow(olExons) > 0){
        transcripts <-
            exonsToTranscripts(exons[exons$gene_id %in%
                                         exons$gene_id[olExons$subjectHits]])
        # find transcripts which contain the cluster region
        olTrans <- as.data.frame(findOverlaps(clusterGRanges.max,
                                              transcripts, type = "within"))
        clusterTranscripts <- transcripts[unique(olTrans$subjectHits)]
        #transcript_exons <- exons[exons$transcript_id %in%
        # clusterTranscripts$transcript_id,]
        clusterExons <- exons[exons$transcript_id %in%
                                  clusterTranscripts$transcript_id,]


        if(length(clusterExons) > 0){

            # add sets to cluster introns
            if(any(clusterGRanges$direction=="-")){
                clusterGRanges.dnre <-
                    addSets(clusterGRanges[clusterGRanges$direction=="-"])
                setlist = unique(clusterGRanges.dnre$set)
                clusterGRanges.dnre$set = match(clusterGRanges.dnre$set, setlist)
            }
            if(any(clusterGRanges$direction=="+")){
                clusterGRanges.upre <-
                    addSets(clusterGRanges[clusterGRanges$direction=="+"])
                setlist = unique(clusterGRanges.upre$set)
                clusterGRanges.upre$set = match(clusterGRanges.upre$set, setlist)
            }
            if(any(clusterGRanges$direction=="+") & any(clusterGRanges$direction=="-")){
                clusterGRanges.upre$set <-
                    clusterGRanges.upre$set + max(clusterGRanges.dnre$set)
                clusterGRanges <- c(clusterGRanges.upre, clusterGRanges.dnre)
            }else if(any(clusterGRanges$direction=="+")){
                clusterGRanges <- clusterGRanges.upre
            }else if(any(clusterGRanges$direction=="-")){
                clusterGRanges <- clusterGRanges.dnre
            }

            clusterGRanges$set <- clusterGRanges$set - min(clusterGRanges$set) + 1

            overlaps <- as.data.frame(findOverlaps(clusterGRanges, clusterExons))
            rmExons <- clusterExons[unique(overlaps$subjectHits)]

            clusterGRanges.intron <- clusterGRanges
            start(clusterGRanges.intron) <- start(clusterGRanges.intron) +1
            end(clusterGRanges.intron) <- end(clusterGRanges.intron) -1

            for(s in seq_along(unique(clusterGRanges$set))){
                verdicts <- clusterGRanges$verdict[clusterGRanges$set==s]

                #make a new intermediate exon for funsies
                if(length(which(verdicts %in% c(c("cryptic_threeprime","cryptic_fiveprime","cryptic_unanchored")))) > 1){
                    clusterGRanges.newExon <- clusterGRanges[clusterGRanges$set==s]
                    clusterGRanges.newExon <- clusterGRanges.newExon[order(start(clusterGRanges.newExon))]

                    newExons <- clusterExons[rep(1, length(clusterGRanges.newExon)-1)]
                    ranges(newExons) <- IRanges(start=end(clusterGRanges.newExon)[-length(clusterGRanges.newExon)],
                                                end=start(clusterGRanges.newExon)[-1])

                    newExons$exon_id <- paste0("newEXON_", c(1:(length(clusterGRanges.newExon)-1)))

                    tids <- unique(clusterExons$transcript_id)
                    newExons <- rep(newExons, length(tids))
                    newExons$transcript_id <- rep(tids, each=length(clusterGRanges.newExon)-1)

                    clusterExons <- c(clusterExons, newExons)
                    clusterExons <- reorderExonNumbers(clusterExons)

                }

            }


            move <- which(clusterGRanges$verdict %in% c("cryptic_threeprime","cryptic_fiveprime","cryptic_unanchored"))
            clusterExons.novel <- NULL
            setTrack <- vector()
            for(m in seq_along(move)){
                if(clusterGRanges$verdict[move[m]] == "cryptic_threeprime"){
                    clusterExons.alt <- makeNewLeafExons(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons, splice.type = 3)
                }else if(clusterGRanges$verdict[move[m]] == "cryptic_fiveprime"){
                    clusterExons.alt <- makeNewLeafExons(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons, splice.type = 5)
                }else if(clusterGRanges$verdict[move[m]] == "cryptic_unanchored"){
                    clusterExons.alt <- makeNewLeafExonsUnanchored(altIntronLocs, clusterGRanges.intron[move[m]], clusterExons)
                }
                setTrack <- c(setTrack, rep(clusterGRanges.intron$set[move[m]], length(clusterExons.alt)))
                clusterExons.novel <- c(clusterExons.novel, clusterExons.alt)
            }
            if(length(move) >1){
                clusterExons.novel <- do.call("c", clusterExons.novel)
            }else if(length(move) == 1){
                clusterExons.novel <- clusterExons.novel[[1]]
            }
            pairs <- which(clusterGRanges$verdict %in% c("novel annotated pair"))
            clusterExons.pairs <- NULL
            for(p in seq_along(pairs)){
                clusterGRanges.start <- clusterGRanges[pairs[p]]
                end(clusterGRanges.start) <- start(clusterGRanges.start)
                ol <- findOverlaps(clusterGRanges.start, clusterExons, type="end")
                clusterExons.start <- clusterExons[ol@to]
                rep.end <- length(clusterExons.start)
                start.ids <- clusterExons.start$transcript_id

                clusterGRanges.end <- clusterGRanges[pairs[p]]
                start(clusterGRanges.end) <- end(clusterGRanges.end)
                ol <- findOverlaps(clusterGRanges.end, clusterExons, type="start")
                clusterExons.end <- clusterExons[ol@to]
                rep.start <- length(clusterExons.end)
                end.ids <- clusterExons.end$transcript_id

                clusterExons.start <- rep(clusterExons.start, rep.start)
                clusterExons.start$transcript_id <- rep(end.ids, each=rep.end)

                clusterExons.end <- rep(clusterExons.end, rep.end)
                clusterExons.end$transcript_id <- rep(start.ids, each=rep.start)

                clusterExons.pairs <- c(clusterExons.pairs, clusterExons.start, clusterExons.end)

                # remember to reorder exon numbers
                setTrack <- c(setTrack, rep(clusterGRanges.intron$set[pairs[p]], length(clusterExons.start)+length(clusterExons.end)))

            }
            if(is.list(clusterExons.pairs) & length(clusterExons.pairs) > 1){
                clusterExons.pairs <- do.call("c", clusterExons.pairs)
            }
            if(!is.null(clusterExons.pairs) & !is.null(clusterExons.novel)){
                clusterExons.novel <- c(clusterExons.novel, clusterExons.pairs)
            }else if(!is.null(clusterExons.pairs) & is.null(clusterExons.novel)){
                clusterExons.novel <- clusterExons.pairs
            }

            #rm(clusterExons.allSets)
            clusterExons.allSets <- NULL

            if(!is.null(junctions)){
                # find junctions which are within the cluster ranges,
                # but aren't the same junctions
                ol <- findOverlaps(clusterGRanges, junctions)
                jnc <- junctions[ol@to]

                newJunc <- which(!(paste0(start(jnc), "-", end(jnc)+1) %in%
                                       paste0(start(clusterGRanges), "-", end(clusterGRanges))))
                jnc <- jnc[newJunc]
            }

            for(i in seq_along(unique(clusterGRanges$set))){

                if(i %in% setTrack){
                    clusterExons.tid <- paste0(clusterExons$exon_id, clusterExons$transcript_id)
                    clusterExons.novel.tid <- paste0(clusterExons.novel$exon_id, clusterExons.novel$transcript_id)[setTrack==i]
                    keep <- which(!(clusterExons.tid %in% clusterExons.novel.tid))
                    clusterExons.replace <- c(clusterExons[keep], clusterExons.novel[setTrack==i])
                }else{
                    clusterExons.replace <- clusterExons
                }

                clusterExons.replace <- reorderExonNumbers(clusterExons.replace)

                clusterGRanges.max <-
                    clusterGRanges.intron[clusterGRanges.intron$set==i]
                start(clusterGRanges.max) <- min(start(clusterGRanges.max))
                end(clusterGRanges.max) <- max(end(clusterGRanges.max))

                overlapsCluster <- findOverlaps(clusterGRanges.max, clusterExons.replace)
                rmExons <- clusterExons.replace[unique(overlapsCluster@to)]
                if(length(rmExons) > 0){
                    clusterExonsBounding <- clusterExons.replace[-unique(overlapsCluster@to)]
                }else{
                    clusterExonsBounding <- clusterExons.replace
                }

                #overlaps start of the intron
                clusterGRanges.start <- clusterGRanges[clusterGRanges$set==i]
                end(clusterGRanges.start) <- start(clusterGRanges.start)
                overlapsStart <-
                    findOverlaps(clusterGRanges.start, rmExons, type = "end")
                #overlaps end of the intron
                clusterGRanges.end <- clusterGRanges[clusterGRanges$set==i]
                start(clusterGRanges.end) <- end(clusterGRanges.end)
                overlapsEnd <-
                    findOverlaps(clusterGRanges.end, rmExons, type = "start")

                exonsStart <- rmExons[overlapsStart@to]
                exonsStart <- removeSameExon(exonsStart)

                exonsEnd <- rmExons[overlapsEnd@to]
                exonsEnd <- removeSameExon(exonsEnd)

                InternalExons <- removeSameExon(c(exonsStart,exonsEnd))

                # make sure exons are within the intronic region
                overlapsIntron <-
                    as.data.frame(findOverlaps(InternalExons,
                                               clusterGRanges.max, type="within"))
                InternalExons <- InternalExons[unique(overlapsIntron$queryHits)]

                overlapsExonStart <- findOverlaps(clusterGRanges.start,
                                                  clusterExonsBounding, type = "end")
                overlapsExonEnd <- findOverlaps(clusterGRanges.end,
                                                clusterExonsBounding, type = "start")

                keepTranscriptIds <-
                    clusterExonsBounding$transcript_id[overlapsExonStart@to]
                keepTranscriptIds <-
                    keepTranscriptIds[clusterExonsBounding$transcript_id[
                        overlapsExonStart@to] %in%
                            clusterExonsBounding$transcript_id[overlapsExonEnd@to]]

                if(length(keepTranscriptIds) > 0){

                    if(any(grepl("_dnre_", keepTranscriptIds) |
                           grepl("_upre_", keepTranscriptIds))){
                        # direction to to remove
                        # if upre, remove dnre isforms
                        removeDirection <- ifelse(clusterGRanges$direction[
                            match(i, clusterGRanges$set)[1]] == "+", "dnre", "upre")
                        keepTranscriptIds <- keepTranscriptIds[
                            !grepl(removeDirection, keepTranscriptIds)]
                    }

                    clusterExonsBounding <- clusterExonsBounding[
                        clusterExonsBounding$transcript_id %in% keepTranscriptIds]

                    if(replaceInternalExons == FALSE){
                        InternalExons.reps <-
                            InternalExons[rep(seq_along(InternalExons),
                                              length(unique(keepTranscriptIds)))]
                        InternalExons.reps$transcript_id <-
                            rep(unique(keepTranscriptIds), each=length(InternalExons))

                        clusterExonsBounding <- c(clusterExonsBounding, InternalExons.reps)
                    }else{
                        if(!is.null(junctions)){
                            replaceExons <- leafcutterIntronsToExons(clusterGRanges[clusterGRanges$set==i])
                            # keep junctions which overlap any exons created from the cluster junctions
                            if(!is.null(replaceExons)){
                                keep <- findOverlaps(jnc, replaceExons, type="within")
                                jnc.set <- jnc[unique(keep@from)]
                                if(length(jnc.set) > 0){
                                    # annotate the same as clusterGRanges
                                    jnc.set$id <- clusterGRanges$id[clusterGRanges$set==i][1]
                                    jnc.set$direction <- clusterGRanges$direction[clusterGRanges$set==i][1]
                                    jnc.set$verdict <- "junction"
                                    jnc.set$deltapsi <- 0
                                    jnc.set$set <- i
                                    jnc.set$count <- NULL

                                    # remove any overlapping junctions
                                    jnc.set <- jnc.set[rev(order(width(jnc.set)))]
                                    ol <- findOverlaps(jnc.set)
                                    rm <- as.data.frame(table(ol@from))
                                    while(any(rm$Freq > 1)){
                                        rm <- rm[rev(order(rm$Freq)),]
                                        jnc.set <- jnc.set[-as.numeric(rm$Var1[1])]
                                        ol <- findOverlaps(jnc.set)
                                        rm <- as.data.frame(table(ol@from))
                                    }
                                    replaceExons <- leafcutterIntronsToExons(c(clusterGRanges[clusterGRanges$set==i], jnc.set))
                                }
                            }

                        }else{
                            replaceExons <- leafcutterIntronsToExons(clusterGRanges[clusterGRanges$set==i])
                        }
                        if(!is.null(replaceExons)){

                            replaceExonsFormat <- clusterExonsBounding[rep(1, length(replaceExons))]
                            ranges(replaceExonsFormat) <- ranges(replaceExons)

                            replaceExonsFormat.reps <-
                                replaceExonsFormat[rep(seq_along(replaceExonsFormat),
                                              length(unique(keepTranscriptIds)))]
                            replaceExonsFormat.reps$transcript_id <-
                                rep(unique(keepTranscriptIds), each=length(replaceExonsFormat))

                            clusterExonsBounding <- c(clusterExonsBounding, replaceExonsFormat.reps)
                        }
                    }
                    clusterExonsBounding <- reorderExonNumbers(clusterExonsBounding)
                    clusterExonsBounding$set <- as.numeric(i)

                    if(!exists("clusterExons.allSets")){
                        clusterExons.allSets <- clusterExonsBounding
                    }else{
                        clusterExons.allSets <- c(clusterExons.allSets,
                                                  clusterExonsBounding)
                    }
                }

            }

            if(!is.null(clusterExons.allSets)){
                if(length(clusterExons.allSets) == 1){
                    clusterExons.allSets <- clusterExons.allSets[[1]]
                }else if(length(clusterExons.allSets) > 1){
                    clusterExons.allSets <- do.call("c", clusterExons.allSets)
                }

                m <- match(clusterExons.allSets$set, clusterGRanges$set)
                #n <- which(!grepl("[+]", clusterExons.allSets$transcript_id))
                clusterExons.allSets$new_transcript_id <- NA
                clusterExons.allSets$new_transcript_id <-
                    paste0(clusterExons.allSets$transcript_id, "+AS ",
                           ifelse(clusterGRanges$direction[m] == "+", "upre","dnre"),
                           "_",
                           gsub("_", "", clusterGRanges$id[m]),
                           "-", clusterExons.allSets$set)

                n <- which(grepl("[+]", clusterExons.allSets$transcript_id))
                clusterExons.allSets$new_transcript_id[n] <-
                    paste0(clusterExons.allSets$transcript_id, ":",
                           gsub("_", "", clusterGRanges$id[m]),
                           "-", clusterExons.allSets$set)[n]

                clusterExons.allSets$transcript_id <-
                    clusterExons.allSets$new_transcript_id
                clusterExons.allSets$new_transcript_id <- NULL
                clusterExons.allSets$set <- NULL

                dnre <- clusterExons.allSets[grep("dnre", clusterExons.allSets$transcript_id)]
                upre <- clusterExons.allSets[grep("upre", clusterExons.allSets$transcript_id)]
                if(length(dnre) > 0 & length(upre) > 0){
                    clusterExons.allSets <- c(
                        removeDuplicateTranscripts(dnre),
                        removeDuplicateTranscripts(upre))
                    return(clusterExons.allSets)
                }else{
                    return(NULL)
                }
            }else{
                return(NULL)
            }
        }else{
            return(NULL)
        }
    }else{
        return(NULL)
    }
}

#' Convert an exon-level gtf annotation to a transcript-level gtf annotation
#'
#' @param exons GRanges object with exons
#' @return GRanges object with transcripts
#' @export
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import
#' @family gtf manipulation
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools"))
#' exons <- gtf[gtf$type=="exon" & gtf$transcript_id=="ENSMUST00000126412.1"]
#' exons
#' transcripts <- exonsToTranscripts(exons)
#' transcripts
exonsToTranscripts <- function(exons){
    transcripts <- exons[!duplicated(exons$transcript_id)]

    minStarts <- aggregate(start ~ transcript_id,
                           as.data.frame(exons), min)
    maxEnds <- aggregate(end ~ transcript_id,
                         as.data.frame(exons), max)


    ranges(transcripts) <-
        IRanges::IRanges(start=as.numeric(
            minStarts$start[match(transcripts$transcript_id,
                                  minStarts$transcript_id)]),
            end=as.numeric(
                maxEnds$end[match(transcripts$transcript_id,
                                  maxEnds$transcript_id)]))


    return(transcripts)
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

#' Create junction ranges from STAR junction files
#'
#' Create junction ranges from STAR junction files
#' @param junctionFiles list of files with STAR junctions
#' @param minReads minimum number of reads a junction requires to be kept
#' @return GRanges object with junctions
#' @export
#' @importFrom IRanges IRanges
#' @import GenomicRanges
#' @family leafcutter splicing isoform creation
#' @author Beth Signal
readLeafcutterJunctions <- function(junctionFiles, minReads=10){

    junctions.all <- NULL
    for(f in seq_along(junctionFiles)){
        junctions <- read.delim(junctionFiles[f], header=F)

        index <- match(paste(junctions$V1,junctions$V2,junctions$V3,junctions$V6, sep="_"),
                       paste(junctions.all$V1,junctions.all$V2,junctions.all$V3,junctions.all$V6, sep="_"))

        junctions.all$V5[index[which(!is.na(index))]] <- junctions.all$V5[index[which(!is.na(index))]] + junctions$V5[which(!is.na(index))]

        junctions.all <- rbind(junctions.all, junctions[which(is.na(index)),])
    }
    junctions.all <- junctions.all[junctions.all$V5 >= minReads,]

    junctionGrange <- GRanges(seqnames=junctions.all$V1,
                              ranges=IRanges(start=junctions.all$V2,
                                             end=junctions.all$V3),
                              strand=junctions.all$V6,
                              count=junctions.all$V5)

    return(junctionGrange)

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
