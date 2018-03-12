#' Read in a list of whippet .jnc.gz files and format as a GRanges object
#' @param files vector of *.jnc.gz file names
#' @return GRanges object with junctions
#' @export
#' @importFrom data.table fread
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' jncFiles <- whippetFiles[grep(".jnc", whippetFiles)]
#' whippetJNC <- readWhippetJNCfiles(jncFiles)
readWhippetJNCfiles <- function(files){

    # check platform for windows/unix specific decompression
    operatingSystem <- .Platform$OS.type
    if(operatingSystem == "unix"){
        ungzip <- "zcat < "
    }else{ #in windows
        ungzip <- "gzip -dc "
    }

    for(f in seq_along(files)){
        whip <- data.table::fread(paste0(ungzip, files[f]), data.table=FALSE)
        colnames(whip) <- c("chromosome","start","end","id","count","strand")

        if(exists("whip.all")){
            m <- match(whip$id, whip.all$id)
            n <- match(whip.all$id, whip$id)
            whip.all <- cbind(whip.all, new_count=whip[n,5])
            addCols <- ncol(whip.all) - 6

            if(any(is.na(m))){
                whip.add <- whip[which(is.na(m)),c(1:4,6)]
                for(c in 1:(addCols)){
                    whip.add <- cbind(whip.add, count=NA)
                }
                whip.add <- cbind(whip.add, count=whip$count[which(is.na(m))])
                colnames(whip.add) <- colnames(whip.all)
                whip.all <- rbind(whip.all,whip.add)
            }
        }else{
            whip.all <- whip[,c(1:4,6,5)]
        }
    }

    colnames(whip.all)[-(1:5)] <- gsub(".jnc.gz","", basename(files))

    geneIds <- unlist(lapply(stringr::str_split(whip.all$id, ":"), "[[", 1))

    jncCoords <- GRanges(seqnames=S4Vectors::Rle(whip.all$chrom),
                         ranges=IRanges::IRanges(start=as.numeric(whip.all$start),
                                                 end=as.numeric(whip.all$end)),
                         strand=whip.all$strand, id=whip.all$id, gene=geneIds)
    return(jncCoords)
}

#' Read in a list of whippet .psi.gz files and format as a data.frame
#' @param files vector of *.psi.gz file names
#' @param attribute which attribute from the PSI files to use (Total_Reads, Psi, CI_width)
#' @param maxNA maximum number of NA values allowed before a site is removed
#' @return data.frame with junction counts for all files
#' @export
#' @importFrom data.table fread
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' psiFiles <- whippetFiles[grep(".psi", whippetFiles)]
#' whippetPSI <- readWhippetPSIfiles(psiFiles)
readWhippetPSIfiles <- function(files, attribute="Total_Reads", maxNA=NA){

    # check platform for windows/unix specific decompression
    operatingSystem <- .Platform$OS.type
    if(operatingSystem == "unix"){
        ungzip <- "zcat < "
    }else{ #in windows
        ungzip <- "gzip -dc "
    }

    for(f in seq_along(files)){
        whip <- data.table::fread(paste0(ungzip, files[f]), data.table=FALSE)
        wantedCol <- which(colnames(whip) == attribute)

        if(exists("whip.all")){
            whip.all <- cbind(whip.all, whip[,wantedCol])
        } else{
            whip.all <- whip[,c(1:5,wantedCol)]
        }
    }

    colnames(whip.all)[-(1:5)] <- gsub(".psi.gz","", basename(files))

    whip.all$na_count <- apply(whip.all[,-c(1:5)], 1, function(x)
        length(which(is.na(x))))
    whip.all$unique_name <- paste(whip.all$Gene,whip.all$Coord,
                                  whip.all$Type,whip.all$Node, sep="_")
    if(!is.na(maxNA)){
        keep <- which(whip.all$na_count <= maxNA)
        whip.all <- whip.all[keep,]
    }
    return(whip.all)
}

#' Read in a list of whippet .diff.gz files and format as a data.frame
#' @param files vector of *.diff.gz file names
#' @return data.frame with junction counts for all files
#' @export
#' @importFrom data.table fread
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' diffFiles <- whippetFiles[grep(".diff", whippetFiles)]
#' whippetDiffSplice <- readWhippetDIFFfiles(diffFiles)
readWhippetDIFFfiles <- function(files){

    # check platform for windows/unix specific decompression
    operatingSystem <- .Platform$OS.type
    if(operatingSystem == "unix"){
        ungzip <- "zcat < "
    }else{ #in windows
        ungzip <- "gzip -dc "
    }

    for(f in seq_along(files)){
        whip <- data.table::fread(paste0(ungzip, files[f]),
                                  data.table=FALSE, skip=1)

        #remove the NA column
        keepcol <- which(apply(whip, 2, function(x) all(!is.na(x))))
        whip <- whip[,keepcol]

        colnames(whip) <- c("gene", "node","coord","strand","type",
                            "psi_a","psi_b","psi_delta",
                            "probability","complexity","entropy")
        whip$unique_name <- paste(whip$gene,whip$coord,whip$type,
                                  whip$node, sep="_")
        whip$comparison <- gsub(".diff.gz","", basename(files[f]))
        conditions <- unique(whip$comparison)
        conditionsSplit <- stringr::str_split(conditions, "_")
        condition1 <- unlist(lapply(conditionsSplit, "[[", 1))
        condition2 <- unlist(lapply(conditionsSplit, function(x) tail(x, n=1)))

        whip$condition_1 <- condition1[match(conditions, whip$comparison)]
        whip$condition_2 <- condition2[match(conditions, whip$comparison)]

        if(exists("whip.all")){
            whip.all <- rbind(whip.all, whip)
        } else{
            whip.all <- whip
        }
    }



    return(whip.all)
}

#' Format Whippet co-ordinates as a GRanges object
#' @param whippet data.frame containing event location information.
#' May be generated by readWhippetDIFFfiles()
#' @return GRanges object with events
#' @export
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @import stringr
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- list.files(system.file("extdata","whippet/",
#' package = "GeneStructureTools"), full.names = TRUE)
#' diffFiles <- whippetFiles[grep(".diff", whippetFiles)]
#' whippetDiffSplice <- readWhippetDIFFfiles(diffFiles)
#' whippetCoords <- formatWhippetEvents(whippetDiffSplice)
formatWhippetEvents <- function(whippet){

    whippet <- whippet[!duplicated(whippet$coord),]

    chromosome <- unlist(lapply(stringr::str_split(whippet$coord,":"), "[[",1))
    range <- unlist(lapply(stringr::str_split(whippet$coord,":"), "[[",2))
    start <- unlist(lapply(stringr::str_split(range,"-"), "[[",1))
    end <- unlist(lapply(stringr::str_split(range,"-"), "[[",2))

    eventCoords <- GRanges(seqnames=S4Vectors::Rle(chromosome),
                           ranges=IRanges::IRanges(start=as.numeric(start),
                                                   end=as.numeric(end)),
                           strand=whippet$strand, id=whippet$coord)
    return(eventCoords)
}

#' Import whippet results files as a whippetDataSet
#' @param filePath path to whippet output files
#' @return whippetDataSet
#' @export
#' @import stringr
#' @import methods
#' @family whippet data processing
#' @author Beth Signal
#' @examples
#' whippetFiles <- system.file("extdata","whippet/",
#' package = "GeneStructureTools")
#' wds <- readWhippetDataSet(whippetFiles)
readWhippetDataSet <- function(filePath="."){

    fileNames <- list.files(filePath)
    wds <- new("whippetDataSet", filePath=filePath)
    filesDiff <- fileNames[grep(".diff.gz", fileNames)]
    if(length(filesDiff) > 0){
        slot(wds, "diffSplicingResults") <- readWhippetDIFFfiles(
            paste0(filePath, "/", filesDiff))
        slot(wds, "comparisons") <- unique(diffSplicingResults(wds)$comparison)
        slot(wds, "coordinates") <-
            formatWhippetEvents(diffSplicingResults(wds))
        m <- match(diffSplicingResults(wds)$coord, coordinates(wds)$id)
        slot(wds, "diffSplicingResults") <-
            cbind(diffSplicingResults(wds), coord_match=m)
    }

    filesJnc <- fileNames[grep(".jnc.gz", fileNames)]
    if(length(filesJnc) > 0){
        slot(wds, "junctions") <-
            readWhippetJNCfiles(paste0(filePath, "/", filesJnc))
    }

    filesPsi <- fileNames[grep(".psi.gz", fileNames)]
    if(length(filesPsi) > 0){
        slot(wds, "readCounts") <-
            readWhippetPSIfiles(paste0(filePath, "/", filesPsi))
    }

    return(wds)
}
