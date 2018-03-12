#' Convert DEXSeq ids to gene ids
#'
#' @param DEXSeqIds vector of DEXSeq group or exon ids
#' @param removeVersion remove the version (.xx) of the gene?
#' @param containsE do the DEXSeq exons ids contain :E00X?
#' @return vector of unique gene ids
#' @export
#' @import stringr
#' @family DEXSeq processing methods
#' @examples
#' # multiple genes in name
#' DEXSeqId <- "ENSMUSG00000027618.17+ENSMUSG00000098950.7+ENSMUSG00000089824.10+ENSMUSG00000074643.12"
#' DEXSeqIdsToGeneIds(DEXSeqId)
#'
#' # exonic part number in id
#' DEXSeqIdsToGeneIds("ENSMUSG00000001017.15:E013", removeVersion=TRUE)
#' @author Beth Signal
DEXSeqIdsToGeneIds <- function(DEXSeqIds, removeVersion=FALSE, containsE=TRUE){

    if(containsE){
        dexSplit <- ":E"
    }else{
        dexSplit <- ":"
    }

    containsExon <- grep(dexSplit, DEXSeqIds)

    if(length(containsExon) >0 ){
        DEXSeqIds[containsExon] <- unlist(lapply(stringr::str_split(
            DEXSeqIds[containsExon], dexSplit), "[[",1))
    }

    geneIds <- unique(unlist(stringr::str_split(DEXSeqIds, "[+]")))

    if(removeVersion==TRUE){
        geneIds <- removeVersion(geneIds)
    }

    return(geneIds)
}

#' Remove version number from ensembl gene/transcript ids
#'
#' @param ids vector of ensembl ids
#' @return vector of ensembl ids without the version number
#' @import stringr
#' @export
#' @examples
#' removeVersion("ENSMUSG00000001017.15")
#' @author Beth Signal
removeVersion <- function(ids){
    return(unlist(lapply(stringr::str_split(ids, "[.]"), "[[",1)))
}

#' Find a DEXSeq exons' biotype
#'
#' @param DEXSeqExonId vector of DEXSeq exon ids
#' @param DEXSeqGtf GRanges object of the DEXSeq formatted gtf
#' @param gtf GRanges object of the GTF annotated with exon biotypes - i.e. exon, CDS, UTR
#' @param set which overlapping set of exon biotypes to return - to, from, and/or overlap
#' @return overlaping types
#' @export
#' @import GenomicRanges
#' @family DEXSeq processing methods
#' @importFrom rtracklayer import
#' @author Beth Signal
#' @examples
#' gtfFile <- system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools")
#' DEXSeqGtfFile <- system.file("extdata","gencode.vM14.dexseq.gtf",
#' package = "GeneStructureTools")
#'
#' gtf <- rtracklayer::import(gtfFile)
#' gtf <- UTR2UTR53(gtf)
#' DEXSeqGtf <- rtracklayer::import(DEXSeqGtfFile)
#'
#' findDEXexonType("ENSMUSG00000032366.15:E028", DEXSeqGtf, gtf)
#'
#' DEXSeqResultsFile <- system.file("extdata","dexseq_results_significant.txt",
#' package = "GeneStructureTools")
#' DEXSeqResults <- read.table(DEXSeqResultsFile, sep="\t")
#'
#' findDEXexonType(rownames(DEXSeqResults), DEXSeqGtf, gtf)
#'
findDEXexonType <- function(DEXSeqExonId, DEXSeqGtf, gtf,set="overlap"){
    DEXSeqGtf$id <- paste0(DEXSeqGtf$gene_id,":E", DEXSeqGtf$exonic_part_number)
    DEXSeqGtf.query <- DEXSeqGtf[match(DEXSeqExonId,DEXSeqGtf$id)]
    overlapTypes <- overlapTypes(DEXSeqGtf.query, gtf, set = set)[,2]
    return(overlapTypes)
}

#' Summarise exon biotypes to broader categories
#' @param types vector of exon biotypes
#' @return vector of broader exon biotypes
#' @export
#' @importFrom rtracklayer import
#' @family DEXSeq processing methods
#' @examples
#' gtfFile <- system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools")
#' DEXSeqGtfFile <- system.file("extdata","gencode.vM14.dexseq.gtf",
#' package = "GeneStructureTools")
#'
#' gtf <- rtracklayer::import(gtfFile)
#' gtf <- UTR2UTR53(gtf)
#' DEXSeqGtf <- rtracklayer::import(DEXSeqGtfFile)
#'
#' findDEXexonType("ENSMUSG00000032366.15:E028", DEXSeqGtf, gtf)
#'
#' DEXSeqResultsFile <- system.file("extdata","dexseq_results_significant.txt",
#' package = "GeneStructureTools")
#' DEXSeqResults <- read.table(DEXSeqResultsFile, sep="\t")
#'
#' types <- findDEXexonType(rownames(DEXSeqResults), DEXSeqGtf, gtf)
#' summarisedTypes <- summariseExonTypes(types)
#' table(types, summarisedTypes)
#' @author Beth Signal
summariseExonTypes <- function(types){

    types <- gsub("protein_coding-CDS:protein_coding-UTR3:protein_coding-UTR5",
                  "protein_coding-CDS", types)
    types <- gsub("protein_coding-CDS:protein_coding-UTR5",
                  "protein_coding-start_codon",
                  types)
    types <- gsub("protein_coding-CDS:protein_coding-UTR3",
                  "protein_coding-stop_codon",
                  types)

    types2 <- types
    types2[grep("protein_coding-start_codon", types2)] <- "start_codon"
    types2[grep("protein_coding-stop_codon", types2)] <- "stop_codon"
    types2[grep("protein_coding-UTR5", types2)] <- "UTR5"
    types2[grep("protein_coding-UTR3", types2)] <- "UTR3"
    types2[grep("protein_coding-UTR", types2)] <- "UTR"
    types2[grep("protein_coding-CDS", types2)] <- "CDS"
    types2[!(types2 %in% c("start_codon",
                           "stop_codon","UTR5","UTR3",
                           "CDS", "UTR")) & !is.na(types2)] <-
        "noncoding_exon"

    return(types2)
}


#' Annotate introns and exonic parts by overlaping exon biotype
#'
#' Annotate introns and exonic parts by overlaping exon biotype
#' @param queryCoords GRanges object of the query regions
#' @param gtf GRanges object of the GTF annotated with exon biotypes - i.e. exon, CDS, UTR
#' @param set which overlapping set of exon biotypes to return - to, from, and/or overlap
#' @return overlaping types in a data.frame
#' @export
#' @import GenomicRanges
#' @importFrom stats aggregate
#' @author Beth Signal
#' @examples
#' gtfFile <- system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools")
#' DEXSeqGtfFile <- system.file("extdata","gencode.vM14.dexseq.gtf",
#' package = "GeneStructureTools")
#'
#' gtf <- rtracklayer::import(gtfFile)
#' gtf <- UTR2UTR53(gtf)
#' DEXSeqGtf <- rtracklayer::import(DEXSeqGtfFile)
#'
#' overlapTypes(DEXSeqGtf[1:10], gtf)
overlapTypes <- function(queryCoords, gtf, set=c("from", "to", "overlap")){
    overlaps <- as.data.frame(GenomicRanges::findOverlaps(queryCoords, gtf))
    gtf.overlap <- gtf[overlaps$subjectHits]
    gtf.overlap$index <- overlaps$queryHits
    gtf.overlap <- gtf.overlap[(gtf.overlap$type %in%
                                    c("exon", "CDS","UTR","UTR3","UTR5"))]

    gtf.overlap <- filterGtfOverlap(gtf.overlap)
    gtf.overlap <- addBroadTypes(gtf.overlap)

    gtf.from <- NULL
    gtf.to <- NULL

    if(any(set=="from")){
        gtf.from <- gtf.overlap[end(gtf.overlap) ==
                                    start(queryCoords[gtf.overlap$index])]
    }
    if(any(set=="to")){
        gtf.to <- gtf.overlap[start(gtf.overlap) ==
                                  end(queryCoords[gtf.overlap$index])]
    }
    #keep only hits with a exon-intron-exon pair
    if(any(set=="to") & any(set=="from") & length(gtf.from) > 0 &
       length(gtf.to) >0){
        tidIndex.from <- paste0(gtf.from$transcript_id, "_",
                                gtf.from$index)
        tidIndex.to <- paste0(gtf.to$transcript_id, "_",
                              gtf.to$index)
        gtf.from <- gtf.from[tidIndex.from %in% tidIndex.to]
        gtf.to <- gtf.to[tidIndex.to %in% tidIndex.from]
    }

    if(any(set=="from") & length(gtf.from) > 0){
        gtf.from$typetype <- paste0(gtf.from$transcript_type_broad,
                                    "-",gtf.from$type)
        #remove nmd/retained introns -- these tend to be isoexons of protein coding exons
        rm <- which(gtf.from$typetype == "retained_intron|exon" |
                        gtf.from$transcript_type_broad == "nmd")

        #not used currently
        fromTypes <- aggregate(type ~ index, mcols(gtf.from),
                               function(x) paste0(sort(unique(x)),
                                                  collapse=":"))
        #not used currently
        fromTranscriptTypes <- aggregate(transcript_type_broad ~ index,
                                         mcols(gtf.from),
                                         function(x) paste0(sort(unique(x)),
                                                            collapse=":"))
        fromTypeTypes <- aggregate(typetype ~ index,
                                   mcols(gtf.from)[-rm,],
                                   function(x) paste0(sort(unique(x)),
                                                      collapse=":"))
    }
    if(any(set=="to") & length(gtf.to) > 0){
        gtf.to$typetype <- paste0(gtf.to$transcript_type_broad,
                                  "-",gtf.to$type)
        #remove nmd/retained introns -- these tend to be isoexons of protein coding exons
        rm <- which(gtf.to$typetype == "retained_intron|exon" |
                        gtf.to$transcript_type_broad == "nmd")

        toTypes <- aggregate(type ~ index, mcols(gtf.to),
                             function(x) paste0(sort(unique(x)),collapse=":"))
        toTranscriptTypes <- aggregate(transcript_type_broad ~ index,
                                       mcols(gtf.to),
                                       function(x) paste0(sort(unique(x)),
                                                          collapse=":"))
        toTypeTypes <- aggregate(typetype ~ index,
                                 mcols(gtf.to)[-rm,],
                                 function(x) paste0(sort(unique(x)),
                                                    collapse=":"))
    }
    if(any(set=="overlap")){
        gtf.overlap$typetype <- paste0(
            gtf.overlap$transcript_type_broad,
            "-",gtf.overlap$type)
        #remove nmd/retained introns -- these tend to be isoexons of protein coding exons
        keep <- which(!(gtf.overlap$typetype == "retained_intron|exon" |
                            gtf.overlap$transcript_type_broad == "nmd"))

        overlapTypes <- aggregate(type ~ index, mcols(gtf.overlap),
                                  function(x) paste0(sort(unique(x)),
                                                     collapse=":"))
        overlapTranscriptTypes <- aggregate(transcript_type_broad ~ index,
                                            mcols(gtf.overlap),
                                            function(x) paste0(sort(unique(x)),
                                                               collapse=":"))
        overlapTypeTypes <- aggregate(typetype ~ index,
                                      mcols(gtf.overlap)[keep,],
                                      function(x) paste0(sort(unique(x)),
                                                         collapse=":"))
    }

    typeTypes <- data.frame(index=1:length(start(ranges(queryCoords))))
    if(any(set=="from") & length(gtf.from) > 0){
        typeTypes$from <- fromTypeTypes$typetype[match(typeTypes$index,
                                                       fromTypeTypes$index)]
    }
    if(any(set=="to") & length(gtf.to) > 0){
        typeTypes$to <- toTypeTypes$typetype[match(typeTypes$index,
                                                   toTypeTypes$index)]
    }
    if(any(set=="overlap")){
        typeTypes$overlap <-
            overlapTypeTypes$typetype[match(typeTypes$index,
                                            overlapTypeTypes$index)]
    }

    return(typeTypes)
}

#' Change transcript biotypes to a broader set
#'
#' Change transcript biotypes to a broader set in a GRanges GTF object
#' @param gtf GRanges object of the GTF
#' @return GRanges object of the GTF with new transcript types
#' @export
#' @import GenomicRanges
#' @importFrom rtracklayer import
#' @family gtf manipulation
#' @author Beth Signal
#' @examples
#' gtfFile <- system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools")
#' gtf <- rtracklayer::import(gtfFile)
#' gtf <- addBroadTypes(gtf)
addBroadTypes <- function(gtf){
    transcriptTypesBroad <- gtf$transcript_type
    transcriptTypesBroad[which(
        transcriptTypesBroad %in% c(
            "3prime_overlapping_ncrna",
            "3prime_overlapping_ncRNA",
            "antisense",
            "bidirectional_promoter_lncRNA",
            "macro_lncRNA",
            "known_ncrna",
            "lincRNA",
            "non_coding",
            "processed_transcript",
            "sense_intronic",
            "sense_overlapping"
        )
    )] <- "lncRNA"
    transcriptTypesBroad[which(
        transcriptTypesBroad %in% c(
            "IG_C_gene",
            "IG_C_pseudogene",
            "IG_D_gene",
            "IG_J_gene",
            "IG_J_pseudogene",
            "IG_V_gene",
            "IG_D_pseudogene",
            "IG_LV_gene",
            "IG_pseudogene",
            "ribozyme",
            "IG_V_pseudogene",
            "miRNA",
            "misc_RNA",
            "Mt_rRNA",
            "Mt_tRNA",
            "rRNA",
            "snoRNA",
            "snRNA",
            "TEC",
            "scaRNA",
            "scRNA",
            "sRNA",
            "TR_C_gene",
            "TR_D_gene",
            "TR_J_gene",
            "TR_J_pseudogene" ,
            "TR_V_gene",
            "TR_V_pseudogene"
        )
    )] <- "short_ncRNA"
    transcriptTypesBroad[which(
        transcriptTypesBroad %in% c(
            "processed_pseudogene",
            " pseudogene",
            "transcribed_processed_pseudogene",
            "transcribed_unitary_pseudogene",
            "transcribed_unprocessed_pseudogene",
            "translated_processed_pseudogene",
            "translated_unprocessed_pseudogene",
            "polymorphic_pseudogene",
            "unitary_pseudogene",
            "unprocessed_pseudogene"
        )
    )] <- "pseudogene"

    transcriptTypesBroad[which(transcriptTypesBroad %in%
                                   c("nonsense_mediated_decay",
                                     "non_stop_decay"))] <-"nmd"


    gtf$transcript_type_broad <-
        transcriptTypesBroad
    return(gtf)
}

#' Filter a GTF overlap to remove exons when exon is annotated as a CDS/UTR
#'
#' Filter a GTF overlap to remove exons when exon is annotated as a CDS/UTR
#' @param gtf.from GRanges object of the GTF produced from an overlap
#' @return GRanges object of the GTF with redundant exons removed
#' @export
#' @import GenomicRanges
#' @importFrom stats aggregate
#' @importFrom rtracklayer import
#' @family gtf manipulation
#' @author Beth Signal
#' @examples
#' gtfFile <- system.file("extdata","example_gtf.gtf",
#' package = "GeneStructureTools")
#' gtf <- rtracklayer::import(gtfFile)
#' overlap <- as.data.frame(GenomicRanges::findOverlaps(gtf[which(gtf$type=="CDS")[1]], gtf))
#' table(gtf$type[overlap$subjectHits])
#' overlapFiltered <- filterGtfOverlap(gtf[overlap$subjectHits])
#' table(overlapFiltered$type[overlap$subjectHits])

#' overlap <- as.data.frame(GenomicRanges::findOverlaps(gtf[which(
#' gtf$transcript_type=="retained_intron")[1]],gtf))
#' table(gtf$type[overlap$subjectHits])
#' overlapFiltered <- filterGtfOverlap(gtf[overlap$subjectHits])
#' table(overlapFiltered$type[overlap$subjectHits])

filterGtfOverlap <- function(gtf.from){
    gtf.fromDF <- as.data.frame(mcols(gtf.from))
    gtf.fromDF$exon_number <- as.numeric(gtf.fromDF$exon_number)
    gtf.fromDF$start_ids <- paste0(start(ranges(gtf.from)),
                                   gtf.from$transcript_id)
    gtf.fromDF$end_ids <- paste0(end(ranges(gtf.from)),
                                 gtf.from$transcript_id)

    rmStart <- which(gtf.fromDF$type == "exon" &
                         gtf.fromDF$start_ids %in%
                         gtf.fromDF$start_ids[gtf.fromDF$type %in%
                                                  c("CDS","UTR","UTR3","UTR5")])
    rmEnd <- which(gtf.fromDF$type == "exon" &
                       gtf.fromDF$end_ids %in%
                       gtf.fromDF$end_ids[gtf.fromDF$type %in%
                                              c("CDS","UTR","UTR3","UTR5")])
    rm <- unique(c(rmEnd, rmStart))
    gtf.from <- gtf.from[-rm]
    return(gtf.from)
}
