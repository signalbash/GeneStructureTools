#' Annotate a GRanges gene model with ORF boundries for visualisation with Gviz
#' @param transcripts GRanges of gene model to be visualised
#' @param orfs ORF predictions. Created by getORFs()
#' @return data.frame of a gene model for visualisation
#' @export
#' @importFrom plyr desc
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "example_gtf.gtf",
#' package="GeneStructureTools"))
#' transcript <- gtf[gtf$type=="exon" & gtf$gene_name=="Neurl1a"]
#' g <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' # longest ORF for each transcripts
#' orfs <- getOrfs(transcript, BSgenome = g, returnLongestOnly = TRUE)
#' geneModelAnnotated <- annotateGeneModel(transcript, orfs)
annotateGeneModel <- function(transcripts, orfs){

    # location to site wise-point
    transcriptsLocation <- as.data.frame(transcripts)
    transcriptsLocation <- transcriptsLocation[,c("start","end","width",
                                                  "strand","exon_id",
                                                  "transcript_id","exon_number")]

    strand <- unique(transcriptsLocation$strand)

    # row for each geneomic location
    if(strand == "+"){
        loctionsAll <- apply(transcriptsLocation, 1, function(x) x[1]:x[2])
    }else{
        loctionsAll <- apply(transcriptsLocation, 1, function(x) x[2]:x[1])
    }
    # add exon/transcript annoations
    loctionsAll.exon <- apply(transcriptsLocation, 1, function(x) rep(x[5], x[3]))
    loctionsAll.transcript <- apply(transcriptsLocation, 1, function(x) rep(x[6], x[3]))
    loctionsAll.exonNum <- apply(transcriptsLocation, 1, function(x) rep(x[7], x[3]))

    transcriptsLocation.bySite <- data.frame(`loc`=unlist(loctionsAll),
                                         `exon`=unlist(loctionsAll.exon),
                                         `transcript`=unlist(loctionsAll.transcript),
                                         `exon_number`=unlist(loctionsAll.exonNum))
    # add relative sites
    if(strand == "+"){
        transcriptsLocation.bySite <- transcriptsLocation.bySite[
            order(transcriptsLocation.bySite$transcript,
                  transcriptsLocation.bySite$loc),]
    }else{
        transcriptsLocation.bySite <- transcriptsLocation.bySite[
            order(transcriptsLocation.bySite$transcript,
                  plyr::desc(transcriptsLocation.bySite$loc)),]
    }
    transcriptsLocation.bySite$site <- unlist(lapply(aggregate(width ~ transcript_id,
                                                               transcriptsLocation, sum)[,2],
                                                     function(x) 1:x))

    # run through each transcript individually
    transcriptIds <- unique(transcriptsLocation.bySite$transcript)
    for(t in seq_along(transcriptIds)){

        utr <- transcripts
        utr <- utr[utr$transcript_id == transcriptIds[t]]
        utr$type <- as.character(utr$transcript_type)

        orfs.t <- orfs[orfs$id == transcriptIds[t],]
        transcriptsLocation.bySite.t <-
            transcriptsLocation.bySite[transcriptsLocation.bySite$transcript == transcriptIds[t],]

        #utr5/3 exon number
        utr5.exon <- as.numeric(as.character(
            transcriptsLocation.bySite.t$exon_number[
                transcriptsLocation.bySite.t$site==
                    orfs.t$start_site_nt[which.max(orfs.t$orf_length)]]))
        utr3.exon <- as.numeric(as.character(
            transcriptsLocation.bySite.t$exon_number[
                transcriptsLocation.bySite.t$site==
                    orfs.t$stop_site_nt[which.max(orfs.t$orf_length)]]))

        #double up the utr containing exons
        utr.cds <- utr[utr$exon_number %in% c(utr3.exon, utr5.exon)]

        # move the UTR5 boundry
        utr$type[as.numeric(utr$exon_number) <= utr5.exon] <- "utr5"
        if(strand == "+"){
            end(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$start_site_nt[which.max(orfs.t$orf_length)]]
        }else{
            start(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$start_site_nt[which.max(orfs.t$orf_length)]]
        }
        # move the UTR3 boundry
        utr$type[as.numeric(utr$exon_number) >= utr3.exon] <- "utr3"
        if(strand == "+"){
            start(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$stop_site_nt[which.max(orfs.t$orf_length)]]
        }else{
            end(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$stop_site_nt[which.max(orfs.t$orf_length)]]
        }
        # move the CDS boundries
        utr <- c(utr, utr.cds)
        utr$type[!(utr$type %in% c("utr3","utr5"))] <- "CDS"

        if(strand == "+"){
            start(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon & utr$type=="CDS"] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$start_site_nt[which.max(orfs.t$orf_length)]] + 1

            end(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon & utr$type=="CDS"] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$stop_site_nt[which.max(orfs.t$orf_length)]] - 1
        }else{
            end(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon & utr$type=="CDS"] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$start_site_nt[which.max(orfs.t$orf_length)]] + 1

            start(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon & utr$type=="CDS"] <-
                transcriptsLocation.bySite.t$loc[
                    transcriptsLocation.bySite.t$site ==
                        orfs.t$stop_site_nt[which.max(orfs.t$orf_length)]] - 1
        }

        if(!exists("geneModel") | t == 1){
            geneModel <- as.data.frame(utr)
        }else{
            geneModel <- rbind(geneModel, as.data.frame(utr))
        }
    }


    n <- match(c('seqnames','start','end',
                 'width','strand','type',
                 'gene_id','exon_id','transcript_id'), colnames(geneModel))
    n.symbol <- which(colnames(geneModel) == "gene_name")
    if(length(n.symbol) == 1){
        geneModel <- geneModel[,c(n, n.symbol)]
    }else{
        geneModel <- geneModel[,n]
        geneModel$gene_symbol <- NA
    }

    colnames(geneModel) <- c("chromosome", "start","end","width","strand",
                             "feature","gene","exon","transcript","symbol")
    return(geneModel)
}
#' Convert GRanges gene model to data.frame for visualisation with Gviz
#' @param transcript GRanges of gene model to be visualised
#' @return data.frame of a gene model for visualisation
#' @export
#' @import GenomicRanges
#' @author Beth Signal
#' @examples
#' gtf <- rtracklayer::import(system.file("extdata", "example_gtf.gtf",
#' package="GeneStructureTools"))
#' transcript <- gtf[gtf$type=="exon" & gtf$gene_name=="Neurl1a"]
#' geneModel <- makeGeneModel(transcript)
makeGeneModel <- function(transcript){
    transcript <- as.data.frame(transcript)
    n <- match(c('seqnames','start','end','width',
                 'strand','type','gene_id',
                 'exon_id','transcript_id'), colnames(transcript))

    replace <- FALSE

    if(is.na(n)[8] & is.na(n)[6]){
        n[8] <- match("transcript_id", colnames(transcript))
        n[6] <- match("transcript_id", colnames(transcript))
        replace <- TRUE
        exonId <- paste0(transcript$transcript_id, "_",transcript$exon_number)
    }else if(is.na(n)[6] & !is.na(n)[8]){
        n[6] <- match("transcript_type", colnames(transcript))
    }
    n.symbol <- which(colnames(transcript) == "gene_name")
    if(length(n.symbol) == 1){
        transcript <- transcript[,c(n, n.symbol)]
    }else{
        transcript <- transcript[,n]
        transcript$gene_symbol <- NA
    }

    colnames(transcript) <- c("chromosome", "start","end","width","strand",
                             "feature","gene","exon","transcript","symbol")
    if(replace){
        transcript$feature <- "protein_coding"
        transcript$exon <- exonId
    }

    return(transcript)

}
