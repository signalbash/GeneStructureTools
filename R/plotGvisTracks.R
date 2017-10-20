#' Annotate a GRanges gene model with ORF boundries for visualisation with Gviz
#' @param transcripts GRanges of gene model to be visualised
#' @param orfs ORF predictions. Created by getORFs()
#' @return data.frame of a gene model for visualisation
#' @export
#' @import plyr
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
annotateGeneModel <- function(transcripts, orfs){

    # location to site wise-point
    transcripts_loc <- as.data.frame(transcripts)
    transcripts_loc <- transcripts_loc[,c("start","end","width","strand","exon_id","transcript_id","exon_number")]

    strand <- unique(transcripts_loc$strand)

    # row for each geneomic location
    if(strand == "+"){
        loc_all <- apply(transcripts_loc, 1, function(x) x[1]:x[2])
    }else{
        loc_all <- apply(transcripts_loc, 1, function(x) x[2]:x[1])
    }
    # add exon/transcript annoations
    loc_all.exon <- apply(transcripts_loc, 1, function(x) rep(x[5], x[3]))
    loc_all.transcript <- apply(transcripts_loc, 1, function(x) rep(x[6], x[3]))
    loc_all.exon_num <- apply(transcripts_loc, 1, function(x) rep(x[7], x[3]))

    transcripts_loc.bySite <- data.frame(loc=unlist(loc_all),
                                         exon=unlist(loc_all.exon),
                                         transcript=unlist(loc_all.transcript),
                                         exon_number=unlist(loc_all.exon_num))
    # add relative sites
    if(strand == "+"){
        transcripts_loc.bySite <- arrange(transcripts_loc.bySite, transcript, loc)
    }else{
        transcripts_loc.bySite <- arrange(transcripts_loc.bySite, transcript, plyr::desc(loc))
    }
    transcripts_loc.bySite$site <- unlist(lapply(aggregate(width ~ transcript_id, transcripts_loc, sum)[,2], function(x) 1:x))

    # run through each transcript individually
    transcript_ids <- unique(transcripts_loc.bySite$transcript)
    for(t in seq_along(transcript_ids)){

        utr <- transcripts
        utr <- utr[utr$transcript_id == transcript_ids[t]]
        utr$type <- as.character(utr$transcript_type)

        orfs_t <- orfs[orfs$id == transcript_ids[t],]
        transcripts_loc.bySite_t <-
            transcripts_loc.bySite[transcripts_loc.bySite$transcript == transcript_ids[t],]

        #utr5/3 exon number
        utr5.exon <- as.numeric(as.character(transcripts_loc.bySite_t$exon_number[transcripts_loc.bySite_t$site==orfs_t$start_site_nt[which.max(orfs_t$orf_length)]]))
        utr3.exon <- as.numeric(as.character(transcripts_loc.bySite_t$exon_number[transcripts_loc.bySite_t$site==orfs_t$stop_site_nt[which.max(orfs_t$orf_length)]]))

        #double up the utr containing exons
        utr.cds <- utr[utr$exon_number %in% c(utr3.exon, utr5.exon)]

        # move the UTR5 boundry
        utr$type[as.numeric(utr$exon_number) <= utr5.exon] <- "utr5"
        if(strand == "+"){
            end(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$start_site_nt[which.max(orfs_t$orf_length)]]
        }else{
            start(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$start_site_nt[which.max(orfs_t$orf_length)]]
        }
        # move the UTR3 boundry
        utr$type[as.numeric(utr$exon_number) >= utr3.exon] <- "utr3"
        if(strand == "+"){
            start(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$stop_site_nt[which.max(orfs_t$orf_length)]]
        }else{
            end(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$stop_site_nt[which.max(orfs_t$orf_length)]]
        }
        # move the CDS boundries
        utr <- c(utr, utr.cds)
        utr$type[!(utr$type %in% c("utr3","utr5"))] <- "CDS"

        if(strand == "+"){
            start(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon & utr$type=="CDS"] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$start_site_nt[which.max(orfs_t$orf_length)]] + 1

            end(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon & utr$type=="CDS"] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$stop_site_nt[which.max(orfs_t$orf_length)]] - 1
        }else{
            end(ranges(utr))[as.numeric(utr$exon_number) == utr5.exon & utr$type=="CDS"] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$start_site_nt[which.max(orfs_t$orf_length)]] + 1

            start(ranges(utr))[as.numeric(utr$exon_number) == utr3.exon & utr$type=="CDS"] <-
                transcripts_loc.bySite_t$loc[transcripts_loc.bySite_t$site == orfs_t$stop_site_nt[which.max(orfs_t$orf_length)]] - 1
        }

        if(!exists("geneModel") | t == 1){
            geneModel <- as.data.frame(utr)
        }else{
            geneModel <- rbind(geneModel, as.data.frame(utr))
        }
    }


    n <- match(c('seqnames','start','end','width','strand','type','gene_id','exon_id','transcript_id'), colnames(geneModel))
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
#' @param transcripts GRanges of gene model to be visualised
#' @return data.frame of a gene model for visualisation
#' @export
#' @import GenomicRanges
#' @examples
#' @author Beth Signal
makeGeneModel <- function(transcript){
    transcript <- as.data.frame(transcript)
    n <- match(c('seqnames','start','end','width','strand','type','gene_id','exon_id','transcript_id'), colnames(transcript))

    replace <- FALSE

    if(is.na(n)[8] & is.na(n)[6]){
        n[8] <- match("transcript_id", colnames(transcript))
        n[6] <- match("transcript_id", colnames(transcript))
        replace <- TRUE
        exon_id <- paste0(transcript$transcript_id, "_",transcript$exon_number)
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
        transcript$exon <- exon_id
    }

    return(transcript)

}
