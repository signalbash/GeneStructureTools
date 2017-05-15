#' Annotate introns and exonic parts by overlaping exon biotype
#'
#' Annotate introns and exonic parts by overlaping exon biotype
#' @param nt_sequence nucleotide sequence(s) to have longest orfs predicted
#' @return data.frame with sequence index, length, and translation frame of the longest ORF predicted
#' @export
#' @import stringr
#' @import Biostrings
#' @importFrom plyr ldply
#' @examples
#' @author Beth Signal
longestORF <- function(nt_sequence){
    for(frame in 1:3){

        normal_seq <- suppressWarnings(as.character(
            Biostrings::translate(
                Biostrings::DNAStringSet(
                    Biostrings::subseq(nt_sequence, start=frame)))))

        start_site <- stringr::str_locate_all(normal_seq, "M")

        start_sites <- plyr::ldply(start_site, data.frame)
        start_sites$index <- unlist(mapply(function(x, y) rep(y,(length(x)/2)),
                                           start_site, 1:length(start_site)))
        start_sites$orf <- normal_seq[start_sites$index]

        start_sites$orf <- str_sub(start_sites$orf, start=start_sites$start)
        start_sites$stop_site <- stringr::str_locate(start_sites$orf, "[*]")[,1]
        start_sites$orf <- str_sub(start_sites$orf, start=1, end=start_sites$stop_site-1)
        start_sites$length <- start_sites$stop_site -1

        longest <- aggregate(length ~ index, start_sites, max)
        longest$frame <- frame

        if(frame == 1){
            longest_all <- longest
        }else{
            m2 <- match(longest$index,longest_all$index)

            if(any(is.na(m2))){
                longest_all <- rbind(longest_all, longest[which(is.na(m2)),])
            }
            m1 <- match(longest_all$index, longest$index)

            replace <- which(longest$length[m1] > longest_all$length)

            if(length(replace) > 0){
                longest_all[replace,] <- longest[m1,][replace,]
            }

        }
    }
    return(longest_all)
}
