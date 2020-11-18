#' Summarise a repeatmasker output by unique element
#'@export 
#'@importFrom dplyr group_by
#'@importFrom dplyr summarise
#'@importFrom dplyr n
summarise_rm_ID <- function(rm_table){
    grouped_rm_table <- group_by(rm_table, ID)
    res <- summarise(grouped_rm_table, 
            score = sum(score),
            p_sub = mean_psub(p_sub, p_del, qstart, qend), 
            p_del =  weighted.mean(p_del, qend-qstart),
            p_ins =  weighted.mean(p_ins, qend-qstart),
            qname = unique(qname),
            qstart = min(qstart), 
            qend=min(qend), 
            qlen = sum(qend + 1- qstart), 
            qextend=min(qextend),
            n_aligned_segments = n(),
            complement=unique(complement),
            tname = unique(tname),
            tclass = unique(tclass),
            tstart = min(tstart), 
            tend=min(tend), 
            tlen = sum(tend + 1- qstart),
            ali_type=unique(ali_type)
    )
    class(res) <- c("rm_summ_id", "tbl_df", "tbl", "data.frame")
    res 

}

#' Summarise a repeatmasker output by family
#'@export 
#'@importFrom dplyr group_by
#'@importFrom dplyr summarise

summarise_rm_family <- function(id_table){
    .summarise_rm(group_by(id_table, tname))
}

#' Summarise a repeatmasker output by TE class
#'@export 

summarise_rm_class <- function(id_table){
    .summarise_rm(group_by(id_table, tclass))
}
    
    
.summarise_rm <- function(grouped_rm_table){
    summarise(grouped_rm_table,
      tclass = unique(tclass),
      n_elements = n(),
      aligned_len = sum(qlen),
      mean_p_sub = mean(p_sub),
      mean_p_del = mean(p_del),
      mean_p_ins = mean(p_ins),
      mean_semgents = mean(n_aligned_segments),
      gt_1_segments = sum(n_aligned_segments > 1)
      )
}

mean_psub <- function(p_sub, p_del, qstart, qend){
    total_qlen <- qend - qstart
    # some query is not aligned, so not part of psub
    # for average
    a_qlen <- total_qlen - (p_del/ 100) * total_qlen
    weighted.mean(p_sub, a_qlen)
}
