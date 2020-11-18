
#'Summarize alignments
#'@export
#'@return result
#'
tidy_rm <- function(rep_table) {
}

#'@export

filter_tclass <- function(rep_table, 
                          to_drop = c("Low_complexity", "Simple_repeat")) {    
    subset(rep_table, !(tclass  %in% to_drop))
}
