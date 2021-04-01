#'@export
filter_by_tclass <- function(rm_table, 
                             exclude=c("Simple_repeat", "Low_complexity", "tRNA",
                                       "ARTEFACT", "scRNA", "rRNA", "snRNA")) {
    subset(rm_table, !(tclass %in% exclude))                            
}


