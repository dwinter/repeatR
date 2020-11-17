#' Read RepeakMasker output
#'@export
#'@importFrom stringr str_split
#'@return A data.frame in which each row represents a (portion of) and alignment
#' between a query sequence (usually a genome) and a reference repeat.
#' Columns are
#' \itemize{
#'    \item score score The Smith-Waterman score for this alignment-section
#'    \item p_sub proportion of mismatches (substitutions)
#'    \item p_del propotion of deletions relative to target
#'    \item p_ins proportions of insertins relative to target
#'    \item qname query (usually chromosome or scaffold) name
#'    \item qstart start position in query sequence
#'    \item qend end position in query sequence
#'    \item qextend query sequence remaining after alignment
#'    \item complement is the alignment reverse and complement ("C") or same strnd ("+")
#'    \item tname target (repeat) name
#'    \item tclass target (repeat) class
#'    \item tstart start position in target sequence
#'    \item tend end position in target sequence
#'    \item textend target sequence remaining after alignment
#'    \item ID alignment ID (can be shared by multiple rows if targets interrupt each other)
#'    \item ali_type alignment type (primary for best alignment for this query region, secondary for others)
#'}

read_rm <- function(file, tibble=FALSE, keep_order=TRUE, include_secondary = FALSE) {
    #read the data in, trim leading white space and split into tokens
    lines <- readLines(file)
    lines <- str_trim(lines[4:length(lines)], "l")
    raw_data <- data.frame( str_split(lines, "\\s+", simplify=TRUE) )
    #we need to deal with complement alignments seperately as the traget
    #information appears in a different order for these alignents. Start by
    #spliting out the comp. alignments, then xhange the order of the target
    #info to match the "+" strand ali (tstart,tend,textend, which is reversed
    #for complement ali)
    by_strand <- split(raw_data, raw_data[,9])
    comp_ali <- by_strand[["C"]][c(1:11,14,13,12,15,16)]
    #rename each data.frame before binding
    col_names <-  c("score", "p_sub", "p_del", "p_ins", "qname", "qstart", 
                    "qend", "qextend", "complement", "tname", "tclass","tstart", 
                    "tend", "textend", "ID", "ali_type")
    names(comp_ali) <- col_names
    names(by_strand[["+"]]) <- col_names
    res <- rbind.data.frame(comp_ali, by_strand[["+"]])
    #stacking the data this way makes all the complement alignments show up
    #first, then all the +ve strand one. Probably usually want to go back to
    #orignal ordering:
    if( keep_order ){
        res <-  res[order(as.numeric(rownames(res))),]
    }
    #repeat masker output includes secondary alignments in a final colums with 
    # * = secondoary; nothing = primiary. Make this explicit (and give user
    # option to pre-filter these
    res$ali_type <-  ifelse(res$ali_type == "*", "secondary", "primary")
    if(!include_secondary){        
      res <- subset(res, ali_type == "primary" )
    }
    #everything is a character at the moment, and the parentheses around 
    # the "extend" is nt helpful in R
    res <- .numerify(res, c(1:4, 6,7, 13))
    res <- .de_paren(res, c("textend", "qextend"))
    if(keep_order){

    }
    class(res) <- c("repeat_table", "data.frame")
    rownames(res) <- NULL
    res
}

#'@export
print.repeat_table <- function(x, ...){
  n_rep <- length(unique(x[["ID"]]))
  cat("RepeatMasker output with ", nrow(x), " entries for ", n_rep, "unique repeat sequences\n")
  #cat( "column names:\n  ", names(x))
  print(head(as.data.frame(x),3))
}

.de_paren <- function(df, column_indices) {
  for(i in column_indices){
   df[[i]] <-  as.numeric(gsub("\\(|\\)", "", df[[i]]))
  }
  df
}


.numerify <- function(df, column_indices){
  for(i in column_indices){
    df[[i]] <- as.numeric(as.character(df[[i]]))
  }
  df
}

