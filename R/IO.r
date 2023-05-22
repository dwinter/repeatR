#' Read RepeakMasker output
#'@export
#'@importFrom stringr str_split
#'@importFrom stringr str_trim
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
    lines <- str_trim(lines[4:length(lines)], "left")
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
    res <- .numerify(res, c(1:4, 6,7, 12:13))
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


#' Convert RM table to BED format
#'
#' @export 
#' @param rm_table a repeat masker outout table, as returend by read_rm
#' @param score which column to use for the bed "score" attribute (one of p_sub
#' (default),  p_del, p_ins or score). 
#' @return A 6 column with one row per alignmened segment (i.e. possibly
#' multiple rows per TE)

rm_to_bed <- function(rm_table, name=c("tname", "tclass", "ID"), 
                      score=c("p_sub", "p_del", "p_ins", "score")                      ){
    name_to_get <- match.arg(name)
    score_to_get <- match.arg(score)
    res <- rm_table[,c("qname", "qstart", "qend", name_to_get, score_to_get)]
    names(res)[1:5] <- c("chrom", "start", "end", "name", "score")
    res$strand <- ifelse(rm_table[["complement"]] == "+", "+", "-")
    res$start <- res$start - 1
    class(res) <- "data.frame"
    res
}

#' Write bed formatted data to file
#'
#'@param bed The data.frame in bed format (note, this function does not check
#'the data is a valid BED file)
#'@param filename name of the file to write to
#'@export
write_bed <- function(bed, filename){
    orig_scipen <- options("scipen")[[1]]
    on.exit(options(scipen = orig_scipen))
    write.table(bed, filename, quote = FALSE, 
                row.names = FALSE, col.names = FALSE, sep = "\t")
}

commafy <- function(n){
    paste(n, collapse=",")
}

#internal fxn used by rm_to_bed12

alignment_to_bed12 <- function(ali, score_to_get) {
    #TODO:
    # If statement to treat n=1 alignments seperately?

    nseg <- nrow(ali)
    total_qstart <- min(ali[["qstart"]])
    total_qend <- max(ali[["qend"]])
    seg_lens <- ali[["qend"]] - ali[["qstart"]] + 1
    block_starts <-  ali[["qstart"]] - total_qstart
    strand <- ifelse(ali[["complement"]][1] == "+", "+", "-")
    ID <- paste0(ali[["tname"]][1], "_", ali[["ID"]][1])
    res <- data.frame(
      qname = ali[["qname"]][1],
      qstart = total_qstart,
      qend = total_qend,
      tname = ID, 
      score = mean(ali[[score_to_get]]), 
      strand = strand,
      thickStart = total_qstart, #Next 3 are graphical things, ignored in this case
      thickEnd = total_qend,
      itemRgb = "0,0,0",
      blockCount = commafy(nseg),
      blockSizes = commafy(seg_lens),
      blockStarts = commafy(block_starts))
    res
}



#' Convert an rm_table to bed12 format
#'
#' @export 
#' @param rm_table a repeat masker outout table, as returend by read_rm
#' @param score which column to use for the bed "score" attribute (one of p_sub
#' (default),  p_del, p_ins or score). 
#' @return A 12 column with one row per alignment (even for split/nested
#' repeats)



rm_to_bed12 <- function(rm_table,
                        score=c("p_sub", "p_del", "p_ins", "score")){
    score_to_get <- match.arg(score)
    by_ID <- split(rm_table, rm_table$ID)
    do.call(rbind.data.frame, lapply(by_ID, alignment_to_bed12, score_to_get))

}
