#' Read RepeakMasker output
#'@export
#'@return A data.frame in which each row represents a (portion of) and alignment
#' between a query sequence (usually a genome) and a reference repeat.
#' Columns are
#' \itemize{
#'    \item \textbf{score} score The Smith-Waterman score for this alignment-section
#'    \item \textbf{p_sub} proportion of mismatches (substitutions)
#'    \item \textbf{p_del} propotion of deletions relative to target
#'    \item \textbf{p_ins} proportions of insertins relative to target
#'    \item \textbf{qname} query (usually chromosome or scaffold) name
#'    \item \textbf{qstart} start position in query sequence
#'    \item \textbf{qend} end position in query sequence
#'    \item \textbf{qextend} query sequence remaining after alignment
#'    \item \textbf{complement} is the alignment reverse and complement ("C") or same strnd ("+")
#'    \item \textbf{tname} target (repeat) name
#'    \item \textbf{tclass} target (repeat) class
#'    \item \textbf{tstart} start position in target sequence
#'    \item \textbf{tend} end position in target sequence
#'    \item \textbf{textend} target sequence remaining after alignment
#'    \item \textbf{ID} alignment ID (can be shared by multiple rows if targets interrupt each other)
#'    \item \textbf{ali_type} alignment type (primary for best alignment for this query region, secondary for others)
#'}
read_rm <- function(file, tibble=TRUE, include_secondary = FALSE, quiet=FALSE){
    lines <- readLines(file)
    lines <- sapply(lines[4:length(lines)], trimws)
    parsed <- strsplit(lines[3:length(lines)], "\\s+")
    res <- do.call(rbind.data.frame, lapply(parsed, .process_row))
    names(res) <- c("score", "p_sub", "p_del", "p_ins", "qname", "qstart", "qend", "qextend", "complement", "tname", "tclass","tstart", "tend", "textend", "ID", "ali_type")
    if(!include_secondary){
      res <- subset(res, ali_type == "primary" )
      res$ID <- droplevels(res$ID)
    }
    res <- .numerify(res, c(1:4, 6,7))
    class(res) <- c("repeat_table", "data.frame")
    rownames(res) <- NULL
    res
}

print.repeat_table <- function(x, ...){
  cat("RepeatMasker output with ", nrow(x), " entries for ", nlevels(x$ID), "unique repeat sequences\n")
  #cat( "column names:\n  ", names(x))
  print(head(as.data.frame(x),3))
}

.process_row <- function(row){
  # for these two rows, the parentheses provide no extra information than the
  # complement field, so remove them from all cases
  row[[8]] <- .de_paren(row[[8]])
  res <- as.list(row[1:11])
  if( row[[9]] == "C"){
    res[[12]] <- as.numeric(row[[13]])
    res[[13]] <- .de_paren(row[[12]])
  } else {
    res[[12]] <- as.numeric(row[[12]])
    res[[13]] <- as.numeric(row[[13]])
  }
  res[[14]] <- .de_paren(row[[14]])
  res[[15]] <- row[[15]]
  res[[16]] <- if(length(row) == 15) "primary" else "secondary"
  res
}

.de_paren <- function(n) {
  as.numeric(gsub("\\(|\\)", "", n))
}


.numerify <- function(df, column_indices){
  for(i in column_indices){
    df[[i]] <- as.numeric(as.character(df[[i]]))
  }
  df
}

