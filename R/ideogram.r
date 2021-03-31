query_sizes <- function(rm_table){
    one_per <- rm_table[ !duplicated(rm_table[["qname"]]),]
    res <- data.frame(qname = one_per[["qname"]], 
                      qsize = one_per[["qextend"]] + one_per[["qend"]]
    )
    res$qname <- reorder(res$qname, res$qsize)
    res

}

order_queries <- function(data){
    data[["qname"]] <- reorder(data[["qname"]], data[["qsize"]])
    data
}

#' plot an ideogram (or genomic regoin) of TE locations
#' @param rm_table A repeat masker output table,
#' @seealso read_rm
rm_ideogram <- function(rm_table){
    qs <- query_sizes(rm_table)
    rm_table$qname <- factor(rm_table$qname, levels=qs$qname)
    ggplot() + geom_rect(data=qs, aes(xmin=0, xmax=qsize, ymax=-0.5, ymin=0.5), colour="black", fill="white") +
        facet_grid(qname ~.) + 
        geom_rect(data=rm_table, aes(xmin=qstart, xmax=qend, ymax=-0.48, ymin=0.48)) +
        theme_ideogram()
}

theme_ideogram <- function() {
    theme <- theme_minimal()
    theme$axis.title.y <- element_blank()
    theme$axis.text.y <- element_blank()
    theme$axis.ticks.y <- element_blank()
    theme$axis.ticks.y <- element_blank()
    theme$strip.background <- element_blank()
    theme$strip.text <- element_blank()
    theme$panel.grid.major = element_blank()
    theme$panel.grid.minor = element_blank()
    theme$panel.grid.minor = element_blank()
    theme$legend.position <- "none"    
    theme

}

#' Number formatters for scales in base pairs
#'
#' For use with \code{\link{ggplot2}} 
#' @param x  The data (in base pairs) to be formatted as Gb, Mb or Kb
#' @return  A character vector with scale labels
#' @examples
#' \dontrun{
#' ali <- read_paf(system.file("extdata", "fungi.paf", package="pafr"))
#' doplot(ali) + scale_x_continuous("Genomic position", label=Mb_lab)
#'}
#' @export
#' @rdname Gb_lab
Gb_lab <- function(x) paste(x / 1e9, "Gb")

#' @export
#' @rdname Gb_lab
Mb_lab <- function(x) paste(x / 1e6, "Mb")

#' @export
#' @rdname Gb_lab
Kb_lab <- function(x) paste(x / 1e3, "Kb")


#' Make a pallette for TE plots
#'@export 
make_TE_pallete <- function(rm_table){
    tc <- unique(rm_table$tclass)
    higher_order <- sapply(strsplit(tc, "/"), "[[", 1)
    #TE_pal in an internal data object shorted in R/sysdata.rda
    combo <- merge(data.frame(tc, type=higher_order), TE_pal, all.x=TRUE)[,2:3]  
    structure(combo$colour, .Names=combo$tc)
}




