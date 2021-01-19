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


rm_ideogram <- function(rm_table){
    qs <- query_sizes(rm_table)
    rm_table$qname <- factor(rm_table$qname, levels=qs$qname)
    ggplot() + geom_rect(data=qs, aes(xmin=0, xmax=qsize, ymax=-0.5, ymin=0.5), colour="black", fill="white") + facet_grid(qname ~.) + 
            geom_rect(data=rm_table, aes(xmin=qstart, xmax=qend, ymax=-0.48, ymin=0.48))
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
