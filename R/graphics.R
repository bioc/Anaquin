#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

.transformPlot <- function(p)
{
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))
    p <- p + theme(panel.border=element_rect(colour="black", fill=NA, size=1))
    p <- p + theme(plot.title=element_text(size=15, face="bold"))
    p <- p + theme(legend.title=element_text(face="bold"))
    p <- p + theme(legend.key=element_blank())

    #
    # Calculating the aspect ratio to make the plot a square.
    #
    
    build <- ggplot_build(p)
    
    minX <- build$layout$panel_ranges[[1]]$x.range[[1]]
    maxX <- build$layout$panel_ranges[[1]]$x.range[[2]]
    minY <- build$layout$panel_ranges[[1]]$y.range[[1]]
    maxY <- build$layout$panel_ranges[[1]]$y.range[[2]]
    
    stopifnot(!is.null(minX))
    stopifnot(!is.null(maxX))
    stopifnot(!is.null(minY))
    stopifnot(!is.null(maxY))    
    
    xrange <- maxX - minX
    yrange <- maxY - minY

    p <- p + coord_fixed(ratio=xrange/yrange)
    
    return (p)    
}