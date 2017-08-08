#
#  Copyright (C) 2017 - Garvan Institute of Medical Research
#

plotConjoint <- function(seqs,
                         units,
                         x,
                         y,
                         title=NULL,
                         xlab=NULL,
                         ylab=NULL)
{
    data <- data.frame(seqs=seqs, units=units, x=x, y=y)

    p <- ggplot(data=data, aes_string(x='data$x', y='data$y')) +
                    xlab(xlab) +
                    ylab(ylab) +
                    ggtitle(title) +
                    geom_point(aes_string(colour='seqs'),
                    size=2.0,
                    alpha=0.5) +
                    guides(colour=FALSE) +
                    ylim(min(data$y), max(data$y)) +
                    geom_smooth(aes_string(colour='seqs'), method='lm', se=FALSE, 
                                linetype='11',
                                size=0.5) +
                    theme_bw() +
                    theme(plot.title = element_text(hjust = 0.5))
        
    p <- .transformPlot(p)
    suppressWarnings(print(p))
}