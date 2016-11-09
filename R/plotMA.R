#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

plotMA <- function(data, ...)
{
    stopifnot(class(data) == 'AnaquinData')
    
    if (analysis(data) != 'PlotMA')
    {
        stop('plotMA requires PlotMA analysis')
    }

    data <- data.frame(row.names=seqs(data),
                       ExpLFC=NA,
                       M=lfc(data),
                       A=mean(data))

    data <- data[!is.na(data$M),]
    data <- data[!is.na(data$A),]    

    # TODO: ... Need to add information for the expected fold-change ...
    data[RnaQuin.genes(row.names(data)),]$ExpLFC <- 0

    p <- ggplot(data, aes(x='A', y='M')) +
             geom_point(colour='grey80', alpha=0.5) +
             geom_point(data=subset(data, !is.na(data$ExpLFC)),
                        colour="orange",size = 2.5)   + 
             geom_hline(aes(yintercept=0), alpha=0.5) +
             theme_bw()

    print(p)
    p
}