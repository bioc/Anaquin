#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute
#

plotLinear <- function(data,
                       title=NULL,
                        xlab=NULL,
                        ylab=NULL,
                      showSD=FALSE,
                     showLOQ=TRUE,
                    showAxis=FALSE, ...)
{
    stopifnot(class(data) == 'AnaquinData')
    
    if (analysis(data) != 'PlotLinear')
    {
        stop('plotLinear requires PlotLinear analysis')
    }

    stopifnot(!is.null(seqs(data)))
    stopifnot(!is.null(input(data)))
    stopifnot(!is.null(measured(data)))
    
    z <- list(...)
    
    if (is.null(z$unitTest)) { z$unitTest <- FALSE }

    tmp <- data.frame(row.names=seqs(data),
                          input=input(data),
                       measured=measured(data))

    if (!is.null(std(data)))
    {
        tmp$sd <- std(data)
    }
    
    # Convert AnaquinData to data.frame
    data <- tmp
    
    if (ncol(data) == 2)
    {
        data <- data[!is.na(data$measured),]
        data <- data[!is.infinite(data$measured),]
    }
    
    data$y <- data[,c(2:ncol(data))]
    data$x <- data$input
    data   <- data[,c(-1,-2)]
    
    data <- data[!is.na(data$x),]
    data$grp <- as.factor(abs(data$x))
        
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) ||
              length(data$x) == nrow((data$y)))
    
    isMultiDF <- is(data$y, 'data.frame')
    
    # Should we show standard deviation?
    isMultiSD <- sum(data$sd) > 0 & showSD

    isMulti <- isMultiDF | isMultiSD
    
    data$ymax <- NULL
    data$ymin <- NULL

    if (isMultiDF)
    {
        if (is.null(data$sd))
        {
            data$sd <- apply(data$y, 1, sd, na.rm=TRUE)
        }
        
        data$y  <- rowMeans(data$y, na.rm=TRUE)
    }
    
    if (isMulti)
    {
        data$ymax <- data$y + data$sd
        data$ymin <- data$y - data$sd
        data <- data[!is.na(data$ymax),]
        data <- data[!is.na(data$ymin),]        
    }
    else
    {
        data$sd <- NULL
    }
    
    data <- data[!is.na(data$y),]
    
    p <- ggplot(data=data, aes_string(x='data$x', y='data$y')) +
                               xlab(xlab) +
                               ylab(ylab) +
                           ggtitle(title) +
                     labs(colour='Ratio') +
                     geom_point(aes_string(colour='grp'),
                                size=2.0,
                               alpha=0.5) +
                geom_smooth(method='lm',
                           formula=y~x,
                          linetype='11',
                             color='black',
                              size=0.5)  +
                theme_bw()

    p <-p + guides(colour=FALSE)
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.7, 0.7)

    if (showAxis)
    {
        p <- p + geom_vline(xintercept=c(0), linetype='solid', size=0.1)
        p <- p + geom_hline(yintercept=c(0), linetype='solid', size=0.1)
    }
    
    overall <- .lm2str(data)
    above   <- NULL

    LOQ <- NULL

    if (showLOQ)
    {
        tryCatch({
            LOQ <- estimateLOQ(data$x, data$y)
        }, error = function(cond)
        {
        })
        
        if (!is.null(LOQ))
        {
            if (LOQ$model$rr > cor(data$x, data$y))
            {
                # Print out the regression above LOQ
                above <- .m2str(LOQ$model$rModel)

                # Assume the break-point is on the log2-scale. Convert it back.
                label <- 2^LOQ$breaks$k
                
                t <- paste('LOQ:', signif(label, 3))
                t <- paste(t, 'attomol/ul')
                
                p <- p + geom_vline(xintercept=c(LOQ$breaks$k),
                                      linetype='33',
                                          size=0.6)
                p <- p + geom_label(aes_string(x='max(LOQ$breaks$k)',
                                               y='min(y)'),
                                           label=t,
                                          colour='black',
                                     show.legend=FALSE,
                                           hjust=0.1,
                                           vjust=0.7)                
            }
            else
            {
                LOQ <- NULL
            }
        }
    }
    
    r <- abs(max(data$y) - min(data$y))
    y_off <- 0.06 * r 

    if (showLOQ)
    {
        a <- paste(c('bold(Overall): ', overall), collapse='')
    }
    else
    {
        a <- overall
    }

    overall <- annotate("text",
                        label=a,
                        x=min(data$x),
                        y=max(data$y)-y_off,
                        size=4.0,
                        colour='grey24',
                        parse=TRUE,
                        hjust=0,
                        vjust=0)
    
    p <- p + overall
    
    if (showLOQ & !is.null(LOQ))
    {
        above <- annotate("text",
                          label=paste(c('bold(Above)~bold(LOQ): ',
                                        above), collapse=''),
                          x=min(data$x),
                          y=max(data$y)-2*y_off,
                          size=4.0,
                          colour='grey24',
                          parse=TRUE,
                          hjust=0,
                          vjust=0)
        p <- p + above
    }

    if (!is.null(data$sd))
    {
        p <- p + geom_errorbar(aes_string(ymax='ymax',
                                          ymin='ymin'),
                                          size=0.2,
                                         alpha=0.5)
    }

    p <- .transformPlot(p)
    print(p)
    
    if (z$unitTest)
    {
        return (LOQ)
    }
}