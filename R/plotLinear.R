#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

plotLinear <- function(data,
                       title=NULL,
                        xlab=NULL,
                        ylab=NULL,
                      showSD=TRUE,
                     showLOQ=TRUE,
                     xBreaks=NULL,
                     yBreaks=NULL,
                    showAxis=FALSE, ...)
{
    stopifnot(class(data) == 'AnaquinData')
    
    if (analysis(data) != 'PlotLinear' &
        analysis(data) != 'plotLinear')
    {
        print(analysis(data))
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
    data <- rename(data, c('input'='x'))
    
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

    # All data points (including all replicates)
    data.all <- NULL
    
    if (isMultiDF)
    {
        for (i in 1:ncol(data$y))
        {
            data.all <- rbind(data.all, data.frame(x=data$x, y=(data$y)[,i]))
        }
        
        if (is.null(data$sd))
        {
            data$sd <- apply(data$y, 1, sd, na.rm=TRUE)
        }
        
        data$y  <- rowMeans(data$y, na.rm=TRUE)
    }
    else
    {
        data.all <- data.frame(x=data$x, y=data$y)
    }

    if (isMulti)
    {
        #
        # There're different types of error bars; SD, SE and CI. Here, we're more intersted in
        # the variation in the data, so we'll implement the SD method.
        #
        
        data$ymax <- data$y + 2*data$sd
        data$ymin <- data$y - 2*data$sd
        
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
                    guides(colour=FALSE) +
                theme_bw()

    if (showSD & !is.null(data$sd))
    {
        p <- p + geom_segment(aes_string(x='data$x',
                                         y='data$ymax',
                                      xend='data$x',
                                      yend='data$ymin'),
                              data=data,
                              size=0.2,
                              alpha=0.5)
    }

    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.7, 0.7)

    if (showAxis)
    {
        p <- p + geom_vline(xintercept=c(0), linetype='solid', size=0.1)
        p <- p + geom_hline(yintercept=c(0), linetype='solid', size=0.1)
    }
    
    if (!is.null(xBreaks))
    {
        p <- p + scale_x_continuous(breaks=xBreaks)
    }

    if (!is.null(yBreaks))
    {
        p <- p + scale_y_continuous(breaks=yBreaks)
    }
    
    overall <- .lm2str(data)
    above   <- NULL

    LOQ <- NULL

    if (showLOQ)
    {
        tryCatch(
        {
            LOQ <- estimateLOQ(data.all$x, data.all$y)
        }, error = function(cond) {})
        
        if (!is.null(LOQ))
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
    }
    
    if (showLOQ & !is.null(LOQ))
    {
        overall <- paste(c('bold(Overall): ', overall), collapse='')
        above   <- paste(c('bold(Above)~bold(LOQ):', above), collapse='')
        label   <- paste(c('atop(', overall, ',', above, ')'), collapse='')
        
        p <- p + annotate("text",
                          label=label,
                          x=min(data$x),
                          y=max(data$y),
                          size=4.0,
                          colour='grey24',
                          parse=TRUE,
                          hjust=0,
                          vjust=1)
    }
    else
    {
        p <- p + annotate("text",
                          label=overall,
                          x=min(data$x),
                          y=max(data$y),
                          size=4.0,
                          colour='grey24',
                          parse=TRUE,
                          hjust=0,
                          vjust=1)
    }

    p <- .transformPlot(p)
    print(p)
    
    if (z$unitTest)
    {
        return (LOQ)
    }
}