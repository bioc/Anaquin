#
#  Copyright (C) 2017 - Garvan Institute of Medical Research
#

plotLinear <- function(seqs,
                       x,
                       y,
                       std=NULL,
                       title=NULL,
                       xlab=NULL,
                       ylab=NULL,
                       showSD=TRUE,
                       showLOQ=TRUE,
                       showStats=TRUE,
                       xBreaks=NULL,
                       yBreaks=NULL,
                       errors=NULL,
                       showLinear=TRUE,
                       showAxis=FALSE)
{
    stopifnot(is.null(errors) | errors == 'SD' | errors == 'Range')
    
    stopifnot(!is.null(seqs))
    stopifnot(!is.null(x))
    stopifnot(!is.null(y))
    
    if (is.factor(x))  { x <- as.numeric(as.character(x)) }
    if (is.factor(y))  { y <- as.numeric(as.character(y)) }    
    
    data <- data.frame(row.names=seqs, x=x, y=y)

    if (!is.null(std)) { data$sd <- std }
    
    if (ncol(data) == 2)
    {
        data <- data[!is.na(data$y),]
        data <- data[!is.infinite(data$y),]
    }
    
    # For mutliple measurments
    data$y <- data[,c(2:ncol(data))]
    
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
        
        data$min <- do.call(pmin, c(as.data.frame(data$y), na.rm=TRUE))
        data$max <- do.call(pmax, c(as.data.frame(data$y), na.rm=TRUE))
        data$y   <- rowMeans(data$y, na.rm=TRUE)
    }
    else
    {
        data.all <- data.frame(x=data$x, y=data$y)
    }
    
    if (isMulti)
    {
        #
        # Quote from: "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2064100"
        #
        #   About two thirds of the data points will lie within the region of
        #   mean ± 1 SD, and ∼95% of the data points will be within 2 SD of
        #   the mean.
        #
        # We want to establish 95% confidence interval.
        #
        if (is.null(errors) || errors == 'SD')
        {
            data$ymax <- data$y + 2*data$sd
            data$ymin <- data$y - 2*data$sd
        }
        
        #
        # Range error encompass the lowest and highest values
        #
        else if (errors == 'Range')
        {
            data$ymax <- data$max
            data$ymin <- data$min
        }
        
        data <- data[!is.na(data$ymax),]
        data <- data[!is.na(data$ymin),]    
    }
    else
    {
        data$sd <- NULL
    }
    
    data <- data[!is.na(data$y),]
    
    p <- ggplot(data=data, aes_string(x='data$x', y='data$y', label='row.names(data)')) +
                    xlab(xlab) +
                    ylab(ylab) +
                    ggtitle(title) +
                    labs(colour='Ratio') +
                    geom_point(aes_string(colour='grp'),
                               size=2.0,
                               alpha=0.5) +
                    guides(colour=FALSE) +
                    #ylim(min(data$y), max(data$y)) +
                    theme_bw() +
                    theme(plot.title = element_text(hjust = 0.5))

    if (showLinear)
    {
       p <- p + geom_smooth(method='lm',
                           formula=y~x,
                          linetype='11',
                             color='black',
                              size=0.5)
    }
        
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

    if (!is.null(xBreaks)) { p <- p + scale_x_continuous(breaks=xBreaks) }
    if (!is.null(yBreaks)) { p <- p + scale_y_continuous(breaks=yBreaks) }
    
    overall <- .lm2str(data)
    above   <- NULL
    
    LOQ <- NULL
    
    if (showLOQ)
    {
        tryCatch(
            {
                LOQ <- estimateLOQ(data.all$x, data.all$y)
            }, error = function(x) {})
        
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
    
    if (showStats)
    {
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
    }
    
    p <- .transformPlot(p)
    suppressWarnings(print(p))
    
    list(LOQ = LOQ)
}