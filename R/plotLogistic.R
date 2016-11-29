#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Garvan Institute of Medical Research
#

plotLogistic <- function(data,
                         title=NULL,
                          xlab=NULL,
                          ylab=NULL,
                       showLOA=TRUE,
                     threshold=0.70, ...)
{
    stopifnot(class(data) == 'AnaquinData')

    if (analysis(data) != 'PlotLogistic' &
        analysis(data) != 'plotLogistic')
    {
        stop('plotLogistic requires PlotLogistic analysis')
    }
      
    stopifnot(!is.null(seqs(data)))
    stopifnot(!is.null(input(data)))
    stopifnot(!is.null(measured(data)))

    z <- list(...)
    
    if (is.null(z$unitTest)) { z$unitTest <- FALSE }

    data <- data.frame(x=input(data),
                       y=measured(data),
                       f=NA,
                       grp=as.factor(round(abs(input(data)))))
    
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) ||
              length(data$x) == nrow((data$y)))
    
    result = tryCatch(
    {
        sigmoid = function(params, x) {
            1 / (1 + exp(-params[1] * (x - params[2])))
        }
        
        x <- data$x
        y <- data$y
        
        fit <- nls(y~1/(1 + exp(-b * (x-c))), start=list(b=1,c=0))
        
        data$f <- sigmoid(coef(fit), x)
    }, error = function(e) {
        showLOA <<- FALSE
    })

    p <- ggplot(data=data, aes_string('x')) +
                        xlab(xlab) +
                        ylab(ylab) +
                    ggtitle(title) +
                        theme_bw() +
                        theme(plot.title = element_text(hjust = 0.5))

    p <- p + geom_point(aes_string(y='y', colour='grp'), size=2.0, alpha=0.5)

    if (!all(is.na(data$f)))
    {
        p <- p + geom_line(aes_string(y='f', colour='"line"'), alpha=1.0)
    }
    
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))

    # Limit of assembly
    LOA <- NULL
   
    if (showLOA)
    {
        t <- data[data$f >= threshold,]
        
        if (nrow(t) == 0)
        {
            warning (paste(c('Failed to estimate LOA.
                              The maximum sensitivity is: ', max(data$f))))
        }
        else
        {
            LOA <- round(min(t$x),2)

            label <- 2^LOA
            label <- paste('LOA:', signif(label, 3))
            label <- paste(label, 'attomol/ul')
            
            p <- p + geom_vline(xintercept=LOA, linetype='33')
            p <- p + geom_label(aes(x=min(data$x),
                                    y=max(data$y)-0.1),
                                label=label,
                               colour='black',
                                vjust='top',
                                hjust='left',
                          show.legend=FALSE)
        }
    }

    p <- p + guides(colour=FALSE)
    
    p <- .transformPlot(p)
    print(p)

    if (z$unitTest)
    {
        return (list(LOA=LOA, fitted=data$f))
    }
}