#
#  Copyright (C) 2017 - Garvan Institute of Medical Research
#

#
# Limit-of-quantification (LOQ) is defined as the level of concentration where
# accurate interpreation becomes questionable.
#

estimateLOQ <- function(x, y, showDetails=FALSE)
{
    percentile <- ecdf(x)

    #
    # LOQ estimation is tricky; inputs could be anything and there's no single
    # model that can possibly model everything. However, we can consider the
    # common scenarios.
    #
    # It's well known that least-square regression is sensitive to outliers. LOQ
    # shouldn't depend on those outliers. The following techniques are used:
    #
    #    - Cook's distance
    #    - Standardized residuals
    #
    # Leverage is unnecessary because the range of the input concentration is
    # fixed.
    #

    data  <- data.frame(x=x, y=y)
    data  <- data[order(x),]
    data  <- data[!is.na(data$y),]
    
    # Compute pecentile for the input concentration
    data$perc <- percentile(data$x)
    
    # Fit a model to all data points
    model.all <- lm(y~x, data=data)

    data$CD <- cooks.distance(model.all)
    data$RS <- rstudent(model.all)
    
    #
    # Attempt 1: Remove outliers after LOQ
    #
    
    data <- data[data$perc < 0.40 | (data$CD < 4/nrow(data) & abs(data$RS) <= 3),]

    #
    # Attempt 2: Remove outliers in the before LOQ. The point of LOQ is to use
    #            stochastic points to estimate the limit, so we should be more
    #            conservative.
    #
    
    data <- data[(data$CD < 8/nrow(data)) & (abs(data$RS) <= 6),]
    
    # Fit a new model withour the outliers
    model <- lm(y~x, data=data)

    #
    # For x=1,2,3,...,n, we would only fit b=3,4,...,n-2.
    # Therefore the length of the frame is n-4.
    #

    r <- data.frame(k=rep(NA,length(x)),
                    sums=rep(NA,length(x)),
                    lR2=rep(NA,length(x)),
                    lSlope=rep(NA,length(x)),                    
                    lInter=rep(NA,length(x)),
                    rR2=rep(NA,length(x)),
                    rSlope=rep(NA,length(x)),                    
                    rInter=rep(NA,length(x)),
                    lr=rep(NA,length(x)),
                    rr=rep(NA,length(x)))

    plm <- function(i)
    {
        #
        # Eg: (1,2,3,4) and i==2
        #
        #   -> (1,2) and (3,4).
        #
        
        d1 <- head(data,i)
        d2 <- tail(data,-i)
        
        stopifnot(nrow(d1) >= 2)
        stopifnot(nrow(d2) >= 2)
        stopifnot(nrow(d1)+nrow(d2) == nrow(data))
        
        m1 <- lm(y~x, data=d1)
        m2 <- lm(y~x, data=d2)
        
        return (list(breaks=data[i,]$x,
                   'lModel'=m1,
                   'rModel'=m2,
                       'lr'=cor(d1$x,d1$y),
                       'rr'=cor(d2$x,d2$y)))
    }

    lapply(2:(nrow(data)-3), function(i)
    {
        options(warn=-1)
        
        # Fit two piecewise linear models
        fit <- plm(i)
        
        options(warn=0)

        # Where this breakpoint occurs
        r$k[i] <<- fit$breaks

        m1 <- fit$lModel # The left regression model
        m2 <- fit$rModel # The right regression model

        sm1 <- summary(m1)
        sm2 <- summary(m2)
        
        r$lr[i]  <<- fit$lr
        r$rr[i]  <<- fit$rr
        r$lR2[i] <<- sm1$r.squared
        r$rR2[i] <<- sm2$r.squared

        if (r$lR2[i] != 0)
        {
            r$lInter[i] <<- sm1$coefficients[1,1]
            r$lSlope[i] <<- sm1$coefficients[2,1]
        }
        
        if (r$rR2[i] != 0)
        {
            r$rInter[i] <<- sm2$coefficients[1,1]
            r$rSlope[i] <<- sm2$coefficients[2,1]
        }

        # Calculate sum of residuals        
        r$sums[i] <<- sum((m1$residuals)^2) + sum((m2$residuals)^2)
    })
    
    r <- r[!is.na(r$k),]
    
    #
    # Construct a plot of breaks vs total SSE
    #
    
    if (showDetails)
    {
        p <- ggplot(data = r, aes_string(x='k', y='sums'))
        p <- p + xlab('Break point')
        p <- p + ylab('Total sum of squares')
        p <- p + geom_line()
        print(p)
    }
    
    #
    # How to define the LOQ breakpoint? There're several measures:
    #
    #   - Quartile
    #   - Total SSE
    #   - Pearson's correlation
    #

    #
    # LOQ is defined be in lower quartiles (eg: we don't want it near
    # the highly expressed genes).
    #
    
    tmp <- r[percentile(r$k) <= 0.40,]
    
    if (nrow(tmp) == 0)
    {
        return (NULL)
    }
    
    # Total SSE
    b1 <- tmp[which.min(tmp$sums),]    

    # Pearson's correlation
    b2 <- tmp[which.max(tmp$rR2),]

    b1q <- percentile(b1$k) 
    b2q <- percentile(b2$k) 

    isReasonableBreak <- function(b)
    {
        if (is.null(b)) { return (NULL) }
        
        #
        # Absolute detection limit. The point of LOQ is to estimate the empirical
        # stochastic detection limit. No point if it's also the absolute limit.
        #

        if (b$k == min(tmp$k)) { return (NULL) }
        
        # Fit the model again
        fit <- plm(as.numeric(row.names(b)))
        
        #
        # Note that the correlation is computed on the filtered data set (not
        # including residuals).
        #
        
        if (fit$rr > cor(data$x, data$y))
        {
            return (fit)
        }
        else
        {
            return (NULL)
        }
    }

    fit <- NULL
    
    if (!is.null(fit <- isReasonableBreak(b2)))
    {
        b <- b2
    }
    else if (!is.null(fit <- isReasonableBreak(b1)))
    {
        b <- b1
    }
    
    if (is.null(b) | is.null(fit))
    {
        return (NULL)
    }

    stopifnot(all.equal(summary(fit$lModel)$r.squared, b$lR2))
    stopifnot(all.equal(summary(fit$rModel)$r.squared, b$rR2))
    stopifnot(all.equal(summary(fit$lModel)$coefficients[1,1], b$lInter))    
    stopifnot(all.equal(summary(fit$rModel)$coefficients[1,1], b$rInter))    
    stopifnot(all.equal(summary(fit$lModel)$coefficients[2,1], b$lSlope))    
    stopifnot(all.equal(summary(fit$rModel)$coefficients[2,1], b$rSlope))    
    
    return (list('model'=fit, 'breaks'=b, 'details'=r))
}