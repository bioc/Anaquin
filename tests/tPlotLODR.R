#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test.PlotLODR_1 <- function()
{
    data(UserGuideData_5.6.3)
    data <- UserGuideData_5.6.3
      
    xlab  <- 'Average Counts'
    ylab  <- 'P-value'
    title <- 'LODR Curves'
      
    anaquin <- AnaquinData(analysis='PlotLODR',
                               seqs=row.names(data),
                           measured=data$Mean,
                              ratio=data$ExpLFC,
                               pval=data$Pval,
                               qval=data$Qval)

    r <- plotLODR(anaquin,
                  xlab=xlab,
                  ylab=ylab,
                  title=title,
                  FDR=0.1,
                  unitTest=TRUE)
    
    checkEquals(0, r[1,1])
    checkTrue(is.infinite(r[1,2]))
    checkEquals(1, r[2,1])
    checkEqualsNumeric(3.617511181, r[2,2])
    checkEquals(2, r[3,1])
    checkEqualsNumeric(5.673070227, r[3,2])
    checkEquals(3, r[4,1])
    checkEqualsNumeric(2.584552757, r[4,2])
    checkEquals(4, r[5,1])
    checkEqualsNumeric(3.896317229, r[5,2])
}

test.PlotLODR_2 <- function()
{
    data(UserGuideData_5.4.6.3)
      
    data <- UserGuideData_5.4.6.3
    data <- AnaquinData(analysis='PlotLinear',
                            seqs=row.names(data),
                           input=log2(data$InputConcent),
                        measured=log2(data$Observed1))
    
    checkException(plotLODR(data,
                            xlab=xlab,
                            ylab=ylab,
                            title=title,
                            FDR=0.1,
                            unitTest=TRUE))
}

test.PlotLODR_1()
test.PlotLODR_2()