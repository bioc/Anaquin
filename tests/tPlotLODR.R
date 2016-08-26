#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotLODR_1 <- function()
{
    data('seqDESeq2')
      
    xlab  <- 'Average Counts'
    ylab  <- 'P-value'
    title <- 'LODR Curves'
      
    anaquin <- AnaquinData(analysis='PlotLODR',
                               seqs=row.names(seqDESeq2),
                           measured=seqDESeq2$Mean,
                              ratio=seqDESeq2$ExpLFC,
                               pval=seqDESeq2$Pval)

    r <- plotLODR(anaquin,
                  xlab=xlab,
                  ylab=ylab,
                  title=title,
                  FDR=0.1,
                  unitTest=TRUE)
    
    checkEquals(0, r[1,1])
    checkTrue(is.infinite(r[1,2]))
    checkEquals(1, r[2,1])
    checkEqualsNumeric(15.446922887, r[2,2])
    checkEquals(2, r[3,1])
    checkEqualsNumeric(5.315685411, r[3,2])
    checkEquals(3, r[4,1])
    checkEqualsNumeric(1.472504127, r[4,2])
    checkEquals(4, r[5,1])
    checkEqualsNumeric(6.890585670, r[5,2])
}

test.PlotLODR_2 <- function()
{
    data('seqCufflinks')
      
    data <- seqCufflinks
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