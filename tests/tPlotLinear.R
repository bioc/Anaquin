#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotLinear_1 <- function()
{
    data('seqCufflinks')
    
    data <- seqCufflinks
    data <- AnaquinData(analysis='PlotLinear',
                            seqs=row.names(data),
                           input=log2(data$InputConcent),
                        measured=log2(data$Observed1))
    
    r <- plotLinear(data, unitTest=TRUE)
    
    checkEqualsNumeric(r$model$breaks, 1.917068888)
}

test.PlotLinear_2 <- function()
{
    data('seqCufflinks')
    
    data <- seqCufflinks
    data <- AnaquinData(analysis='PlotLinear',
                            seqs=row.names(data),
                           input=log2(data$InputConcent),
                        measured=log2(data[,c(2:4)]))

    r <- plotLinear(data, unitTest=TRUE)
    
    checkEqualsNumeric(r$model$breaks, 1.917068888)
}

test.PlotLinear_3 <- function()
{
    data('seqDESeq2')
      
    data <- AnaquinData(analysis='PlotLODR',
                            seqs=row.names(seqDESeq2),
                        measured=seqDESeq2$Mean,
                           ratio=seqDESeq2$ExpLFC,
                            pval=seqDESeq2$Pval)
  
    checkException(plotLinear(data, unitTest=TRUE))
}

test.PlotLinear_1()
test.PlotLinear_2()
test.PlotLinear_3()