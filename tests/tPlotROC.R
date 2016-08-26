#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.PlotROC_1 <- function()
{
    data('seqDESeq2')
    data <- seqDESeq2

    data$label <- ifelse(abs(data$ExpLFC) <= 0, 'FP', 'TP')

    data <- AnaquinData(analysis='PlotROC',
                            seqs=row.names(data),
                           ratio=data$ExpLFC,
                           score=1-data$Pval,
                           label=data$label)

    r <- plotROC(data, refRats=0, unitTest=TRUE)
    
    checkEquals(r$AUC[1,]$Ratio, 1)
    checkEqualsNumeric(r$AUC[1,]$AUC, 0.7692)
    checkEquals(r$AUC[2,]$Ratio, 2)
    checkEqualsNumeric(r$AUC[2,]$AUC, 0.8223)
    checkEquals(r$AUC[3,]$Ratio, 3)
    checkEqualsNumeric(r$AUC[3,]$AUC, 0.9485)
    checkEquals(r$AUC[4,]$Ratio, 4)
    checkEqualsNumeric(r$AUC[4,]$AUC, 0.8608)
}

test.PlotROC_2 <- function()
{
    data('seqCufflinks')
      
    data <- seqCufflinks
    data <- AnaquinData(analysis='PlotLinear',
                            seqs=row.names(data),
                           input=log2(data$InputConcent),
                        measured=log2(data$Observed1))
    
    checkException(plotROC(data, refRats=0, unitTest=TRUE))
}

test.PlotROC_1()
test.PlotROC_2()