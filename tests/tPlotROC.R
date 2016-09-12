#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test.PlotROC_1 <- function()
{
    data(UserGuideData_5.6.3)
    data <- UserGuideData_5.6.3

    data$label <- ifelse(abs(data$ExpLFC) <= 0, 'FP', 'TP')

    data <- AnaquinData(analysis='PlotROC',
                            seqs=row.names(data),
                           ratio=data$ExpLFC,
                           score=1-data$Pval,
                           label=data$label)

    r <- plotROC(data, refRats=0, unitTest=TRUE)

    checkEquals(as.character(r$AUC[1,]$Ratio), '1')
    checkEqualsNumeric(r$AUC[1,]$AUC, 0.6713)
    checkEquals(as.character(r$AUC[2,]$Ratio), '2')
    checkEqualsNumeric(r$AUC[2,]$AUC, 0.7955)
    checkEquals(as.character(r$AUC[3,]$Ratio), '3')
    checkEqualsNumeric(r$AUC[3,]$AUC, 0.8939)
    checkEquals(as.character(r$AUC[4,]$Ratio), '4')
    checkEqualsNumeric(r$AUC[4,]$AUC, 0.9062)
}

test.PlotROC_2 <- function()
{
    data(UserGuideData_5.4.6.3)
      
    data <- UserGuideData_5.4.6.3
    data <- AnaquinData(analysis='PlotLinear',
                            seqs=row.names(data),
                           input=log2(data$InputConcent),
                        measured=log2(data$Observed1))
    
    checkException(plotROC(data, refRats=0, unitTest=TRUE))
}

test.PlotROC_1()
test.PlotROC_2()