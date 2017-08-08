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
    x <- UserGuideData_5.6.3
    x$label <- ifelse(abs(x$ExpLFC) <= 0, 'FP', 'TP')

    r <- plotROC(row.names(x),
                 1-x$Pval,
                 abs(x$ExpLFC),
                 x$label,
                 refGroup=0)

    checkEquals(as.character(r$AUC[1,][1]), '4')
    checkEqualsNumeric(r$AUC[1,]$AUC, 0.9062)
    checkEquals(as.character(r$AUC[2,][1]), '3')
    checkEqualsNumeric(r$AUC[2,]$AUC, 0.8939)
    checkEquals(as.character(r$AUC[3,][1]), '2')
    checkEqualsNumeric(r$AUC[3,]$AUC, 0.7955)
    checkEquals(as.character(r$AUC[4,][1]), '1')
    checkEqualsNumeric(r$AUC[4,]$AUC, 0.6713)
}

test.PlotROC_1()