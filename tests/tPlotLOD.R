#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test.PlotLOD_1 <- function()
{
    data(UserGuideData_5.6.3)
    data <- UserGuideData_5.6.3
      
    xlab  <- 'Average Counts'
    ylab  <- 'P-value'
    title <- 'LODR Curves'
      
    r <- plotLOD(measured = data$Mean,
                    ratio = data$ExpLFC,
                     pval = data$Pval,
                     qval = data$Qval,
                     xlab = xlab,
                     ylab = ylab,
                    title = title,
                      FDR = 0.1)
    
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

test.PlotLOD_1()