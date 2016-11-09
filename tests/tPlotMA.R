#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(DESeq2)
library(Anaquin)

test.PlotMA.1 <- function()
{
    data('UserGuideData_5.6.3')
    data <- AnaquinData(analysis='PlotMA',
                        seqs=row.names(UserGuideData_5.6.3),
                        mean=log2(UserGuideData_5.6.3$Mean),
                        lfc=UserGuideData_5.6.3$ObsLFC)
    r <- plotMA(data)
    
    checkTrue(inherits(r, 'ggplot'))
}

test.PlotMA.1()