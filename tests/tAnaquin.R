#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.AnaquinData_PlotLinear_1 <- function()
{
    data('seqCufflinks')
      
    data <- seqCufflinks      
    anaquin <- new("AnaquinData", analysis = 'PlotLinear',
                                      seqs = row.names(data),
                                     input = log2(data$InputConcent),
                                  measured = log2(data[,c(2:4)]))
 
    checkEquals('PlotLinear', analysis(anaquin))
    checkEquals(row.names(data), seqs(anaquin))
    checkEquals(log2(data[,c(2:4)]), measured(anaquin))
    checkEquals(log2(data$InputConcent), input(anaquin))    
}

test.AnaquinData_PlotLinear_2 <- function()
{
    data('seqCufflinks')
      
    data <- seqCufflinks
    checkException(new("AnaquinData", analysis = 'PlotLinear??',
                                          seqs = row.names(data),
                                         input = log2(data$InputConcent),
                                      measured = log2(data[,c(2:4)])))
}

test.AnaquinData_PlotROC_1 <- function()
{
    data('seqDESeq2')

    data <- seqDESeq2
    data$label <- ifelse(abs(data$ExpLFC) <= 0, 'FP', 'TP')

    anaquin <- new("AnaquinData", analysis = 'PlotROC',
                                      seqs = row.names(data),
                                     ratio = data$ExpLFC,
                                     score = 1-data$Pval,
                                     label = data$label)

    checkEquals('PlotROC', analysis(anaquin))
    checkEquals(data$ExpLFC, ratio(anaquin))
    checkEquals(data$label, label(anaquin))      
    checkEquals(1-data$Pval, score(anaquin))
    checkEquals(row.names(data), seqs(anaquin))
}

test.AnaquinData_PlotLogistic_1 <- function()
{
    data('seqCuffcompare')

    data <- seqCuffcompare
    anaquin <- new("AnaquinData", analysis = 'PlotLogistic',
                                      seqs = row.names(data),
                                     input = log2(data$InputConcent),
                                  measured = data$Sn)
      
    checkEquals('PlotLogistic', analysis(anaquin))
    checkEquals(row.names(data), seqs(anaquin))
    checkEquals(log2(data$InputConcent), input(anaquin))
    checkEquals(data$Sn, measured(anaquin))
}

test.AnaquinData_PlotLogistic_2 <- function()
{
      data('seqCuffcompare')

      data <- seqCuffcompare
      checkException(new("AnaquinData", analysis = 'PlotLogistic',
                                            seqs = row.names(data),
                                           input = log2(data$InputConcent)))
}

test.AnaquinData_PlotLODR_1 <- function()
{
    data('seqDESeq2')

    data <- seqDESeq2
    anaquin <- new("AnaquinData", analysis = 'PlotLODR',
                                      seqs = row.names(data),
                                     ratio = data$ExpLFC,
                                  measured = data$Mean,
                                      pval = data$Pval,
                                      qval = data$Qval)
    
    checkEquals('PlotLODR', analysis(anaquin))
    checkEquals(row.names(data), seqs(anaquin))
    checkEquals(data$ExpLFC, ratio(anaquin))
    checkEquals(data$Mean, measured(anaquin)) 
    checkEquals(data$Pval, pval(anaquin))    
    checkEquals(data$Qval, qval(anaquin))
}

test.AnaquinData_PlotLinear_1()
#test.AnaquinData_PlotLinear_2()
test.AnaquinData_PlotROC_1()
test.AnaquinData_PlotLogistic_1()
#test.AnaquinData_PlotLogistic_2()
test.AnaquinData_PlotLODR_1()
